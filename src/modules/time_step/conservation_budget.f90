      module conservation_budget
      ! THE TERMS OF THE FINITE-VOLUME CONSERVATION ROWS OF ONE EVALUATION,
      ! written out at the point where the hydrodynamic rows are assembled,
      ! so that the discrete balance of each cell can be rebuilt by a reader
      ! that shares no arithmetic with this code.
      !
      ! WHAT IT IS FOR.  The stationary residual of a cell is a difference of
      ! terms of very different size, and the files a run leaves carry the
      ! residual and the mass face flux alone: neither the momentum and
      ! energy faces, nor the geometry, nor the sources, nor the
      ! gravitational work the energy flux difference carries inside it.  A
      ! budget cannot be rebuilt from them.  This writer exports every term
      ! of the three rows of every cell at the precision of the double the
      ! assembly holds, together with the geometry, the potential and the
      ! statement of which assembly produced them, and nothing else.  It
      ! decides nothing and changes nothing: dF, S and R are written as the
      ! assembly formed them.
      !
      ! DEFAULT OFF.  It writes only when EXHALE_CONSERVATION_BUDGET is set
      ! to a positive integer, which is the number of assemblies to export;
      ! `1` exports the first one.  Absent, empty, `0` or unreadable, nothing
      ! is written, no file is opened and no branch of the assembly changes,
      ! so a production run is the run an unarmed build gives.
      !
      ! THE FOUR MOMENTUM BRANCHES.  The row the assembly builds is not the
      ! same expression in the four combinations of reconstruction and
      ! balance option, and a reader that assumes one of them is wrong in
      ! the other three (A+ / A- the face areas, V the cell volume, F the
      ! stored numerical momentum flux, p the face pressure, q the face
      ! pressure measured from the cell's own hydrostatic equilibrium):
      !
      !   PLM,   ordinary       dF = (A+ F+ - A- F-)/V
      !                         S  = the weight + (A+ - A-) p_c/V
      !                         and F carries the face pressure (Phys_flux)
      !   WENO3, ordinary       dF = (A+ F+ - A- F-)/V + (p+ - p-)/dr
      !                         S  = the weight
      !   PLM,   well balanced  dF = (A+ (F+ + q_up+) - A- (F- + q_dn-))/V
      !                         S  = 0
      !   WENO3, well balanced  dF = (A+ F+ - A- F-)/V
      !                              + (q_up+ - q_dn-)/dr
      !                         S  = 0
      !
      ! The branch in force is named in the header and every quantity the
      ! four need is exported, so no reader has to infer it.
      !
      ! THE GRAVITATIONAL WORK OF THE ENERGY ROW IS INSIDE THE FLUX
      ! DIFFERENCE, not in the source: S(3) = 0 in every branch (Source.f90)
      ! and the work rides on the MASS face flux inside dF(3) (RK_rhs),
      !
      !   dF3p     = A+ F_mass+ (phi_i+ - phi_c) - A- F_mass- (phi_i- - phi_c)
      !   dF(3,j)  = (A+ F_E+ - A- F_E- + dF3p)/V
      !
      ! so the column grav_work_over_volume below is dF3p/V and it enters
      ! dF_E with a POSITIVE sign.  A reader that adds a second
      ! gravitational work term to the energy budget has counted it twice.
      !
      ! SPHERICAL NORMALIZATION.  The assembly carries the face area as
      ! r_edg^2 and the cell volume as (r+^3 - r-^3)/3, both with the common
      ! solid-angle factor 4 pi OMITTED, and the export carries the same two
      ! numbers.  The factor cancels in every row, so a reader must omit it
      ! as well rather than restore it on one side of a balance.
      !
      ! GHOST CELLS ARE INCLUDED, and are marked.  The lowest ghost cell has
      ! no lower face, carries no equation and its row is zero by
      ! construction, so it is not exported at all; the remaining cells of
      ! the padded range are, with a flag saying which are physical.  A
      ! budget is taken over the physical cells; the ghosts are there so
      ! that the face a physical cell shares with a ghost can be read from
      ! both sides.
      !
      ! PRECISION.  Every real is written ES25.16E3, seventeen significant
      ! decimal digits, which round-trips a binary64.  That is the precision
      ! of the DOUBLE the assembly returns, and it is not the precision the
      ! assembly worked in: the kind-generic rows may form the difference in
      ! a wider kind and convert each returned array separately
      ! (hydrodynamic_rows), so a budget rebuilt from these faces can differ
      ! from the exported dF by the rounding of that conversion.  The header
      ! states the assembly, so a reader can tell the two apart.

      use global_parameters
      use grid_construction, only: spherical_cell_volume
      use RK_integration, only: face_flux, face_p, face_q_up, face_q_dn, &
                                momentum_ram_divergence,                 &
                                momentum_pressure_gradient,              &
                                momentum_gravity,                        &
                                equilibrium_pressure_force
      use utils, only: write_provenance_header,                          &
                       write_coupling_state_header,                      &
                       file_rolling_checksum
      use viscous_conduction, only: transport_active
      use hydrodynamic_rows, only: ROWS_PRODUCTION, ROWS_QUADRUPLE,      &
                                   ROWS_GENERIC_DOUBLE
      use ionization_equilibrium, only: ieq_sweep_state_kind,            &
                                        ieq_state_marching,             &
                                        ieq_state_steady_iterate,       &
                                        ieq_state_steady_candidate

      implicit none
      private
      public :: conservation_budget_exports_left,                        &
                write_conservation_budget_terms

      ! How many further assemblies are to be exported. Negative means the
      ! environment has not been read yet; zero means nothing more is
      ! written. It is read once and then counted down, so the cost of an
      ! unarmed run is one getenv for the whole run.
      integer, save :: exports_left = -1
      ! Index of the export about to be written, so that the file name and
      ! the header of each one say which assembly of the process it is.
      integer, save :: export_index = 0

      contains

      ! ------------------------------------------------------!

      logical function conservation_budget_exports_left() result(more)
      ! Whether another assembly is to be exported.
      !
      ! EXHALE_CONSERVATION_BUDGET = <n>, a positive integer: export the
      ! next n assemblies of the stationary residual and then stop. Anything
      ! else, including an absent or empty value, exports nothing. The value
      ! is read at the first call alone: an environment that changes inside
      ! a process would otherwise make two assemblies of one run export
      ! under different settings with nothing saying so.
      character(len=64) :: env
      integer :: ios, want
      if (exports_left .lt. 0) then
         env = ''
         call get_environment_variable('EXHALE_CONSERVATION_BUDGET', env)
         exports_left = 0
         if (len_trim(env) .gt. 0) then
            read(env, *, iostat=ios) want
            if (ios .eq. 0 .and. want .gt. 0) exports_left = want
         endif
      endif
      more = (exports_left .gt. 0)
      end function conservation_budget_exports_left

      ! ------------------------------------------------------!

      subroutine write_conservation_budget_terms(u, dF, S, Rrow,         &
                                                 heat, cool,             &
                                                 Smom, Sene, rows_kind)
      ! Write one assembly's complete row terms to
      ! output/conservation_budget_<nnnn>.txt.
      !
      ! The arguments are the arrays the assembly has just formed, in the
      ! assembly's own convention (steady_residual):
      !
      !   Rrow(1) = dF(1) - S(1)
      !   Rrow(2) = dF(2) - S(2) - Smom
      !   Rrow(3) = dF(3) - S(3) - (heat - cool) - Sene
      !
      ! and the face data is the module state the same assembly stored.
      ! Nothing is recomputed here except the geometry, which is a function
      ! of the grid alone.
      !
      ! The residual dummy is NOT named R: the grid radius r of
      ! global_parameters is in scope here and Fortran matches names
      ! without regard to case, so an R here hides it.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u, dF, S, Rrow
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: heat, cool, Smom, Sene
      integer, intent(in) :: rows_kind

      character(len=64) :: fname
      character(len=32) :: kind_name, state_name, branch_name, recon_name
      integer :: uu, ios, j
      real*8  :: dr, rp, rm, dAp, dAm, dV, dF3p
      real*8  :: q_up_hi, q_dn_lo, ram, pgr, grv, epf
      logical :: is_plm, wb, tr, have_terms, have_epf
      real*8  :: not_applicable

      if (exports_left .le. 0) return
      exports_left = exports_left - 1
      export_index = export_index + 1

      ! A quantity this branch does not use is written as a NaN rather than
      ! as a zero: the face departures hold whatever an earlier call left
      ! when the well-balanced option is off, and a zero there would read as
      ! a measured zero.
      not_applicable = quiet_nan_value()

      is_plm = assembled_reconstruction_is_plm()
      wb     = well_balanced
      tr     = transport_active()
      have_terms = allocated(momentum_ram_divergence)
      have_epf   = wb .and. allocated(equilibrium_pressure_force)

      if (is_plm) then
         recon_name = 'PLM'
      else
         recon_name = 'WENO3'
      endif
      if (wb) then
         if (is_plm) then
            branch_name = 'plm_well_balanced'
         else
            branch_name = 'weno3_well_balanced'
         endif
      else
         if (is_plm) then
            branch_name = 'plm_ordinary'
         else
            branch_name = 'weno3_ordinary'
         endif
      endif

      select case (rows_kind)
      case (ROWS_PRODUCTION)
         kind_name = 'production_double'
      case (ROWS_QUADRUPLE)
         kind_name = 'kind_generic_quadruple'
      case (ROWS_GENERIC_DOUBLE)
         kind_name = 'kind_generic_double'
      case default
         kind_name = 'unrecognized'
      end select

      select case (ieq_sweep_state_kind)
      case (ieq_state_marching)
         state_name = 'marching'
      case (ieq_state_steady_iterate)
         state_name = 'steady_iterate'
      case (ieq_state_steady_candidate)
         state_name = 'steady_candidate'
      case default
         state_name = 'unrecognized'
      end select

      write(fname,'(A,I4.4,A)') 'output/conservation_budget_',           &
                                export_index, '.txt'
      open(newunit=uu, file=trim(fname), status='replace',               &
           action='write', iostat=ios)
      if (ios .ne. 0) then
         write(*,'(A)') ' (conservation_budget) cannot open '//          &
              trim(fname)//'; nothing exported'
         return
      endif

      write(uu,'(A)') '# EXHALE conservation_budget schema 1'
      write(uu,'(A,I0)') '# export index in this process: ', export_index
      call write_provenance_header(uu)
      call write_coupling_state_header(uu)
      write(uu,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '# rows ', N+2*Ng-1,        &
           ': ', Ng-1, ' lower ghost cells and ', Ng,                    &
           ' upper ghost cells; physical cells are rows ', Ng,           &
           ' to ', Ng+N-1
      write(uu,'(A)') '# the lowest ghost cell has no lower face,'//     &
           ' carries no equation and is not exported'

      ! WHICH ASSEMBLY PRODUCED THESE ROWS.
      write(uu,'(A)') '# assembly rows_kind='//trim(kind_name)//         &
           ' ieq_sweep_state_kind='//trim(state_name)//                  &
           ' reconstruction='//trim(recon_name)//                        &
           ' momentum_branch='//trim(branch_name)
      write(uu,'(A,L1,A,L1,A,L1,A,L1,A,L1)')                             &
           '# flags well_balanced=', wb,                                 &
           ' transport_active=', tr,                                     &
           ' use_plm=', use_plm, ' use_weno3=', use_weno3,               &
           ' recon_lambda_on=', recon_lambda_on
      write(uu,'(A,ES25.16E3)') '# recon_lambda ', recon_lambda
      write(uu,'(A)') '# rec_method '//trim(rec_method)//                &
           ' (the flag pair does not name the branch while the'//        &
           ' continuation is armed; the branch above does)'

      ! THE STATE THESE ROWS BELONG TO.  The conserved variables are
      ! exported column by column, and the restart pair this run loaded is
      ! identified by the rolling checksum of its bytes.
      write(uu,'(A,I0,A,I0)') '# state restart_ck_hydro=',               &
           file_rolling_checksum('output/Hydro_ioniz_IC.txt'),           &
           ' restart_ck_species=',                                       &
           file_rolling_checksum('output/Ion_species_IC.txt')

      ! THE NORMALIZATION.  Every column is in code units; these are the
      ! factors that take them to cgs.
      write(uu,'(A)') '# units: all columns are in code units.'//        &
           ' Lengths R0, mass density n0*mu, velocity v0,'//             &
           ' time t_s = R0/v0, pressure and energy density'//            &
           ' p0 = n0*mu*v0^2, potential v0^2, row rate p0/t_s'//         &
           ' for the energy row and n0*mu/t_s, n0*mu*v0/t_s for'//       &
           ' the mass and momentum rows'
      write(uu,'(A,5(1X,ES25.16E3))') '# normalization R0[cm] n0[cm-3]'//&
           ' mu[g] v0[cm/s] T0[K]', R0, n0, mu, v0, T0
      write(uu,'(A,3(1X,ES25.16E3))') '# normalization t_s[s] p0[cgs]'// &
           ' b0', t_s, p0, b0
      write(uu,'(A)') '# geometry: area = r_edg^2 and'//                 &
           ' volume = (r+^3 - r-^3)/3, the common solid-angle factor'//  &
           ' 4*pi OMITTED from both, as the assembly omits it'

      ! WHAT EACH QUANTITY IS.
      write(uu,'(A)') '# kinds: face_* are FLUX DENSITIES'//             &
           ' (transport per unit area per unit time);'//                 &
           ' area_*, volume, dr, r_* are GEOMETRY;'//                    &
           ' phi_* is a POTENTIAL;'//                                    &
           ' dF_*, S_*, R_*, heat, cool, Smom, Sene,'//                  &
           ' grav_work_over_volume and the three momentum terms are'//   &
           ' CONTRIBUTIONS ALREADY DIVIDED BY THE CELL VOLUME'//         &
           ' (rates of change of a conserved density)'
      write(uu,'(A)') '# grav_work_over_volume is dF3p/V and it is'//    &
           ' ALREADY INSIDE dF_energy with a POSITIVE sign;'//           &
           ' S_energy is zero in every branch, so no second'//           &
           ' gravitational work term exists to subtract'
      write(uu,'(A)') '# face_q_up_hi and face_q_dn_lo are the face'//   &
           ' pressures measured from the cell equilibria of the two'//   &
           ' sides; they are written only under the well-balanced'//     &
           ' option and are NaN otherwise, because the arrays then'//    &
           ' hold whatever an earlier call left'
      write(uu,'(A)') '# the two values of one face differ by the'//     &
           ' equilibrium pressure jump of its two cells by'//            &
           ' construction, so they do not cancel across an interior'//   &
           ' face and their failure to cancel is not a defect'
      write(uu,'(A)') '# Sene is the COMBINED transport energy source'// &
           ' w*F_mu + q_mu + conduction, as viscous_conduction_sources'//&
           ' returns it; the three are not separated at this point'
      write(uu,'(A)') '# momentum_ram, momentum_pressure and'//          &
           ' momentum_gravity are the production attribution of the'//   &
           ' momentum row, not inputs to the identity; their sum is'//   &
           ' dF_momentum - S_momentum'

      write(uu,'(A)') '# columns j physical r_cell r_face_lo r_face_hi'//&
           ' area_lo area_hi volume dr'//                                &
           ' phi_face_lo phi_face_hi phi_cell'//                         &
           ' face_mass_lo face_mass_hi'//                                &
           ' face_momentum_lo face_momentum_hi'//                        &
           ' face_energy_lo face_energy_hi'//                            &
           ' face_p_lo face_p_hi face_q_dn_lo face_q_up_hi'//            &
           ' grav_work_over_volume'//                                    &
           ' rho momentum_density energy_density'//                      &
           ' dF_mass dF_momentum dF_energy'//                            &
           ' S_mass S_momentum S_energy heat cool Smom Sene'//           &
           ' R_mass R_momentum R_energy'//                               &
           ' momentum_ram momentum_pressure momentum_gravity'//          &
           ' equilibrium_pressure_force'

      do j = 2-Ng, N+Ng
         dr  = dr_j(j)
         rp  = r_edg(j)
         rm  = r_edg(j-1)
         dAp = rp*rp
         dAm = rm*rm
         dV  = spherical_cell_volume(j)
         dF3p = dAp*face_flux(1,j)*(Gphi_i(j)   - Gphi_c(j))             &
              - dAm*face_flux(1,j-1)*(Gphi_i(j-1) - Gphi_c(j))
         if (wb) then
            q_up_hi = face_q_up(j)
            q_dn_lo = face_q_dn(j-1)
         else
            q_up_hi = not_applicable
            q_dn_lo = not_applicable
         endif
         if (have_terms) then
            ram = momentum_ram_divergence(j)
            pgr = momentum_pressure_gradient(j)
            grv = momentum_gravity(j)
         else
            ram = not_applicable
            pgr = not_applicable
            grv = not_applicable
         endif
         if (have_epf) then
            epf = equilibrium_pressure_force(j)
         else
            epf = not_applicable
         endif
         write(uu,'(1X,I6,1X,I2,41(1X,ES25.16E3))')                      &
              j, merge(1, 0, j .ge. 1 .and. j .le. N),                   &
              r(j), rm, rp, dAm, dAp, dV, dr,                            &
              Gphi_i(j-1), Gphi_i(j), Gphi_c(j),                         &
              face_flux(1,j-1), face_flux(1,j),                          &
              face_flux(2,j-1), face_flux(2,j),                          &
              face_flux(3,j-1), face_flux(3,j),                          &
              face_p(j-1), face_p(j), q_dn_lo, q_up_hi,                  &
              dF3p/dV,                                                   &
              u(1,j), u(2,j), u(3,j),                                    &
              dF(1,j), dF(2,j), dF(3,j),                                 &
              S(1,j), S(2,j), S(3,j),                                    &
              heat(j), cool(j), Smom(j), Sene(j),                        &
              Rrow(1,j), Rrow(2,j), Rrow(3,j),                           &
              ram, pgr, grv, epf
      enddo

      close(uu)
      write(*,'(A)') ' (conservation_budget) wrote '//trim(fname)//      &
           ', assembly '//trim(kind_name)//', momentum branch '//        &
           trim(branch_name)
      if (exports_left .eq. 0)                                           &
         write(*,'(A)') ' (conservation_budget) the requested number'//  &
              ' of exports has been written; no further assembly'//      &
              ' will be exported'

      end subroutine write_conservation_budget_terms

      ! ------------------------------------------------------!

      real*8 function quiet_nan_value() result(x)
      ! A quiet NaN, for a column this branch does not define. It is formed
      ! from the bit pattern rather than from 0/0 so that no floating-point
      ! exception is raised in a build that traps them, and written out
      ! rather than taken from ieee_arithmetic so that the generated module
      ! dependency graph stays over the source tree, as finite_real
      ! (ionization_equilibrium.f90) is.
      integer*8, parameter :: qnan_bits = 9221120237041090560_8
      real*8 :: mold
      mold = 0.0d0
      x = transfer(qnan_bits, mold)
      end function quiet_nan_value

      end module conservation_budget
