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
      !
      ! THE FACE ATTRIBUTION OF THE MASS FLUX (schema 3).  Each row also
      ! carries, for the face at its r_face_hi, the two face states the
      ! Riemann problem of that face was solved on, the well-balanced
      ! pressure data (wb_dp_eq, wb_dev_L, wb_dev_R and the jump
      !     dp_WB = wb_dp_eq + wb_dev_R - wb_dev_L
      ! the Roe flux uses in place of p_R - p_L), the central transport
      ! 0.5 (F_L + F_R) of the mass flux, and the part of the numerical mass
      ! flux the pressure jump carries.  At a Mach number M the Roe mass
      ! flux holds, besides the transport, about -dp_WB/(2 c): a pressure
      ! mismatch of order M of the pressure is as large as the wind itself
      ! (md/atomic_heh97_reduced_xuv_20260929_review.md, section 3.3).  The
      ! part is formed as the DIFFERENCE of two calls of the production
      ! Num_flux on the same two states, one with the production jump and
      ! one with the three well-balanced arguments set to zero.  It is an
      ! attribution of the flux the assembly used, not a flux: the
      ! zero-jump call is no discretization of anything, since removing the
      ! jump removes legitimate pressure-driven transport as well.  The flux
      ! function is the production one; nothing of it is written again here.
      ! The two calls are measurements: the counters of the flux branches
      ! Num_flux advances are restored after them, and they write nothing
      ! the assembly or any later evaluation reads.
      !
      ! THE BASE BLOCK.  The header carries the lower boundary as the
      ! assembly saw it: the reservoir at the level, the interior state
      ! the outgoing characteristic read at the base face (the Wi of
      ! characteristic_base_face_state: cell 1's own constant-density
      ! equilibrium at the face with the well-balanced reconstruction, its
      ! hydrostatic isentrope without it; formed by the same routine,
      ! interior_state_at_level_face, on the inputs the boundary was derived
      ! from and checked bitwise against the velocity and temperature the
      ! boundary stored), the characteristic face velocity, the two states
      ! of face 0, and the numerical base mass flow.
      !
      ! THE CONTINUITY ROW (schema 3).  Each row also carries |R_mass|, the
      ! scale the stationary solve divides it by (residual_row_scale, passed
      ! in by the assembly), their ratio, and the signed running sum of
      ! V R_mass over the physical cells, which telescopes to
      ! A_{j+1/2} F_{j+1/2} - A_{1/2} F_{1/2}: the base inflow and every
      ! face after it can be read off it.

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
      use Numerical_Fluxes, only: Num_flux, Phys_flux,                   &
                                  n_faces_roe_hlle, n_faces_llf
      use Reconstruction_step, only: wb_dp_eq, wb_dev_L, wb_dev_R
      use BC_Apply, only: base_face_W, bc_W_hold, bc_npart1_hold,        &
                          base_boundary_cache_is_current
      use base_boundary, only: interior_state_at_level_face,             &
                               reservoir_state_at_level,                 &
                               wind_window_mass_flux,                    &
                               base_face_vi_last, base_face_vb_last,     &
                               base_face_T_i_last
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
                write_conservation_budget_terms,                         &
                conservation_budget_request_next

      ! How many further assemblies are to be exported. Negative means the
      ! environment has not been read yet; zero means nothing more is
      ! written. It is read once and then counted down, so the cost of an
      ! unarmed run is one getenv for the whole run.
      integer, save :: exports_left = -1
      ! Index of the export about to be written, so that the file name and
      ! the header of each one say which assembly of the process it is.
      integer, save :: export_index = 0
      ! The name of the state the next export belongs to, set by a
      ! stage-labeled request and cleared by the export that consumes it.
      character(len=96), save :: pending_label = ''

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

      subroutine conservation_budget_request_next(label)
      ! EXPORT THE NEXT ASSEMBLY UNDER A NAME, whatever the environment
      ! asked for.  The stage export of the stationary alternation
      ! (EXHALE_STAGE_EXPORT in EXHALE_main.f90) calls this at the state it
      ! names, so that the export is identified by the state it was
      ! assembled on and not by a count of assemblies, which finite-
      ! difference probes and trial evaluations also consume.  The next
      ! assembly after the request is the one exported; the label goes into
      ! its header.  It adds to any count EXHALE_CONSERVATION_BUDGET set.
      character(len=*), intent(in) :: label
      logical :: read_once
      read_once = conservation_budget_exports_left()
      exports_left = exports_left + 1
      pending_label = label
      end subroutine conservation_budget_request_next

      ! ------------------------------------------------------!

      subroutine write_conservation_budget_terms(u, WL, WR, dF, S, Rrow, &
                                                 heat, cool,             &
                                                 Smom, Sene, Sidf,       &
                                                 continuity_scale,       &
                                                 rows_kind)
      ! Write one assembly's complete row terms to
      ! output/conservation_budget_<nnnn>.txt.
      !
      ! The arguments are the arrays the assembly has just formed, in the
      ! assembly's own convention (steady_residual):
      !
      !   Rrow(1) = dF(1) - S(1)
      !   Rrow(2) = dF(2) - S(2) - Smom
      !   Rrow(3) = dF(3) - S(3) - (heat - cool) - Sene + Sidf
      !
      ! with Sidf the divergence of the enthalpy flux of the element and
      ! carrier fluxes (zero unless He_diffusion moves the elements or a
      ! molecular carrier is transported).
      !
      ! and the face data is the module state the same assembly stored.
      ! WL(:,f), WR(:,f) are the left and right states of the face f at
      ! r_edg(f) that the assembly's Riemann problems were solved on, and
      ! continuity_scale(j) is residual_row_scale(1, j, u) of the same
      ! state.  Nothing is recomputed here except the geometry, which is a
      ! function of the grid alone, and the attribution of the face mass
      ! flux and the continued base state described in the module header.
      !
      ! The residual dummy is NOT named R: the grid radius r of
      ! global_parameters is in scope here and Fortran matches names
      ! without regard to case, so an R here hides it.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u, dF, S, Rrow
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: WL, WR
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: heat, cool, Smom, Sene
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: Sidf
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: continuity_scale
      integer, intent(in) :: rows_kind

      character(len=64) :: fname
      character(len=32) :: kind_name, state_name, branch_name, recon_name
      integer :: uu, ios, j
      real*8  :: dr, rp, rm, dAp, dAm, dV, dF3p
      real*8  :: q_up_hi, q_dn_lo, ram, pgr, grv, epf
      logical :: is_plm, wb, tr, have_terms, have_epf
      real*8  :: not_applicable
      ! The attribution of the face mass flux, one column set for the face
      ! at r_face_hi of every exported row (the module header).
      integer, parameter :: n_face_cols = 13
      real*8, dimension(n_face_cols,1-Ng:N+Ng) :: face_cols
      real*8  :: flux_left(3), flux_right(3), flux_jump(3), flux_nojump(3)
      real*8  :: p_face_unused, q_up_unused, q_dn_unused
      integer :: jr_cell, n_roe_hlle_saved, n_llf_saved
      logical :: rows_blended, traces_defined, production_cols_defined
      logical :: wb_cols_defined, jump_col_defined
      ! The continuity row's running budget over the physical cells.
      real*8, dimension(1-Ng:N+Ng) :: cum_mass
      real*8  :: cum_sum
      ! The base block.
      real*8  :: rho_res_b, t_res_b, p_res_b
      real*8  :: nhat_cont, t_cell1, rho_cont, t_cont, p_cont, v_cont
      real*8  :: w_cont(3)
      real*8  :: f_wind_b, d_window_b, area0, base_mass_flow
      logical :: have_wind_b, cont_defined, cont_matches, cache_current

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

      ! WHICH FACE COLUMNS THIS ASSEMBLY DEFINES.  On the open interval of
      ! the reconstruction continuation, WL and WR are the lambda-weighted
      ! blend of the PLM and WENO3 face states and no Riemann problem was
      ! solved on them, so no face column is defined.  The kind-generic rows
      ! return their face states rounded to double but solve the Riemann
      ! problem with their own text, not with Num_flux, and do not fill the
      ! well-balanced face data of Reconstruction_step; there the traces are
      ! written and every column formed with Num_flux or from that face data
      ! is not.  The pressure-jump part is a statement about the Roe flux
      ! with the well-balanced jump: the HLLC flux takes the jump into its
      ! contact speed and the LLF flux does not read it at all, and without
      ! the well-balanced option there is no such jump.
      rows_blended = recon_lambda_on .and. recon_lambda .gt. 0.0d0       &
                     .and. recon_lambda .lt. 1.0d0
      traces_defined          = .not. rows_blended
      production_cols_defined = traces_defined .and.                     &
                                rows_kind .eq. ROWS_PRODUCTION
      wb_cols_defined  = production_cols_defined .and. wb .and.          &
                         allocated(wb_dp_eq)
      jump_col_defined = wb_cols_defined .and. trim(flux) .eq. 'ROE'

      face_cols = not_applicable
      n_roe_hlle_saved = n_faces_roe_hlle
      n_llf_saved      = n_faces_llf
      do j = 2-Ng, N+Ng
         jr_cell = min(j+1, N+Ng)
         if (traces_defined) then
            face_cols(1:3,j) = WL(:,j)
            face_cols(4:6,j) = WR(:,j)
         endif
         if (wb_cols_defined) then
            face_cols(7,j)  = wb_dp_eq(j)
            face_cols(8,j)  = wb_dev_L(j)
            face_cols(9,j)  = wb_dev_R(j)
            ! The jump exactly as Num_flux forms it (dp_wb there).
            face_cols(10,j) = wb_dp_eq(j) + wb_dev_R(j) - wb_dev_L(j)
         endif
         if (production_cols_defined) then
            ! The same two cells RK_rhs names for this face.
            call Phys_flux(WL(:,j), flux_left,  j)
            call Phys_flux(WR(:,j), flux_right, jr_cell)
            face_cols(11,j) = 0.5d0*(flux_left(1) + flux_right(1))
            if (wb) then
               call Num_flux(WL(:,j), WR(:,j), flux_jump, p_face_unused, &
                             j, jr_cell, wb_dev_L(j), wb_dev_R(j),       &
                             wb_dp_eq(j), q_up_unused, q_dn_unused)
            else
               call Num_flux(WL(:,j), WR(:,j), flux_jump, p_face_unused, &
                             j, jr_cell)
            endif
            face_cols(13,j) = flux_jump(1)
            if (jump_col_defined) then
               call Num_flux(WL(:,j), WR(:,j), flux_nojump,              &
                             p_face_unused, j, jr_cell,                  &
                             0.0d0, 0.0d0, 0.0d0,                        &
                             q_up_unused, q_dn_unused)
               face_cols(12,j) = flux_jump(1) - flux_nojump(1)
            endif
         endif
      enddo
      n_faces_roe_hlle = n_roe_hlle_saved
      n_faces_llf      = n_llf_saved

      ! THE SIGNED RUNNING MASS BUDGET, sum over physical k <= j of V_k
      ! R_mass,k in the order of the cells.
      cum_mass = not_applicable
      cum_sum  = 0.0d0
      do j = 1, N
         cum_sum     = cum_sum + spherical_cell_volume(j)*Rrow(1,j)
         cum_mass(j) = cum_sum
      enddo

      ! THE BASE BLOCK.  The interior state at the face is formed from the
      ! primitive state and the cell-1 particle count the boundary was last
      ! derived from (BC_Apply), by the routine the boundary itself calls;
      ! the boundary does the same only where cell 1 is an admissible gas
      ! state, and so does this.
      call reservoir_state_at_level(rho_res_b, t_res_b, p_res_b)
      cache_current = base_boundary_cache_is_current(u)
      cont_defined  = .false.
      cont_matches  = .false.
      rho_cont = not_applicable;  t_cont = not_applicable
      p_cont   = not_applicable;  v_cont = not_applicable
      f_wind_b = not_applicable;  d_window_b = not_applicable
      have_wind_b = .false.
      if (allocated(bc_W_hold)) then
         if (bc_W_hold(1,1) .gt. 0.0d0 .and. bc_W_hold(3,1) .gt. 0.0d0   &
             .and. bc_npart1_hold .gt. 0.0d0) then
            nhat_cont = bc_npart1_hold/bc_W_hold(1,1)
            t_cell1   = bc_W_hold(3,1)/bc_npart1_hold
            call interior_state_at_level_face(bc_W_hold(:,1), nhat_cont, &
                                              t_cell1, w_cont, t_cont)
            rho_cont = w_cont(1)
            v_cont   = w_cont(2)
            p_cont   = w_cont(3)
            cont_defined = .true.
            cont_matches = (t_cont .eq. base_face_T_i_last) .and.        &
                           (v_cont .eq. base_face_vi_last)
         endif
         call wind_window_mass_flux(bc_W_hold, f_wind_b, d_window_b,     &
                                    have_wind_b)
      endif
      area0          = r_edg(0)*r_edg(0)
      base_mass_flow = area0*face_flux(1,0)

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

      write(uu,'(A)') '# EXHALE conservation_budget schema 3'
      write(uu,'(A,I0)') '# export index in this process: ', export_index
      if (len_trim(pending_label) .gt. 0) then
         write(uu,'(A,A)') '# stage: ', trim(pending_label)
         pending_label = ''
      endif
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
           ' dF_*, S_*, R_*, heat, cool, Smom, Sene, Sidf,'//            &
           ' grav_work_over_volume and the three momentum terms are'//   &
           ' CONTRIBUTIONS ALREADY DIVIDED BY THE CELL VOLUME'//         &
           ' (rates of change of a conserved density);'//                &
           ' trace_* are FACE STATES (density, velocity, pressure);'//   &
           ' wb_dp_eq_hi, wb_dev_*_hi and dp_WB_hi are PRESSURES;'//     &
           ' central_mass_hi, pressure_jump_mass_hi and'//               &
           ' num_flux_mass_recomputed_hi are FLUX DENSITIES like'//      &
           ' face_*; R_mass_abs and mass_row_scale are rates like'//     &
           ' R_mass; R_mass_normalized is DIMENSIONLESS;'//              &
           ' cumulative_mass_budget is a MASS FLOW (area times flux'//   &
           ' density, the 4*pi omitted as for the geometry)'
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
      write(uu,'(A)') '# Sidf is the divergence of the interdiffusion'// &
           ' enthalpy flux q_d = sum_s h_s J_s of the element and'//  &
           ' carrier fluxes, ADDED to R_energy: R_energy = dF_energy'//  &
           ' - S_energy - (heat - cool) - Sene + Sidf; zero unless'//    &
           ' He_diffusion moves the elements or a molecular carrier'//   &
           ' is transported (schema 2 adds this column)'
      write(uu,'(A)') '# Sidf of cell 1 is zero for the ELEMENT fluxes:'//&
           ' cell 1 is the element reservoir, whose steady composition'//&
           ' makes the face-0 supply equal the face-1 flux (A_0 q_0 ='//  &
           ' A_1 q_1), so its sum over cells 1..N is A_N q_N - A_1 q_1,'//&
           ' the outflow minus the reservoir energy supply A_1 q_1'//     &
           ' (binary_element_diffusion, interdiffusion_enthalpy_divergence)'
      write(uu,'(A)') '# momentum_ram, momentum_pressure and'//          &
           ' momentum_gravity are the production attribution of the'//   &
           ' momentum row, not inputs to the identity; their sum is'//   &
           ' dF_momentum - S_momentum'

      ! THE FACE ATTRIBUTION AND THE CONTINUITY ROW (schema 3).
      write(uu,'(A)') '# schema 3 appends, for the face at r_face_hi'//  &
           ' of each row (face j, between cells j and j+1; the lowest'// &
           ' face, below the first exported row, has none):'//           &
           ' trace_{rho,v,p}_{L,R}_hi, the left and right face'//        &
           ' states the Riemann problem of that face was solved on'//    &
           ' (face 0: L is the lower boundary state, R the'//            &
           ' reconstruction of cell 1); wb_dp_eq_hi, wb_dev_L_hi,'//     &
           ' wb_dev_R_hi and dp_WB_hi = wb_dp_eq + wb_dev_R - wb_dev_L,'//&
           ' the well-balanced pressure data and the jump the flux'//    &
           ' uses in place of p_R - p_L; central_mass_hi ='//            &
           ' 0.5 (F_L + F_R), the mass component of the physical'//      &
           ' flux (Phys_flux) at the two states; pressure_jump_mass_hi'//&
           ' = Num_flux(production jump) - Num_flux(jump set to'//       &
           ' zero) at the same two states, mass component;'//            &
           ' num_flux_mass_recomputed_hi, the production Num_flux at'//  &
           ' the two states, mass component'
      write(uu,'(A)') '# pressure_jump_mass_hi is an ATTRIBUTION of'//   &
           ' the flux the assembly used, NOT a flux: the zero-jump'//    &
           ' call is no discretization, since removing the jump'//       &
           ' removes legitimate pressure-driven transport as well;'//     &
           ' face_mass_hi - central_mass_hi - pressure_jump_mass_hi is'//&
           ' the remaining Roe dissipation of the mass row; on a face'// &
           ' where the Roe flux fell back to HLLE the part is exactly'// &
           ' zero'
      write(uu,'(A)') '# num_flux_mass_recomputed_hi equals'//           &
           ' face_mass_hi bitwise unless the stagnant-layer contact'//   &
           ' dissipation (low_mach_dissipation) is added to the stored'//&
           ' flux, which Num_flux does not carry'
      write(uu,'(A,L1,A,L1,A,L1,A,L1)') '# face columns defined:'//     &
           ' traces=', traces_defined,                                   &
           ' central_and_recomputed=', production_cols_defined,          &
           ' well_balanced_data=', wb_cols_defined,                      &
           ' pressure_jump_mass=', jump_col_defined
      write(uu,'(A)') '# undefined face columns are NaN: all of them'//  &
           ' on the open interval 0 < recon_lambda < 1, whose face'//    &
           ' states are a blend no Riemann problem was solved on;'//     &
           ' all but the traces under the kind-generic rows, which'//    &
           ' solve the face with their own text and not with'//          &
           ' Num_flux; the well-balanced data and the jump part'//       &
           ' without the well-balanced option, which forms no such'//    &
           ' jump; and the jump part under HLLC (the jump enters its'//  &
           ' contact speed, not a Roe wave strength) and LLF (the'//     &
           ' flux does not read it). This run: numerical flux '//        &
           trim(flux)
      write(uu,'(A)') '# R_mass_abs = |R_mass|; mass_row_scale ='//      &
           ' residual_row_scale(1, j, u), the scale the stationary'//    &
           ' solve divides the continuity row by (max of |A F| of the'// &
           ' two faces over V); R_mass_normalized = R_mass_abs /'//      &
           ' mass_row_scale; cumulative_mass_budget = sum over'//        &
           ' physical cells k <= j of volume_k R_mass,k, which'//        &
           ' telescopes to area_hi face_mass_hi - A_0 F_0 (NaN on'//     &
           ' ghost rows)'

      ! THE BASE BLOCK.
      write(uu,'(A)') '# base block: the lower boundary of this'//       &
           ' assembly. rho,v,p in code units; the face is r_edg(0)'
      write(uu,'(A,ES25.16E3)') '# base r_face_0 ', r_edg(0)
      write(uu,'(A,3(1X,ES25.16E3))') '# base reservoir rho_res'//       &
           ' T_res p_res', rho_res_b, t_res_b, p_res_b
      write(uu,'(A,L1,A,L1,A,L1)') '# base continuation_defined=',      &
           cont_defined, ' continuation_matches_boundary=',              &
           cont_matches, ' boundary_cache_belongs_to_this_state=',       &
           cache_current
      write(uu,'(A)') '# base the interior continuation Wi is cell 1'//  &
           ' carried to the face: its own constant-density'//            &
           ' equilibrium there with the well-balanced'//                 &
           ' reconstruction, its hydrostatic isentrope without it'//     &
           ' (interior_state_at_level_face on the primitive state'//     &
           ' and particle count the boundary was derived from);'//       &
           ' continuation_matches_boundary says its T and v equal'//     &
           ' bitwise the ones characteristic_base_face_state stored'
      write(uu,'(A,4(1X,ES25.16E3))') '# base interior_continuation'//   &
           ' rho_i v_i p_i T_i', rho_cont, v_cont, p_cont, t_cont
      write(uu,'(A,ES25.16E3)') '# base characteristic_face_velocity'//  &
           ' v_b ', base_face_vb_last
      write(uu,'(A,3(1X,ES25.16E3))') '# base face_0_left_trace'//       &
           ' rho v p (the boundary state)', WL(:,0)
      write(uu,'(A,3(1X,ES25.16E3))') '# base face_0_right_trace'//      &
           ' rho v p (cell 1 reconstructed down to the face)', WR(:,0)
      write(uu,'(A,3(1X,ES25.16E3))') '# base cached_boundary_face_state'//&
           ' rho v p (BC_Apply base_face_W)', base_face_W
      write(uu,'(A,3(1X,ES25.16E3))') '# base area_0 face_mass_0'//      &
           ' area_0*face_mass_0', area0, face_flux(1,0), base_mass_flow
      write(uu,'(A,L1,2(1X,ES25.16E3))') '# base wind_window'//          &
           ' have_F=', have_wind_b, f_wind_b, d_window_b
      write(uu,'(A)') '# base wind_window is the mean r^2 rho v of the'//&
           ' cells from j_flux outward and its relative spread, as'//    &
           ' the boundary reads them (wind_window_mass_flux), for'//     &
           ' comparison with area_0*face_mass_0'

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
           ' S_mass S_momentum S_energy heat cool Smom Sene Sidf'//      &
           ' R_mass R_momentum R_energy'//                               &
           ' momentum_ram momentum_pressure momentum_gravity'//          &
           ' equilibrium_pressure_force'//                               &
           ' trace_rho_L_hi trace_v_L_hi trace_p_L_hi'//                 &
           ' trace_rho_R_hi trace_v_R_hi trace_p_R_hi'//                 &
           ' wb_dp_eq_hi wb_dev_L_hi wb_dev_R_hi dp_WB_hi'//             &
           ' central_mass_hi pressure_jump_mass_hi'//                    &
           ' num_flux_mass_recomputed_hi'//                              &
           ' R_mass_abs mass_row_scale R_mass_normalized'//              &
           ' cumulative_mass_budget'

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
         write(uu,'(1X,I6,1X,I2,59(1X,ES25.16E3))')                      &
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
              heat(j), cool(j), Smom(j), Sene(j), Sidf(j),               &
              Rrow(1,j), Rrow(2,j), Rrow(3,j),                           &
              ram, pgr, grv, epf,                                        &
              face_cols(:,j),                                            &
              abs(Rrow(1,j)), continuity_scale(j),                       &
              abs(Rrow(1,j))/continuity_scale(j),                        &
              cum_mass(j)
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
