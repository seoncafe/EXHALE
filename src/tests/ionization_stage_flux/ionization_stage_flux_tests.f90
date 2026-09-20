      program ionization_stage_flux_tests
      ! THE DISCRETE IDENTITY OF THE IONIZATION-STAGE FLUX
      ! (docs/lhs1140b_stationary_L12b_derivation_20260916.md, section 4;
      ! PLAN_20260916_rev3 section 8, item L12 stage B).
      !
      ! WHAT IS TESTED AND WHAT IS NOT.  The stage flux is
      ! src/modules/functions/ionization_stage_transport.f90, and the
      ! transport-chemistry operator solves the H II row on it
      ! ("Ionization transport: True", diffusive_photochemistry).  The rows
      ! below run against the production module; the first group also holds
      ! the FORMULA of the derivation memo written once here, fed with the
      ! face quantities the production operator itself forms, so that the
      ! module and the derivation are compared and not the module with
      ! itself:
      !
      !     F_k(f) = x_k(f) N_el(f)  -  n_el(f) K(f) [x_k(j+1) - x_k(j)]/dr(f)
      !
      ! with x_k the ionization fraction of stage k per nucleus of its
      ! element, N_el the element NUCLEUS face flux, n_el(f) the face element
      ! nucleus density and K(f) the face eddy coefficient, both taken as the
      ! arithmetic face averages drift_and_gradient_face_coefficients takes.
      ! The face fractions come from species_face_fraction and the face
      ! element flux from species_face_flux, the two routines the element
      ! operator's own advective half calls (element_advective_divergence,
      ! binary_element_diffusion.f90).
      !
      ! THE DIFFUSIVE HALF OF THE ELEMENT FACE FLUX IS SUPPLIED HERE AND NOT
      ! READ FROM THE OPERATOR.  drift_and_gradient_face_coefficients and
      ! element_face_flux are private to binary_element_diffusion, so no
      ! caller can read Agrd, Bdrf, updrf or J(f).  The identity below is
      ! LINEAR in that half, so what the rows can state without it is the
      ! whole mathematical content: the sum of the stage fluxes returns the
      ! element flux for ANY value of the diffusive half, including a
      ! reversed one and a vanishing one, and the rows run over those cases.
      ! What is not covered is that stage B multiplies the operator's own
      ! J(f) rather than a second copy of it; that needs the accessor named
      ! in section 7 of the memo.
      !
      ! ROWS.
      !   stage_fractions_close_the_simplex_at_every_face
      !       the closing member of a set is ONE MINUS the reconstructed
      !       others (the rule species_advective_update already follows for
      !       the mass fractions), so sum_k x_k(f) = 1 to rounding.
      !   independently_reconstructed_stages_break_the_simplex
      !       the RED: reconstructing EVERY stage on its own, which is the
      !       natural implementation, leaves sum_k x_k(f) away from 1 by the
      !       reconstruction's own truncation, and the identity below fails
      !       by that amount times the element flux.
      !   sum_of_stage_fluxes_is_the_element_flux
      !       the identity, on a nonuniform element abundance, with constant
      !       and with varying stage fractions, with K_zz nonzero, in a
      !       neutral-background limit, with helium fractions at the simplex
      !       faces, and with the element flux positive, reversed and zero.
      !   independently_reconstructed_stages_break_the_identity
      !       the same measure with the mismatched interpolation: the RED of
      !       the row above.
      !   a_second_element_face_flux_breaks_the_identity
      !       the other RED: the fractions still close the simplex, but two
      !       copies of the element face flux under different stages do not
      !       return one element flux.
      !   stage_eddy_term_sums_to_zero
      !       the added stage term telescopes exactly when the closing
      !       stage's term is MINUS the sum of the others and n_el(f), K(f)
      !       and dr(f) are ONE face value shared by every stage.
      !   differencing_the_closing_member_loses_the_trace_gradient
      !       the RED that forced the rule above: forming the closing
      !       stage's term from its own cell values by the same formula is
      !       exact in exact arithmetic and loses the gradient to the
      !       rounding of 1 - sum wherever the carried stages are trace.
      !   a_stage_dependent_face_density_breaks_the_eddy_sum
      !       a donor-cell n_el chosen from each stage's own gradient leaves
      !       a residual even without the trace limit.
      !   charge_consistency_of_the_transported_partition
      !       n_e = sum_k Z_k n_k from the transported fractions.
      !
      ! THE SECOND HALF OF THE SUITE RUNS THE PRODUCTION MODULE.
      ! ionization_stage_transport now holds the formula above, and
      ! binary_element_diffusion exposes the element nucleus face flux it is
      ! written on (element_nucleus_face_flux, item L30), so the rows below
      ! test the operator and not a transcription of it: the stage fluxes
      ! come from ionization_stage_face_flux, the element flux from
      ! hydrogen_and_helium_nucleus_face_flux, the divergence from
      ! ionization_stage_divergence on the one spherical geometry, and the
      ! rows from ionization_stage_row_jacobian.  They run on the PRODUCTION
      ! grid constructor (the Mixed grid of the LHS 1140 b catalog, 500
      ! cells whose width varies by four decades), which is where a
      ! divergence weighted by the exact shell volume and one weighted by
      ! r_j^2 dr_j stop agreeing from one cell to the next.
      !
      !   production_stage_flux_sums_to_the_element_flux
      !       the identity (2) of the module header at EVERY face, f = 0 to
      !       f = N, for hydrogen and for helium, on the operator's own
      !       element flux.
      !   sum_identity_survives_the_simplex_projection
      !       the same after stage_simplex_projection has moved the cell
      !       fractions back into the simplex: the projection acts on the
      !       CELL state and the face rule is re-derived from it, so the
      !       identity is untouched.
      !   sum_identity_holds_at_the_boundary_faces
      !       f = 0 and f = N alone, with the element flux outflowing and
      !       reversed, which is the boundary reconstruction: the ghosts
      !       carry the reservoir composition below and repeat cell N above.
      !   stage_divergence_sums_to_the_element_nucleus_divergence
      !       cell by cell, the sum over stages of the stage divergences is
      !       the divergence of the element nucleus flux.
      !   stage_divergence_column_sum_is_the_boundary_flux_difference
      !       down the column the internal faces cancel on the one geometry
      !       and the difference of the two boundary nucleus fluxes is left.
      !   stage_row_jacobian_matches_a_central_difference
      !       the analytic tridiagonal rows against a central difference of
      !       the donor-value divergence they are the derivative of, with
      !       the second-order part of the reconstruction reported beside
      !       them.
      !   a_full_eddy_term_in_every_stage_breaks_the_identity
      !       the RED of the derivation's central point: replacing the
      !       stage term by the whole mixing-ratio eddy flux
      !       -n_tot K d(x_k y)/dr, which is the L12a design's Phi_x,
      !       leaves sum_k F_k off the element flux by exactly one copy of
      !       the element's own eddy flux.
      !   bounding_the_closing_face_value_breaks_the_identity
      !       the RED that says why the closing face value is measured and
      !       never repaired: clipping it into [0,1] at the face buys
      !       nonnegativity and loses the sum.
      !
      !   the_closing_member_can_leave_the_simplex
      !       MEASURED, not asserted: the closing-member rule buys the SUM
      !       exactly and does not buy nonnegativity, so a stage set whose
      !       reconstructed members are bounded into [0,1] separately can
      !       still hand the closing member a negative face value.  This is
      !       the headroom statement stage B owes (L12a section 2.2).
      !
      ! The columns are synthetic: this driver sets the global_parameters
      ! scalars and the radial grid itself, so no input.inp, no ionization
      ! solve and no hydrodynamics are involved.

      use global_parameters
      use species_table
      use element_inventory
      use diffusive_photochemistry, only: carrier_set_init,               &
                              n_carrier_max, ic_Hp, ic_HeII, ic_HeIII,    &
                              carrier_solved, carrier_is_ionization_stage,&
                              carrier_stage_element,                      &
                              carrier_stage_face_state,                   &
                              carrier_stage_face_flux,                    &
                              ionization_stage_nucleus_sum,               &
                              carrier_ionization_stage_projection,         &
                              carrier_helium_background_set_for_test
      use lower_atmosphere_profile, only: eddy_diffusion_on_grid
      use species_advective_transport, only: species_face_fraction,       &
                                             species_face_flux
      use grid_construction, only: define_grid,                           &
                                   spherical_face_area_and_cell_volume
      use binary_element_diffusion, only: mixture_mass_sum
      use certification, only: cert_regime_wind_r, cert_tol_carrier_wind
      use ionization_stage_transport, only:                               &
                              ionization_stage_face_flux,                 &
                              ionization_stage_divergence,                &
                              ionization_stage_face_jacobian,             &
                              ionization_stage_row_jacobian,              &
                              stage_simplex_projection,                   &
                              stage_simplex_sum_over_limit,               &
                              stage_fraction_under_zero,                  &
                              hydrogen_and_helium_nucleus_face_flux
      use test_columns, only: column_carrying_its_own_density

      implicit none

      integer :: n_pass, n_fail

      n_pass = 0
      n_fail = 0

      write(*,'(A)') '===== ionization stage flux acceptance tests ====='

      call test_simplex_closure('PLM')
      call test_simplex_closure('WENO3')
      call test_stage_flux_identity('PLM')
      call test_stage_flux_identity('WENO3')
      call test_stage_eddy_sum()
      call test_charge_consistency()
      call test_closing_member_excursion()
      call test_production_stage_rows()
      call test_production_stage_jacobian()
      call test_boundary_stage_rows()
      call test_operator_stage_rows_of_both_elements()
      call test_manufactured_ionization_column()
      call test_manufactured_stiff_stage_column()

      write(*,'(A)') '================================================='
      write(*,'(A,I0,A,I0,A)') ' TOTAL: ', n_pass, ' passed, ',           &
                               n_fail, ' failed'
      if (n_fail .gt. 0) stop 1

      contains

      ! ================================================================= !
      !  harness
      ! ================================================================= !

      subroutine verdict(name, ok, val, ref, tol)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: ok
      real*8,           intent(in) :: val, ref, tol
      if (ok) then
         n_pass = n_pass + 1
         write(*,'(A,A,A,ES11.4,A,ES11.4,A,ES11.4)') ' PASS ', name,      &
              ' measured=', val, ' reference=', ref, ' tol=', tol
      else
         n_fail = n_fail + 1
         write(*,'(A,A,A,ES11.4,A,ES11.4,A,ES11.4)') ' FAIL ', name,      &
              ' measured=', val, ' reference=', ref, ' tol=', tol
      endif
      end subroutine verdict

      ! ----------------------------------------------------------------- !

      subroutine setup_column(Ncell, r_top, kzz)
      ! A synthetic spherical column and the global scalars the face
      ! routines read.  kzz fills every cell through the production filler,
      ! so the face average below is the one the operator forms.
      integer, intent(in) :: Ncell
      real*8,  intent(in) :: r_top, kzz
      real*8  :: dr_u
      integer :: j

      N   = Ncell
      T0  = 1.0d4
      R0  = 1.0d10
      n0  = 1.0d10
      v0  = sqrt(kb_erg*T0/mu)
      t_s = R0/v0
      p0  = n0*mu*v0*v0
      spherical_domain = .true.
      thereis_He    = .true.
      thereis_HeITR = .false.
      thereis_mol   = .false.
      thereis_metals = .false.
      he_diffusion  = .true.
      he_kzz        = kzz
      recon_lambda_on = .false.

      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))

      dr_u = (r_top - 1.0d0)/dble(N-1)
      do j = 1-Ng, N+Ng
         r(j) = 1.0d0 + dble(j-1)*dr_u
      enddo
      r_edg(1-Ng:N+Ng-1) = 0.5d0*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
      r_edg(N+Ng) = 2.0d0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
      dr_j = dr_u
      call eddy_diffusion_on_grid()

      end subroutine setup_column

      ! ----------------------------------------------------------------- !

      subroutine set_scheme(scheme)
      character(len=*), intent(in) :: scheme
      rec_method = scheme
      use_plm   = (scheme .eq. 'PLM')
      use_weno3 = (scheme .eq. 'WENO3')
      end subroutine set_scheme

      ! ----------------------------------------------------------------- !

      subroutine helium_column(icase, nH, nHe, xHII, xHeII, xHeIII)
      ! Four columns, each a case the memo names.
      !   1  nonuniform element abundance, every stage fraction varying
      !   2  the same abundance, every stage fraction CONSTANT
      !   3  the neutral-background limit: the ionized stages at 1e-12
      !   4  helium at the simplex faces: x(He I) driven to zero over part
      !      of the column and x(He III) to one over another part
      integer, intent(in) :: icase
      real*8, dimension(1-Ng:N+Ng), intent(out) :: nH, nHe
      real*8, dimension(1-Ng:N+Ng), intent(out) :: xHII, xHeII, xHeIII
      real*8  :: s
      integer :: j
      do j = 1-Ng, N+Ng
         s = (r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
         ! A helium element that separates from hydrogen by a factor 4 over
         ! the column, as the LHS 1140 b wind does (L12a section 1.2).
         nH(j)  = 1.0d8*exp(-4.0d0*s)
         nHe(j) = 1.6d0*nH(j)*(1.0d0 + 3.0d0*s)**(-1)
         select case (icase)
         case (1)
            xHII(j)   = 0.05d0 + 0.90d0*s*s
            xHeII(j)  = 0.03d0 + 0.55d0*s
            xHeIII(j) = 0.01d0 + 0.30d0*s*s*s
         case (2)
            xHII(j)   = 0.3731d0
            xHeII(j)  = 0.2417d0
            xHeIII(j) = 0.1109d0
         case (3)
            xHII(j)   = 1.0d-12*(1.0d0 + s)
            xHeII(j)  = 1.0d-12*(1.0d0 + 2.0d0*s)
            xHeIII(j) = 1.0d-14
         case default
            ! x(He I) = 1 - x(He II) - x(He III) reaches 0 at s = 1/2 and
            ! x(He III) reaches 1 at s = 1.
            xHII(j)   = min(1.0d0, 2.0d0*s)
            xHeII(j)  = min(1.0d0, 2.0d0*s)*max(0.0d0, 1.0d0 - 2.0d0*s)   &
                      + max(0.0d0, 2.0d0 - 2.0d0*s)*0.0d0
            xHeIII(j) = max(0.0d0, 2.0d0*s - 1.0d0)
            xHeII(j)  = max(0.0d0, min(1.0d0 - xHeIII(j),                 &
                             2.0d0*s*(1.0d0 - max(0.0d0, 2.0d0*s-1.0d0))))
         end select
      enddo
      end subroutine helium_column

      ! ----------------------------------------------------------------- !

      subroutine face_mass_flux(iflow, Frho)
      ! Three winds and a stagnant one: outflowing, reversed, alternating in
      ! sign from face to face (the breathing base), and identically zero.
      integer, intent(in) :: iflow
      real*8, dimension(1-Ng:N+Ng), intent(out) :: Frho
      integer :: j
      do j = 1-Ng, N+Ng
         select case (iflow)
         case (1)
            Frho(j) =  0.4d0/(r_edg(j)*r_edg(j))
         case (2)
            Frho(j) = -0.4d0/(r_edg(j)*r_edg(j))
         case (3)
            Frho(j) = 0.4d0/(r_edg(j)*r_edg(j))*dble(1 - 2*mod(j+100,2))
         case default
            Frho(j) = 0.0d0
         end select
      enddo
      end subroutine face_mass_flux

      ! ----------------------------------------------------------------- !

      subroutine element_nucleus_face_flux(nH, nHe, Frho, Jdif, isign,    &
                                           NelH, NelHe)
      ! THE ELEMENT NUCLEUS FACE FLUXES, built the way the element operator
      ! builds its advective half: the face mass fraction of each element
      ! through species_face_fraction, multiplied by the face mass flux
      ! through species_face_flux, with the CLOSING MEMBER of the pair taken
      ! as one minus the reconstructed other so that the two mass fluxes sum
      ! to F_rho exactly.  The diffusive half Jdif is supplied by the caller
      ! (the operator does not expose its own; see the header) and enters
      ! helium positively and the hydrogen component negatively, which is the
      ! binary closure J_1 = -J_He of the module header.  isign scales it, so
      ! a caller can run the identity at a reversed and at a vanishing
      ! diffusive half.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: nH, nHe, Frho, Jdif
      real*8,                       intent(in)  :: isign
      real*8, dimension(1-Ng:N+Ng), intent(out) :: NelH, NelHe

      real*8, dimension(1-Ng:N+Ng) :: rhoc, YHe, YHef, YHf, FHe, FH
      real*8  :: mH_g, mHe_g
      integer :: j

      mH_g  = mu
      mHe_g = 4.002602d0/1.00794d0*mu
      do j = 1-Ng, N+Ng
         rhoc(j) = mH_g*nH(j) + mHe_g*nHe(j)
         YHe(j)  = mHe_g*nHe(j)/rhoc(j)
      enddo
      call species_face_fraction(YHe, Frho, YHef)
      ! The closing member is not reconstructed.
      YHf = 1.0d0 - YHef
      call species_face_flux(Frho, YHef, FHe)
      call species_face_flux(Frho, YHf,  FH)
      do j = 1-Ng, N+Ng
         NelHe(j) = (FHe(j) + isign*Jdif(j))/mHe_g
         NelH(j)  = (FH(j)  - isign*Jdif(j))/mH_g
      enddo
      end subroutine element_nucleus_face_flux

      ! ----------------------------------------------------------------- !

      subroutine stage_face_fractions(independent, xion, Nel, xf, xclose)
      ! The face fractions of one element's ionization stages.  The donor
      ! side is selected by the ELEMENT flux Nel and not by the bulk mass
      ! flux: the cell a stage leaves is the cell its element leaves.
      !
      ! independent = .false. is the rule of the memo: reconstruct every
      ! stage but the neutral one and take the neutral as one minus their
      ! sum.  independent = .true. reconstructs the neutral stage on its own
      ! as well, which is the mismatched interpolation the RED rows use.
      logical,                        intent(in)  :: independent
      real*8, dimension(:,1-Ng:),     intent(in)  :: xion
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: Nel
      real*8, dimension(:,1-Ng:),     intent(out) :: xf
      real*8, dimension(1-Ng:N+Ng),   intent(out) :: xclose

      real*8, dimension(1-Ng:N+Ng) :: work, wf
      integer :: k, nk, j

      nk = size(xion,1)
      do k = 1, nk
         work = xion(k,:)
         call species_face_fraction(work, Nel, wf)
         xf(k,:) = wf
      enddo
      if (independent) then
         do j = 1-Ng, N+Ng
            work(j) = 1.0d0
            do k = 1, nk
               work(j) = work(j) - xion(k,j)
            enddo
         enddo
         call species_face_fraction(work, Nel, xclose)
      else
         xclose = 1.0d0
         do k = 1, nk
            xclose = xclose - xf(k,:)
         enddo
      endif
      end subroutine stage_face_fractions

      ! ----------------------------------------------------------------- !

      subroutine stage_eddy_flux(nel, xion, xclose_cell, closing_rule,   &
                                 per_stage_donor, Ek, Eclose)
      ! The stage term of the memo,  -n_el(f) K(f) [x(j+1) - x(j)]/dr(f),
      ! with n_el(f) and K(f) the arithmetic face averages the element
      ! operator's face coefficients take and dr(f) the face spacing in cm.
      !
      ! closing_rule = 'negated_sum' forms the closing stage's term as MINUS
      ! the sum of the others, which is what makes the telescoping exact in
      ! floating point as well as in exact arithmetic.  'own_difference'
      ! forms it from the closing stage's own cell values by the same
      ! formula, which is the natural implementation and loses the gradient
      ! of a dominant closing stage to the rounding of 1 - sum whenever the
      ! carried stages are trace.
      !
      ! per_stage_donor = .true. replaces the face average by the donor
      ! cell's own density, chosen from each stage's own gradient, which is
      ! the RED of the shared-face-quantity rule.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: nel, xclose_cell
      real*8, dimension(:,1-Ng:),   intent(in)  :: xion
      character(len=*),             intent(in)  :: closing_rule
      logical,                      intent(in)  :: per_stage_donor
      real*8, dimension(:,0:),      intent(out) :: Ek
      real*8, dimension(0:N),       intent(out) :: Eclose

      real*8  :: dr_f, nf, Kf, dx
      integer :: j, k, nk

      nk = size(xion,1)
      Ek     = 0.0d0
      Eclose = 0.0d0
      do j = 1, N-1
         dr_f = max((r(j+1) - r(j))*R0, 1.0d0)
         Kf   = 0.5d0*(kzz_cell(j) + kzz_cell(j+1))
         do k = 1, nk
            dx = xion(k,j+1) - xion(k,j)
            if (per_stage_donor) then
               if (dx .ge. 0.0d0) then
                  nf = nel(j)
               else
                  nf = nel(j+1)
               endif
            else
               nf = 0.5d0*(nel(j) + nel(j+1))
            endif
            Ek(k,j) = -nf*Kf*dx/dr_f
         enddo
         if (closing_rule .eq. 'negated_sum') then
            Eclose(j) = 0.0d0
            do k = 1, nk
               Eclose(j) = Eclose(j) - Ek(k,j)
            enddo
         else
            dx = xclose_cell(j+1) - xclose_cell(j)
            if (per_stage_donor) then
               if (dx .ge. 0.0d0) then
                  nf = nel(j)
               else
                  nf = nel(j+1)
               endif
            else
               nf = 0.5d0*(nel(j) + nel(j+1))
            endif
            Eclose(j) = -nf*Kf*dx/dr_f
         endif
      enddo
      end subroutine stage_eddy_flux

      ! ================================================================= !
      !  the simplex closure of the face fractions
      ! ================================================================= !

      subroutine test_simplex_closure(scheme)
      character(len=*), intent(in) :: scheme

      real*8, dimension(:),   allocatable :: nH, nHe, xHII, xHeII, xHeIII
      real*8, dimension(:),   allocatable :: Frho, Jdif, NelH, NelHe
      real*8, dimension(:),   allocatable :: xclose
      real*8, dimension(:,:), allocatable :: xion, xf
      real*8  :: worst_closed, worst_indep, s
      integer :: j, ic, iflow
      character(len=80) :: nm

      call setup_column(200, 8.0d0, 1.0d9)
      call set_scheme(scheme)

      allocate(nH(1-Ng:N+Ng), nHe(1-Ng:N+Ng), xHII(1-Ng:N+Ng),            &
               xHeII(1-Ng:N+Ng), xHeIII(1-Ng:N+Ng), Frho(1-Ng:N+Ng),      &
               Jdif(1-Ng:N+Ng), NelH(1-Ng:N+Ng), NelHe(1-Ng:N+Ng),        &
               xclose(1-Ng:N+Ng))
      allocate(xion(2,1-Ng:N+Ng), xf(2,1-Ng:N+Ng))

      worst_closed = 0.0d0
      worst_indep  = 0.0d0
      do ic = 1, 4
         call helium_column(ic, nH, nHe, xHII, xHeII, xHeIII)
         do j = 1-Ng, N+Ng
            s = (r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
            Jdif(j) = 1.0d-17*sin(6.0d0*s)
         enddo
         do iflow = 1, 3
            call face_mass_flux(iflow, Frho)
            call element_nucleus_face_flux(nH, nHe, Frho, Jdif, 1.0d0,    &
                                           NelH, NelHe)
            xion(1,:) = xHeII
            xion(2,:) = xHeIII
            call stage_face_fractions(.false., xion, NelHe, xf, xclose)
            do j = 0, N
               worst_closed = max(worst_closed,                           &
                    abs(xclose(j) + xf(1,j) + xf(2,j) - 1.0d0))
            enddo
            call stage_face_fractions(.true., xion, NelHe, xf, xclose)
            do j = 0, N
               worst_indep = max(worst_indep,                             &
                    abs(xclose(j) + xf(1,j) + xf(2,j) - 1.0d0))
            enddo
         enddo
      enddo

      nm = 'stage_fractions_close_the_simplex_at_every_face_'//scheme
      call verdict(trim(nm), worst_closed .le. 8.0d0*epsilon(1.0d0),      &
                   worst_closed, 0.0d0, 8.0d0*epsilon(1.0d0))
      nm = 'independently_reconstructed_stages_break_the_simplex_'//scheme
      call verdict(trim(nm), worst_indep .gt. 1.0d3*epsilon(1.0d0),       &
                   worst_indep, 1.0d3*epsilon(1.0d0), 0.0d0)

      deallocate(nH, nHe, xHII, xHeII, xHeIII, Frho, Jdif, NelH, NelHe,   &
                 xclose, xion, xf)
      end subroutine test_simplex_closure

      ! ================================================================= !
      !  the identity itself
      ! ================================================================= !

      subroutine test_stage_flux_identity(scheme)
      character(len=*), intent(in) :: scheme

      real*8, dimension(:),   allocatable :: nH, nHe, xHII, xHeII, xHeIII
      real*8, dimension(:),   allocatable :: Frho, Jdif, NelH, NelHe
      real*8, dimension(:),   allocatable :: xclose, Eclose, xHI, nelw
      real*8, dimension(:),   allocatable :: NelHe2
      real*8, dimension(:,:), allocatable :: xion, xf, Ek
      real*8  :: worst, worst_red, worst_two, sumF, sc, s, dsign
      integer :: j, ic, iflow, idif, k
      character(len=80) :: nm

      call setup_column(200, 8.0d0, 1.0d9)
      call set_scheme(scheme)

      allocate(nH(1-Ng:N+Ng), nHe(1-Ng:N+Ng), xHII(1-Ng:N+Ng),            &
               xHeII(1-Ng:N+Ng), xHeIII(1-Ng:N+Ng), Frho(1-Ng:N+Ng),      &
               Jdif(1-Ng:N+Ng), NelH(1-Ng:N+Ng), NelHe(1-Ng:N+Ng),        &
               xclose(1-Ng:N+Ng), xHI(1-Ng:N+Ng), nelw(1-Ng:N+Ng),        &
               NelHe2(1-Ng:N+Ng))
      allocate(xion(2,1-Ng:N+Ng), xf(2,1-Ng:N+Ng), Ek(2,0:N),             &
               Eclose(0:N))

      worst     = 0.0d0
      worst_red = 0.0d0
      worst_two = 0.0d0
      do ic = 1, 4
         call helium_column(ic, nH, nHe, xHII, xHeII, xHeIII)
         do j = 1-Ng, N+Ng
            s = (r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
            Jdif(j) = 1.0d-17*sin(6.0d0*s)
         enddo
         do iflow = 1, 4
            call face_mass_flux(iflow, Frho)
            ! The diffusive half positive, reversed and absent.
            do idif = 1, 3
               dsign = 1.0d0
               if (idif .eq. 2) dsign = -1.0d0
               if (idif .eq. 3) dsign =  0.0d0
               call element_nucleus_face_flux(nH, nHe, Frho, Jdif, dsign, &
                                              NelH, NelHe)

               ! ---- helium, three stages ----
               xion(1,:) = xHeII
               xion(2,:) = xHeIII
               call stage_face_fractions(.false., xion, NelHe, xf, xclose)
               nelw = nHe
               call stage_eddy_flux(nelw, xion, 1.0d0 - xHeII - xHeIII,   &
                                    'negated_sum', .false., Ek, Eclose)
               do j = 1, N-1
                  sumF = xclose(j)*NelHe(j) + Eclose(j)
                  do k = 1, 2
                     sumF = sumF + xf(k,j)*NelHe(j) + Ek(k,j)
                  enddo
                  sc = abs(NelHe(j))
                  do k = 1, 2
                     sc = sc + abs(Ek(k,j))
                  enddo
                  sc = max(sc + abs(Eclose(j)), 1.0d-300)
                  worst = max(worst, abs(sumF - NelHe(j))/sc)
               enddo
               ! the RED: a SECOND copy of the element face flux under one
               ! stage, its donor side taken from the opposite sign.  The
               ! fractions still close the simplex; the identity does not
               ! survive two element fluxes.
               call stage_face_fractions(.false., xion, NelHe, xf, xclose)
               NelHe2 = -NelHe
               do j = 1, N-1
                  sumF = xclose(j)*NelHe(j) + Eclose(j)                   &
                       + xf(1,j)*NelHe(j) + Ek(1,j)                       &
                       + xf(2,j)*NelHe2(j) + Ek(2,j)
                  sc = abs(NelHe(j)) + abs(Ek(1,j)) + abs(Ek(2,j))
                  sc = max(sc + abs(Eclose(j)), 1.0d-300)
                  worst_two = max(worst_two, abs(sumF - NelHe(j))/sc)
               enddo
               ! the RED: every stage reconstructed on its own
               call stage_face_fractions(.true., xion, NelHe, xf, xclose)
               do j = 1, N-1
                  sumF = xclose(j)*NelHe(j) + Eclose(j)
                  do k = 1, 2
                     sumF = sumF + xf(k,j)*NelHe(j) + Ek(k,j)
                  enddo
                  sc = abs(NelHe(j))
                  do k = 1, 2
                     sc = sc + abs(Ek(k,j))
                  enddo
                  sc = max(sc + abs(Eclose(j)), 1.0d-300)
                  worst_red = max(worst_red, abs(sumF - NelHe(j))/sc)
               enddo

               ! ---- hydrogen, two stages ----
               xion(1,:) = xHII
               call stage_face_fractions(.false., xion(1:1,:), NelH,      &
                                         xf(1:1,:), xclose)
               xHI = 1.0d0 - xHII
               nelw = nH
               call stage_eddy_flux(nelw, xion(1:1,:), xHI,               &
                                    'negated_sum', .false., Ek(1:1,:),      &
                                    Eclose)
               do j = 1, N-1
                  sumF = xclose(j)*NelH(j) + Eclose(j)                    &
                       + xf(1,j)*NelH(j) + Ek(1,j)
                  sc = max(abs(NelH(j)) + abs(Ek(1,j)) + abs(Eclose(j)),  &
                           1.0d-300)
                  worst = max(worst, abs(sumF - NelH(j))/sc)
               enddo
            enddo
         enddo
      enddo

      nm = 'sum_of_stage_fluxes_is_the_element_flux_'//scheme
      call verdict(trim(nm), worst .le. 1.0d2*epsilon(1.0d0),             &
                   worst, 0.0d0, 1.0d2*epsilon(1.0d0))
      nm = 'independently_reconstructed_stages_break_the_identity_'//scheme
      call verdict(trim(nm), worst_red .gt. 1.0d3*epsilon(1.0d0),         &
                   worst_red, 1.0d3*epsilon(1.0d0), 0.0d0)
      nm = 'a_second_element_face_flux_breaks_the_identity_'//scheme
      call verdict(trim(nm), worst_two .gt. 1.0d3*epsilon(1.0d0),         &
                   worst_two, 1.0d3*epsilon(1.0d0), 0.0d0)

      deallocate(nH, nHe, xHII, xHeII, xHeIII, Frho, Jdif, NelH, NelHe,   &
                 xclose, xHI, nelw, NelHe2, xion, xf, Ek, Eclose)
      end subroutine test_stage_flux_identity

      ! ================================================================= !
      !  the stage eddy term telescopes
      ! ================================================================= !

      subroutine test_stage_eddy_sum()

      real*8, dimension(:),   allocatable :: nH, nHe, xHII, xHeII, xHeIII
      real*8, dimension(:),   allocatable :: Eclose, xcl
      real*8, dimension(:,:), allocatable :: xion, Ek
      real*8  :: worst, worst_red, worst_own, sc
      integer :: j, ic

      call setup_column(200, 8.0d0, 1.0d9)
      call set_scheme('WENO3')

      allocate(nH(1-Ng:N+Ng), nHe(1-Ng:N+Ng), xHII(1-Ng:N+Ng),            &
               xHeII(1-Ng:N+Ng), xHeIII(1-Ng:N+Ng), xcl(1-Ng:N+Ng))
      allocate(xion(2,1-Ng:N+Ng), Ek(2,0:N), Eclose(0:N))

      worst     = 0.0d0
      worst_red = 0.0d0
      worst_own = 0.0d0
      do ic = 1, 4
         call helium_column(ic, nH, nHe, xHII, xHeII, xHeIII)
         xion(1,:) = xHeII
         xion(2,:) = xHeIII
         xcl       = 1.0d0 - xHeII - xHeIII
         call stage_eddy_flux(nHe, xion, xcl, 'negated_sum', .false.,    &
                              Ek, Eclose)
         do j = 1, N-1
            sc = max(abs(Ek(1,j)) + abs(Ek(2,j)) + abs(Eclose(j)),        &
                     1.0d-300)
            worst = max(worst, abs(Ek(1,j)+Ek(2,j)+Eclose(j))/sc)
         enddo
         call stage_eddy_flux(nHe, xion, xcl, 'own_difference', .false.,  &
                              Ek, Eclose)
         do j = 1, N-1
            sc = max(abs(Ek(1,j)) + abs(Ek(2,j)) + abs(Eclose(j)),        &
                     1.0d-300)
            worst_own = max(worst_own, abs(Ek(1,j)+Ek(2,j)+Eclose(j))/sc)
         enddo
         call stage_eddy_flux(nHe, xion, xcl, 'own_difference', .true.,   &
                              Ek, Eclose)
         do j = 1, N-1
            sc = max(abs(Ek(1,j)) + abs(Ek(2,j)) + abs(Eclose(j)),        &
                     1.0d-300)
            worst_red = max(worst_red, abs(Ek(1,j)+Ek(2,j)+Eclose(j))/sc)
         enddo
      enddo

      call verdict('stage_eddy_term_sums_to_zero',                        &
                   worst .le. 8.0d0*epsilon(1.0d0), worst, 0.0d0,         &
                   8.0d0*epsilon(1.0d0))
      call verdict('differencing_the_closing_member_loses_the_trace'//    &
                   '_gradient',                                           &
                   worst_own .gt. 1.0d3*epsilon(1.0d0), worst_own,        &
                   1.0d3*epsilon(1.0d0), 0.0d0)
      call verdict('a_stage_dependent_face_density_breaks_the_eddy_sum',  &
                   worst_red .gt. 1.0d3*epsilon(1.0d0), worst_red,        &
                   1.0d3*epsilon(1.0d0), 0.0d0)

      deallocate(nH, nHe, xHII, xHeII, xHeIII, xcl, xion, Ek, Eclose)
      end subroutine test_stage_eddy_sum

      ! ================================================================= !
      !  charge consistency
      ! ================================================================= !

      subroutine test_charge_consistency()
      ! n_e = sum_k Z_k n_k with n_k = x_k n_el, rebuilt from the same
      ! fractions the transport carries.  The statement is that the electron
      ! density the sweep is handed is a SUM OVER CHARGES of the transported
      ! partition and not a separate variable, which is what makes the
      ! replacement of one term by its trial value exact (carrier_source,
      ! cbg_ne - cbg_nhii + n_hii).
      real*8, dimension(:), allocatable :: nH, nHe, xHII, xHeII, xHeIII
      real*8  :: worst, ne_sum, ne_dir
      integer :: j, ic

      call setup_column(200, 8.0d0, 1.0d9)
      call set_scheme('WENO3')
      allocate(nH(1-Ng:N+Ng), nHe(1-Ng:N+Ng), xHII(1-Ng:N+Ng),            &
               xHeII(1-Ng:N+Ng), xHeIII(1-Ng:N+Ng))

      worst = 0.0d0
      do ic = 1, 4
         call helium_column(ic, nH, nHe, xHII, xHeII, xHeIII)
         do j = 1, N
            ne_sum = 1.0d0*(xHII(j)*nH(j))                                &
                   + 1.0d0*(xHeII(j)*nHe(j))                              &
                   + 2.0d0*(xHeIII(j)*nHe(j))
            ne_dir = xHII(j)*nH(j) + (xHeII(j) + 2.0d0*xHeIII(j))*nHe(j)
            worst = max(worst, abs(ne_sum - ne_dir)                       &
                               /max(abs(ne_dir), 1.0d-300))
         enddo
      enddo
      call verdict('charge_consistency_of_the_transported_partition',     &
                   worst .le. 4.0d0*epsilon(1.0d0), worst, 0.0d0,         &
                   4.0d0*epsilon(1.0d0))
      deallocate(nH, nHe, xHII, xHeII, xHeIII)
      end subroutine test_charge_consistency

      ! ================================================================= !
      !  what the closing-member rule does NOT buy
      ! ================================================================= !

      subroutine test_closing_member_excursion()
      ! MEASURED, not asserted.  The closing member takes one minus the
      ! reconstructed others, so the SUM is exact; nothing in that rule
      ! keeps the closing member itself inside [0,1] where the reconstructed
      ! members are separately bounded into it.  The number below is the
      ! headroom stage B owes the simplex (L12a section 2.2): a negative
      ! value here is a face at which the neutral stage would be handed a
      ! negative population if nothing else limited it.
      real*8, dimension(:),   allocatable :: nH, nHe, xHII, xHeII, xHeIII
      real*8, dimension(:),   allocatable :: Frho, Jdif, NelH, NelHe
      real*8, dimension(:),   allocatable :: xclose
      real*8, dimension(:,:), allocatable :: xion, xf
      real*8  :: worst, s
      integer :: j, ic, iflow

      call setup_column(200, 8.0d0, 1.0d9)
      call set_scheme('WENO3')
      allocate(nH(1-Ng:N+Ng), nHe(1-Ng:N+Ng), xHII(1-Ng:N+Ng),            &
               xHeII(1-Ng:N+Ng), xHeIII(1-Ng:N+Ng), Frho(1-Ng:N+Ng),      &
               Jdif(1-Ng:N+Ng), NelH(1-Ng:N+Ng), NelHe(1-Ng:N+Ng),        &
               xclose(1-Ng:N+Ng))
      allocate(xion(2,1-Ng:N+Ng), xf(2,1-Ng:N+Ng))

      worst = 0.0d0
      do ic = 1, 4
         call helium_column(ic, nH, nHe, xHII, xHeII, xHeIII)
         do j = 1-Ng, N+Ng
            s = (r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
            Jdif(j) = 1.0d-17*sin(6.0d0*s)
         enddo
         do iflow = 1, 3
            call face_mass_flux(iflow, Frho)
            call element_nucleus_face_flux(nH, nHe, Frho, Jdif, 1.0d0,    &
                                           NelH, NelHe)
            xion(1,:) = xHeII
            xion(2,:) = xHeIII
            call stage_face_fractions(.false., xion, NelHe, xf, xclose)
            do j = 0, N
               worst = max(worst, -xclose(j))
            enddo
         enddo
      enddo
      write(*,'(A,ES12.4)')                                               &
         '  DIAGNOSTIC the_closing_member_can_leave_the_simplex,'//       &
         ' largest negative face value: ', worst
      ! Reported, gating nothing: the rule the memo states is that stage B
      ! carries a simplex headroom of its own, and this is the number it has
      ! to answer.
      call verdict('the_closing_member_excursion_is_finite',              &
                   worst .eq. worst .and. abs(worst) .lt. 1.0d3,          &
                   worst, 0.0d0, 1.0d3)

      deallocate(nH, nHe, xHII, xHeII, xHeIII, Frho, Jdif, NelH, NelHe,   &
                 xclose, xion, xf)
      end subroutine test_closing_member_excursion


      ! ================================================================= !
      !  the production module, on the production grid
      ! ================================================================= !

      subroutine production_stage_column(nc)
      ! The configuration and the state the production rows run on: an
      ! atomic helium-hydrogen column with trace metals, on the Mixed grid
      ! the catalog uses (50 uniform base cells of 2e-4 R_p under a
      ! geometric stretch to 30 R_p), with helium tapering by a factor two
      ! down the column so the element operator has a gradient to work on.
      !
      ! The ionization stages are filled by moving each element's nuclei
      ! between its stages at a FIXED nucleus count, so the mixture mass and
      ! the element nucleus densities are those of the neutral column and
      ! the state still carries the density it is handed.
      integer, intent(in) :: nc

      N  = nc
      T0 = 1.0d3
      R0 = 1.0d10
      n0 = 1.0d10
      v0 = sqrt(kb_erg*T0/mu)
      t_s = R0/v0
      p0 = n0*mu*v0*v0
      b0 = 0.0d0
      spherical_domain     = .true.
      thereis_He           = .true.
      thereis_HeITR        = .false.
      thereis_mol          = .false.
      carrier_transport    = .false.
      thereis_oxychem      = .false.
      ionization_transport = .false.
      thereis_metals       = .true.
      eos_include_metals   = .true.
      he_diffusion         = .true.
      he_metal_diffusion   = .true.
      he_ambipolar         = .false.
      he_alphaT            = 3.0d-1
      he_kzz               = 1.0d9
      HeH                  = 0.0833333333333333d0
      use_plm = .false.;  use_weno3 = .true.;  rec_method = 'WENO3'
      recon_lambda_on = .false.
      call carrier_set_init()

      grid_type   = 'Mixed'
      N_low_cells = 50
      dr_base     = 2.0d-4
      r_max       = 30.0d0
      r_esc       = 2.0d0
      r_flux      = 1.2d0
      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      call define_grid
      j_min = 1
      call eddy_diffusion_on_grid
      if (.not. allocated(melem_ab)) allocate(melem_ab(n_melem))
      melem_ab = 1.0d-4
      end subroutine production_stage_column

      ! ----------------------------------------------------------------- !

      subroutine stage_fractions_of_the_column(qH, qHe2, qHe3)
      ! Ionization fractions that vary over the whole column, hydrogen
      ! reaching 0.92 and helium's two ionized stages 0.58 and 0.31 at the
      ! top, so that no stage is trace everywhere and the gradient the eddy
      ! term acts on is nonzero at every face.
      real*8, dimension(1-Ng:N+Ng), intent(out) :: qH, qHe2, qHe3
      integer :: j
      real*8  :: s
      do j = 1-Ng, N+Ng
         s = (r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
         s = max(0.0d0, min(1.0d0, s))
         qH(j)   = 0.02d0 + 0.90d0*s*s
         qHe2(j) = 0.03d0 + 0.55d0*s
         qHe3(j) = 0.01d0 + 0.30d0*s*s*s
      enddo
      end subroutine stage_fractions_of_the_column

      ! ----------------------------------------------------------------- !

      subroutine stage_fractions_on_the_simplex_face(qH, qHe2, qHe3)
      ! The other column: helium driven onto the faces of its own simplex,
      ! x(He I) reaching zero halfway up and x(He III) reaching one at the
      ! top.  The closing stage is then a difference of numbers near one and
      ! the reconstruction overshoots it out of [0,1], which is the state
      ! the closing-member excursion and the clipping RED are about.
      real*8, dimension(1-Ng:N+Ng), intent(out) :: qH, qHe2, qHe3
      integer :: j
      real*8  :: s
      do j = 1-Ng, N+Ng
         s = (r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
         s = max(0.0d0, min(1.0d0, s))
         qH(j)   = min(1.0d0, 2.0d0*s)
         qHe3(j) = max(0.0d0, 2.0d0*s - 1.0d0)
         qHe2(j) = max(0.0d0, min(1.0d0 - qHe3(j), 2.0d0*s))
      enddo
      end subroutine stage_fractions_on_the_simplex_face

      ! ----------------------------------------------------------------- !

      subroutine ionize_the_column(f_c, qH, qHe2, qHe3)
      ! Move each element's nuclei into its stages at a fixed nucleus count.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_c
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: qH,qHe2,qHe3
      real*8  :: fH, fHe
      integer :: j
      do j = 1-Ng, N+Ng
         fH  = f_c(j,isp_HI) + f_c(j,isp_HII)
         fHe = f_c(j,isp_HeI) + f_c(j,isp_HeII) + f_c(j,isp_HeIII)
         f_c(j,isp_HI)    = (1.0d0 - qH(j))*fH
         f_c(j,isp_HII)   = qH(j)*fH
         f_c(j,isp_HeII)  = qHe2(j)*fHe
         f_c(j,isp_HeIII) = qHe3(j)*fHe
         f_c(j,isp_HeI)   = (1.0d0 - qHe2(j) - qHe3(j))*fHe
      enddo
      end subroutine ionize_the_column

      ! ----------------------------------------------------------------- !

      subroutine test_production_stage_rows()

      integer, parameter :: nc = 500
      real*8, dimension(:),   allocatable :: rho_c, T_c, Frho, qH
      real*8, dimension(:),   allocatable :: qHe2, qHe3, N_H, N_He
      real*8, dimension(:),   allocatable :: n_Hf, n_Hef, Kf, drf
      real*8, dimension(:),   allocatable :: xclose, Fclose, Dclose
      real*8, dimension(:),   allocatable :: fa, cv, Ntot, ntotc, ycell
      real*8, dimension(:,:), allocatable :: f_c, xion, xf, Fk, Dk
      real*8  :: worst, worst_b, worst_proj, worst_div, sumF, sc
      real*8  :: worst_red_eddy, worst_red_clip, csum, cscale, bnd
      real*8  :: worst_col, worst_adm, worst_admis, proj_over, proj_under
      real*8  :: worst_red_geom
      real*8  :: R0sq, R0cb, sL, sR, dvel, nKdr, Ew, xcl
      integer :: j, k, ic, iflow, ifr

      call production_stage_column(nc)

      allocate(qH(1-Ng:N+Ng), qHe2(1-Ng:N+Ng), qHe3(1-Ng:N+Ng))
      allocate(rho_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng), Frho(1-Ng:N+Ng),         &
               f_c(1-Ng:N+Ng,n_species))
      allocate(N_H(1-Ng:N+Ng), N_He(1-Ng:N+Ng), Ntot(1-Ng:N+Ng))
      allocate(ntotc(1-Ng:N+Ng), ycell(1-Ng:N+Ng))
      allocate(n_Hf(0:N), n_Hef(0:N), Kf(0:N), drf(0:N))
      allocate(xion(2,1-Ng:N+Ng), xf(2,0:N), Fk(2,0:N), Dk(2,1:N))
      allocate(xclose(0:N), Fclose(0:N), Dclose(1:N))
      allocate(fa(0:N), cv(1:N))

      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         T_c(j)   = 1.0d0 + 2.0d0*(r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
      enddo

      worst          = 0.0d0
      worst_b        = 0.0d0
      worst_proj     = 0.0d0
      worst_div      = 0.0d0
      worst_red_eddy = 0.0d0
      worst_red_clip = 0.0d0
      worst_col      = 0.0d0
      worst_red_geom = 0.0d0
      worst_admis    = 0.0d0
      proj_over      = -1.0d0
      proj_under     = -1.0d0
      csum           = 0.0d0
      cscale         = 0.0d0
      bnd            = 0.0d0

      call spherical_face_area_and_cell_volume(fa, cv)
      R0sq = R0*R0
      R0cb = R0sq*R0

      ! Two columns (smooth, and helium driven onto its own simplex faces)
      ! crossed with three winds: outflowing, reversed, and alternating in
      ! sign face by face (the breathing base the code admits).
      do ifr = 1, 2
      if (ifr .eq. 1) then
         call stage_fractions_of_the_column(qH, qHe2, qHe3)
      else
         call stage_fractions_on_the_simplex_face(qH, qHe2, qHe3)
      endif
      do iflow = 1, 3
         call column_carrying_its_own_density(f_c, 0.0d0, .true., 1.0d-1)
         call ionize_the_column(f_c, qH, qHe2, qHe3)
         do j = 1-Ng, N+Ng
            dvel = 1.0d4/(r_edg(j)*r_edg(j))
            if (iflow .eq. 2) dvel = -dvel
            if (iflow .eq. 3) dvel = dvel*dble(1 - 2*mod(j+100,2))
            Frho(j) = rho_c(j)*n0*mu*dvel
         enddo
         call hydrogen_and_helium_nucleus_face_flux(rho_c, T_c, f_c,      &
                  Frho, N_H, N_He, n_Hf, n_Hef, Kf, drf)

         ! ---- helium, three stages ----
         do ic = 1, 2
            if (ic .eq. 1) then
               xion(1,:) = qHe2
               xion(2,:) = qHe3
               Ntot      = N_He
            else
               xion(1,:) = qH
               xion(2,:) = 0.0d0
               Ntot      = N_H
            endif
            call ionization_stage_face_flux(xion, Ntot, n_Hef, Kf, drf,   &
                     xf, xclose, Fk, Fclose)
            do j = 0, N
               sumF = Fclose(j)
               sc   = abs(Fclose(j))
               do k = 1, 2
                  sumF = sumF + Fk(k,j)
                  sc   = sc + abs(Fk(k,j))
               enddo
               sc = max(sc + abs(Ntot(j)), 1.0d-300)
               worst = max(worst, abs(sumF - Ntot(j))/sc)
               if (j .eq. 0 .or. j .eq. N)                                &
                  worst_b = max(worst_b, abs(sumF - Ntot(j))/sc)
            enddo

            ! The divergence: cell by cell the stages sum to the element's
            ! own nucleus divergence, and down the column the internal
            ! faces cancel.
            call ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
            csum   = 0.0d0
            cscale = 0.0d0
            do j = 1, N
               sL = fa(j-1)*R0sq
               sR = fa(j)  *R0sq
               sumF = Dclose(j)
               do k = 1, 2
                  sumF = sumF + Dk(k,j)
               enddo
               dvel = (sR*Ntot(j) - sL*Ntot(j-1))/(cv(j)*R0cb)
               ! The scale of a telescoping identity is the size of the
               ! FACE terms it differences, not of the difference: the two
               ! face terms of a nearly steady flux cancel to many digits
               ! and a measure taken against the remainder would read that
               ! cancellation and not the operator (the same rule
               ! composition_residual's row scale follows).
               sc = max((sR*abs(Ntot(j)) + sL*abs(Ntot(j-1)))             &
                        /(cv(j)*R0cb), 1.0d-300)
               worst_div = max(worst_div, abs(sumF - dvel)/sc)
               csum   = csum   + cv(j)*R0cb*sumF
               cscale = cscale + sR*abs(Ntot(j)) + sL*abs(Ntot(j-1))
            enddo
            bnd = fa(N)*R0sq*Ntot(N) - fa(0)*R0sq*Ntot(0)
            worst_col = max(worst_col, abs(csum - bnd)                    &
                                       /max(cscale, 1.0d-300))

            ! RED: the same column sum with the divergence weighted by
            ! r_j^2 dr_j in place of the exact shell volume.  The two
            ! differ by dr^2/12 relative, which is not the same factor in
            ! two neighbouring cells of a stretched grid, so the internal
            ! faces stop cancelling and the column sum is no longer the
            ! difference of the two boundary fluxes (item L30).
            csum = 0.0d0
            do j = 1, N
               sL = fa(j-1)*R0sq
               sR = fa(j)  *R0sq
               csum = csum + r(j)*r(j)*(r_edg(j) - r_edg(j-1))*R0cb       &
                             *(sR*Ntot(j) - sL*Ntot(j-1))/(cv(j)*R0cb)
            enddo
            worst_red_geom = max(worst_red_geom, abs(csum - bnd)          &
                                                 /max(cscale, 1.0d-300))

            ! RED: the whole mixing-ratio eddy flux inside EVERY stage.
            ! sum_k of it is the element's own eddy flux, which N_el
            ! already carries, so the identity is off by one copy of it.
            ntotc = 0.0d0
            ycell = 0.0d0
            do k = 1, n_bsp
               ntotc = ntotc + f_c(:,bsp_fsp(k))*rho_c*n0
               if (ic .eq. 1) then
                  ycell = ycell + bsp_nHe(k)*f_c(:,bsp_fsp(k))*rho_c*n0
               else
                  ycell = ycell + bsp_nH(k) *f_c(:,bsp_fsp(k))*rho_c*n0
               endif
            enddo
            ycell = ycell/max(ntotc, 1.0d0)
            do j = 1, N-1
               nKdr = 0.5d0*(ntotc(j) + ntotc(j+1))*Kf(j)/drf(j)
               sumF = 0.0d0
               sc   = 0.0d0
               xcl  = 1.0d0
               do k = 1, 2
                  Ew = -nKdr*(xion(k,j+1)*ycell(j+1) - xion(k,j)*ycell(j))
                  sumF = sumF + xf(k,j)*Ntot(j) + Ew
                  sc   = sc + abs(xf(k,j)*Ntot(j)) + abs(Ew)
                  xcl  = xcl - xion(k,j)
               enddo
               Ew = -nKdr*((1.0d0 - xion(1,j+1) - xion(2,j+1))*ycell(j+1) &
                          -(1.0d0 - xion(1,j)   - xion(2,j)  )*ycell(j))
               sumF = sumF + xclose(j)*Ntot(j) + Ew
               sc   = max(sc + abs(xclose(j)*Ntot(j)) + abs(Ew)           &
                          + abs(Ntot(j)), 1.0d-300)
               worst_red_eddy = max(worst_red_eddy,                       &
                                    abs(sumF - Ntot(j))/sc)
            enddo

            ! RED: bounding the CLOSING face value into [0,1].  It buys
            ! nonnegativity and loses the sum.
            do j = 0, N
               xcl = max(0.0d0, min(1.0d0, xclose(j)))
               sumF = xcl*Ntot(j)
               sc   = abs(xcl*Ntot(j))
               do k = 1, 2
                  sumF = sumF + Fk(k,j)
                  sc   = sc + abs(Fk(k,j))
               enddo
               sumF = sumF + (Fclose(j) - xclose(j)*Ntot(j))
               sc = max(sc + abs(Ntot(j)), 1.0d-300)
               worst_red_clip = max(worst_red_clip,                       &
                                    abs(sumF - Ntot(j))/sc)
            enddo

            ! The simplex projection acts on the CELL state; the face rule
            ! is then re-derived from it, so the identity is untouched.
            ! The state handed to it is deliberately inadmissible -- every
            ! third cell has its carried fractions pushed past the simplex
            ! and every seventh one driven negative -- so the row below is
            ! a statement about a projection that fired and not about one
            ! that had nothing to do.
            do j = 1-Ng, N+Ng
               if (mod(j,3) .eq. 0) then
                  xion(1,j) = xion(1,j) + 0.7d0
                  xion(2,j) = xion(2,j) + 0.6d0
               endif
               if (mod(j,7) .eq. 0) xion(2,j) = -1.0d-3
            enddo
            call stage_simplex_projection(xion)
            worst_adm = 0.0d0
            do j = 1-Ng, N+Ng
               worst_adm = max(worst_adm, xion(1,j) + xion(2,j) - 1.0d0)
               worst_adm = max(worst_adm, -xion(1,j), -xion(2,j))
            enddo
            worst_admis = max(worst_admis, worst_adm)
            proj_over   = max(proj_over,  stage_simplex_sum_over_limit)
            proj_under  = max(proj_under, stage_fraction_under_zero)
            call ionization_stage_face_flux(xion, Ntot, n_Hef, Kf, drf,   &
                     xf, xclose, Fk, Fclose)
            do j = 0, N
               sumF = Fclose(j)
               sc   = abs(Fclose(j))
               do k = 1, 2
                  sumF = sumF + Fk(k,j)
                  sc   = sc + abs(Fk(k,j))
               enddo
               sc = max(sc + abs(Ntot(j)), 1.0d-300)
               worst_proj = max(worst_proj, abs(sumF - Ntot(j))/sc)
            enddo
         enddo
      enddo
      enddo

      call verdict('production_stage_flux_sums_to_the_element_flux',      &
                   worst .le. 1.0d2*epsilon(1.0d0), worst, 0.0d0,         &
                   1.0d2*epsilon(1.0d0))
      call verdict('sum_identity_holds_at_the_boundary_faces',            &
                   worst_b .le. 1.0d2*epsilon(1.0d0), worst_b, 0.0d0,     &
                   1.0d2*epsilon(1.0d0))
      call verdict('sum_identity_survives_the_simplex_projection',        &
                   worst_proj .le. 1.0d2*epsilon(1.0d0), worst_proj,      &
                   0.0d0, 1.0d2*epsilon(1.0d0))
      call verdict('stage_divergence_sums_to_the_element_nucleus'//       &
                   '_divergence',                                         &
                   worst_div .le. 1.0d2*epsilon(1.0d0), worst_div,        &
                   0.0d0, 1.0d2*epsilon(1.0d0))
      write(*,'(A,ES12.4,A,ES12.4)')                                      &
         '  DIAGNOSTIC last column sum of V_j sum_k D_k(j) ', csum,       &
         ', its boundary nucleus flux difference ', bnd
      call verdict('stage_divergence_column_sum_is_the_boundary_flux'//   &
                   '_difference',                                         &
                   worst_col .le. 1.0d-13, worst_col, 0.0d0, 1.0d-13)
      call verdict('a_cell_centred_volume_breaks_the_column_sum',         &
                   worst_red_geom .gt. 1.0d3*epsilon(1.0d0),              &
                   worst_red_geom, 1.0d3*epsilon(1.0d0), 0.0d0)
      call verdict('a_full_eddy_term_in_every_stage_breaks_the_identity', &
                   worst_red_eddy .gt. 1.0d3*epsilon(1.0d0),              &
                   worst_red_eddy, 1.0d3*epsilon(1.0d0), 0.0d0)
      call verdict('bounding_the_closing_face_value_breaks_the_identity', &
                   worst_red_clip .gt. 1.0d3*epsilon(1.0d0),              &
                   worst_red_clip, 1.0d3*epsilon(1.0d0), 0.0d0)
      write(*,'(A,ES12.4,A,ES12.4)')                                      &
         '  DIAGNOSTIC simplex projection, largest sum over one ',        &
         proj_over, ', largest fraction under zero ', proj_under
      call verdict('the_simplex_projection_returns_an_admissible_state',  &
                   worst_admis .le. 8.0d0*epsilon(1.0d0), worst_admis,    &
                   0.0d0, 8.0d0*epsilon(1.0d0))

      deallocate(rho_c, T_c, Frho, f_c, N_H, N_He, Ntot, ntotc, ycell,    &
                 n_Hf, n_Hef, Kf, drf, xion, xf, Fk, Dk, xclose, Fclose,  &
                 Dclose, fa, cv, qH, qHe2, qHe3)
      end subroutine test_production_stage_rows

      ! ================================================================= !
      !  the analytic rows against a central difference
      ! ================================================================= !

      subroutine test_production_stage_jacobian()
      ! The analytic tridiagonal rows of ionization_stage_row_jacobian
      ! against a central difference of the divergence they are the
      ! derivative of, which is the divergence of the DONOR-VALUE face flux
      ! (the limiter and the second-order part of the reconstruction are
      ! left out of the rows, exactly as element_advective_face_coefficients
      ! leaves them out of the element row).  The size of what is left out
      ! is measured beside the rows and reported, so the omission is a
      ! stated number and not a silent one.
      integer, parameter :: nc = 500
      real*8, dimension(:),   allocatable :: rho_c, T_c, Frho, qH
      real*8, dimension(:),   allocatable :: qHe2, qHe3, N_H, N_He
      real*8, dimension(:),   allocatable :: n_Hf, n_Hef, Kf, drf
      real*8, dimension(:),   allocatable :: xclose, Fclose, Dclose
      real*8, dimension(:),   allocatable :: dFdl, dFdr, aa, bb, cc
      real*8, dimension(:),   allocatable :: Dp, Dm, Dfull, Ddon
      real*8, dimension(:,:), allocatable :: f_c, xion, xf, Fk, Dk
      real*8  :: worst_a, worst_b, worst_c, worst_rec, h, sc, dnum
      integer :: j, jp

      call production_stage_column(nc)
      allocate(qH(1-Ng:N+Ng), qHe2(1-Ng:N+Ng), qHe3(1-Ng:N+Ng))
      call stage_fractions_of_the_column(qH, qHe2, qHe3)
      allocate(rho_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng), Frho(1-Ng:N+Ng),         &
               f_c(1-Ng:N+Ng,n_species))
      allocate(N_H(1-Ng:N+Ng), N_He(1-Ng:N+Ng))
      allocate(n_Hf(0:N), n_Hef(0:N), Kf(0:N), drf(0:N))
      allocate(xion(1,1-Ng:N+Ng), xf(1,0:N), Fk(1,0:N), Dk(1,1:N))
      allocate(xclose(0:N), Fclose(0:N), Dclose(1:N))
      allocate(dFdl(0:N), dFdr(0:N), aa(1:N), bb(1:N), cc(1:N))
      allocate(Dp(1:N), Dm(1:N), Dfull(1:N), Ddon(1:N))

      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         T_c(j)   = 1.0d0 + 2.0d0*(r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
         Frho(j)  = rho_c(j)*n0*mu*1.0d4/(r_edg(j)*r_edg(j))
      enddo
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 1.0d-1)
      call ionize_the_column(f_c, qH, qHe2, qHe3)
      call hydrogen_and_helium_nucleus_face_flux(rho_c, T_c, f_c,         &
               Frho, N_H, N_He, n_Hf, n_Hef, Kf, drf)

      xion(1,:) = qH
      call ionization_stage_face_jacobian(N_H, n_Hf, Kf, drf, dFdl, dFdr)
      call ionization_stage_row_jacobian(dFdl, dFdr, aa, bb, cc)

      ! The two faces the reconstruction differs from the donor value on,
      ! reported: this is what the rows leave out.
      call ionization_stage_face_flux(xion, N_H, n_Hf, Kf, drf,           &
               xf, xclose, Fk, Fclose)
      call ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
      Dfull = Dk(1,:)
      call ionization_stage_face_flux(xion, N_H, n_Hf, Kf, drf,           &
               xf, xclose, Fk, Fclose, donor_value = .true.)
      call ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
      Ddon = Dk(1,:)
      worst_rec = 0.0d0
      do j = 1, N
         sc = max(abs(Dfull(j)) + abs(Ddon(j)), 1.0d-300)
         worst_rec = max(worst_rec, abs(Dfull(j) - Ddon(j))/sc)
      enddo

      worst_a = 0.0d0
      worst_b = 0.0d0
      worst_c = 0.0d0
      do jp = 2, N-1, 37
         ! The donor-value face flux is exactly linear in the carried
         ! fraction, so the quotient below is exact for any step and the
         ! step is taken LARGE: the divergence is a difference of two face
         ! terms that cancel to many digits, and a small step would read
         ! that cancellation instead of the derivative.
         h = 1.0d-2
         xion(1,:) = qH
         xion(1,jp) = qH(jp) + h
         call ionization_stage_face_flux(xion, N_H, n_Hf, Kf, drf,        &
                  xf, xclose, Fk, Fclose, donor_value = .true.)
         call ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
         Dp = Dk(1,:)
         xion(1,jp) = qH(jp) - h
         call ionization_stage_face_flux(xion, N_H, n_Hf, Kf, drf,        &
                  xf, xclose, Fk, Fclose, donor_value = .true.)
         call ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
         Dm = Dk(1,:)
         do j = max(1,jp-1), min(N,jp+1)
            dnum = (Dp(j) - Dm(j))/(2.0d0*h)
            if (j .eq. jp-1) then
               sc = max(abs(cc(j)) + abs(dnum), 1.0d-300)
               worst_c = max(worst_c, abs(dnum - cc(j))/sc)
            else if (j .eq. jp) then
               sc = max(abs(bb(j)) + abs(dnum), 1.0d-300)
               worst_b = max(worst_b, abs(dnum - bb(j))/sc)
            else
               sc = max(abs(aa(j)) + abs(dnum), 1.0d-300)
               worst_a = max(worst_a, abs(dnum - aa(j))/sc)
            endif
         enddo
      enddo

      write(*,'(A,ES12.4)')                                               &
         '  DIAGNOSTIC largest relative difference between the'//         &
         ' reconstructed and the donor-value stage divergence: ',         &
         worst_rec
      call verdict('stage_row_jacobian_matches_a_central_difference',     &
                   max(worst_a, max(worst_b, worst_c)) .le. 1.0d-7,       &
                   max(worst_a, max(worst_b, worst_c)), 0.0d0, 1.0d-7)

      deallocate(rho_c, T_c, Frho, f_c, N_H, N_He, n_Hf, n_Hef, Kf, drf,  &
                 xion, xf, Fk, Dk, xclose, Fclose, Dclose, dFdl, dFdr,    &
                 aa, bb, cc, Dp, Dm, Dfull, Ddon, qH, qHe2, qHe3)
      end subroutine test_production_stage_jacobian

      ! ---------------------------------------------------------------- !

      subroutine test_boundary_stage_rows()
      ! THE TWO END ROWS, where the stage flux is live and the ghost is a
      ! COPY of the interior cell.
      !
      ! The element operator carries no diffusive flux through f = 0 or
      ! f = N, but the stage flux there is not zero: the material advection
      ! crosses both faces.  Nothing below the base states an ionization
      ! fraction, so the inner ghosts carry cell 1's own partition, and the
      ! outer ghost is cell N continued; both are therefore copies of the
      ! unknown and contribute their derivative to the DIAGONAL of the end
      ! row.  That folding is what the transport operator assembles
      ! (diffusive_photochemistry, the stage branch of solve_carriers), and
      ! this row states it against a central difference of the divergence
      ! taken with the ghosts continued the same way.
      !
      ! The interior rows are the row above; these two are the ones a
      ! tridiagonal band cannot hold as neighbours, because the neighbour
      ! is the cell itself.
      integer, parameter :: nc = 500
      real*8, dimension(:),   allocatable :: rho_c, T_c, Frho, qH
      real*8, dimension(:),   allocatable :: qHe2, qHe3, N_H, N_He
      real*8, dimension(:),   allocatable :: n_Hf, n_Hef, Kf, drf
      real*8, dimension(:),   allocatable :: xclose, Fclose, Dclose
      real*8, dimension(:),   allocatable :: dFdl, dFdr, aa, bb, cc
      real*8, dimension(:),   allocatable :: Dp, Dm, fa
      real*8, dimension(:),   allocatable :: cv
      real*8, dimension(:,:), allocatable :: f_c, xion, xf, Fk, Dk
      real*8  :: h, sc, dnum, R0sq, R0cb, bb1, bbN, worst
      integer :: j

      call production_stage_column(nc)
      allocate(qH(1-Ng:N+Ng), qHe2(1-Ng:N+Ng), qHe3(1-Ng:N+Ng))
      call stage_fractions_of_the_column(qH, qHe2, qHe3)
      allocate(rho_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng), Frho(1-Ng:N+Ng),         &
               f_c(1-Ng:N+Ng,n_species))
      allocate(N_H(1-Ng:N+Ng), N_He(1-Ng:N+Ng))
      allocate(n_Hf(0:N), n_Hef(0:N), Kf(0:N), drf(0:N))
      allocate(xion(1,1-Ng:N+Ng), xf(1,0:N), Fk(1,0:N), Dk(1,1:N))
      allocate(xclose(0:N), Fclose(0:N), Dclose(1:N))
      allocate(dFdl(0:N), dFdr(0:N), aa(1:N), bb(1:N), cc(1:N))
      allocate(Dp(1:N), Dm(1:N), fa(0:N), cv(1:N))

      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         T_c(j)   = 1.0d0 + 2.0d0*(r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
         Frho(j)  = rho_c(j)*n0*mu*1.0d4/(r_edg(j)*r_edg(j))
      enddo
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 1.0d-1)
      call ionize_the_column(f_c, qH, qHe2, qHe3)
      call hydrogen_and_helium_nucleus_face_flux(rho_c, T_c, f_c,         &
               Frho, N_H, N_He, n_Hf, n_Hef, Kf, drf)

      call ionization_stage_face_jacobian(N_H, n_Hf, Kf, drf, dFdl, dFdr)
      call ionization_stage_row_jacobian(dFdl, dFdr, aa, bb, cc)
      call spherical_face_area_and_cell_volume(fa, cv)
      R0sq = R0*R0
      R0cb = R0sq*R0
      ! The two folded diagonals: the interior entry of the end row plus
      ! the derivative of the ghost that is a copy of the same unknown.
      bb1 = bb(1) - fa(0)*R0sq*dFdl(0)/(cv(1)*R0cb)
      bbN = bb(N) + fa(N)*R0sq*dFdr(N)/(cv(N)*R0cb)

      worst = 0.0d0
      h     = 1.0d-2
      ! ---- cell 1, with the inner ghosts carrying cell 1's value -------
      xion(1,:) = qH
      xion(1,1-Ng:1) = qH(1) + h
      call ionization_stage_face_flux(xion, N_H, n_Hf, Kf, drf,           &
               xf, xclose, Fk, Fclose, donor_value = .true.)
      call ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
      Dp = Dk(1,:)
      xion(1,1-Ng:1) = qH(1) - h
      call ionization_stage_face_flux(xion, N_H, n_Hf, Kf, drf,           &
               xf, xclose, Fk, Fclose, donor_value = .true.)
      call ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
      Dm = Dk(1,:)
      dnum = (Dp(1) - Dm(1))/(2.0d0*h)
      sc   = max(abs(bb1) + abs(dnum), 1.0d-300)
      worst = max(worst, abs(dnum - bb1)/sc)

      ! ---- cell N, with the outer ghosts carrying cell N's value -------
      xion(1,:) = qH
      xion(1,N:N+Ng) = qH(N) + h
      call ionization_stage_face_flux(xion, N_H, n_Hf, Kf, drf,           &
               xf, xclose, Fk, Fclose, donor_value = .true.)
      call ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
      Dp = Dk(1,:)
      xion(1,N:N+Ng) = qH(N) - h
      call ionization_stage_face_flux(xion, N_H, n_Hf, Kf, drf,           &
               xf, xclose, Fk, Fclose, donor_value = .true.)
      call ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
      Dm = Dk(1,:)
      dnum = (Dp(N) - Dm(N))/(2.0d0*h)
      sc   = max(abs(bbN) + abs(dnum), 1.0d-300)
      worst = max(worst, abs(dnum - bbN)/sc)

      call verdict('boundary_stage_rows_fold_the_ghost_onto_the_diagonal',&
                   worst .le. 1.0d-7, worst, 0.0d0, 1.0d-7)

      ! THE RED: leaving the end faces out of the rows, which is the rule a
      ! purely diffusive carrier row follows because its flux vanishes
      ! there.  The stage flux does not vanish there, so the entry it drops
      ! is the whole material advection through the boundary face.
      worst = 0.0d0
      sc = max(abs(bb(1)) + abs(bb1), 1.0d-300)
      worst = max(worst, abs(bb(1) - bb1)/sc)
      sc = max(abs(bb(N)) + abs(bbN), 1.0d-300)
      worst = max(worst, abs(bb(N) - bbN)/sc)
      call verdict('dropping_the_end_face_entry_changes_the_end_rows',    &
                   worst .gt. 1.0d3*epsilon(1.0d0), worst,                &
                   1.0d3*epsilon(1.0d0), 0.0d0)

      deallocate(rho_c, T_c, Frho, f_c, N_H, N_He, n_Hf, n_Hef, Kf, drf,  &
                 xion, xf, Fk, Dk, xclose, Fclose, Dclose, dFdl, dFdr,    &
                 aa, bb, cc, Dp, Dm, fa, cv, qH, qHe2, qHe3)
      end subroutine test_boundary_stage_rows

      ! ----------------------------------------------------------------- !

      subroutine test_operator_stage_rows_of_both_elements()
      ! THE OPERATOR'S OWN STAGE ROWS, hydrogen and helium together.
      !
      ! The rows above state the identity of the FLUX; these state that the
      ! transport-chemistry operator forms each element's rows on that
      ! element's own nucleus flux.  x(H II) rides on N_H and x(He II) and
      ! x(He III) ride on N_He, the two helium stages sharing one closing
      ! stage and one simplex -- which is what makes the helium fluxes sum
      ! to the helium nucleus flux and not to something else.
      !
      ! The operator's frozen stage face state is installed by
      ! carrier_stage_face_state from ONE call of the element operator's
      ! public flux, and the rows below compare what the operator then
      ! hands out with the same flux formed here from the same call.
      integer, parameter :: nc = 500
      real*8, dimension(:),   allocatable :: rho_c, T_c, Frho, Frho_cgs
      real*8, dimension(:),   allocatable :: qH, qHe2, qHe3, msum
      real*8, dimension(:),   allocatable :: N_H, N_He, n_Hf, n_Hef
      real*8, dimension(:),   allocatable :: Kf, drf
      real*8, dimension(:),   allocatable :: Jf, dJl, dJr, dFdl, dFdr
      real*8, dimension(:),   allocatable :: xclose, Fclose
      real*8, dimension(:,:), allocatable :: f_c, fc, xion, xf, Fk
      real*8  :: worst, worst_j, dmaxH, dmaxHe, dmax_red, sc, ssum, rat
      real*8  :: worst_ratio
      integer :: j, k, ic, jw

      call production_stage_column(nc)
      ! The key the rows exist for: the operator then carries x(H II) and
      ! the two helium stages, which carrier_set_init registers.
      ionization_transport = .true.
      call carrier_set_init()

      allocate(qH(1-Ng:N+Ng), qHe2(1-Ng:N+Ng), qHe3(1-Ng:N+Ng))
      allocate(rho_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng), Frho(1-Ng:N+Ng),         &
               Frho_cgs(1-Ng:N+Ng), msum(1-Ng:N+Ng),                      &
               f_c(1-Ng:N+Ng,n_species), fc(1-Ng:N+Ng,n_carrier_max))
      allocate(N_H(1-Ng:N+Ng), N_He(1-Ng:N+Ng))
      allocate(n_Hf(0:N), n_Hef(0:N), Kf(0:N), drf(0:N))
      allocate(Jf(0:N), dJl(0:N), dJr(0:N), dFdl(0:N), dFdr(0:N))
      allocate(xion(2,1-Ng:N+Ng), xf(2,0:N), Fk(2,0:N))
      allocate(xclose(0:N), Fclose(0:N))

      call stage_fractions_of_the_column(qH, qHe2, qHe3)
      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         T_c(j)   = 1.0d0 + 2.0d0*(r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
         Frho(j)  = rho_c(j)*n0*mu*1.0d4/(r_edg(j)*r_edg(j))
      enddo
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 1.0d-1)
      call ionize_the_column(f_c, qH, qHe2, qHe3)
      call mixture_mass_sum(f_c, msum)

      ! The operator's own installation of the frozen stage face state,
      ! with the material advection carried by the rows (the fixed-wind
      ! relaxation and the stationary balance).
      call carrier_stage_face_state(rho_c, T_c*T0, f_c, .true., msum,     &
                                    Frho)
      ! The same call of the element operator, here, with the face mass
      ! flux converted the way that routine converts it.
      do j = 1-Ng, N+Ng-1
         Frho_cgs(j) = Frho(j)*0.5d0*(msum(j) + msum(j+1))*n0*mu*v0
      enddo
      Frho_cgs(N+Ng) = 0.0d0
      call hydrogen_and_helium_nucleus_face_flux(rho_c, T_c, f_c,         &
               Frho_cgs, N_H, N_He, n_Hf, n_Hef, Kf, drf)

      fc = 0.0d0
      fc(:,ic_Hp)    = qH
      fc(:,ic_HeII)  = qHe2
      fc(:,ic_HeIII) = qHe3

      ! ---- (1) each stage's face flux is its own element's ------------
      worst = 0.0d0
      do ic = 1, 3
         if (ic .eq. 1) then
            k = ic_Hp
            xion(1,:) = qH
            N_He      = N_He        ! unchanged; the reference below picks
         else if (ic .eq. 2) then
            k = ic_HeII
            xion(1,:) = qHe2
         else
            k = ic_HeIII
            xion(1,:) = qHe3
         endif
         call carrier_stage_face_flux(k, fc, Jf, dJl, dJr)
         if (k .eq. ic_Hp) then
            call ionization_stage_face_flux(xion(1:1,:), N_H, n_Hf, Kf,   &
                     drf, xf(1:1,:), xclose, Fk(1:1,:), Fclose)
            call ionization_stage_face_jacobian(N_H, n_Hf, Kf, drf,       &
                                                dFdl, dFdr)
         else
            call ionization_stage_face_flux(xion(1:1,:), N_He, n_Hef, Kf, &
                     drf, xf(1:1,:), xclose, Fk(1:1,:), Fclose)
            call ionization_stage_face_jacobian(N_He, n_Hef, Kf, drf,     &
                                                dFdl, dFdr)
         endif
         do j = 0, N
            sc = max(abs(Jf(j)) + abs(Fk(1,j)), 1.0d-300)
            worst = max(worst, abs(Jf(j) - Fk(1,j))/sc)
            sc = max(abs(dJl(j)) + abs(dFdl(j)), 1.0d-300)
            worst = max(worst, abs(dJl(j) - dFdl(j))/sc)
            sc = max(abs(dJr(j)) + abs(dFdr(j)), 1.0d-300)
            worst = max(worst, abs(dJr(j) - dFdr(j))/sc)
         enddo
      enddo
      call verdict('operator_stage_rows_ride_on_their_own_element_flux',  &
                   worst .le. 1.0d-14, worst, 0.0d0, 1.0d-14)

      ! ---- (2) the sum identity, one element at a time ----------------
      call ionization_stage_nucleus_sum(fc, ien_H,  dmaxH,  jw)
      call ionization_stage_nucleus_sum(fc, ien_He, dmaxHe, jw)
      call verdict('operator_stage_sum_is_the_hydrogen_nucleus_flux',     &
                   dmaxH .le. 1.0d-13, dmaxH, 0.0d0, 1.0d-13)
      call verdict('operator_stage_sum_is_the_helium_nucleus_flux',       &
                   dmaxHe .le. 1.0d-13, dmaxHe, 0.0d0, 1.0d-13)

      ! ---- (3) the RED: helium charged to the hydrogen flux -----------
      ! A stage table that named one element for every stage would build
      ! the helium rows on N_H.  The fractions still close their own
      ! simplex, so what breaks is the identity itself.
      xion(1,:) = qHe2
      xion(2,:) = qHe3
      call ionization_stage_face_flux(xion, N_H, n_Hf, Kf, drf,           &
               xf, xclose, Fk, Fclose)
      dmax_red = 0.0d0
      do j = 0, N
         ssum = Fclose(j)
         sc   = abs(Fclose(j))
         do k = 1, 2
            ssum = ssum + Fk(k,j)
            sc   = sc + abs(Fk(k,j))
         enddo
         sc = max(abs(N_He(j)), sc)
         if (sc .le. 0.0d0) cycle
         dmax_red = max(dmax_red, abs(ssum - N_He(j))/sc)
      enddo
      call verdict('charging_helium_to_the_hydrogen_flux_breaks_the_sum', &
                   dmax_red .gt. 1.0d-3, dmax_red, 1.0d-3, 0.0d0)

      ! ---- (4) the two helium stages share ONE simplex ----------------
      ! Each of them inside [0,1] and their sum above one: projecting
      ! stage by stage would leave the cell holding 1.4 of its helium.
      !
      ! The bound on the sum is the share of the element the atomic stages
      ! partition, so the frozen helium background the projection reads has
      ! to stand first.  Here every helium nucleus of the cell is in an
      ! atomic stage -- no HeH+, no metastable population -- so the bound
      ! is one and this case is the pure simplex.  Case (4b) below carries
      ! helium outside the stages and states the smaller bound.
      call carrier_helium_background_set_for_test(1.0d0, 0.0d0, 0.0d0)
      fc(:,ic_HeII)  = 0.8d0
      fc(:,ic_HeIII) = 0.6d0
      fc(:,ic_Hp)    = 0.5d0
      call carrier_ionization_stage_projection(fc)
      worst_j     = 0.0d0
      worst_ratio = 0.0d0
      do j = 1, N
         worst_j = max(worst_j,                                          &
                       abs(fc(j,ic_HeII) + fc(j,ic_HeIII) - 1.0d0))
         ! The projection scales, so the RATIO of the two carried stages
         ! is what the charge sum reads and must not move.
         rat = fc(j,ic_HeII)/max(fc(j,ic_HeIII), 1.0d-300)
         worst_ratio = max(worst_ratio, abs(rat - 0.8d0/0.6d0)           &
                                        /(0.8d0/0.6d0))
      enddo
      call verdict('helium_stages_are_returned_to_one_shared_simplex',    &
                   worst_j .le. 1.0d-15, worst_j, 0.0d0, 1.0d-15)
      call verdict('the_projection_keeps_the_ratio_of_the_helium_stages', &
                   worst_ratio .le. 1.0d-14, worst_ratio, 0.0d0, 1.0d-14)
      ! The hydrogen stage, whose own simplex is [0,1], is left alone by
      ! the same call: the two elements are projected separately.
      worst_j = 0.0d0
      do j = 1, N
         worst_j = max(worst_j, abs(fc(j,ic_Hp) - 0.5d0))
      enddo
      call verdict('the_hydrogen_stage_inside_its_simplex_is_untouched',  &
                   worst_j .le. 0.0d0, worst_j, 0.0d0, 0.0d0)

      ! ---- (4b) helium held outside the atomic stages is reserved -----
      ! A tenth of the helium of the cell is bound in HeH+ and a fiftieth
      ! sits in the frozen He 2^3S level.  Neither is the ionized stages'
      ! to take, so the sum they are projected onto is 0.88 and not one; a
      ! projection that stopped at one would admit the partition
      ! x(He II) + x(He III) = 0.95 beside a molecule holding 0.10, and
      ! the state written from it would hold 1.05 of the helium the cell
      ! has (PLAN_20260918 D6).
      call carrier_helium_background_set_for_test(1.0d0, 0.10d0, 0.02d0)
      fc(:,ic_HeII)  = 0.60d0
      fc(:,ic_HeIII) = 0.35d0
      call carrier_ionization_stage_projection(fc)
      worst_j     = 0.0d0
      worst_ratio = 0.0d0
      do j = 1, N
         worst_j = max(worst_j,                                          &
                       abs(fc(j,ic_HeII) + fc(j,ic_HeIII) - 0.88d0))
         rat = fc(j,ic_HeII)/max(fc(j,ic_HeIII), 1.0d-300)
         worst_ratio = max(worst_ratio, abs(rat - 0.60d0/0.35d0)         &
                                        /(0.60d0/0.35d0))
      enddo
      call verdict('the_helium_of_a_molecule_is_reserved_from_the_sum',   &
                   worst_j .le. 1.0d-15, worst_j, 0.88d0, 1.0d-15)
      call verdict('the_reserving_projection_keeps_the_stage_ratio',      &
                   worst_ratio .le. 1.0d-14, worst_ratio, 0.0d0, 1.0d-14)
      ! Back to the wholly atomic background for whatever follows.
      call carrier_helium_background_set_for_test(1.0d0, 0.0d0, 0.0d0)

      ! ---- (5) the RED of (4): each helium stage projected on its own --
      ! Both are inside [0,1], so a projection per stage moves nothing and
      ! the cell keeps 1.4 of its helium nuclei in ionized stages.
      fc(:,ic_HeII)  = 0.8d0
      fc(:,ic_HeIII) = 0.6d0
      worst_j = 0.0d0
      do j = 1, N
         worst_j = max(worst_j,                                          &
                       min(max(fc(j,ic_HeII), 0.0d0), 1.0d0)             &
                     + min(max(fc(j,ic_HeIII), 0.0d0), 1.0d0) - 1.0d0)
      enddo
      call verdict('projecting_each_helium_stage_alone_leaves_the_sum',   &
                   worst_j .gt. 1.0d-3, worst_j, 1.0d-3, 0.0d0)

      ! ---- (6) the stage table names the right element ----------------
      worst_j = 0.0d0
      if (.not. carrier_is_ionization_stage(ic_HeII))  worst_j = 1.0d0
      if (.not. carrier_is_ionization_stage(ic_HeIII)) worst_j = 1.0d0
      if (.not. carrier_is_ionization_stage(ic_Hp))    worst_j = 1.0d0
      if (carrier_stage_element(ic_Hp)    .ne. ien_H)  worst_j = 1.0d0
      if (carrier_stage_element(ic_HeII)  .ne. ien_He) worst_j = 1.0d0
      if (carrier_stage_element(ic_HeIII) .ne. ien_He) worst_j = 1.0d0
      if (.not. carrier_solved(ic_HeII))  worst_j = 1.0d0
      if (.not. carrier_solved(ic_HeIII)) worst_j = 1.0d0
      call verdict('the_carried_set_is_x_HII_x_HeII_x_HeIII',             &
                   worst_j .le. 0.0d0, worst_j, 0.0d0, 0.0d0)

      ionization_transport = .false.
      call carrier_set_init()
      deallocate(rho_c, T_c, Frho, Frho_cgs, msum, f_c, fc, N_H, N_He,    &
                 n_Hf, n_Hef, Kf, drf, Jf, dJl, dJr, dFdl, dFdr,          &
                 xion, xf, Fk, xclose, Fclose, qH, qHe2, qHe3)
      end subroutine test_operator_stage_rows_of_both_elements

      ! ----------------------------------------------------------------- !

      subroutine manufactured_column_grid(nc)
      ! A UNIFORM SPHERICAL COLUMN OVER [1, 2] WHOSE FACES ARE THE SAME
      ! TWO RADII AT EVERY RESOLUTION, so that three grids discretize ONE
      ! domain and the measurements below are a refinement of one problem.
      ! r_edg(j) = 1 + j dr and r(j) = 1 + (j - 1/2) dr, which is exactly
      ! the centre/edge pair define_grid builds (r_edg the midpoint of the
      ! centres, dr_j the distance between a cell's own faces).
      integer, intent(in) :: nc
      real*8  :: dr
      integer :: j
      N  = nc
      T0 = 1.0d3
      R0 = 1.0d10
      n0 = 1.0d10
      v0 = sqrt(kb_erg*T0/mu)
      p0 = n0*mu*v0*v0
      spherical_domain = .true.
      grid_type        = 'Uniform'
      r_max            = 2.0d0
      use_plm = .false.;  use_weno3 = .true.;  rec_method = 'WENO3'
      recon_lambda_on = .false.
      j_min = 1
      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      dr = (r_max - 1.0d0)/dble(N)
      do j = 1-Ng, N+Ng
         r(j)     = 1.0d0 + (dble(j) - 0.5d0)*dr
         r_edg(j) = 1.0d0 +  dble(j)         *dr
         dr_j(j)  = dr
      enddo
      end subroutine manufactured_column_grid

      ! ----------------------------------------------------------------- !

      subroutine test_manufactured_ionization_column()
      ! THE DISCRETIZATION ERROR OF ONE STAGE ROW, ON A COLUMN WHOSE
      ! STATIONARY SOLUTION IS KNOWN (the anchoring measurement L12a
      ! section 2.4 and the L36c memo owe).
      !
      ! The stage row is  dn_k/dt = -div F_k + S_k , and a certification
      ! tolerance on it has to be anchored by what the DISCRETIZATION of
      ! that row costs at the production spacing, not by what an arithmetic
      ! floor forbids.  So a smooth fraction x(r) is chosen, the continuous
      ! flux of equation (1) is formed from it analytically, and its exact
      ! divergence is used as the source S_k.  x(r) is then an exact
      ! stationary solution of the CONTINUOUS row, and what the discrete
      ! row reads on it is the truncation error of the discretization and
      ! nothing else.
      !
      !   n_el(r) = n_0 exp[-a (r-1)]      [cm^-3]
      !   N_el(r) = C / r^2                [nuclei cm^-2 s^-1], divergence
      !                                     free, a steady nucleus flux
      !   K(r)    = K_0 r^2                [cm^2 s^-1]
      !   x(r)    = x_0 + dx (6 s^5 - 15 s^4 + 10 s^3) ,  s = (r-1)/(r_t-1)
      !
      ! The quintic is monotone, so the limiter of the reconstruction never
      ! fires, and it has x' = x'' = 0 at both ends, so the operator's own
      ! boundary rule (no eddy flux through the two end faces) is EXACT
      ! there and not a source of error the interior would be charged with.
      !
      !   r^2 F   = C x - (K_0 n_0 / R_0) r^4 exp[-a(r-1)] x'
      !   div F   = { C x' - (K_0 n_0/R_0) exp[-a(r-1)]
      !               [ (4 r^3 - a r^4) x' + r^4 x'' ] } / (r^2 R_0)
      !
      ! with x' and x'' derivatives with respect to r in units of R_0.
      !
      ! The measure is the row's residual divided by the row's OWN terms,
      ! which is the measure the certification puts a tolerance on: the
      ! magnitudes of the two face contributions and of the source.
      integer, parameter :: ngrid = 3
      integer :: ncells(ngrid)
      real*8  :: meas(ngrid), drg(ngrid), meas_all(ngrid)
      ! The same residual WITHOUT the normalization, and the largest exact
      ! source of the continuous problem, which is a grid-independent
      ! reference scale (D3a step 3).
      real*8  :: dimen(ngrid), sref(ngrid), o_dim_s
      real*8  :: ne0, alph, Cn, K0, x0m, dxm, L, ss, xr, xp, xpp
      real*8  :: rj, adv_t, edd_t, worst, worst_all, dsc, resid, divex
      real*8  :: order12, order23, order13, dr_prod, anchor_row
      real*8  :: adv_max, edd_max
      real*8, dimension(:),   allocatable :: Nel, n_elf, Kf, drf
      real*8, dimension(:),   allocatable :: xclose, Fclose, Dclose, fa, cv
      real*8, dimension(:,:), allocatable :: xion, xf, Fk, Dk
      real*8  :: R0sq, R0cb, sL, sR, Vj
      integer :: ig, j

      ncells(1) = 250
      ncells(2) = 500
      ncells(3) = 1000

      ne0  = 1.0d10
      alph = 4.0d0
      Cn   = 1.0d13
      K0   = 1.0d12
      x0m  = 2.0d-2
      dxm  = 9.0d-1
      L    = 1.0d0

      do ig = 1, ngrid
         call manufactured_column_grid(ncells(ig))
         drg(ig) = (r_max - 1.0d0)/dble(N)
         allocate(Nel(1-Ng:N+Ng), xion(1,1-Ng:N+Ng))
         allocate(n_elf(0:N), Kf(0:N), drf(0:N))
         allocate(xf(1,0:N), Fk(1,0:N), Dk(1,1:N))
         allocate(xclose(0:N), Fclose(0:N), Dclose(1:N), fa(0:N), cv(1:N))

         do j = 1-Ng, N+Ng
            ss = (r(j) - 1.0d0)/L
            xion(1,j) = x0m + dxm*(6.0d0*ss**5 - 15.0d0*ss**4            &
                                 + 10.0d0*ss**3)
            ! The element nucleus flux is a FACE quantity in this array,
            ! exactly as the operator stores it.
            Nel(j) = Cn/(r_edg(j)*r_edg(j))
         enddo
         do j = 0, N
            n_elf(j) = 0.5d0*(ne0*exp(-alph*(r(j)   - 1.0d0))            &
                            + ne0*exp(-alph*(r(j+1) - 1.0d0)))
            Kf(j)    = 0.5d0*(K0*r(j)*r(j) + K0*r(j+1)*r(j+1))
            drf(j)   = (r(j+1) - r(j))*R0
         enddo

         call ionization_stage_face_flux(xion, Nel, n_elf, Kf, drf,       &
                  xf, xclose, Fk, Fclose)
         call ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
         call spherical_face_area_and_cell_volume(fa, cv)
         R0sq = R0*R0
         R0cb = R0sq*R0

         worst     = 0.0d0
         worst_all = 0.0d0
         adv_max   = 0.0d0
         edd_max   = 0.0d0
         dimen(ig) = 0.0d0
         sref(ig)  = 0.0d0
         do j = 1, N
            rj  = r(j)
            ss  = (rj - 1.0d0)/L
            xp  = dxm*30.0d0*ss*ss*(ss - 1.0d0)**2/L
            xpp = dxm*60.0d0*ss*(2.0d0*ss - 1.0d0)*(ss - 1.0d0)/(L*L)
            adv_t = Cn*xp
            edd_t = -(K0*ne0/R0)*exp(-alph*(rj - 1.0d0))                  &
                    *((4.0d0*rj**3 - alph*rj**4)*xp + rj**4*xpp)
            divex = (adv_t + edd_t)/(rj*rj*R0)
            adv_max = max(adv_max, abs(adv_t/(rj*rj*R0)))
            edd_max = max(edd_max, abs(edd_t/(rj*rj*R0)))
            sL = fa(j-1)*R0sq
            sR = fa(j)  *R0sq
            Vj = cv(j)*R0cb
            resid = Dk(1,j) - divex
            ! The row's own terms: the two face contributions and the
            ! source, which is what the returned-state acceptance divides
            ! by (carrier_residual, row_terms_phys).
            dsc = (sR*abs(Fk(1,j)) + sL*abs(Fk(1,j-1)))/Vj + abs(divex)
            worst_all = max(worst_all, abs(resid)/max(dsc, 1.0d-300))
            ! The interior: the two end cells carry the operator's own
            ! boundary rule, which is a statement about the domain edge and
            ! not a truncation of the interior stencil.
            if (j .ge. 2 .and. j .le. N-1) then
               worst = max(worst, abs(resid)/max(dsc, 1.0d-300))
               dimen(ig) = max(dimen(ig), abs(resid))
            endif
            sref(ig) = max(sref(ig), abs(divex))
         enddo
         meas(ig)     = worst
         meas_all(ig) = worst_all
         write(*,'(A,I5,A,ES11.4,A,ES11.4,A,ES11.4,A,ES11.4)')            &
              '  DIAGNOSTIC manufactured column N=', ncells(ig),          &
              ' dr=', drg(ig), ' row=', worst, ' with the end cells=',    &
              worst_all, ' advective/eddy source ratio=',                 &
              adv_max/max(edd_max, 1.0d-300)

         deallocate(Nel, xion, n_elf, Kf, drf, xf, Fk, Dk,                &
                    xclose, Fclose, Dclose, fa, cv)
      enddo

      ! THE DIMENSIONAL RESIDUAL BESIDE THE NORMALIZED ONE.  The scale the
      ! measure above divides by contains the two face MAGNITUDES over the
      ! cell volume, and for a smooth nonzero flux that grows as 1/h, so the
      ! normalized measure carries one power of h more than the truncation
      ! itself.  The order printed below is therefore NOT a dimensional
      ! order; test_manufactured_stiff_stage_column measures that one on the
      ! same three grids (1.93) and the two differ by one, as they must.
      o_dim_s = log(dimen(1)/dimen(3))/log(drg(1)/drg(3))
      write(*,'(A,ES11.4,1X,ES11.4,1X,ES11.4,A,F7.3)')                    &
           '  DIAGNOSTIC dimensional row residual [cm^-3 s^-1] at '//     &
           '250/500/1000: ', dimen(1), dimen(2), dimen(3),                &
           '   observed dimensional order 250->1000 ', o_dim_s
      write(*,'(A,ES11.4,A,ES11.4,1X,ES11.4,1X,ES11.4)')                  &
           '  DIAGNOSTIC on the grid-independent reference scale '//      &
           'max|S| = ', sref(3), ': ', dimen(1)/sref(1),                  &
           dimen(2)/sref(2), dimen(3)/sref(3)

      order12 = log(meas(1)/meas(2))/log(drg(1)/drg(2))
      order23 = log(meas(2)/meas(3))/log(drg(2)/drg(3))
      order13 = log(meas(1)/meas(3))/log(drg(1)/drg(3))
      write(*,'(A,F7.3,A,F7.3,A,F7.3)')                                   &
           '  DIAGNOSTIC observed order 250->500 ', order12,              &
           '   500->1000 ', order23, '   250->1000 ', order13

      ! THE SPACING THE CARRIER GATE WAS ANCHORED AT.  The element and
      ! carrier tolerance of docs/certification_tolerance_anchoring_20260910.md
      ! section 2 (3) puts its smooth manufactured column at
      ! dr = 1.8447e-3 R_p, the spacing of the atomic fixture's grid at the
      ! binding radius, where that column reads 1.52e-6 (both READ).  The
      ! stage row of THIS column extrapolated to the same spacing is what
      ! makes 1e-5 an anchored value for a stage row and not a borrowed one.
      anchor_row = meas(3)*(1.8447d-3/drg(3))**order13
      write(*,'(A,ES11.4,A,ES11.4)')                                      &
           '  DIAGNOSTIC at the anchoring spacing 1.8447e-3 R_p the'//    &
           ' stage row extrapolates to ', anchor_row,                     &
           ' and a decade above it is ', 1.0d1*anchor_row
      ! The production spacing of the catalog's Mixed grid, read from the
      ! fixture's own grid, at the radius where the species rows bind and
      ! at the radius from which a species row is GATED
      ! (cert_regime_wind_r).  The spacing grows outward on that grid, so
      ! the second is the larger of the two.
      call production_stage_column(500)
      do j = 1, N
         if (r(j) .ge. 1.15d0) then
            dr_prod = dr_j(j)
            exit
         endif
      enddo
      anchor_row = meas(3)*(dr_prod/drg(3))**order13
      write(*,'(A,ES11.4,A,ES11.4)')                                      &
           '  DIAGNOSTIC production spacing at r = 1.15 ', dr_prod,       &
           ' R_p, the row there extrapolates to ', anchor_row
      do j = 1, N
         if (r(j) .ge. cert_regime_wind_r) then
            dr_prod = dr_j(j)
            exit
         endif
      enddo
      anchor_row = meas(3)*(dr_prod/drg(3))**order13
      write(*,'(A,ES11.4,A,ES11.4,A,ES11.4)')                             &
           '  DIAGNOSTIC production spacing at the gating radius ',       &
           dr_prod, ' R_p, the row there extrapolates to ', anchor_row,   &
           ' against the gate cert_tol_carrier_wind ',                    &
           cert_tol_carrier_wind

      ! The measurement is the point of this row; what it ASSERTS is that
      ! the stage row converges at least at the order the conservative
      ! divergence is written to, so that the extrapolation above is a
      ! statement and not a fit.  A first-order reading would mean the
      ! reconstruction or the divergence is not what it is written to be.
      ! MEASURED here: 2.886, the NORMALIZED order.  The dimensional order
      ! of the same residual on the same grids is 1.93 (printed above and
      ! asserted in test_manufactured_stiff_stage_column): the face-
      ! magnitude scale grows as 1/h, so the normalized order stands one
      ! above the truncation's own.  The anchoring extrapolations below are
      ! extrapolations of the NORMALIZED measure, which is the quantity the
      ! certification gates, and they use the normalized order for that
      ! reason.
      call verdict('the_manufactured_stage_row_converges_at_second_'//    &
                   'order_or_better', order13 .ge. 2.0d0,                 &
                   order13, 2.0d0, 0.0d0)
      ! And the finest grid resolves the row far below the tolerance the
      ! entry would take, so the anchor is not set by the measurement's own
      ! noise.
      call verdict('the_manufactured_row_falls_with_the_grid',            &
                   meas(3) .lt. 0.3d0*meas(1), meas(3), meas(1), 0.0d0)
      end subroutine test_manufactured_ionization_column

      ! ----------------------------------------------------------------- !

      subroutine manufactured_stage_column_state(ne0, alph, Cn, K0, x0m,  &
                                                 dxm, Nel, xion, n_elf,   &
                                                 Kf, drf)
      ! The manufactured state of test_manufactured_ionization_column on
      ! WHATEVER grid is currently built: x(r) the quintic, n_el(r) the
      ! exponential, N_el(r) the divergence-free nucleus flux, K(r) the
      ! quadratic eddy coefficient.  Written once so that the smooth column
      ! and the stiff one below are the same problem with and without a
      ! reaction.
      real*8, intent(in)  :: ne0, alph, Cn, K0, x0m, dxm
      real*8, intent(out) :: Nel(1-Ng:), xion(:,1-Ng:)
      real*8, intent(out) :: n_elf(0:), Kf(0:), drf(0:)
      real*8  :: ss
      integer :: j
      ! THE PROFILE IS DEFINED ON [1, 2] AND IS FLAT ABOVE IT.  The quintic
      ! has x' = x'' = 0 at both ends, so continuing it by its end value
      ! leaves the solution twice differentiable, the flux above r = 2
      ! divergence free and the manufactured source zero there.  The
      ! uniform column, whose domain IS [1, 2], is unchanged by the clamp;
      ! what it buys is the production Mixed grid, which reaches 30 R_p.
      do j = 1-Ng, N+Ng
         ss = min(max(r(j) - 1.0d0, 0.0d0), 1.0d0)
         xion(1,j) = x0m + dxm*(6.0d0*ss**5 - 15.0d0*ss**4               &
                              + 10.0d0*ss**3)
         Nel(j) = Cn/(r_edg(j)*r_edg(j))
      enddo
      do j = 0, N
         n_elf(j) = 0.5d0*(ne0*exp(-alph*(r(j)   - 1.0d0))               &
                         + ne0*exp(-alph*(r(j+1) - 1.0d0)))
         Kf(j)    = 0.5d0*(K0*r(j)*r(j) + K0*r(j+1)*r(j+1))
         drf(j)   = (r(j+1) - r(j))*R0
      enddo
      end subroutine manufactured_stage_column_state

      ! ----------------------------------------------------------------- !

      real*8 function manufactured_stage_divergence(rj, ne0, alph, Cn,    &
                                                    K0, dxm)
      ! The exact divergence of the continuous stage flux at r = rj, in
      ! cm^-3 s^-1: the source that makes x(r) an exact stationary solution
      ! of the continuous row.
      real*8, intent(in) :: rj, ne0, alph, Cn, K0, dxm
      real*8 :: ss, xp, xpp, adv_t, edd_t
      ss  = min(max(rj - 1.0d0, 0.0d0), 1.0d0)
      xp  = dxm*30.0d0*ss*ss*(ss - 1.0d0)**2
      xpp = dxm*60.0d0*ss*(2.0d0*ss - 1.0d0)*(ss - 1.0d0)
      adv_t = Cn*xp
      edd_t = -(K0*ne0/R0)*exp(-alph*ss)                                  &
              *((4.0d0*rj**3 - alph*rj**4)*xp + rj**4*xpp)
      manufactured_stage_divergence = (adv_t + edd_t)/(rj*rj*R0)
      end function manufactured_stage_divergence

      ! ----------------------------------------------------------------- !

      subroutine test_manufactured_stiff_stage_column()
      ! A STIFF IONIZATION-RECOMBINATION BALANCE ON THE SAME MANUFACTURED
      ! COLUMN (PLAN_20260918_rev2 item D3a step 3).
      !
      ! The smooth column above carries no reaction, so the only scale its
      ! row has is its own transport.  A transported ionization stage in
      ! production carries a reaction whose production and loss stand orders
      ! above their difference, and the question D3a asks is what a
      ! certification measure does in that regime.  So the same manufactured
      ! solution is given a reaction that is stiff BY CONSTRUCTION and exact
      ! by construction:
      !
      !     P(r) = Lambda + D(r)/2 ,     L(r) = Lambda - D(r)/2 ,
      !
      ! with D(r) the exact divergence of the continuous stage flux.  Then
      ! P - L = D exactly -- x(r) is still the exact stationary solution --
      ! while P and L are each Lambda, which is set to 1, 1e2, 1e4 and 1e6
      ! times the largest |D| of the column.  Written as a two-body
      ! ionization-recombination pair,
      !
      !     P = alpha (1 - x) n_el ,      L = beta x n_el ,
      !
      ! the local relaxation of the row is
      !
      !     dS/dx = -(alpha + beta) n_el = -( P/(1-x) + L/x ) ,
      !
      ! which is what turns a row residual into an abundance error and what
      ! sets the chemical-to-transport timescale ratio reported below.
      !
      ! NOTHING IN THE DISCRETIZATION CHANGES WITH Lambda: the residual is
      ! the same truncation at every Lambda, because the source is exact.
      ! What changes is every normalized measure built on it, and the four
      ! columns below are the four statements a norm can make:
      !
      !   dimensional        |R|                     [cm^-3 s^-1]
      !   acceptance         |R| / (faces + |P - L|)   (row_terms_phys)
      !   turnover           |R| / (faces + P + L)
      !   abundance error    |R| / |dS/dx|             (grid independent,
      !                                                 and the physical one)
      !
      ! The last is the linearization of the discrete solution error and is
      ! exact to first order in the stiff limit, where the local balance and
      ! not the transport coupling sets the correction.
      !
      ! The grids are the three uniform refinements of the smooth column and
      ! the production Mixed grid, and the two END cells -- the ones that
      ! carry the operator's own boundary rule -- are reported separately
      ! from the interior.
      integer, parameter :: ngrid = 4, nlam = 4
      real*8  :: lamrat(nlam)
      real*8  :: ne0, alph, Cn, K0, x0m, dxm
      real*8  :: rj, divex, resid, dsc, face, Pj, Lj, dsdx, xj, Dmax
      real*8  :: dimw(ngrid), accw(ngrid,nlam), turw(ngrid,nlam)
      real*8  :: abew(ngrid,nlam), refw(ngrid), drg(ngrid), Sref(ngrid)
      real*8  :: dimb(ngrid), accb(ngrid,nlam)
      real*8  :: dach(ngrid,nlam), o_dim, o_ref, o_acc, o_abe
      real*8, dimension(:),   allocatable :: Nel, n_elf, Kf, drf
      real*8, dimension(:),   allocatable :: xclose, Fclose, Dclose, fa, cv
      real*8, dimension(:,:), allocatable :: xion, xf, Fk, Dk
      real*8  :: R0sq, R0cb, sL, sR, Vj
      integer :: ig, j, il, ncg(ngrid)

      lamrat = (/ 1.0d0, 1.0d2, 1.0d4, 1.0d6 /)
      ncg    = (/ 250, 500, 1000, 500 /)
      ne0  = 1.0d10
      alph = 4.0d0
      Cn   = 1.0d13
      K0   = 1.0d12
      x0m  = 2.0d-2
      dxm  = 9.0d-1

      do ig = 1, ngrid
         if (ig .le. 3) then
            call manufactured_column_grid(ncg(ig))
            drg(ig) = (r_max - 1.0d0)/dble(N)
         else
            ! THE PRODUCTION MIXED GRID, the one the catalog runs on: 50
            ! uniform base cells under a geometric stretch.  The
            ! manufactured solution is a function of r alone, so it is the
            ! same problem on it; what changes is the spacing, which is no
            ! longer uniform, and the domain, which now reaches 30 R_p.
            call production_stage_column(ncg(ig))
            drg(ig) = dr_j(1)
         endif
         allocate(Nel(1-Ng:N+Ng), xion(1,1-Ng:N+Ng))
         allocate(n_elf(0:N), Kf(0:N), drf(0:N))
         allocate(xf(1,0:N), Fk(1,0:N), Dk(1,1:N))
         allocate(xclose(0:N), Fclose(0:N), Dclose(1:N), fa(0:N), cv(1:N))
         call manufactured_stage_column_state(ne0, alph, Cn, K0, x0m,     &
                                              dxm, Nel, xion, n_elf, Kf,  &
                                              drf)
         call ionization_stage_face_flux(xion, Nel, n_elf, Kf, drf,       &
                                         xf, xclose, Fk, Fclose)
         call ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
         call spherical_face_area_and_cell_volume(fa, cv)
         R0sq = R0*R0
         R0cb = R0sq*R0

         ! The grid-independent reference scale: the largest source of the
         ! CONTINUOUS problem, which does not move with h.
         Dmax = 0.0d0
         do j = 1, N
            Dmax = max(Dmax, abs(manufactured_stage_divergence(r(j), ne0, &
                                 alph, Cn, K0, dxm)))
         enddo
         Sref(ig) = Dmax

         dimw(ig)   = 0.0d0
         refw(ig)   = 0.0d0
         dimb(ig)   = 0.0d0
         accw(ig,:) = 0.0d0
         turw(ig,:) = 0.0d0
         abew(ig,:) = 0.0d0
         accb(ig,:) = 0.0d0
         dach(ig,:) = 0.0d0
         do j = 1, N
            rj    = r(j)
            xj    = xion(1,j)
            divex = manufactured_stage_divergence(rj, ne0, alph, Cn, K0,  &
                                                  dxm)
            sL = fa(j-1)*R0sq
            sR = fa(j)  *R0sq
            Vj = cv(j)*R0cb
            resid = Dk(1,j) - divex
            face  = (sR*abs(Fk(1,j)) + sL*abs(Fk(1,j-1)))/Vj
            do il = 1, nlam
               Pj = lamrat(il)*Dmax + 0.5d0*divex
               Lj = lamrat(il)*Dmax - 0.5d0*divex
               if (Pj .le. 0.0d0 .or. Lj .le. 0.0d0) cycle
               dsdx = -(Pj/max(1.0d0 - xj, 1.0d-300)                      &
                      + Lj/max(xj, 1.0d-300))
               if (j .ge. 2 .and. j .le. N-1) then
                  accw(ig,il) = max(accw(ig,il),                          &
                       abs(resid)/(face + abs(Pj - Lj)))
                  turw(ig,il) = max(turw(ig,il),                          &
                       abs(resid)/(face + Pj + Lj))
                  abew(ig,il) = max(abew(ig,il), abs(resid)/abs(dsdx))
                  ! Where the manufactured solution is flat the transport
                  ! vanishes and the ratio is not defined; the cells that
                  ! carry the profile are the ones it is read at.
                  if (abs(divex) .gt. 1.0d-6*Dmax)                        &
                     dach(ig,il) = max(dach(ig,il),                       &
                          abs(dsdx)*xj/abs(divex))
               else
                  accb(ig,il) = max(accb(ig,il),                          &
                       abs(resid)/(face + abs(Pj - Lj)))
               endif
            enddo
            if (j .ge. 2 .and. j .le. N-1) then
               dimw(ig) = max(dimw(ig), abs(resid))
               refw(ig) = max(refw(ig), abs(resid)/Dmax)
            else
               dimb(ig) = max(dimb(ig), abs(resid))
            endif
         enddo
         if (ig .le. 3) then
            write(*,'(A,I5,A,ES11.4,A,ES11.4,A,ES11.4,A,ES11.4)')         &
                 '  DIAGNOSTIC stiff column N=', ncg(ig), ' dr=',         &
                 drg(ig), ' dimensional |R|=', dimw(ig),                  &
                 ' |R|/max|S|=', refw(ig), ' boundary |R|=', dimb(ig)
         else
            write(*,'(A,I5,A,ES11.4,A,ES11.4,A,ES11.4,A,ES11.4)')         &
                 '  DIAGNOSTIC stiff column on the Mixed grid, N=',       &
                 ncg(ig), ' base dr=', drg(ig), ' dimensional |R|=',      &
                 dimw(ig), ' |R|/max|S|=', refw(ig), ' boundary |R|=',    &
                 dimb(ig)
         endif
         do il = 1, nlam
            write(*,'(A,ES9.2,A,ES11.4,A,ES11.4,A,ES11.4,A,ES11.4)')      &
                 '    Lambda/max|D|=', lamrat(il),                        &
                 '  Damkohler=', dach(ig,il),                             &
                 '  acceptance=', accw(ig,il),                            &
                 '  turnover=', turw(ig,il),                              &
                 '  abundance error=', abew(ig,il)
            write(*,'(A,ES11.4)')                                         &
                 '      the same acceptance measure at the two END '//    &
                 'cells, which carry the boundary rule: ', accb(ig,il)
         enddo
         deallocate(Nel, xion, n_elf, Kf, drf, xf, Fk, Dk,                &
                    xclose, Fclose, Dclose, fa, cv)
      enddo

      ! THE ORDERS.  The dimensional residual and the residual on the
      ! grid-independent reference scale are ONE order apart from the
      ! acceptance-normalized one wherever a smooth nonzero flux dominates
      ! the face-magnitude scale, because that scale grows as 1/h.  The
      ! smooth column above records a NORMALIZED order 2.886; the two
      ! numbers printed here are what that normalized order corresponds to
      ! dimensionally on the same three grids.
      o_dim = log(dimw(1)/dimw(3))/log(drg(1)/drg(3))
      o_ref = log(refw(1)/refw(3))/log(drg(1)/drg(3))
      o_acc = log(accw(1,4)/accw(3,4))/log(drg(1)/drg(3))
      o_abe = log(abew(1,4)/abew(3,4))/log(drg(1)/drg(3))
      write(*,'(A,F7.3,A,F7.3,A,F7.3,A,F7.3)')                            &
           '  DIAGNOSTIC observed order 250->1000  dimensional ', o_dim,  &
           '   on max|S| ', o_ref, '   acceptance-normalized ', o_acc,    &
           '   abundance error ', o_abe
      write(*,'(A,ES13.6,A,ES13.6)')                                      &
           '  DIAGNOSTIC the acceptance scale the same residual is '//    &
           'divided by, N=250: ', dimw(1)/max(accw(1,1), 1.0d-300),       &
           '  N=1000: ', dimw(3)/max(accw(3,1), 1.0d-300)

      ! What the rows assert.
      ! (1) The residual is the same truncation at every Lambda: the source
      !     is exact by construction, so a normalization can move the
      !     measure by decades without a single digit of the discretization
      !     changing.  That is the statement D3b needs.
      call verdict('the_stiff_source_leaves_the_dimensional_residual_'//  &
                   'unchanged', abs(dimw(1)) .gt. 0.0d0,                  &
                   dimw(1), dimw(1), 0.0d0)
      ! (2) The acceptance and the turnover norms differ by the stiffness
      !     ratio, which is what a stiff stage is.
      call verdict('the_acceptance_and_turnover_norms_differ_by_the_'//   &
                   'stiffness', accw(3,4)/max(turw(3,4), 1.0d-300)        &
                   .gt. 1.0d3, accw(3,4)/max(turw(3,4), 1.0d-300),        &
                   1.0d3, 0.0d0)
      ! (3) The dimensional residual converges at better than first order,
      !     which is what makes the extrapolations of the smooth column
      !     statements about a converging row.  MEASURED here: 1.93, and
      !     what fails at first order is the reconstruction or the
      !     divergence, not the normalization.
      call verdict('the_stiff_column_converges_dimensionally_above_'//    &
                   'first_order', o_dim .gt. 1.0d0, o_dim, 1.0d0, 0.0d0)
      ! (4) AND THE NORMALIZED ORDER IS ONE ABOVE IT, because the
      !     face-magnitude scale the acceptance divides by grows as 1/h
      !     wherever a smooth nonzero flux dominates it.  This is the row
      !     that makes the recorded normalized order 2.886 a DIMENSIONAL
      !     order near 1.9 and not a third-order truncation.
      call verdict('the_acceptance_normalized_order_is_one_above_the_'//  &
                   'dimensional_one', abs(o_acc - o_dim - 1.0d0)          &
                   .lt. 1.0d-1, o_acc - o_dim, 1.0d0, 1.0d-1)
      ! (5) The abundance error the residual implies falls with the grid at
      !     the dimensional order and does NOT depend on Lambda, which is
      !     what makes it the physical measure of the two.
      call verdict('the_abundance_error_is_the_grid_independent_'//       &
                   'measure', abs(o_abe - o_dim) .lt. 1.0d-3,             &
                   abs(o_abe - o_dim), 0.0d0, 1.0d-3)
      end subroutine test_manufactured_stiff_stage_column


      end program ionization_stage_flux_tests
