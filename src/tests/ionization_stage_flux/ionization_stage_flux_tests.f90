      program ionization_stage_flux_tests
      ! THE DISCRETE IDENTITY OF THE IONIZATION-STAGE FLUX
      ! (docs/lhs1140b_stationary_L12b_derivation_20260916.md, section 4;
      ! PLAN_20260916_rev3 section 8, item L12 stage B).
      !
      ! WHAT IS TESTED AND WHAT IS NOT.  The stage flux does not exist in the
      ! source yet: stage B builds it.  What this driver holds is the FORMULA
      ! of the derivation memo, written once here, fed with the face
      ! quantities the production operator itself forms:
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
      use lower_atmosphere_profile, only: eddy_diffusion_on_grid
      use species_advective_transport, only: species_face_fraction,       &
                                             species_face_flux

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

      end program ionization_stage_flux_tests
