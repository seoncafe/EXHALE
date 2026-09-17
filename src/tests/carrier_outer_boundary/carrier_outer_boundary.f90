      program carrier_outer_boundary
      ! THE OUTER BOUNDARY OF THE CARRIER COLUMN AND THE DEFERRED
      ! RECONSTRUCTION TERMS OF ITS JACOBIAN, stated as measurements.
      !
      ! Three properties of the operator this suite exists for:
      !
      !   1. The diffusive, eddy and settling face coefficients are formed
      !      for the interior faces 1 .. N-1 alone, so neither the face at
      !      r_{1/2} nor the face at r_{N+1/2} carries such a flux: the
      !      carrier column is a CLOSED BOX for diffusion while the wind
      !      still advects across both faces.  At the outer end the reason
      !      is that the gas there is not collisional -- the bulk Knudsen
      !      number of these domains reaches 1.4 to 2.4 at the top -- so a
      !      Chapman-Enskog coefficient is not a fluid quantity; at the
      !      inner end the composition the layer below hands over enters on
      !      the inflowing base face and a diffusive flux would state the
      !      same handoff twice.
      !
      !   2. ONE OUTER GHOST RULE.  The outflow continuation X_ghost = X_N
      !      is applied by carrier_face_mass_fraction, which every
      !      evaluation of the advective term goes through, so the row the
      !      certification measures and the row the fixed-wind relaxation
      !      drives to zero are the same function of the interior state.
      !      The row below asserts that BITWISE: a column whose outer ghost
      !      is DATA, as the ionization sweep leaves it in the species
      !      vector the stationary evaluation reads, and the same column
      !      with the ghost the relaxation's line search writes beside every
      !      trial, give the identical outward face term of cell N.
      !
      !   3. The advective entries of the block-tridiagonal Jacobian are the
      !      FIRST-ORDER donor-cell linearization; the limited slopes of
      !      species_face_fraction are not differentiated.  On a profile
      !      that is constant or linear in the mass fraction the two agree
      !      exactly, which is what src/tests/carrier_boundary_jacobian
      !      measures.  On a CURVED profile they do not, and the rows below
      !      measure the difference and show that the banded central
      !      difference of the operator itself removes it.
      !
      ! EXHALE_TEST_CURVED=1 gives the carrier profile a quadratic part,
      ! which is what makes a ghost that is not a copy reach the outflow
      ! face and what separates the first-order entries from the deferred
      ! ones.
      use global_parameters
      use species_table
      use diffusive_photochemistry, only:                                 &
           carrier_set_init, carrier_advective_divergence,                &
           carrier_advective_face_coefficients,                           &
           carrier_face_coefficients, l22b_setup,                         &
           carrier_outflow_ghost, l22b_advective_band_derivative,         &
           carrier_mass_amu, n_carrier_max, n_carrier, ic_H2
      use assertion_report, only: check_relative, check_absolute,         &
                                  check_positive, check_at_least,         &
                                  assertion_failures
      implicit none

      real*8, allocatable :: fc(:,:), msum(:), Frho(:), wfac(:)
      real*8, allocatable :: ntot(:), TK(:), mbar(:), gphys(:), rp(:)
      real*8, allocatable :: Dco(:,:)
      real*8, allocatable :: Agrd(:,:), Bdrf(:,:)
      integer, allocatable :: updrf(:,:)
      real*8, allocatable :: advp(:,:), advmn(:,:)
      real*8, allocatable :: advj(:), advm(:), dadv(:,:,:)
      real*8  :: h, dnum, dana, dband, curv
      real*8  :: yfN, yN, yf_data, yf_copy
      character(len=32) :: curved
      integer :: ic, jtest

      ic = ic_H2
      call get_environment_variable('EXHALE_TEST_CURVED', curved)
      curv = 0.0d0
      if (trim(curved) .eq. '1') curv = 1.0d0

      call setup_globals()
      call l22b_setup()
      allocate(fc(1-Ng:N+Ng, n_carrier_max), msum(1-Ng:N+Ng),             &
               Frho(1-Ng:N+Ng), wfac(1-Ng:N+Ng), ntot(1-Ng:N+Ng),         &
               TK(1-Ng:N+Ng), mbar(1-Ng:N+Ng), gphys(1-Ng:N+Ng),          &
               rp(1-Ng:N+Ng), Dco(1-Ng:N+Ng, n_carrier_max),              &
               Agrd(0:N, n_carrier_max), Bdrf(0:N, n_carrier_max),        &
               updrf(0:N, n_carrier_max),                                 &
               advp(1:N, n_carrier_max), advmn(1:N, n_carrier_max),       &
               advj(1:N), advm(1:N),                                      &
               dadv(1:N, n_carrier_max, -1:1))

      call fill_state(fc, msum, Frho, wfac, ntot, TK, mbar, gphys, rp,    &
                      Dco, curv, ic)
      call carrier_face_coefficients(ntot, TK, mbar, gphys, Dco, rp,      &
                                     Agrd, Bdrf, updrf)

      ! ---- the two end faces of the column --------------------------- !
      call check_absolute('the_outer_diffusive_face_carries_no'//         &
           '_gradient_flux', Agrd(N,ic), 0.0d0, 0.0d0)
      call check_absolute('the_outer_diffusive_face_carries_no'//         &
           '_settling_drift', Bdrf(N,ic), 0.0d0, 0.0d0)
      call check_absolute('the_inner_diffusive_face_carries_no'//         &
           '_gradient_flux', Agrd(0,ic), 0.0d0, 0.0d0)
      call check_absolute('the_inner_diffusive_face_carries_no'//         &
           '_settling_drift', Bdrf(0,ic), 0.0d0, 0.0d0)
      ! and the face inside them does, so the zeros are the boundary and
      ! not an empty state.
      call check_positive('the_face_inside_it_carries_a_gradient'//       &
           '_flux', Agrd(N-1,ic))

      ! ---- one outer ghost rule, one outflow face -------------------- !
      ! The outward face term of cell N, formed once from a column whose
      ! outer ghost is DATA and once from the same column with the ghost
      ! the line search writes.  The two are the same number BITWISE, or
      ! the certification and the relaxation are measuring different rows
      ! of one state.
      if (curv .gt. 0.0d0) then
         call continue_outer_ghost(fc, ic)
         yf_data = face_fraction_outer(fc, msum, Frho, ic)
         call carrier_outflow_ghost(fc)
         yf_copy = face_fraction_outer(fc, msum, Frho, ic)
         call check_absolute('one_outer_ghost_rule_gives_one_outflow'//   &
              '_face', yf_data - yf_copy, 0.0d0, 0.0d0)
         ! and what that one face carries is the donor cell average, which
         ! is what the outflow continuation leaves the MC limiter of PLM.
         yN  = carrier_mass_amu(ic)*fc(N,ic)/msum(N)
         yfN = face_fraction_outer(fc, msum, Frho, ic)
         call check_absolute('the_outflow_face_is_the_donor_cell'//       &
              '_average', abs(yfN - yN)/yN, 0.0d0, 1.0d-13)
      endif

      ! ---- the deferred reconstruction terms ------------------------- !
      ! At an interior cell of a CURVED profile, the first-order donor-cell
      ! derivative the assembly writes against a central difference of the
      ! operator itself, and then the banded central difference against the
      ! same.  The first row is the size of what the Jacobian defers; the
      ! second says the banded derivative is the operator's own.
      jtest = N/2
      call carrier_advective_face_coefficients(Frho, advj, advm)
      h = 1.0d-7*fc(jtest,ic)
      call perturbed_divergence(fc, msum, Frho, ic, jtest,  h, advp)
      call perturbed_divergence(fc, msum, Frho, ic, jtest, -h, advmn)
      dnum = (advp(jtest,ic) - advmn(jtest,ic))/(2.0d0*h)
      dana = 0.0d0
      if (Frho(jtest) .ge. 0.0d0)                                         &
         dana = dana + advj(jtest)*msum(jtest)/msum(jtest)
      if (Frho(jtest-1) .lt. 0.0d0)                                       &
         dana = dana + advm(jtest)*msum(jtest)/msum(jtest)
      call l22b_advective_band_derivative(fc, msum, Frho, dadv)
      dband = dadv(jtest,ic,0)
      if (curv .le. 0.0d0) then
         call check_relative('the_first_order_entry_is_the'//             &
              '_derivative_on_a_constant_profile', dana, dnum, 1.0d-6)
      else
         call check_at_least('the_first_order_entry_misses_the'//         &
              '_reconstruction_on_a_curved_profile',                      &
              abs(dana - dnum)/max(abs(dnum), 1.0d-300), 1.0d-3)
      endif
      call check_relative('the_banded_derivative_is_the_operator'//       &
           's_own', dband, dnum, 1.0d-5)

      if (assertion_failures .gt. 0) stop 1

      contains

      subroutine setup_globals()
      ! A short column with the molecular chemistry and the carrier
      ! transport on, and no handoff partition, so that the base ghosts are
      ! a copy of cell 1.
      integer :: j
      N   = 12
      n0  = 1.0d10
      R0  = 1.0d10
      v0  = 1.0d6
      T0  = 1000.0d0
      HeH = 0.0793d0
      p_base_bar         = 1.0d-6
      thereis_He         = .true.
      thereis_HeITR      = .false.
      thereis_metals     = .false.
      thereis_mol        = .true.
      thereis_oxychem    = .false.
      eos_include_metals = .false.
      he_diffusion       = .false.
      carrier_transport  = .true.
      ionization_transport = .false.
      rec_method = 'PLM'
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         r(j)     = 1.0d0 + 0.01d0*dble(j)
         r_edg(j) = 1.0d0 + 0.01d0*(dble(j) + 0.5d0)
         dr_j(j)  = 0.01d0
      enddo
      allocate(kzz_cell(1-Ng:N+Ng))
      kzz_cell = 1.0d8
      allocate(melem_ab(n_melem))
      melem_ab = 0.0d0
      call carrier_set_init()
      end subroutine setup_globals

      subroutine fill_state(fc, msum, Frho, wfac, ntot, TK, mbar, gphys,  &
                            rp, Dco, curv, ic)
      ! A column whose carrier mass fraction is either exactly constant or
      ! curved.  The limited reconstruction returns the donor value exactly
      ! in the first case and does not in the second, which is what
      ! separates the first-order entries from the deferred ones and what
      ! makes the outflow face carry a slope.
      real*8, intent(out) :: fc(1-Ng:N+Ng, n_carrier_max)
      real*8, intent(out) :: msum(1-Ng:N+Ng), Frho(1-Ng:N+Ng)
      real*8, intent(out) :: wfac(1-Ng:N+Ng), ntot(1-Ng:N+Ng)
      real*8, intent(out) :: TK(1-Ng:N+Ng), mbar(1-Ng:N+Ng)
      real*8, intent(out) :: gphys(1-Ng:N+Ng), rp(1-Ng:N+Ng)
      real*8, intent(out) :: Dco(1-Ng:N+Ng, n_carrier_max)
      real*8,  intent(in) :: curv
      integer, intent(in) :: ic
      real*8  :: x, Yc
      integer :: j
      fc   = 0.0d0
      do j = 1-Ng, N+Ng
         x       = dble(j)/dble(N)
         msum(j) = 2.0d0 + 1.0d-2*x
         Yc      = 0.30d0*(1.0d0 + curv*(0.20d0*x + 0.80d0*x*x))
         fc(j,ic) = Yc*msum(j)/carrier_mass_amu(ic)
         Frho(j) = 1.0d-3
         ntot(j) = 1.0d10*exp(-2.0d0*x)
         TK(j)   = 500.0d0
         mbar(j) = 2.0d0*mu
         gphys(j)= 1.0d3
         rp(j)   = r(j)*R0
         wfac(j) = 1.0d0
         Dco(j,:) = 1.0d6
      enddo
      ! The outer ghosts follow cell N, the rule the residual uses.
      do j = N+1, N+Ng
         msum(j)  = msum(N)
         fc(j,ic) = fc(N,ic)
      enddo
      end subroutine fill_state

      double precision function face_fraction_outer(fc0, msum, Frho, ic)  &
                                                    result(Yf)
      ! The face mass fraction the OPERATOR put at r_{N+1/2}, recovered
      ! from the two faces of cell N that carrier_advective_divergence
      ! reports.
      real*8, intent(in)  :: fc0(1-Ng:N+Ng, n_carrier_max)
      real*8, intent(in)  :: msum(1-Ng:N+Ng), Frho(1-Ng:N+Ng)
      integer, intent(in) :: ic
      real*8, allocatable :: a1(:,:), m1(:,:), ain(:,:), aout(:,:)
      real*8  :: rpf, rmf, dV, tconv
      allocate(a1(1:N,n_carrier_max), m1(1:N,n_carrier_max),              &
               ain(1:N,n_carrier_max), aout(1:N,n_carrier_max))
      call carrier_advective_divergence(fc0, msum, Frho, a1, m1,          &
                                        advin = ain, advout = aout)
      tconv = n0*v0/R0
      rpf   = r_edg(N)
      rmf   = r_edg(N-1)
      dV    = (rpf*rpf*rpf - rmf*rmf*rmf)/3.0d0
      Yf = aout(N,ic)*dV*carrier_mass_amu(ic)                             &
           /(rpf*rpf*Frho(N)*msum(N)*tconv)
      deallocate(a1, m1, ain, aout)
      end function face_fraction_outer

      subroutine continue_outer_ghost(fc0, ic)
      ! An outer ghost that is DATA and not a copy: a state the ionization
      ! sweep could have left in the species vector, which is what the
      ! stationary evaluation used to read there.
      real*8, intent(inout) :: fc0(1-Ng:N+Ng, n_carrier_max)
      integer, intent(in)   :: ic
      integer :: j
      do j = N+1, N+Ng
         fc0(j,ic) = fc0(N,ic)                                            &
                   + dble(j - N)*(fc0(N,ic) - fc0(N-1,ic))
      enddo
      end subroutine continue_outer_ghost

      subroutine perturbed_divergence(fc0, msum, Frho, ic, k, h, adv)
      real*8, intent(in)  :: fc0(1-Ng:N+Ng, n_carrier_max)
      real*8, intent(in)  :: msum(1-Ng:N+Ng), Frho(1-Ng:N+Ng)
      integer, intent(in) :: ic, k
      real*8, intent(in)  :: h
      real*8, intent(out) :: adv(1:N, n_carrier_max)
      real*8, allocatable :: fw(:,:), advmag(:,:)
      allocate(fw(1-Ng:N+Ng, n_carrier_max), advmag(1:N, n_carrier_max))
      fw = fc0
      fw(k,ic) = fw(k,ic) + h
      call carrier_outflow_ghost(fw)
      call carrier_advective_divergence(fw, msum, Frho, adv, advmag)
      deallocate(fw, advmag)
      end subroutine perturbed_divergence

      end program carrier_outer_boundary
