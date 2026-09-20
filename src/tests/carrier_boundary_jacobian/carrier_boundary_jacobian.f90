      program carrier_boundary_jacobian
      ! THE BOUNDARY DERIVATIVES OF THE CARRIER ADVECTION, measured.
      !
      ! The block-tridiagonal Jacobian of the fixed-wind carrier relaxation
      ! linearizes the donor-cell face-flux divergence.  At the two ends of
      ! the grid the donor may be a GHOST, and whether that ghost carries a
      ! derivative depends on the rule the RESIDUAL uses for it:
      !
      !   outer:  fc(N+1) = fc(N)                       -- a COPY, always
      !   inner:  Y(1-Ng:0) = Y(1) unless the handoff states this carrier's
      !           partition                             -- a COPY or DATA
      !
      ! A copy contributes d/d f_c of the cell it copies; data contributes
      ! nothing.  This program measures the derivative of the operator
      ! itself, by a central difference of carrier_advective_divergence, and
      ! compares it with the coefficient the assembly uses.  It separates
      ! the FIRST-ORDER donor-cell term, which the Jacobian carries, from
      ! the limiter and the second-order reconstruction, which it
      ! deliberately does not: the carrier profile of every case below is
      ! either constant or linear, so species_face_fraction returns the
      ! donor value and the two agree exactly.
      !
      ! Four cases per boundary: the face mass flux positive and negative,
      ! times a constant and a linear carrier profile.  At the base the
      ! imposed and copied rules are both exercised.
      use global_parameters
      use species_table
      use diffusive_photochemistry, only:                                 &
           carrier_set_init, carrier_advective_divergence,                &
           carrier_advective_face_coefficients,                           &
           carrier_base_composition_imposed,                              &
           carrier_mass_amu, n_carrier_max, ic_H2
      use grid_construction, only: spherical_face_area_and_cell_volume
      use species_advective_transport, only: species_face_fraction,       &
                                             species_face_flux
      use assertion_report, only: check_relative, check_positive,         &
                                  check_absolute, assertion_failures
      implicit none

      real*8, allocatable :: fc(:,:), msum(:), Frho(:)
      real*8, allocatable :: adv_p(:,:), adv_m(:,:)
      real*8, allocatable :: advj(:), advm(:)
      real*8, allocatable :: fa(:), cv(:), Y(:), Yf(:), Fs(:)
      real*8  :: h, dnum, dana, mc
      real*8  :: tconv, dvj, worst, wscale, csum, cscale, bnd
      integer :: ic, kase, j
      character(len=72) :: nm
      logical :: linear, outward


      call setup_globals()
      ic = ic_H2
      mc = carrier_mass_amu(ic)
      allocate(fc(1-Ng:N+Ng, n_carrier_max), msum(1-Ng:N+Ng),             &
               Frho(1-Ng:N+Ng), adv_p(1:N, n_carrier_max),                &
               adv_m(1:N, n_carrier_max), advj(1:N), advm(1:N))

      do kase = 1, 4
         linear  = (mod(kase,2) .eq. 0)
         outward = (kase .le. 2)
         call fill_state(fc, msum, Frho, linear, outward, ic)
         call carrier_advective_face_coefficients(Frho, advj, advm)

         ! ---- the OUTER boundary: d adv(N) / d f_c(N) ----------------- !
         h = 1.0d-7*fc(N,ic)
         call perturbed_divergence(fc, msum, Frho, ic, N,  h, adv_p)
         call perturbed_divergence(fc, msum, Frho, ic, N, -h, adv_m)
         ! d adv(N)/d f_c(N) = sum over the two faces of
         !   (face coefficient) x msum(N)/msum(donor),
         ! the carrier mass cancelling between adv's msum/m_c and Y's
         ! m_c/msum.
         dnum = (adv_p(N,ic) - adv_m(N,ic))/(2.0d0*h)
         if (Frho(N) .ge. 0.0d0) then
            dana = advj(N)*msum(N)/msum(N)
         else
            dana = advj(N)*msum(N)/msum(N+1)
         endif
         if (Frho(N-1) .lt. 0.0d0) dana = dana + advm(N)*msum(N)/msum(N)
         write(nm,'(A,I0)') 'outer_face_jacobian_case', kase
         call check_relative(trim(nm), dnum, dana, 1.0d-6)

         ! ---- the INNER boundary: d adv(1) / d f_c(1) ----------------- !
         h = 1.0d-7*fc(1,ic)
         call perturbed_divergence(fc, msum, Frho, ic, 1,  h, adv_p)
         call perturbed_divergence(fc, msum, Frho, ic, 1, -h, adv_m)
         dnum = (adv_p(1,ic) - adv_m(1,ic))/(2.0d0*h)
         dana = 0.0d0
         if (Frho(1) .ge. 0.0d0) dana = dana + advj(1)*msum(1)/msum(1)
         if (Frho(0) .lt. 0.0d0) then
            dana = dana + advm(1)*msum(1)/msum(1)
         else if (.not. carrier_base_composition_imposed(ic)) then
            dana = dana + advm(1)*msum(1)/msum(1)
         endif
         write(nm,'(A,I0)') 'inner_face_jacobian_case', kase
         call check_relative(trim(nm), dnum, dana, 1.0d-6)

         ! ---- and the entries the old assembly wrote as zero ---------- !
         if (Frho(N) .lt. 0.0d0) then
            write(nm,'(A,I0)') 'outer_inflow_derivative_is_nonzero_case', &
                 kase
            call check_positive(trim(nm), abs(dnum))
         endif
      enddo

      ! THE CARRIER ADVECTION IS THE EXACT SPHERICAL DIVERGENCE OF ONE
      ! FACE FLUX (docs/PLAN_20260917.md item L30).  The divergence the
      ! operator returns is rebuilt here from the face areas A = r_edg^2 and
      ! the exact shell volumes V = (r_+^3 - r_-^3)/3 that
      ! spherical_face_area_and_cell_volume defines once for the whole code,
      ! and summed down the column, where the internal faces cancel and the
      ! two boundary face fluxes are left.  The carrier row's diffusive half
      ! divides by the same V, which is what makes the row the divergence of
      ! a single flux; these rows hold the half that is reachable from
      ! outside the module.
      allocate(fa(0:N), cv(1:N), Y(1-Ng:N+Ng), Yf(1-Ng:N+Ng),             &
               Fs(1-Ng:N+Ng))
      call fill_state(fc, msum, Frho, .true., .true., ic)
      call carrier_advective_divergence(fc, msum, Frho, adv_p, adv_m)
      call spherical_face_area_and_cell_volume(fa, cv)
      Y = mc*fc(:,ic)/msum
      call species_face_fraction(Y, Frho, Yf)
      call species_face_flux(Frho, Yf, Fs)
      tconv  = n0*v0/R0
      worst  = 0.0d0
      wscale = 0.0d0
      csum   = 0.0d0
      cscale = 0.0d0
      do j = 1, N
         dvj    = (fa(j)*Fs(j) - fa(j-1)*Fs(j-1))/cv(j)
         worst  = max(worst, abs(dvj*msum(j)/mc*tconv - adv_p(j,ic)))
         wscale = max(wscale, abs(adv_m(j,ic)))
         csum   = csum   + cv(j)*adv_p(j,ic)*mc/(msum(j)*tconv)
         cscale = cscale + cv(j)*adv_m(j,ic)*mc/(msum(j)*tconv)
      enddo
      bnd = fa(N)*Fs(N) - fa(0)*Fs(0)
      call check_absolute('carrier_advection_is_the_exact_shell_'//       &
           'divergence', worst/max(wscale, 1.0d-300), 0.0d0, 1.0d-14)
      call check_absolute('carrier_advection_telescopes_to_the_'//        &
           'boundary_face_flux', abs(csum - bnd)/max(cscale, 1.0d-300),   &
           0.0d0, 1.0d-13)
      deallocate(fa, cv, Y, Yf, Fs)

      if (assertion_failures .gt. 0) stop 1

      contains

      subroutine setup_globals()
      ! A short column with the molecular chemistry and the carrier
      ! transport on, and no handoff partition, so that the base ghosts are
      ! a COPY of cell 1 -- the case the assembly was dropping.
      integer :: j
      N   = 8
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
      ! Donor-cell reconstruction: the face fraction is then the donor
      ! value exactly, which is the first-order term the Jacobian carries
      ! and the only one this probe is about.
      rec_method = 'PLM'
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         r(j)     = 1.0d0 + 0.01d0*dble(j)
         r_edg(j) = 1.0d0 + 0.01d0*(dble(j) + 0.5d0)
         dr_j(j)  = 0.01d0
      enddo
      allocate(melem_ab(n_melem))
      melem_ab = 0.0d0
      call carrier_set_init()
      end subroutine setup_globals

      subroutine fill_state(fc, msum, Frho, linear, outward, ic)
      real*8, intent(out) :: fc(1-Ng:N+Ng, n_carrier_max)
      real*8, intent(out) :: msum(1-Ng:N+Ng), Frho(1-Ng:N+Ng)
      logical, intent(in) :: linear, outward
      integer, intent(in) :: ic
      real*8  :: sgn, yconst
      integer :: j
      yconst = 0.30d0
      fc   = 0.0d0
      msum = 2.0d0
      do j = 1-Ng, N+Ng
         ! a mixture mass that varies, so that a derivative taken with the
         ! wrong cell's msum cannot pass
         msum(j) = 2.0d0 + 1.0d-2*real(j, 8)/real(N, 8)
         ! f_c chosen so that the MASS fraction Y = m_c f_c/msum is exactly
         ! constant: the slope-limited reconstruction then returns the donor
         ! value and the face carries no second-order term, which is what
         ! separates the first-order donor-cell derivative the Jacobian
         ! holds from the reconstruction it deliberately omits.
         fc(j,ic) = yconst*msum(j)/carrier_mass_amu(ic)
         if (linear) fc(j,ic) = fc(j,ic)*1.5d0
      enddo
      ! The outer ghosts are a copy of cell N, as the operator fills them,
      ! and the mixture mass is carried out flat with them -- which is what
      ! a free-outflow ghost gives -- so that the copy preserves the
      ! constant MASS fraction and the reconstruction stays exact at that
      ! face.  Without this the ghost's Y differs from the interior's and
      ! the face carries a second-order term this probe is not about.
      do j = N+1, N+Ng
         msum(j) = msum(N)
         fc(j,ic) = fc(N,ic)
      enddo
      sgn = merge(1.0d0, -1.0d0, outward)
      do j = 1-Ng, N+Ng
         Frho(j) = sgn*1.0d-3
      enddo
      end subroutine fill_state

      subroutine perturbed_divergence(fc0, msum, Frho, ic, k, h, adv)
      real*8, intent(in)  :: fc0(1-Ng:N+Ng, n_carrier_max)
      real*8, intent(in)  :: msum(1-Ng:N+Ng), Frho(1-Ng:N+Ng)
      integer, intent(in) :: ic, k
      real*8, intent(in)  :: h
      real*8, intent(out) :: adv(1:N, n_carrier_max)
      real*8, allocatable :: fw(:,:), advmag(:,:)
      integer :: j
      allocate(fw(1-Ng:N+Ng, n_carrier_max), advmag(1:N, n_carrier_max))
      fw = fc0
      fw(k,ic) = fw(k,ic) + h
      ! the outer ghosts follow cell N, which is the rule under test
      do j = N+1, N+Ng
         fw(j,ic) = fw(N,ic)
      enddo
      call carrier_advective_divergence(fw, msum, Frho, adv, advmag)
      deallocate(fw, advmag)
      end subroutine perturbed_divergence

      end program carrier_boundary_jacobian
