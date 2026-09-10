      program xuv_cell_mean_attenuation
      ! QUANTITY UNDER TEST
      !   the attenuation of the stellar XUV beam that the photoionization
      !   rate of a CELL is evaluated with, in PH_heat_H (util_ion_eq.f90).
      !
      !   calc_column_dens accumulates from the top down, so N1(j) holds the
      !   WHOLE of cell j and the depth built from it is the depth at that
      !   cell's INNER (planet-ward) face.  The rate that belongs to the cell
      !   is the mean over the cell of the field, which for a uniform
      !   absorber density -- the rectangle rule the column integration
      !   itself uses -- is
      !
      !     <exp(-tau)> = exp(-tau_out) (1 - exp(-dtau))/dtau ,
      !
      !   with tau_out the depth at the cell's star-ward face and dtau the
      !   cell's own optical depth.  Two statements are asserted:
      !
      !     (1) RATE.  P_HI(j) equals the analytic cell mean, cell by cell,
      !         at dtau = 0.01, 0.3, 1 and 3.  The reference is an
      !         INDEPENDENT quadrature of the same integral (composite
      !         three-point Gauss-Legendre in the depth across the cell,
      !         128 segments, truncation below 1e-16 relative), not the
      !         closed form the code uses.
      !
      !     (2) PHOTON NUMBER.  With one absorber, every photon the beam
      !         loses is one H I photoionization, so the ionizations of the
      !         whole column,
      !           sum_j P_HI(j) n_HI dr ,
      !         must equal the photons the beam loses between its two ends,
      !           Phi_in [1 - exp(-tau_total)] .
      !         The cell mean satisfies this identity exactly, because
      !         <exp(-tau)> dtau = exp(-tau_out) - exp(-tau_in) telescopes
      !         over the column.  This is the physical statement: the
      !         photons the cells absorb are the photons the beam loses.
      !         Tolerance 1e-12: the columns are accumulated cell by cell,
      !         so a last-digit perturbation of the accumulated depth is
      !         amplified by tau_total, which reaches 36 in the dtau = 3
      !         case, and twelve such terms are summed.
      !
      !   The INNER-FACE rule satisfies neither.  Its departure is reported
      !   at every dtau and asserted against its own closed form,
      !     face/mean = dtau exp(-dtau)/(1 - exp(-dtau)) ,
      !   which is 1 - dtau/2 + O(dtau^2) and 0.582 at dtau = 1, i.e. a rate
      !   low by 41.8 per cent there.  That assertion is on the reference
      !   arithmetic alone and keeps this test discriminating: it fails if
      !   the two forms ever stop differing.
      !
      ! CONFIGURATION.  A twelve-cell uniform slab (N = 8 plus two ghost
      ! cells at each end) of H I alone, thereis_He false, no metals, no
      ! secondary ionization, opa_pf = 1, a_tau = 0 so that the field is the
      ! pure exponential.  The spectrum is MONOCHROMATIC at 20 eV: with one
      ! photon energy the cell's optical depth is one number rather than a
      ! spectrum of them, so "dtau = 1" means dtau = 1 and the 41.8 per cent
      ! above is the error of the rate itself.  n_HI is set from the target
      ! dtau; nothing else changes between the four cases.

      use global_parameters
      use energy_vectors_construct, only: set_energy_vectors
      use utils_ion_eq,           only: PH_heat_H
      use assertion_report

      implicit none

      ! Slab geometry.  dr_geo is the geometric width of every cell [cm].
      real*8, parameter :: dr_geo  = 1.0d4
      integer, parameter :: n_case = 4
      real*8, parameter :: dtau_case(n_case) = (/ 1.0d-2, 3.0d-1,          &
                                                  1.0d0,  3.0d0 /)
      ! Segments of the composite Gauss-Legendre reference.
      integer, parameter :: nseg_ref = 128
      ! The cross sections are carried in units of 1e-18 cm^2 and the
      ! production integrands restore the unit with the double-precision
      ! literal 1.0d-18 (util_ion_eq.f90); the reference here uses the same
      ! value so the two agree to round-off.
      real*8, parameter :: sigma_unit = 1.0d-18

      real*8, dimension(:), allocatable :: nhi,xion
      real*8, dimension(:), allocatable :: P_HI,heat,q
      real*8  :: n_HI, dtau, tau_out, dtau_bin
      real*8  :: p_mean, p_face, phi_in, ions, beam_loss, face_sum
      integer :: ic, j
      character(len=16) :: tag

      ! --- the run-wide input, as input_read resolves it ------------------
      is_PL_sed    = .false.
      do_read_sed  = .false.
      is_monochr   = .true.
      e_low        = 20.0d0
      LEUV         = 30.42d0
      LX           = -3.0d2
      a_orb        = 0.02544d0*AU
      appx_mth     = 'Rate/2 + Mdot/2'
      ! a_tau = 0 makes the field the pure exponential whose cell mean has
      ! the closed form asserted here.
      a_tau        = 0.0d0
      thereis_He   = .false.
      thereis_HeITR       = .false.
      thereis_Xray        = .false.
      thereis_lowIP_metal = .false.
      use_sec_ion    = .false.
      sec_ion_active = .false.

      call set_energy_vectors

      ! --- the grid -------------------------------------------------------
      N  = 8
      R0 = 1.0d10
      call allocate_grid_arrays
      dr_j   = dr_geo/R0
      opa_pf = 1.0d0

      allocate(nhi(1-Ng:N+Ng), xion(1-Ng:N+Ng))
      allocate(P_HI(1-Ng:N+Ng), heat(1-Ng:N+Ng), q(1-Ng:N+Ng))
      xion = 0.0d0

      ! Photons entering the slab per unit area and time [cm^-2 s^-1].
      phi_in = F_XUV(1)*erg2eV/e_v(1)*de_v(1)

      write(*,'(a,es12.5,a,es12.5)')                                       &
           '  DIAGNOSTIC monochromatic e_v [eV] = ', e_v(1),               &
           '   s_hi [1e-18 cm^2] = ', s_hi(1)

      do ic = 1,n_case

         dtau = dtau_case(ic)
         ! n_HI from the target optical depth of one cell.  The cross
         ! sections carry the code's internal unit of 1e-18 cm^2.
         n_HI = dtau/(s_hi(1)*dr_geo*sigma_unit)
         nhi  = n_HI

         call PH_heat_H(nhi, xion, P_HI, heat, q)

         write(tag,'(a,f0.2)') 'dtau=', dtau

         ! (1) the rate of every cell against the independent quadrature of
         !     its cell mean.  tau_out is the depth of everything OUTSIDE
         !     the cell; nothing sits above the outermost cell.
         do j = 1-Ng,N+Ng
            tau_out = dble(N+Ng-j)*dtau
            p_mean  = rate_from_attenuation(mean_attenuation(tau_out,dtau))
            call check_relative('rate_is_cell_mean['//trim(tag)//']',      &
                                P_HI(j), p_mean, 1.0d-12)
         enddo

         ! (2) photon number over the whole column.
         ions      = 0.0d0
         do j = 1-Ng,N+Ng
            ions = ions + P_HI(j)*n_HI*dr_geo
         enddo
         beam_loss = phi_in*(1.0d0 - exp(-dble(N+2*Ng)*dtau))
         call check_relative('column_ionizations_equal_beam_loss['         &
                             //trim(tag)//']', ions, beam_loss, 1.0d-12)

         ! The inner-face rule, on the reference arithmetic alone: its
         ! ratio of one cell's rate to the cell mean, and the column sum
         ! it would give.  Reported at every dtau; the ratio is asserted so that the
         ! test fails if the two discretizations ever coincide.
         tau_out  = 0.0d0
         p_face   = rate_from_attenuation(exp(-dtau))
         p_mean   = rate_from_attenuation(mean_attenuation(0.0d0,dtau))
         face_sum = 0.0d0
         do j = 1-Ng,N+Ng
            dtau_bin = dble(N+Ng-j)*dtau
            face_sum = face_sum                                            &
                     + rate_from_attenuation(exp(-(dtau_bin+dtau)))        &
                       *n_HI*dr_geo
         enddo
         write(*,'(a,f6.2,a,f9.6,a,f9.6)')                                 &
              '  DIAGNOSTIC inner-face rule at dtau = ', dtau,             &
              ' : rate ratio to the cell mean = ', p_face/p_mean,          &
              ' ; column photon budget ratio = ', face_sum/beam_loss
         call check_relative('inner_face_rule_ratio['//trim(tag)//']',     &
                             p_face/p_mean,                                &
                             dtau*exp(-dtau)/(1.0d0-exp(-dtau)), 1.0d-12)

      enddo

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'xuv_cell_mean_attenuation: ',                &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'xuv_cell_mean_attenuation: all assertions passed'

      contains

      !--------------!

      double precision function rate_from_attenuation(tr) result(p)
      ! The H I photoionization rate [s^-1] of a cell whose beam
      ! attenuation is tr, on the one-bin grid: the same frequency integral
      ! PH_heat_H forms, with the attenuation supplied instead of built.
      real*8, intent(in) :: tr
      p = F_XUV(1)*tr*s_hi(1)/e_v(1)*de_v(1)*sigma_unit*erg2eV
      end function rate_from_attenuation

      !--------------!

      double precision function mean_attenuation(tau_lo, d) result(tr)
      ! int_0^1 exp(-(tau_lo + s d)) ds by composite three-point
      ! Gauss-Legendre on nseg_ref equal segments.  An INDEPENDENT
      ! quadrature of the integral the code evaluates in closed form: the
      ! rule is exact to degree five, so the truncation of one segment is
      ! below (d/nseg_ref)^6/2016000 relative, i.e. under 1e-16 for every
      ! d used here.
      real*8, intent(in) :: tau_lo, d
      real*8, parameter :: gl_s(3) = (/ 0.1127016653792583d0,              &
                                        0.5d0,                             &
                                        0.8872983346207417d0 /)
      real*8, parameter :: gl_w(3) = (/ 5.0d0/18.0d0, 8.0d0/18.0d0,        &
                                        5.0d0/18.0d0 /)
      real*8  :: ds, s_lo, acc
      integer :: i, g
      ds  = 1.0d0/dble(nseg_ref)
      acc = 0.0d0
      do i = 0,nseg_ref-1
         s_lo = dble(i)*ds
         do g = 1,3
            acc = acc + gl_w(g)*ds*exp(-(tau_lo + (s_lo + gl_s(g)*ds)*d))
         enddo
      enddo
      tr = acc
      end function mean_attenuation

      end program xuv_cell_mean_attenuation
