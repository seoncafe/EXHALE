      program roe_low_mach_velocity_jump
      ! The low-Mach scaling of the normal velocity jump on the ROE branch:
      ! Rieper (2011, J. Comput. Phys. 230, 5263) eq. 3.15 with the local
      ! Mach number 3.16 of the Roe averages,
      !     dvel -> min(Ma_Roe, 1) dvel,   Ma_Roe = (|U| + |V|)/a,
      ! which in one dimension is |U_Roe|/a_Roe.  The key is
      ! "Low Mach velocity jump:", the module variable
      ! Numerical_Fluxes::low_mach_velocity_jump, default .false.
      !
      ! Everything exercised here is the production routine Num_flux of
      ! src/modules/flux/Num_Fluxes.f90, taking its ROE branch because
      ! global_parameters' flux is set to 'ROE'.
      !
      ! The reference of every assertion is the analytic Roe flux of the
      ! three-wave linearization (Toro 2009, eq. 11.29), formed here from the
      ! Roe averages of the pair:
      !     NF(f) = (F_L + F_R)/2
      !           - (a1(f)|l1|K1 + a2|l2|K2 + a3(f)|l3|K3)/2
      !     a1(f) = (dp - rho_avg a_avg f dvel)/(2 a_avg^2)
      !     a2    = drho - dp/a_avg^2
      !     a3(f) = (dp + rho_avg a_avg f dvel)/(2 a_avg^2)
      ! so that NF is AFFINE in f and only a1 and a3 carry the velocity
      ! jump.  The state pairs below are chosen so that neither eigenvalue
      ! is transonic and the rarefaction entropy correction of the branch
      ! does not fire; that the analytic reference reproduces the returned
      ! flux to rounding is itself the check of that choice.
      !
      ! The four properties measured:
      !   (1) with the key off the returned flux is the uncorrected Roe
      !       flux, and a call with the key on leaves nothing behind that
      !       moves a later call with the key off by a single bit;
      !   (2) at a supersonic face the factor is exactly 1 and the flux is
      !       bitwise the uncorrected one;
      !   (3) at Ma_Roe = 1e-3 the velocity-jump part of the dissipation is
      !       scaled by that Mach number while the pressure-jump and
      !       entropy parts are untouched;
      !   (4) where the two velocities cancel in the Roe average the factor
      !       is 0: the velocity-jump dissipation is removed entirely, so
      !       the 2 dr alternating velocity mode is left undamped.  The
      !       mass row of that pair is zero with and without the factor;
      !       the damping that disappears is the momentum row's
      !       rho_avg a_avg dvel / 2.
      use global_parameters, only: gamma_ad, T0, use_plm, flux,          &
                                   well_balanced
      use Numerical_Fluxes, only: Num_flux, Phys_flux,                   &
                                  low_mach_velocity_jump
      use caloric_eos, only: caloric_mixture_active
      use assertion_report
      implicit none

      real*8 :: wl(3), wr(3)
      real*8 :: nf_off(3), nf_on(3), nf_off_again(3), pf
      real*8 :: ref_off(3), ref_on(3), ref_nojump(3), implied(3)
      real*8 :: mach, rho_avg, a_avg, v_avg, dvel
      real*8 :: central(3), fl(3), fr(3)
      real*8 :: target_mach, u_alt, rho_alt, p_alt
      integer :: k
      character(len=1) :: row(3) = (/ 'm', 'p', 'e' /)

      ! The atomic, single-gamma configuration: one gamma on both sides of
      ! every face, which is what the analytic reference above assumes.
      ! use_plm with the well-balanced option off puts the pressure inside
      ! the momentum flux, the form Phys_flux then returns.
      T0 = 1000.0d0
      use_plm = .true.
      well_balanced = .false.
      caloric_mixture_active = .false.
      flux = 'ROE'

      ! --- (1) the key off is the uncorrected Roe flux ---------------- !
      ! A subsonic pair with jumps in all three primitive variables, so
      ! a1, a2 and a3 are all nonzero.
      wl = (/ 1.0d0, 0.30d0, 1.0d0 /)
      wr = (/ 1.2d0, 0.34d0, 1.1d0 /)
      low_mach_velocity_jump = .false.
      call Num_flux(wl, wr, nf_off, pf, 1, 2)
      call roe_reference(wl, wr, 1.0d0, ref_off)
      do k = 1, 3
         call check_relative('roe_key_off_row_'//row(k), nf_off(k),      &
                             ref_off(k), 1.0d-12)
      enddo

      ! A corrected call between two uncorrected ones: the branch keeps no
      ! state, so the second uncorrected flux is bitwise the first.
      low_mach_velocity_jump = .true.
      call Num_flux(wl, wr, nf_on, pf, 1, 2)
      low_mach_velocity_jump = .false.
      call Num_flux(wl, wr, nf_off_again, pf, 1, 2)
      do k = 1, 3
         call check_absolute('roe_key_off_repeatable_row_'//row(k),      &
                             nf_off_again(k) - nf_off(k), 0.0d0, 0.0d0)
      enddo

      ! --- (2) a supersonic pair: the factor is exactly 1 -------------- !
      ! a_avg is near sqrt(gamma p / rho) ~ 1.3 here, so Ma_Roe ~ 4.
      wl = (/ 1.0d0, 5.0d0, 1.0d0 /)
      wr = (/ 1.2d0, 5.2d0, 1.1d0 /)
      call roe_average(wl, wr, rho_avg, v_avg, a_avg)
      mach = abs(v_avg)/a_avg
      call check_at_least('roe_supersonic_face_mach', mach, 1.0d0)
      low_mach_velocity_jump = .false.
      call Num_flux(wl, wr, nf_off, pf, 1, 2)
      low_mach_velocity_jump = .true.
      call Num_flux(wl, wr, nf_on, pf, 1, 2)
      do k = 1, 3
         call check_absolute('roe_supersonic_unchanged_row_'//row(k),    &
                             nf_on(k) - nf_off(k), 0.0d0, 0.0d0)
      enddo

      ! --- (3) Ma_Roe = 1e-3: the velocity jump scaled, the pressure --- !
      !         jump untouched
      ! Equal densities put the Roe weight s at 1/2 exactly, so v_avg is
      ! the arithmetic mean of the two velocities and the face Mach number
      ! is placed where it is wanted; the pressure jump is what makes a2
      ! and the pressure part of a1 and a3 nonzero, and the velocity jump
      ! is 4e-3, four times the mean.
      target_mach = 1.0d-3
      wl = (/ 1.0d0, 0.0d0, 1.0d0 /)
      wr = (/ 1.0d0, 0.0d0, 1.3d0 /)
      call roe_average(wl, wr, rho_avg, v_avg, a_avg)
      ! a_avg of the pair at rest; the velocities below move it only in
      ! the 7th digit, and the Mach number compared below is the one built
      ! from the Roe average the flux itself forms.
      wl(2) = target_mach*a_avg - 2.0d-3
      wr(2) = target_mach*a_avg + 2.0d-3
      call roe_average(wl, wr, rho_avg, v_avg, a_avg)
      mach = abs(v_avg)/a_avg
      dvel = wr(2) - wl(2)
      call check_relative('roe_low_mach_face_mach', mach, target_mach,   &
                          5.0d-3)

      low_mach_velocity_jump = .false.
      call Num_flux(wl, wr, nf_off, pf, 1, 2)
      low_mach_velocity_jump = .true.
      call Num_flux(wl, wr, nf_on, pf, 1, 2)
      call roe_reference(wl, wr, mach, ref_on)
      do k = 1, 3
         call check_relative('roe_low_mach_scaled_row_'//row(k),         &
                             nf_on(k), ref_on(k), 1.0d-12)
      enddo

      ! The flux is affine in the factor, so the two calls determine the
      ! flux at factor 0, the one with no velocity-jump dissipation at all.
      ! That it equals the analytic reference at factor 0 says the
      ! pressure-jump and entropy parts of the dissipation did not move.
      call roe_reference(wl, wr, 0.0d0, ref_nojump)
      do k = 1, 3
         implied(k) = (nf_on(k) - mach*nf_off(k))/(1.0d0 - mach)
         call check_relative('roe_pressure_jump_part_unchanged_row_'//   &
                             row(k), implied(k), ref_nojump(k), 1.0d-9)
      enddo

      ! --- (4) cancelling velocities: the factor is 0 ------------------ !
      ! v_L = -v_R at equal densities gives s = 1/2 and v_avg = 0 exactly,
      ! the 2 dr alternating mode of a base at rest.  With equal densities
      ! and pressures drho and dp vanish too, so the corrected flux is the
      ! central flux alone.
      rho_alt = 1.0d0
      p_alt   = 1.0d0
      u_alt   = 1.0d-4
      wl = (/ rho_alt, -u_alt, p_alt /)
      wr = (/ rho_alt,  u_alt, p_alt /)
      call roe_average(wl, wr, rho_avg, v_avg, a_avg)
      call check_absolute('roe_alternating_face_velocity', v_avg,        &
                          0.0d0, 0.0d0)
      low_mach_velocity_jump = .false.
      call Num_flux(wl, wr, nf_off, pf, 1, 2)
      low_mach_velocity_jump = .true.
      call Num_flux(wl, wr, nf_on, pf, 1, 2)

      call Phys_flux(wl, fl, 1)
      call Phys_flux(wr, fr, 2)
      central = 0.5d0*(fr + fl)
      do k = 1, 3
         call check_absolute('roe_alternating_central_flux_row_'//row(k),&
                             nf_on(k) - central(k), 0.0d0, 1.0d-15)
      enddo

      ! The mass row of this pair is zero with and without the factor:
      ! the central mass flux is rho(-u + u)/2 = 0 and the acoustic
      ! coefficients a1 and a3 are equal and opposite there, so their
      ! contributions to the mass row cancel.  The dissipation the factor
      ! removes is the momentum row's rho_avg a_avg dvel / 2.
      call check_absolute('roe_alternating_mass_flux_uncorrected',       &
                          nf_off(1), 0.0d0, 1.0d-18)
      call check_absolute('roe_alternating_mass_flux_corrected',         &
                          nf_on(1), 0.0d0, 1.0d-18)
      dvel = wr(2) - wl(2)
      call check_relative('roe_alternating_momentum_damping_removed',    &
                          nf_off(2) - nf_on(2),                          &
                          -0.5d0*rho_avg*a_avg*dvel, 1.0d-10)

      if (assertion_failures .gt. 0) then
         write(*,*) 'roe_low_mach_velocity_jump: FAILED'
         error stop 1
      endif
      write(*,*) 'roe_low_mach_velocity_jump: PASSED'

      contains

      !--------------!

      subroutine roe_average(WL, WR, rho_avg, v_avg, a_avg, H_out)
      ! The Roe averages of the pair for an ideal gas of one gamma:
      ! the geometric mean density, the square-root-density weighted
      ! velocity and total enthalpy, and the sound speed built from them
      ! (Roe 1981, J. Comput. Phys. 43, 357).
      real*8, intent(in)  :: WL(3), WR(3)
      real*8, intent(out) :: rho_avg, v_avg, a_avg
      real*8, intent(out), optional :: H_out
      real*8 :: s, HL, HR, H_avg, EL, ER
      s = sqrt(WL(1))/(sqrt(WL(1)) + sqrt(WR(1)))
      rho_avg = sqrt(WL(1)*WR(1))
      EL = 0.5d0*WL(1)*WL(2)*WL(2) + WL(3)/(gamma_ad - 1.0d0)
      ER = 0.5d0*WR(1)*WR(2)*WR(2) + WR(3)/(gamma_ad - 1.0d0)
      HL = (EL + WL(3))/WL(1)
      HR = (ER + WR(3))/WR(1)
      v_avg = s*WL(2) + (1.0d0 - s)*WR(2)
      H_avg = s*HL + (1.0d0 - s)*HR
      a_avg = sqrt((gamma_ad - 1.0d0)*(H_avg - 0.5d0*v_avg*v_avg))
      if (present(H_out)) H_out = H_avg
      end subroutine roe_average

      !--------------!

      subroutine roe_reference(WL, WR, factor, NF)
      ! The three-wave Roe flux with the normal velocity jump multiplied
      ! by `factor`, without the transonic rarefaction correction (the
      ! pairs used here have none).
      real*8, intent(in)  :: WL(3), WR(3), factor
      real*8, intent(out) :: NF(3)
      real*8 :: rho_avg, v_avg, a_avg
      real*8 :: drho, dvel, dp, a1, a2, a3
      real*8 :: l1, l2, l3, K1(3), K2(3), K3(3), H_avg
      real*8 :: FL(3), FR(3)
      call roe_average(WL, WR, rho_avg, v_avg, a_avg, H_avg)
      drho = WR(1) - WL(1)
      dvel = factor*(WR(2) - WL(2))
      dp   = WR(3) - WL(3)
      a1 = 0.5d0/a_avg**2*(dp - rho_avg*a_avg*dvel)
      a2 = drho - dp/a_avg**2
      a3 = 0.5d0/a_avg**2*(dp + rho_avg*a_avg*dvel)
      l1 = v_avg - a_avg
      l2 = v_avg
      l3 = v_avg + a_avg
      K1 = (/ 1.0d0, v_avg - a_avg, H_avg - v_avg*a_avg /)
      K2 = (/ 1.0d0, v_avg, 0.5d0*v_avg*v_avg /)
      K3 = (/ 1.0d0, v_avg + a_avg, H_avg + v_avg*a_avg /)
      call Phys_flux(WL, FL, 1)
      call Phys_flux(WR, FR, 2)
      NF = 0.5d0*(FR + FL)                                              &
         - 0.5d0*(a1*abs(l1)*K1 + a2*abs(l2)*K2 + a3*abs(l3)*K3)
      end subroutine roe_reference

      end program roe_low_mach_velocity_jump
