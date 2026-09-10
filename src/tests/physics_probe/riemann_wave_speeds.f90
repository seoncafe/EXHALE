      program riemann_wave_speeds
      ! Signal speeds at a cell face, and the state one conservative update
      ! with the resulting flux leaves behind.
      !
      ! Origin: docs/audit_20260905/audit_probe.f90 (its LLF and ROE
      ! normalization sections), roe_equal_pressure.f90 and
      ! roe_flux_vacuum_probe.f90.  Those print numbers; this program
      ! compares them with the analytic reference and exits nonzero.
      !
      ! Every routine exercised below is the production one, compiled from
      ! src/modules: lax_friedrichs_flux and Num_flux of Num_Fluxes.f90 and
      ! speed_estimate_ROE of speed_estimate_ROE.f90.  Num_flux takes the
      ! ROE branch because global_parameters' flux is set to 'ROE', which
      ! is what the input key "Numerical flux:" writes.
      !
      ! References:
      !   Lax-Friedrichs/Rusanov signal speed: the flux is monotone and
      !   bound-preserving only if its viscosity coefficient is not smaller
      !   than the spectral radius of the flux Jacobian on either side.
      !   The eigenvalues are v-a, v and v+a, so that radius is |v| + a and
      !   the coefficient must be
      !       alpha = max(|v_L| + a_L, |v_R| + a_R).
      !   The positivity lemma Num_Fluxes.f90 cites for this flux (Perthame
      !   & Shu 1996; Zhang & Shu 2010) is proved under that bound.  With
      !   WL=(1,-2,0.6), WR=(2,-2,1.2) and gamma = 5/3 both sides have
      !   a = 1, so alpha = 3 and the mass flux is
      !   0.5*(-2-4-3*(2-1)) = -4.5.
      !   Mirror invariance: (rho,v,p) -> (rho,-v,p) with left and right
      !   exchanged is an exact symmetry of the Euler equations, so the mass
      !   flux must change sign.  It does so for any viscosity coefficient
      !   that is itself symmetric under the map; max(|v+a|) is not.
      !   Star state of a symmetric expansion WL=(1,-u,1), WR=(1,u,1):
      !   the exact two-rarefaction solution has
      !   c* = sqrt(gamma) - (gamma-1)u/2, which is 0.957661... at u = 1.
      !   Normalization: the Euler equations are invariant under
      !   (rho,p) -> (s rho, s p) at fixed velocity, so every sound speed
      !   of the star state is unchanged.
      !   Vacuum: at u = 10 the two rarefactions separate
      !   (u > 2 sqrt(gamma)/(gamma-1)) and the exact interface state is
      !   vacuum with zero flux; whatever the estimate returns, one Euler
      !   update of the left cell must not produce a negative density or a
      !   negative thermal energy.
      use global_parameters, only: gamma_ad, T0, use_plm, flux
      use Numerical_Fluxes, only: Num_flux, Phys_flux, lax_friedrichs_flux
      use Conversion, only: W_to_U_comp
      use caloric_eos, only: caloric_mixture_active
      use S_estimate_ROE, only: speed_estimate_ROE, roe_star_state,   &
                                ROE_STAR_OK, ROE_STAR_VACUUM
      use assertion_report
      implicit none

      real*8 :: wl(3), wr(3), ml(3), mr(3)
      real*8 :: nf(3), nf_mirror(3), fl(3), fr(3), pf
      real*8 :: u(3), unew(3)
      real*8 :: aL, aR, alpha_required, alpha_inferred
      real*8 :: c_star_exact, dt_dx, thermal
      real*8, parameter :: scale = 10.0d0
      type(roe_star_state) :: star, star_scaled
      integer :: roe_status, roe_status_scaled

      ! The atomic, single-gamma configuration: no molecular mixture, so
      ! both faces carry gamma_ad and the analytic references above apply.
      T0 = 1000.0d0
      use_plm = .true.
      caloric_mixture_active = .false.

      ! --- Lax-Friedrichs viscosity coefficient ----------------------- !
      wl = (/ 1.0d0, -2.0d0, 0.6d0 /)
      wr = (/ 2.0d0, -2.0d0, 1.2d0 /)
      aL = sqrt(gamma_ad*wl(3)/wl(1))
      aR = sqrt(gamma_ad*wr(3)/wr(1))
      alpha_required = max(abs(wl(2))+aL, abs(wr(2))+aR)

      call lax_friedrichs_flux(wl, wr, nf, pf, 1, 2)
      call Phys_flux(wl, fl, 1)
      call Phys_flux(wr, fr, 2)
      alpha_inferred = (fl(1) + fr(1) - 2.0d0*nf(1))/(wr(1)-wl(1))

      call check_relative('llf_signal_speed', alpha_inferred,             &
                          alpha_required, 1.0d-12)
      call check_relative('llf_mass_flux', nf(1), -4.5d0, 1.0d-12)

      ! Mirror the pair: (rho,v,p) -> (rho,-v,p) and exchange the sides.
      ml = (/ wr(1), -wr(2), wr(3) /)
      mr = (/ wl(1), -wl(2), wl(3) /)
      call lax_friedrichs_flux(ml, mr, nf_mirror, pf, 1, 2)
      call check_relative('llf_mirror_mass_flux', nf_mirror(1),           &
                          -nf(1), 1.0d-12)

      ! --- ROE star state: invariance under the density/pressure scale - !
      ! Two-rarefaction states (p_max/p_min = 10 selects that branch).
      wl = (/ 1.0d0, -1.0d0, 1.0d0 /)
      wr = (/ 1.0d0,  1.0d0, 0.1d0 /)
      call speed_estimate_ROE(wl, wr, gamma_ad, star, roe_status)
      wl(1) = wl(1)*scale; wl(3) = wl(3)*scale
      wr(1) = wr(1)*scale; wr(3) = wr(3)*scale
      call speed_estimate_ROE(wl, wr, gamma_ad, star_scaled,              &
                              roe_status_scaled)
      call check_relative('roe_rarefaction_scale_invariance_left',        &
                          star_scaled%c_L_star, star%c_L_star, 1.0d-12)
      call check_relative('roe_rarefaction_scale_invariance_right',       &
                          star_scaled%c_R_star, star%c_R_star, 1.0d-12)

      ! Two-shock states.
      wl = (/ 2.0d0,  5.0d0, 2.0d0 /)
      wr = (/ 8.0d0, -5.0d0, 0.2d0 /)
      call speed_estimate_ROE(wl, wr, gamma_ad, star, roe_status)
      wl(1) = wl(1)*scale; wl(3) = wl(3)*scale
      wr(1) = wr(1)*scale; wr(3) = wr(3)*scale
      call speed_estimate_ROE(wl, wr, gamma_ad, star_scaled,              &
                              roe_status_scaled)
      call check_relative('roe_shock_scale_invariance_left',              &
                          star_scaled%c_L_star, star%c_L_star, 1.0d-12)
      call check_relative('roe_shock_scale_invariance_right',             &
                          star_scaled%c_R_star, star%c_R_star, 1.0d-12)

      ! --- ROE star state of the equal-pressure expansion ------------- !
      ! An approximate estimator: compared at the accuracy of its
      ! approximation, 5 per cent of the exact two-rarefaction value.
      wl = (/ 1.0d0, -1.0d0, 1.0d0 /)
      wr = (/ 1.0d0,  1.0d0, 1.0d0 /)
      call speed_estimate_ROE(wl, wr, gamma_ad, star, roe_status)
      c_star_exact = sqrt(gamma_ad) - 0.5d0*(gamma_ad-1.0d0)
      call check_absolute('roe_expansion_status', dble(roe_status),       &
                          dble(ROE_STAR_OK), 0.0d0)
      call check_positive('roe_expansion_star_pressure', star%p_star)
      call check_positive('roe_expansion_star_sound_speed_positive_left', &
                          star%c_L_star)
      call check_positive('roe_expansion_star_sound_speed_positive_right',&
                          star%c_R_star)
      call check_relative('roe_expansion_star_sound_speed_left',          &
                          star%c_L_star, c_star_exact, 5.0d-2)
      call check_relative('roe_expansion_star_sound_speed_right',         &
                          star%c_R_star, c_star_exact, 5.0d-2)

      ! --- Stationary contact ----------------------------------------- !
      ! Equal pressures and zero velocity on both sides: the exact solution
      ! is the data itself, a contact at rest, with p_* = p and u_* = 0 and
      ! each star density equal to its own data density.  The estimator has
      ! to reproduce it exactly, not approximately: the linearized (PVRS)
      ! branch is exact on a state with no velocity jump and no pressure
      ! jump, so the tolerance is round-off and not the estimator's 5 per
      ! cent.
      wl = (/ 1.0d0,   0.0d0, 1.0d0 /)
      wr = (/ 0.125d0, 0.0d0, 1.0d0 /)
      call speed_estimate_ROE(wl, wr, gamma_ad, star, roe_status)
      call check_absolute('roe_contact_status', dble(roe_status),         &
                          dble(ROE_STAR_OK), 0.0d0)
      call check_relative('roe_contact_star_pressure', star%p_star,       &
                          1.0d0, 1.0d-12)
      call check_absolute('roe_contact_star_velocity', star%u_star,       &
                          0.0d0, 1.0d-12)
      call check_relative('roe_contact_star_sound_speed_left',            &
                          star%c_L_star,                                  &
                          sqrt(gamma_ad*wl(3)/wl(1)), 1.0d-12)
      call check_relative('roe_contact_star_sound_speed_right',           &
                          star%c_R_star,                                  &
                          sqrt(gamma_ad*wr(3)/wr(1)), 1.0d-12)

      ! --- Asymmetric two-shock collision ----------------------------- !
      ! WL=(1,0.5,1), WR=(3,-0.5,1.5), gamma = 5/3.  Both nonlinear waves
      ! are shocks.  The exact star state, from the root of Toro (2009)
      ! eq. 4.5 with the shock relations 4.20 and the density 4.50, is
      ! p_* = 2.3768663335, u_* = -0.2357045218,
      ! c_*L = 1.5505344803, c_*R = 1.0024512384.
      ! Compared at 5 per cent, the accuracy of a one-step estimator.
      wl = (/ 1.0d0,  0.5d0, 1.0d0 /)
      wr = (/ 3.0d0, -0.5d0, 1.5d0 /)
      call speed_estimate_ROE(wl, wr, gamma_ad, star, roe_status)
      call check_absolute('roe_shock_status', dble(roe_status),           &
                          dble(ROE_STAR_OK), 0.0d0)
      call check_relative('roe_shock_star_pressure', star%p_star,         &
                          2.3768663335d0, 5.0d-2)
      call check_relative('roe_shock_star_velocity', star%u_star,         &
                          -0.2357045218d0, 5.0d-2)
      call check_relative('roe_shock_star_sound_speed_left',              &
                          star%c_L_star, 1.5505344803d0, 5.0d-2)
      call check_relative('roe_shock_star_sound_speed_right',             &
                          star%c_R_star, 1.0024512384d0, 5.0d-2)

      ! --- Asymmetric two-rarefaction expansion ----------------------- !
      ! WL=(1,-2,1), WR=(0.5,1,0.2), gamma = 5/3: both nonlinear waves are
      ! rarefactions and the two gas edges still meet, so a star state
      ! exists.  Exact values from the same root, with the isentropic
      ! relations 4.53:
      ! p_* = 0.0201759337, u_* = 0.0987433798,
      ! c_*L = 0.5914133221, c_*R = 0.5160777075.
      wl = (/ 1.0d0, -2.0d0, 1.0d0 /)
      wr = (/ 0.5d0,  1.0d0, 0.2d0 /)
      call speed_estimate_ROE(wl, wr, gamma_ad, star, roe_status)
      call check_absolute('roe_rarefaction_status', dble(roe_status),     &
                          dble(ROE_STAR_OK), 0.0d0)
      call check_relative('roe_rarefaction_star_pressure', star%p_star,   &
                          0.0201759337d0, 5.0d-2)
      call check_relative('roe_rarefaction_star_velocity', star%u_star,   &
                          0.0987433798d0, 5.0d-2)
      call check_relative('roe_rarefaction_star_sound_speed_left',        &
                          star%c_L_star, 0.5914133221d0, 5.0d-2)
      call check_relative('roe_rarefaction_star_sound_speed_right',       &
                          star%c_R_star, 0.5160777075d0, 5.0d-2)

      ! --- Vacuum-producing expansion: the status and the gas edges --- !
      ! WL=(1,-10,1), WR=(1,10,1): the two rarefaction fans separate, no
      ! star state exists, and the gas edges are u -/+ 2c/(gamma-1) with
      ! c = sqrt(gamma), that is -6.127017 and +6.127017 (Toro 2009,
      ! eq. 4.76).
      wl = (/ 1.0d0, -10.0d0, 1.0d0 /)
      wr = (/ 1.0d0,  10.0d0, 1.0d0 /)
      call speed_estimate_ROE(wl, wr, gamma_ad, star, roe_status)
      call check_absolute('roe_vacuum_status', dble(roe_status),          &
                          dble(ROE_STAR_VACUUM), 0.0d0)
      call check_relative('roe_vacuum_edge_speed_left', star%v_edge_L,    &
                          -6.127017d0, 1.0d-6)
      call check_relative('roe_vacuum_edge_speed_right', star%v_edge_R,   &
                          6.127017d0, 1.0d-6)

      ! --- One conservative update with the production ROE flux ------- !
      ! The vacuum-producing expansion of roe_flux_vacuum_probe.f90.
      flux = 'ROE'
      wl = (/ 1.0d0, -10.0d0, 1.0d0 /)
      wr = (/ 1.0d0,  10.0d0, 1.0d0 /)
      call Num_flux(wl, wr, nf, pf, 1, 2)
      call W_to_U_comp(wl, u, 1)
      call Phys_flux(wl, fl, 1)
      dt_dx = 0.01d0
      unew = u - dt_dx*(nf - fl)
      thermal = unew(3) - 0.5d0*unew(2)**2/unew(1)
      call check_positive('roe_vacuum_update_density', unew(1))
      call check_positive('roe_vacuum_update_thermal_energy', thermal)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'riemann_wave_speeds: ', assertion_failures, &
              ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'riemann_wave_speeds: all assertions passed'

      end program riemann_wave_speeds
