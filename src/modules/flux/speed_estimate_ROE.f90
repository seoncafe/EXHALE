      module S_estimate_ROE
      ! Star state of the Riemann problem at a cell face, as the
      ! approximate-state solvers of Toro (2009) chapter 9 build it: the
      ! primitive-variable estimate (eq. 9.20), the two-rarefaction estimate
      ! (eq. 9.32) and the two-shock estimate (eq. 9.42), selected by the
      ! pressure ratio of the two sides, with the vacuum condition of
      ! eq. 4.76 tested before any of them.
      !
      ! VALIDITY: one ideal gas with a CONSTANT adiabatic index on both
      ! sides.  gam is that index; the caller passes the mean of the two
      ! face gammas, which is that value exactly wherever the two sides
      ! agree (every face of an atomic gas).  A run whose caloric EOS makes
      ! the index vary across the front is refused the Roe flux at input
      ! (input_read.f90; docs/a2_roe_interface.md section 3).
      !
      ! The estimate carries a status because the star state need not exist.
      ! Two rarefactions that separate fast enough leave vacuum between the
      ! gas edges (ROE_STAR_VACUUM), and no branch is entitled to a star
      ! pressure or a star density that is not positive
      ! (ROE_STAR_INADMISSIBLE).  Sound speeds are formed only from a state
      ! that passed that test, so no square root of a negative ratio is
      ! taken and no NaN leaves this routine.
      !
      ! Interface: docs/a2_roe_interface.md, approved 2026-09-05.

      implicit none

      ! Status of the returned star state.
      integer, parameter :: ROE_STAR_OK           = 0   ! admissible star state
      integer, parameter :: ROE_STAR_VACUUM       = 1   ! separating flow, vacuum between the gas edges
      integer, parameter :: ROE_STAR_INADMISSIBLE = 2   ! no branch gave p_star > 0 and rho_star > 0

      type roe_star_state
         real*8 :: p_star                 ! star pressure (0 without an admissible star state)
         real*8 :: rho_L_star, rho_R_star ! star densities (0 without an admissible star state)
         real*8 :: u_star                 ! contact velocity (0 and undefined without one)
         real*8 :: c_L_star, c_R_star     ! star sound speeds (0 without one)
         real*8 :: v_edge_L, v_edge_R     ! gas-edge speeds u_L + 2 c_L/(gam-1), u_R - 2 c_R/(gam-1)
         integer :: branch                ! 1 PVRS, 2 two-rarefaction, 3 two-shock, 0 none
      end type roe_star_state

      contains

      subroutine speed_estimate_ROE(WL,WR,gam,star,status)
      real*8, intent(in) :: WL(3), WR(3)
      real*8, intent(in) :: gam
      type(roe_star_state), intent(out) :: star
      integer, intent(out) :: status
      real*8 :: rhoL,uL,pL,cL
      real*8 :: rhoR,uR,pR,cR
      real*8 :: rho_bar,c_bar
      real*8 :: p_min,p_max,Q
      real*8 :: p_pv,u_pv,rhoL_pv,rhoR_pv
      real*8 :: z,pLR
      real*8 :: AL,BL,AR,BR,gL,gR,p_guess
      logical :: pvrs_admissible
      ! Toro's switch between the estimates: the primitive-variable guess is
      ! kept only while the two pressures are within this ratio of each
      ! other (Toro 2009, eq. 9.43 and the discussion of Q_user).
      real*8, parameter :: Q_user = 2.0d0

      ! Extract left state
      rhoL = WL(1)
      uL   = WL(2)
      pL   = WL(3)

      ! Extract right state
      rhoR = WR(1)
      uR   = WR(2)
      pR   = WR(3)

      ! Nothing is defined until a branch defines it.
      star%p_star     = 0.0d0
      star%rho_L_star = 0.0d0
      star%rho_R_star = 0.0d0
      star%u_star     = 0.0d0
      star%c_L_star   = 0.0d0
      star%c_R_star   = 0.0d0
      star%v_edge_L   = uL
      star%v_edge_R   = uR
      star%branch     = 0
      status = ROE_STAR_INADMISSIBLE

      ! A face state without a positive density and a positive pressure has
      ! no sound speed, so none of the estimates below is defined on it, and
      ! the gas-edge speeds keep the bare velocities.  Written as the
      ! negation of the admissibility test, so that a NaN takes this exit
      ! instead of passing it.
      if (.not. (rhoL .gt. 0.0d0 .and. rhoR .gt. 0.0d0 .and.              &
                 pL   .gt. 0.0d0 .and. pR   .gt. 0.0d0)) return

      ! Sound speeds of the two data states
      cL = sqrt(gam*pL/rhoL)
      cR = sqrt(gam*pR/rhoR)

      ! Speeds of the leading characteristics of the two rarefaction fans,
      ! the edges the gas would reach if it expanded into vacuum (Toro 2009,
      ! eq. 4.76).  Meaningful in every case, and the only meaningful
      ! quantities in the vacuum case.
      star%v_edge_L = uL + 2.0d0*cL/(gam-1.0d0)
      star%v_edge_R = uR - 2.0d0*cR/(gam-1.0d0)

      ! Vacuum test first: when the right gas edge does not reach the left
      ! one the two states never meet, there is no star state, and the exact
      ! solution carries vacuum across the face (Toro 2009, eq. 4.76).
      if (uR - uL .ge. 2.0d0*(cL+cR)/(gam-1.0d0)) then
         status = ROE_STAR_VACUUM
         return
      endif

      !-------------------------------------------------------------!

      p_min = min(pL,pR)
      p_max = max(pL,pR)
      Q     = p_max/p_min          ! p_min > 0 by the test above
      rho_bar = 0.5d0*(rhoL+rhoR)
      c_bar   = 0.5d0*(cL+cR)

      ! Primitive-variable (PVRS) estimate, Toro (2009) eq. 9.20: the star
      ! state of the linearized system about the mean of the two data
      ! states.  It is the accurate one while the two pressures are close,
      ! and it carries no bound of its own, hence the admissibility test.
      p_pv    = 0.5d0*(pL+pR) - 0.5d0*(uR-uL)*rho_bar*c_bar
      u_pv    = 0.5d0*(uL+uR) - 0.5d0*(pR-pL)/(rho_bar*c_bar)
      rhoL_pv = rhoL + (p_pv-pL)/(cL*cL)
      rhoR_pv = rhoR + (p_pv-pR)/(cR*cR)

      pvrs_admissible = (Q .le. Q_user)      .and.                        &
                        (p_pv .ge. p_min)    .and.                        &
                        (p_pv .le. p_max)    .and.                        &
                        (p_pv .gt. 0.0d0)    .and.                        &
                        (rhoL_pv .gt. 0.0d0) .and.                        &
                        (rhoR_pv .gt. 0.0d0)

      if (pvrs_admissible) then

         star%p_star     = p_pv
         star%u_star     = u_pv
         star%rho_L_star = rhoL_pv
         star%rho_R_star = rhoR_pv
         star%branch     = 1

      elseif (p_pv .lt. p_min) then

         ! Two-rarefaction (TRRS) estimate, Toro (2009) eq. 9.32, exact when
         ! both nonlinear waves are rarefactions.  The exponent of the
         ! pressure ratio is 1/z with z = (gam-1)/(2 gam); the base is
         ! positive because the vacuum test above already removed the states
         ! for which it is not.  The star densities follow the isentrope
         ! through each data state (Toro eq. 4.53).
         z   = 0.5d0*(gam-1.0d0)/gam
         pLR = (pL/pR)**z

         star%p_star = ((cL + cR - 0.5d0*(gam-1.0d0)*(uR-uL))/            &
                        (cL/pL**z + cR/pR**z))**(1.0d0/z)
         star%u_star = (pLR*uL/cL + uR/cR                                 &
                        + 2.0d0*(pLR-1.0d0)/(gam-1.0d0))/                 &
                       (pLR/cL + 1.0d0/cR)
         star%rho_L_star = rhoL*(star%p_star/pL)**(1.0d0/gam)
         star%rho_R_star = rhoR*(star%p_star/pR)**(1.0d0/gam)
         star%branch     = 2

      else

         ! Two-shock (TSRS) estimate, Toro (2009) eq. 9.42, exact when both
         ! nonlinear waves are shocks: one step of the shock relations from
         ! the PVRS guess.  The star densities are the Rankine-Hugoniot
         ! compression ratios TIMES the upstream densities (Toro eq. 4.50),
         ! which is what makes them densities and not ratios.
         AL = 2.0d0/((gam+1.0d0)*rhoL)
         BL = (gam-1.0d0)/(gam+1.0d0)*pL
         AR = 2.0d0/((gam+1.0d0)*rhoR)
         BR = (gam-1.0d0)/(gam+1.0d0)*pR

         p_guess = max(0.0d0, p_pv)
         gL = sqrt(AL/(p_guess + BL))
         gR = sqrt(AR/(p_guess + BR))

         star%p_star = (gL*pL + gR*pR - (uR-uL))/(gL+gR)
         star%u_star = 0.5d0*(uL+uR)                                      &
                     + 0.5d0*((star%p_star-pR)*gR - (star%p_star-pL)*gL)
         star%rho_L_star = rhoL*(star%p_star/pL + (gam-1.0d0)/(gam+1.0d0))&
                           /((gam-1.0d0)/(gam+1.0d0)*star%p_star/pL       &
                             + 1.0d0)
         star%rho_R_star = rhoR*(star%p_star/pR + (gam-1.0d0)/(gam+1.0d0))&
                           /((gam-1.0d0)/(gam+1.0d0)*star%p_star/pR       &
                             + 1.0d0)
         star%branch     = 3

      endif

      !-------------------------------------------------------------!

      ! Sound speeds are formed only from a state with a positive pressure
      ! and positive densities.  The test is the negation of admissibility,
      ! so a NaN in any branch above ends here and not in a sqrt.
      if (.not. (star%p_star     .gt. 0.0d0 .and.                         &
                 star%rho_L_star .gt. 0.0d0 .and.                         &
                 star%rho_R_star .gt. 0.0d0)) then
         star%p_star     = 0.0d0
         star%rho_L_star = 0.0d0
         star%rho_R_star = 0.0d0
         star%u_star     = 0.0d0
         star%branch     = 0
         status = ROE_STAR_INADMISSIBLE
         return
      endif

      star%c_L_star = sqrt(gam*star%p_star/star%rho_L_star)
      star%c_R_star = sqrt(gam*star%p_star/star%rho_R_star)
      status = ROE_STAR_OK

      ! End of subroutine
      end subroutine speed_estimate_ROE

      ! End of module
      end module S_estimate_ROE
