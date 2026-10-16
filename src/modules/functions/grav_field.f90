      module grav_func
      ! Expression for the gravitational field and for its derivative,
      ! and the steady isothermal flow speed in that field (the Bernoulli
      ! root shared by the transonic initial condition and the outer ghost).

      use global_parameters
      
      implicit none
      
      contains
      
      !-------------------------------------------------------!
      
      ! Gravitational field strenght (function)
      double precision function phi(r_in)
      real*8, intent(in)  :: r_in

      if (spherical_domain) then
         ! Pure planetary potential (tidal + centrifugal dropped) so the
         ! domain can extend far past the Hill/L1 radius (Huang Case A-like).
         phi = -b0/r_in
      else
         ! Full Roche potential: planet + stellar-tidal + centrifugal.
         phi = -b0/r_in                                &
         	-b0*Mrapp/(atilde-r_in)                    &
         	-b0*(1.0+Mrapp)/(2.0*atilde**3.0)*         &
         	 (atilde*Mrapp/(1.0+Mrapp)-r_in)**2.0
      endif

      ! End of function
      end function phi
      
      !-------------------------------------------------------!
      
      ! Gravitational field gradient
      double precision function Dphi(r_in)
      real*8, intent(in)  :: r_in

      if (spherical_domain) then
         ! d/dr (-b0/r) = b0/r^2
         DPhi = b0/r_in**2.0
      else
         DPhi = b0/r_in**2.0                       &
         	 -b0*Mrapp/(atilde-r_in)**2.0          &
         	 +b0*(1.0 + Mrapp)/atilde**3.0*        &
         	 (atilde*Mrapp/(1.0+Mrapp)-r_in)
      endif

      ! End of function
      end function Dphi

      !-------------------------------------------------------!

      double precision function isothermal_wind_speed(Kr, c2, cs,      &
                                                      supersonic)
      ! THE SPEED OF A STEADY ISOTHERMAL FLOW AT A GIVEN BERNOULLI LEVEL.
      ! Steady, isothermal (p/rho = c2), spherical flow in the potential phi
      ! conserves  1/2 v^2 - c2 ln v - 2 c2 ln r + phi  along the flow, the
      ! Bernoulli integral of the momentum equation with rho v r^2 constant.
      ! At a radius where that integral leaves h(v) = 1/2 v^2 - c2 ln v = Kr,
      ! this returns the root v on the requested branch. h decreases on
      ! (0, cs) and increases on (cs, inf), cs = sqrt(c2), so each branch has
      ! one root when Kr >= h(cs); bisection is unconditionally robust.
      real*8, intent(in)  :: Kr, c2, cs
      logical, intent(in) :: supersonic
      real*8 :: vlo, vhi, vmid, hmid
      integer :: it

      if (supersonic) then
         vlo = cs
         vhi = max(2.0d0*cs, sqrt(2.0d0*max(Kr, c2)))
         ! grow the upper bracket until h(vhi) > Kr
         do it = 1, 200
            if (0.5d0*vhi*vhi - c2*log(vhi) .gt. Kr) exit
            vhi = 2.0d0*vhi
         enddo
      else
         vhi = cs
         vlo = cs*1.0d-12                  ! h -> +inf as v -> 0, so a root exists
      endif

      do it = 1, 200
         vmid = 0.5d0*(vlo + vhi)
         hmid = 0.5d0*vmid*vmid - c2*log(vmid)
         if (supersonic) then
            ! h increasing: if h<Kr the root is to the right
            if (hmid .lt. Kr) then
               vlo = vmid
            else
               vhi = vmid
            endif
         else
            ! h decreasing: if h<Kr the root is to the left (smaller v)
            if (hmid .lt. Kr) then
               vhi = vmid
            else
               vlo = vmid
            endif
         endif
         if (abs(vhi - vlo) .lt. 1.0d-12*max(vmid,1.0d-30)) exit
      enddo
      isothermal_wind_speed = 0.5d0*(vlo + vhi)

      ! End of function
      end function isothermal_wind_speed
      
      ! End of module
      end module grav_func
