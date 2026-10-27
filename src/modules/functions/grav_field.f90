      module grav_func
      ! Expression for the gravitational field and for its derivative,
      ! and the steady isothermal flow speed in that field (the Bernoulli
      ! root shared by the transonic initial condition and the outer ghost).

      use global_parameters
      
      implicit none
      
      contains
      
      
      ! Gravitational potential [code-unit energy per mass].
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

      end function phi
      
      
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

      end function Dphi


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
      real*8 :: vlo, vhi, vmid, hmid, log_lo, log_hi, log_mid, level
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
         ! For q=ln(v/cs), solve exp(2q)/2-q=level on q<=0.
         ! level>=1/2 is the stated Bernoulli-domain condition. The lower
         ! endpoint -level-1 has residual >1, even for exponentially small
         ! speeds; logarithmic bisection resolves them without a velocity floor.
         level = Kr/c2 + log(cs)
         if (level .eq. 0.5d0) then
            isothermal_wind_speed = cs
            return
         endif
         log_lo = -level - 1.0d0
         log_hi = 0.0d0
         do it = 1, 200
            log_mid = 0.5d0*(log_lo + log_hi)
            ! Near the sonic minimum, subtracting 1/2 loses q^2.
            ! The series through q^8 avoids that cancellation for |q|<1e-3.
            if (abs(log_mid) .lt. 1.0d-3) then
               hmid = log_mid**2*(1.0d0+log_mid*(2.0d0/3.0d0+log_mid*   &
                  (1.0d0/3.0d0+log_mid*(2.0d0/15.0d0+log_mid*          &
                  (2.0d0/45.0d0+log_mid*(4.0d0/315.0d0+                &
                  log_mid/315.0d0))))))
            else
               hmid = 0.5d0*exp(2.0d0*log_mid)-log_mid-0.5d0
            endif
            if (hmid .lt. level-0.5d0) then
               log_hi = log_mid
            else
               log_lo = log_mid
            endif
            if (abs(log_hi-log_lo) .le. 4.0d0*epsilon(log_mid)*           &
                max(1.0d0,abs(log_mid))) exit
         enddo
         ! A speed below the representable exponential range underflows to
         ! zero; no finite velocity floor replaces that physical limit.
         isothermal_wind_speed = exp(log(cs)+0.5d0*(log_lo+log_hi))
         return
      endif

      do it = 1, 200
         vmid = 0.5d0*(vlo + vhi)
         hmid = 0.5d0*vmid*vmid - c2*log(vmid)
         ! The subsonic branch has returned; h increases on this branch.
         if (hmid .lt. Kr) then
            vlo = vmid
         else
            vhi = vmid
         endif
         if (abs(vhi - vlo) .lt. 1.0d-12*max(vmid,1.0d-30)) exit
      enddo
      isothermal_wind_speed = 0.5d0*(vlo + vhi)

      end function isothermal_wind_speed
      
      end module grav_func
