      module S_estimate_ROE
      ! Estimate for left and right velocities for the Roe method
      !     based on Toro (2009) [ch. 11.3-4]
      
      use global_parameters
            
      implicit none
      
      contains
      
      ! gam is the SINGLE adiabatic index this estimate is entitled to.
      ! Toro's two-rarefaction and two-shock p_star estimates below are
      ! derived for one ideal gas on both sides, so the caller passes the
      ! mean of the two face gammas rather than two of them; see the note
      ! in Num_flux.  For a molecule-free face that mean is the monatomic
      ! constant exactly.
      subroutine speed_estimate_ROE(WL,WR,u_star,cL_star,cR_star,gam)
      real*8, intent(in) :: WL(3), WR(3)
      real*8, intent(in) :: gam
      real*8 :: rhoL,uL,pL,cL            
      real*8 :: rhoR,uR,pR,cR            
      real*8 :: rho_bar,c_bar            
      real*8 :: Q,Q_user                
      real*8 :: p_min,p_max,p_star      
      real*8 :: z,pLR                    
      real*8 :: AL,BL,AR,BR,gL,gR        
      real*8 :: rhoL_star,rhoR_star
      real*8, intent(out) :: u_star,cL_star,cR_star
                
      ! Exctract left state
      rhoL = WL(1)
      uL = WL(2)
      pL = WL(3)
      
      ! Exctract right state
      rhoR = WR(1)
      uR = WR(2)
      pR = WR(3)
      
      ! Sound speed
      cL = sqrt(gam*pL/rhoL)
      cR = sqrt(gam*pR/rhoR)
      
      ! Average density and sound speed
      rho_bar = 0.5*(rhoL+rhoR)
      c_bar   = 0.5*(cL+cR)     
      
      !-------------------------------------------------------------!
      
      ! Substitutions
      p_min = min(pL,pR)
      p_max = max(pL,pR)
      ! Floor: guards Q = p_max/p_min against a pathological zero/negative
      ! left/right pressure (code review 2026-07-02). No-op for healthy
      ! states (any positive p_min above 1e-30*p_max is unchanged).
      p_min = max(p_min, 1.0d-30*max(p_max,1.0d-300))

      Q = p_max/p_min
      Q_user = 2.0

      ! Initial PVRS guess
      p_star = 0.5*(pL+pR)-0.5*(uR-uL)*rho_bar*c_bar
      u_star = 0.5*(uL+uR)-0.5*(pR-pL)/(c_bar*rho_bar)
      rhoL_star = rhoL + (p_star-pL)/(cL*cL)
      rhoR_star = rhoR + (p_star-pR)/(cR*cR)
      
      ! Construct middle state approximation
      if(Q.gt.Q_user) then      ! PVRS
      
         if(p_star.le.p_min) then      ! TRRS
               
            ! Subs
            z = 0.5*(gam-1.0)/gam
            pLR = (pL/pR)**z
            
            p_star = sqrt((cL+cR-0.5*(gam-1.0)*(uR-uL))/ &
                           (cL/pL**z + cR/pR**z)) 
            u_star = (pLR*uL/cL+uR/cR                      &
                        +2.0*(pLR-1.0)/(gam-1.0))/     &
                     (pLR/cL+1.0/cR)
            
            rhoL_star = rhoL*(p_star/pL)**(1.0/gam)
            rhoR_star = rhoR*(p_star/pR)**(1.0/gam)
               
         elseif(p_star.ge.p_max) then  ! TSRS
               
            ! Subs
            AL = 2.0/((gam+1.0)*rhoL)
            BL = (gam-1.0)/(gam+1.0)*pL
            AR = 2.0/((gam+1.0)*rhoR)
            BR = (gam-1.0)/(gam+1.0)*pR
            
            p_star = max(0.0, p_star)
            gL = sqrt(AL/(p_star + BL))
            gR = sqrt(AR/(p_star + BR))
            
            p_star = (gL*pL+gR*pR-uR+uL)/(gL+gR)
            u_star = 0.5*(uL+uR)+ &
                     0.5*((p_star-pR)*gR-(p_star-pL)*gL)
            
            rhoL_star = (p_star/pL + (gam-1.0)/(gam+1.0))/      &
                        ((gam-1.0)/(gam+1.0)*p_star/pL+1.0)
            rhoR_star = (p_star/pR + (gam-1.0)/(gam+1.0))/      &
                        ((gam-1.0)/(gam+1.0)*p_star/pR+1.0)
               
         endif
      endif
      
      !-------------------------------------------------------------!
      
      ! Construct sound speed estimates      
      cL_star = sqrt(gam*p_star/rhoL_star)
      cR_star = sqrt(gam*p_star/rhoR_star)
      
      ! End of subroutine
      end subroutine speed_estimate_ROE
      
      ! End of module
      end module S_estimate_ROE
