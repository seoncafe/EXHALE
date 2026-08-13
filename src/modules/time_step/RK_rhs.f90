      module RK_integration
      ! Evaluate RK right hand side (convection + source)
      
      use global_parameters
      use Numerical_Fluxes
      use source_func
      use low_mach_dissipation, only: low_mach_damping_active,          &
                                      contact_mode_dissipation_flux

      implicit none

      contains

      subroutine RK_rhs(u_in,WL,WR,alpha,dF,S)
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u_in
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: WL,WR
      real*8, intent(in) :: alpha
      integer :: j
      real*8 :: dr
      real*8 :: rp,rm
      real*8 :: dAp,dAm
      real*8 :: dV
      real*8, dimension(3) ::  Fp,Fm
      real*8 :: dF3p
      real*8 :: pL,pR
      ! Gated fourth-difference dissipation of the stagnant-layer contact
      ! mode, added to the numerical flux below so that the marching RHS and
      ! the steady residual (which reaches this routine through
      ! assemble_residual) solve the same equation. Off by default: the array
      ! is then never filled and never read.
      logical :: damp_lowmach
      real*8, dimension(3,1-Ng:N+Ng) :: Ddis
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: dF
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: S

      damp_lowmach = low_mach_damping_active()
      if (damp_lowmach) call contact_mode_dissipation_flux(u_in,Ddis)

      do j = 2-Ng,N+Ng
      
         ! Substitutions
         dr = dr_j(j)
         rp = r_edg(j)
         rm = r_edg(j-1)
         dAp = rp*rp
         dAm = rm*rm
         dV = (dAp*rp - dAm*rm)/3.0
            
         ! Evaluate numerical fluxes
         if (j.eq.(2-Ng)) then
         
               ! Use flux from previous step
               call Num_flux(WL(:,j-1),WR(:,j-1),Fm,alpha,pL)
               if (damp_lowmach) Fm = Fm + Ddis(:,j-1)
         else

               Fm = Fp
               pL = pR
         endif

         ! Evaluate flux at the right interface
         call Num_flux(WL(:,j),WR(:,j),Fp,alpha,pR)
         if (damp_lowmach) Fp = Fp + Ddis(:,j)

         ! Evaluate source
         call source(j,dr,dAp,dAm,dV,    &
                     u_in(:,j),WR(:,j-1),WL(:,j),S(:,j))
      
         ! Evaluate flux differences
         dF(1,j) = (dAp*Fp(1) - dAm*Fm(1))/dV
         dF(2,j) = (dAp*Fp(2) - dAm*Fm(2))/dV 
         
         ! Correct for WENO3 discretization
         if (use_weno3)  dF(2,j) = dF(2,j) + (pR - pL)/dr
         
         dF3p  = dAp*Fp(1)*(Gphi_i(j) - Gphi_c(j))         &
               - dAm*Fm(1)*(Gphi_i(j-1) - Gphi_c(j))
         dF(3,j) = (dAp*Fp(3) - dAm*Fm(3) + dF3p)/dV 
      
      enddo
      
      ! End of subroutine
      end subroutine RK_rhs
      
      ! End of module
      end module RK_integration
