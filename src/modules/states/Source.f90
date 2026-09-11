   module source_func
   ! Subroutine to evaluate the source function
   !     (gravitational + geometrical)
   
   use global_parameters
   use Conversion
   use grav_func
   use caloric_eos, only: pressure_from_energy_density
   
   implicit none
   
   contains
   
   subroutine source(j,dr,dAp,dAm,dV,u,WR,WL,S)
   integer, intent(in) :: j
   real*8, intent(in) :: dr,dAp,dAm,dV
   real*8, intent(in) :: u(3),WR(3),WL(3) 
   real*8 :: rhoL,rhoR
   real*8 :: pC
   real*8, intent(out) :: S(3)
   
   !--- Extract physical quantities ---!

   ! Density
   rhoL = WR(1)
   rhoR = WL(1)       
   
   ! Central ressure
   if (.not. well_balanced)                                              &
      pC = pressure_from_energy_density(j, u(1), u(3)-0.5*u(2)*u(2)/u(1))
   
   !--- Evaluare source term ---!

   ! WELL-BALANCED OPTION ("Well balanced:").  The gravitational source and the
   ! geometric pressure term are then not evaluated here at all: they cancel
   ! the equilibrium part of the momentum flux difference IDENTICALLY, and
   ! RK_rhs assembles the momentum row from what is left, the departure of
   ! the face pressure from the cell's own hydrostatic equilibrium (the
   ! identities are written out at the cell loop of RK_rhs; the source is the
   ! face-pressure difference of the cell's own equilibrium, Kaeppeli and
   ! Mishra 2014, J. Comput. Phys. 259, 199, eq. 2.26).  Returning zero here
   ! is therefore not the absence of gravity: it is gravity carried in the
   ! flux difference in a form that vanishes on the discrete equilibrium.
   if (well_balanced) then
      S = 0.0d0
      return
   endif

   S(1) = 0.0
   S(2) = - 0.5*(rhoL + rhoR)*(Gphi_i(j) - Gphi_i(j-1))/dr
   S(3) = 0.0
   
   ! Add source term explicitly if PLM is used
   if (use_plm) S(2) = S(2) + (dAp-dAm)/dV*pC

   ! End of subroutine
   end subroutine source
   
   ! End of module
   end module source_func
