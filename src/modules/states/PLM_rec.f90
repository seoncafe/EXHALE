   module PLM_reconstruction
   ! Piecewise linear reconstruction
   
   use global_parameters
   use Conversion
   
   implicit none            
   
   contains
   
   subroutine PLM_rec(u_in,WL_rec,WR_rec,P_up,P_dn,dev_L,dev_R)
   ! WELL-BALANCED OPTION.  With the four optional arrays present (the caller
   ! passes them when "Well balanced:" is set) the PRESSURE is reconstructed
   ! in the coordinate of the departure from the cell's own local hydrostatic
   ! equilibrium: the stencil data becomes the neighboring cell's pressure
   ! measured against that equilibrium continued through the face between
   ! them (this cell's density up to the face, the neighbor's beyond it),
   ! with the value at the cell's own centre identically zero, the limiter
   ! acts on that unchanged, and the equilibrium face value P_up / P_dn is
   ! added back at the end.  The departure itself is returned so that the
   ! Riemann jump and the pressure force can be formed from small numbers.
   ! Density and velocity are reconstructed as they are without it
   ! (Kaeppeli and Mishra 2016, A&A 587, A94, section 2.1.3).
   integer :: j,k
   real*8, dimension(3,1-Ng:N+Ng),intent(in) :: u_in
   real*8 :: x(-1:1)
   real*8 :: W(3,-1:1)
   real*8, dimension(3) :: sp,sm,sd,sc
   real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL_rec,WR_rec
   real*8, dimension(1-Ng:N+Ng), intent(in),  optional :: P_up,P_dn
   real*8, dimension(1-Ng:N+Ng), intent(out), optional :: dev_L,dev_R
   logical :: wb

   wb = present(P_up) .and. present(P_dn) .and.                          &
        present(dev_L) .and. present(dev_R)

   ! Cell-local: the limited slope of cell j is built from the three cell
   ! averages j-1, j, j+1 and written into that cell's own two face states.
   ! No reduction, so a cell's arithmetic is unchanged and the result is
   ! bitwise identical at any number of threads.
   !$omp parallel do default(shared) schedule(static)                   &
   !$omp   private(j,k,x,W,sp,sm,sd,sc)
   do j = 2-Ng,N+Ng-1

      ! Extract stencil grid
      x = r(j-1:j+1)

      ! Convert to local primitive variables (2nd order conversion)
      do k = j-1,j+1
            call U_to_W_comp(u_in(:,k),W(:,k-j),k)
      enddo

      ! WELL-BALANCED OPTION: the pressure stencil in departure coordinates.
      ! The departure at a neighboring cell centre is measured against the
      ! equilibrium continued to it through the face, that is with THIS
      ! cell's density up to the face and the NEIGHBOUR's beyond it, which is
      ! what makes the data vanish exactly on the discrete equilibrium the
      ! scheme preserves (the two neighboring face extrapolations agree
      ! there; Kaeppeli and Mishra 2016, A&A 587, A94, their eqs. 18 and 26,
      ! where the same construction appears as the average density of a
      ! uniform mesh).  Continuing the cell's own constant density over the
      ! whole gap instead would leave data of the size (rho_j - rho_j+1) x
      ! (potential difference), and the reconstruction would be second order
      ! and not exact.
      ! The difference of the two cell pressures is taken FIRST, while it is
      ! still an exact floating-point operation, and the hydrostatic terms
      ! are added to it afterwards; the reverse order would bury the small
      ! terms in the rounding of an O(1) sum.
      if (wb) then
         W(3,1)  = (W(3,1)  - W(3,0))                                    &
                 + W(1,0)*(Gphi_i(j)   - Gphi_c(j))                      &
                 + W(1,1)*(Gphi_c(j+1) - Gphi_i(j))
         W(3,-1) = (W(3,-1) - W(3,0))                                    &
                 - W(1,0) *(Gphi_c(j)   - Gphi_i(j-1))                   &
                 - W(1,-1)*(Gphi_i(j-1) - Gphi_c(j-1))
         W(3,0)  = 0.0d0
      endif

      ! Compute derivative approximations
      sp = (W(:,1) - W(:,0))/(x(1)-x(0))
      sm = (W(:,0) - W(:,-1))/(x(0)-x(-1))
      sd = (W(:,1) - W(:,-1))/(x(1)-x(-1))
      
      ! Compute limited slope with MC limiter
      call Minmod_MC(sp,sm,sd,sc)
      
      ! Compute recontructed boundary values
      WL_rec(:,j)   = W(:,0) + 0.5*sc*(x(1)-x(0))
      WR_rec(:,j-1) = W(:,0) - 0.5*sc*(x(0)-x(-1))

      ! The departure the limited slope returned, and the face pressure it
      ! belongs to
      if (wb) then
         dev_L(j)        = WL_rec(3,j)
         dev_R(j-1)      = WR_rec(3,j-1)
         WL_rec(3,j)     = P_up(j) + dev_L(j)
         WR_rec(3,j-1)   = P_dn(j) + dev_R(j-1)
      endif

   enddo
   !$omp end parallel do
   
   ! End of subroutine
   end subroutine PLM_rec
   
   !--------------------------------------------
   
   ! Minmod generalized slope limiter
   
   subroutine Minmod_MC(a,b,c,d)
   real*8, dimension(1,3), intent(in) :: a,b,c
   integer :: i
   
   real*8, dimension(1,3), intent(out) :: d

   ! Generalized Minmod limiter (Kappeli 2016), one component at a time.
   ! The limiter itself is minmod_mc_slope below and exists only there, so
   ! the scalar reconstruction of a mass fraction and the reconstruction of
   ! the primitive vector cannot drift apart.
            
   do i = 1,3
      d(1,i) = minmod_mc_slope(a(1,i),b(1,i),c(1,i))
   enddo
   
   ! End of subroutine
   end subroutine Minmod_MC

   !--------------------------------------------

   real*8 function minmod_mc_slope(a,b,c) result(d)
   ! Generalized minmod (monotonized central) limited slope of one quantity
   ! from its right, left and central difference quotients a, b, c
   ! (Kappeli 2016).  theta = 2 is the least diffusive value for which the
   ! scheme stays total-variation diminishing.  The slope is zero whenever
   ! the three arguments do not all share a sign, which is what keeps a
   ! reconstructed face value between the neighboring cell averages.
   real*8, intent(in) :: a,b,c
   real*8 :: v_arg(3)
   real*8 :: theta = 2.0

   v_arg(1) = theta*a
   v_arg(2) = theta*b
   v_arg(3) = c

   if (maxval(v_arg).lt.(0.0)) then
      d =maxval(v_arg)
   elseif (minval(v_arg).gt.(0.0)) then
      d = minval(v_arg)
   else
      d = 0.0
   endif     
   
   ! End of function
   end function minmod_mc_slope

   !--------------------------------------------

   subroutine PLM_rec_scalar(q,qL,qR)
   ! Piecewise linear reconstruction of ONE cell-centered scalar, with the
   ! same limited slope minmod_mc_slope the primitive variables are
   ! reconstructed with, so that a mass fraction rides on the same
   ! discretization as the density whose face flux carries it.
   !
   ! qL(j) is the state on the LEFT of face r_edg(j), extrapolated from cell
   ! j; qR(j) the state on its right, extrapolated from cell j+1.  That is
   ! the pairing WL/WR use.  The two outermost faces have no complete
   ! three-cell stencil and keep the donor cell average, which is the
   ! first-order value the zero-gradient outer ghosts carry in any case.
   integer :: j
   real*8, dimension(1-Ng:N+Ng), intent(in)  :: q
   real*8 :: x(-1:1)
   real*8 :: sp,sm,sd,sc
   real*8, dimension(1-Ng:N+Ng), intent(out) :: qL,qR

   do j = 1-Ng,N+Ng
      qL(j) = q(j)
      qR(j) = q(min(j+1,N+Ng))
   enddo

   !$omp parallel do default(shared) schedule(static)                   &
   !$omp   private(j,x,sp,sm,sd,sc)
   do j = 2-Ng,N+Ng-1

      x = r(j-1:j+1)

      sp = (q(j+1) - q(j))/(x(1)-x(0))
      sm = (q(j) - q(j-1))/(x(0)-x(-1))
      sd = (q(j+1) - q(j-1))/(x(1)-x(-1))

      sc = minmod_mc_slope(sp,sm,sd)

      qL(j)   = q(j) + 0.5*sc*(x(1)-x(0))
      qR(j-1) = q(j) - 0.5*sc*(x(0)-x(-1))

   enddo
   !$omp end parallel do

   ! End of subroutine
   end subroutine PLM_rec_scalar

   ! End of module
   end module PLM_reconstruction
