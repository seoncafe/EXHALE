   module PLM_reconstruction
   ! Piecewise linear reconstruction
   
   use global_parameters
   use Conversion
   
   implicit none            
   
   contains
   
   subroutine PLM_rec(u_in,WL_rec,WR_rec) 
   integer :: j,k
   real*8, dimension(3,1-Ng:N+Ng),intent(in) :: u_in
   real*8 :: x(-1:1)
   real*8 :: W(3,-1:1)
   real*8, dimension(3) :: sp,sm,sd,sc
   real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL_rec,WR_rec
   
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

      ! Compute derivative approximations
      sp = (W(:,1) - W(:,0))/(x(1)-x(0))
      sm = (W(:,0) - W(:,-1))/(x(0)-x(-1))
      sd = (W(:,1) - W(:,-1))/(x(1)-x(-1))
      
      ! Compute limited slope with MC limiter
      call Minmod_MC(sp,sm,sd,sc)
      
      ! Compute recontructed boundary values
      WL_rec(:,j)   = W(:,0) + 0.5*sc*(x(1)-x(0))
      WR_rec(:,j-1) = W(:,0) - 0.5*sc*(x(0)-x(-1))
   
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
   ! reconstructed face value between the neighbouring cell averages.
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
