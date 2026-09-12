! Extracted production routines for direct comparison; bodies unchanged.
module scalar_reconstruction_probe
implicit none
integer :: N=500
integer, parameter :: Ng=2
real*8, allocatable :: r(:)
contains
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
end module
program probe
use scalar_reconstruction_probe
use omp_lib
implicit none
real*8, allocatable :: q(:),ql(:),qr(:),left_reference(:),right_reference(:)
real*8 :: start_time,elapsed,difference,elapsed_best
integer :: j,k,rep,nt,team,case_id
integer, parameter :: repeats=3000
call omp_set_dynamic(.false.)
do case_id=1,2
N=500
if (case_id==2) N=8000
allocate(r(1-Ng:N+Ng),q(1-Ng:N+Ng),ql(1-Ng:N+Ng),qr(1-Ng:N+Ng), &
left_reference(1-Ng:N+Ng),right_reference(1-Ng:N+Ng))
do j=1-Ng,N+Ng
r(j)=exp(2d0*dble(j)/dble(N))
q(j)=1d0+0.2d0*sin(10d0*r(j))
enddo
call omp_set_num_threads(1)
call PLM_rec_scalar(q,left_reference,right_reference)
do k=0,4
nt=2**k
call omp_set_num_threads(nt)
team=0
!$omp parallel
!$omp single
team=omp_get_num_threads()
!$omp end single
!$omp end parallel
elapsed_best=huge(1d0)
difference=0d0
do rep=1,3
start_time=omp_get_wtime()
do j=1,repeats
call PLM_rec_scalar(q,ql,qr)
enddo
elapsed=omp_get_wtime()-start_time
elapsed_best=min(elapsed_best,elapsed)
difference=max(difference,maxval(abs(ql-left_reference)),maxval(abs(qr-right_reference)))
enddo
write(*,'(I6,1X,I3,1X,I3,1X,ES15.7,1X,ES15.7)') N,nt,team,elapsed_best/repeats,difference
enddo
deallocate(r,q,ql,qr,left_reference,right_reference)
enddo
end program
