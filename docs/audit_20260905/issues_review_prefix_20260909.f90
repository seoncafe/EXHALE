! Isolated numerical review harness. The shell driver inserts the current
! production pgmres and dogleg_step routines without editing them.
! This is not a hydrodynamic integration or a planetary convergence test.
module numerical_review
  implicit none
  integer, parameter :: nvar_jac=1, N=1, Ng=0, n_species=1
  integer, parameter :: kl_jac=0, ku_jac=0, nspec_row=1
  real*8 :: operator_value=1d0
contains
  subroutine jv_product(Y,F0,f_sp_base,v,Jv,ok)
    real*8, intent(in) :: Y(1),F0(1),f_sp_base(1,1),v(1)
    real*8, intent(out) :: Jv(1)
    logical, intent(out) :: ok
    Jv=operator_value*v
    ok=.true.
  end subroutine

  subroutine dgbtrs(trans,nrow,kl,ku,nrhs,abf,ldab,ipiv,z,ldz,info)
    ! Identity preconditioner: isolate Arnoldi and triangular solution.
    character, intent(in) :: trans
    integer, intent(in) :: nrow,kl,ku,nrhs,ldab,ldz,ipiv(nrow)
    real*8, intent(in) :: abf(ldab,nrow)
    real*8, intent(inout) :: z(ldz)
    integer, intent(out) :: info
    info=0
  end subroutine
