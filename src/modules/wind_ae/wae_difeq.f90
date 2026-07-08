      module wae_difeq
      ! Finite-difference Jacobian assembly, ported verbatim from
      ! wind-ae difeq.c. Builds the solvde s-matrix at the base BC
      ! (k==k1), the sonic-point BC (k>k2), and the interior points
      ! (numerical Jacobian of the residual equations via wae_eval_eqn).
      ! ymod(0:ne,0:1): col 0 perturbs cell k-1, col 1 perturbs cell k.
      use wae_config, only: nsp => wae_nspecies, wae_m, DERIVDIV=>wae_derivdiv
      use wae_types,  only: wae_i_eqnvars, wae_varlist
      use wae_params, only: par => wae_par
      use wae_grid,   only: x => wae_x
      use wae_soe,    only: wae_get_rad, wae_linearize_dvdr_crit
      use wae_eqns,   only: wae_eval_eqn, VEQN, RHOEQN, IONEQN, NCOLEQN,  &
                            TEQN, SPVEQN, SPCRITEQN
      implicit none
      contains

      subroutine wae_difeq_eval(k, k1, k2, jsf, is1, isf, indexv, ne, s, y)
      integer, intent(in)    :: k, k1, k2, jsf, is1, isf, ne
      integer, intent(in)    :: indexv(ne)
      real*8,  intent(inout) :: s(ne, 2*ne+1)
      real*8,  intent(in)    :: y(ne, wae_m)
      if (k .eq. k1) then
         call set_bcs_base(s, y, ne, indexv, jsf)
         call wae_linearize_dvdr_crit(x, y, wae_m, ne)  ! once per iter
      else if (k .gt. k2) then
         call set_bcs_sp(s, y, ne, indexv, jsf)
      else
         call set_interior_eqns(s, y, ne, indexv, jsf, k)
      end if
      end subroutine wae_difeq_eval

      !---------------------------------------------------------------!
      subroutine zero_ymod(ymod, ne)
      integer, intent(in) :: ne
      real*8, intent(out) :: ymod(0:ne,0:1)
      ymod = 0.0d0
      end subroutine zero_ymod

      !---------------------------------------------------------------!
      subroutine set_bcs_base(s, y, ne, indexv, jsf)
      integer, intent(in)    :: ne, indexv(ne), jsf
      real*8,  intent(inout) :: s(ne, 2*ne+1)
      real*8,  intent(in)    :: y(ne, wae_m)
      integer :: j, i, last
      last = 2*nsp + 4
      ! density fixed at the base
      s(last, ne+indexv(1)) = 0.0d0
      s(last, ne+indexv(2)) = 0.0d0
      s(last, ne+indexv(3)) = 1.0d0
      s(last, ne+indexv(4)) = 0.0d0
      do j = 0, nsp-1
         s(last, ne+indexv(j+5))     = 0.0d0
         s(last, ne+indexv(j+5+nsp)) = 0.0d0
      end do
      s(last, jsf) = y(3,1) - par%rho_rmin
      ! ionization fraction fixed at the base, all species
      do j = 0, nsp-1
         s(last-j-1, ne+indexv(1)) = 0.0d0
         s(last-j-1, ne+indexv(2)) = 0.0d0
         s(last-j-1, ne+indexv(3)) = 0.0d0
         s(last-j-1, ne+indexv(4)) = 0.0d0
         do i = 0, nsp-1
            s(last-j-1, ne+indexv(i+5+nsp)) = 0.0d0
            if (i .eq. j) then
               s(last-j-1, ne+indexv(i+5)) = 1.0d0
            else
               s(last-j-1, ne+indexv(i+5)) = 0.0d0
            end if
         end do
         s(last-j-1, jsf) = y(j+5,1) - par%Ys_rmin(j+1)
      end do
      ! temperature fixed at the base
      s(last-nsp-1, ne+indexv(1)) = 0.0d0
      s(last-nsp-1, ne+indexv(2)) = 0.0d0
      s(last-nsp-1, ne+indexv(3)) = 0.0d0
      s(last-nsp-1, ne+indexv(4)) = 1.0d0
      do j = 0, nsp-1
         s(last-nsp-1, ne+indexv(j+5))     = 0.0d0
         s(last-nsp-1, ne+indexv(j+5+nsp)) = 0.0d0
      end do
      s(last-nsp-1, jsf) = y(4,1) - par%T_rmin
      end subroutine set_bcs_base

      !---------------------------------------------------------------!
      subroutine set_bcs_sp(s, y, ne, indexv, jsf)
      integer, intent(in)    :: ne, indexv(ne), jsf
      real*8,  intent(inout) :: s(ne, 2*ne+1)
      real*8,  intent(in)    :: y(ne, wae_m)
      real*8 :: eqn
      type(wae_varlist) :: eqn_pl, eqn_mi, varmod
      real*8 :: ymod(0:ne,0:1)
      integer :: j, i, M
      M = wae_m

      ! Numerator (Bernoulli) condition at the sonic point
      call get_bc_evals(varmod, eqn, eqn_pl, eqn_mi, M, y, ymod, SPCRITEQN, ne)
      s(1, ne+indexv(3)) = (eqn_pl%rho-eqn_mi%rho)/varmod%rho
      s(1, ne+indexv(1)) = (eqn_pl%v-eqn_mi%v)/varmod%v
      s(1, ne+indexv(2)) = (eqn_pl%z-eqn_mi%z)/varmod%z
      s(1, ne+indexv(4)) = (eqn_pl%T-eqn_mi%T)/varmod%T
      do j = 0, nsp-1
         s(1, ne+indexv(j+5))     = (eqn_pl%Ys(j+1)-eqn_mi%Ys(j+1))/varmod%Ys(j+1)
         s(1, ne+indexv(j+5+nsp)) = (eqn_pl%Ncol(j+1)-eqn_mi%Ncol(j+1))/varmod%Ncol(j+1)
      end do
      s(1, jsf) = eqn

      ! velocity (Mach-1) condition at the sonic point
      call get_bc_evals(varmod, eqn, eqn_pl, eqn_mi, M, y, ymod, SPVEQN, ne)
      s(2, ne+indexv(3)) = (eqn_pl%rho-eqn_mi%rho)/varmod%rho
      s(2, ne+indexv(1)) = (eqn_pl%v-eqn_mi%v)/varmod%v
      s(2, ne+indexv(2)) = (eqn_pl%z-eqn_mi%z)/varmod%z
      s(2, ne+indexv(4)) = (eqn_pl%T-eqn_mi%T)/varmod%T
      do j = 0, nsp-1
         s(2, ne+indexv(j+5))     = (eqn_pl%Ys(j+1)-eqn_mi%Ys(j+1))/varmod%Ys(j+1)
         s(2, ne+indexv(j+5+nsp)) = (eqn_pl%Ncol(j+1)-eqn_mi%Ncol(j+1))/varmod%Ncol(j+1)
      end do
      s(2, jsf) = eqn

      ! column density fixed at the sonic point per species
      do j = 0, nsp-1
         s(3+j, ne+indexv(1)) = 0.0d0
         s(3+j, ne+indexv(2)) = 0.0d0
         s(3+j, ne+indexv(3)) = 0.0d0
         s(3+j, ne+indexv(4)) = 0.0d0
         do i = 0, nsp-1
            s(3+j, ne+indexv(i+5)) = 0.0d0
            if (i .eq. j) then
               s(3+j, ne+indexv(i+5+nsp)) = 1.0d0
            else
               s(3+j, ne+indexv(i+5+nsp)) = 0.0d0
            end if
         end do
         s(3+j, jsf) = y(j+5+nsp, M) - par%Ncol_sp(j+1)
      end do
      end subroutine set_bcs_sp

      !---------------------------------------------------------------!
      subroutine set_interior_eqns(s, y, ne, indexv, jsf, k)
      integer, intent(in)    :: ne, indexv(ne), jsf, k
      real*8,  intent(inout) :: s(ne, 2*ne+1)
      real*8,  intent(in)    :: y(ne, wae_m)
      type(wae_i_eqnvars) :: avg
      type(wae_varlist)   :: km1_pl, km1_mi, k_pl, k_mi, varmod
      real*8 :: eqn, ravg, ymod(0:ne,0:1)
      integer :: i, j, last
      last = 4 + 2*nsp

      avg%q   = 0.5d0*(x(k)+x(k-1))
      avg%v   = 0.5d0*(y(1,k)+y(1,k-1))
      avg%z   = 0.5d0*(y(2,k)+y(2,k-1))
      avg%rho = 0.5d0*(y(3,k)+y(3,k-1))
      avg%T   = 0.5d0*(y(4,k)+y(4,k-1))
      do j = 0, nsp-1
         avg%Ys(j+1)   = 0.5d0*(y(j+5,k)+y(j+5,k-1))
         avg%Ncol(j+1) = 0.5d0*(y(j+5+nsp,k)+y(j+5+nsp,k-1))
      end do
      call wae_get_rad(ravg, avg)
      avg%r = ravg

      ! Velocity equation (row 1)
      call get_eqn_evals(varmod, eqn, km1_pl, km1_mi, k_pl, k_mi, avg, k, y, ymod, VEQN, 1, ne)
      call fill_row(s, 1, indexv, ne, km1_pl, km1_mi, k_pl, k_mi, varmod, eqn, jsf)

      ! Density equation (row 2)
      call get_eqn_evals(varmod, eqn, km1_pl, km1_mi, k_pl, k_mi, avg, k, y, ymod, RHOEQN, 1, ne)
      call fill_row(s, 2, indexv, ne, km1_pl, km1_mi, k_pl, k_mi, varmod, eqn, jsf)

      ! z-parameter equation (row 3): direct, z_k - z_{k-1}
      s(3, indexv(1)) = 0.0d0; s(3, indexv(2)) = -1.0d0
      s(3, indexv(3)) = 0.0d0; s(3, indexv(4)) = 0.0d0
      do j = 0, nsp-1
         s(3, indexv(j+5)) = 0.0d0; s(3, indexv(j+5+nsp)) = 0.0d0
      end do
      s(3, ne+indexv(1)) = 0.0d0; s(3, ne+indexv(2)) = 1.0d0
      s(3, ne+indexv(3)) = 0.0d0; s(3, ne+indexv(4)) = 0.0d0
      do j = 0, nsp-1
         s(3, ne+indexv(j+5)) = 0.0d0; s(3, ne+indexv(j+5+nsp)) = 0.0d0
      end do
      s(3, jsf) = y(2,k) - y(2,k-1)

      ! Column density equations (rows 4+j)
      do j = 0, nsp-1
         call get_eqn_evals(varmod, eqn, km1_pl, km1_mi, k_pl, k_mi, avg, k, y, ymod, NCOLEQN, j+1, ne)
         call fill_row(s, 4+j, indexv, ne, km1_pl, km1_mi, k_pl, k_mi, varmod, eqn, jsf)
      end do

      ! Ionization equations (rows 4+j+nsp)
      do j = 0, nsp-1
         call get_eqn_evals(varmod, eqn, km1_pl, km1_mi, k_pl, k_mi, avg, k, y, ymod, IONEQN, j+1, ne)
         call fill_row(s, 4+j+nsp, indexv, ne, km1_pl, km1_mi, k_pl, k_mi, varmod, eqn, jsf)
      end do

      ! Temperature equation (row last)
      call get_eqn_evals(varmod, eqn, km1_pl, km1_mi, k_pl, k_mi, avg, k, y, ymod, TEQN, 1, ne)
      call fill_row(s, last, indexv, ne, km1_pl, km1_mi, k_pl, k_mi, varmod, eqn, jsf)
      end subroutine set_interior_eqns

      !---------------------------------------------------------------!
      subroutine fill_row(s, row, indexv, ne, km1_pl, km1_mi, k_pl, k_mi, varmod, eqn, jsf)
      ! Common pattern: s[row][indexv[var]]    = d(eqn)/d(var at k-1)
      !                 s[row][ne+indexv[var]] = d(eqn)/d(var at k)
      integer, intent(in)    :: row, ne, indexv(ne), jsf
      real*8,  intent(inout) :: s(ne, 2*ne+1)
      type(wae_varlist), intent(in) :: km1_pl, km1_mi, k_pl, k_mi, varmod
      real*8, intent(in) :: eqn
      integer :: j
      s(row, indexv(1)) = (km1_pl%v-km1_mi%v)/varmod%v
      s(row, indexv(2)) = (km1_pl%z-km1_mi%z)/varmod%z
      s(row, indexv(3)) = (km1_pl%rho-km1_mi%rho)/varmod%rho
      s(row, indexv(4)) = (km1_pl%T-km1_mi%T)/varmod%T
      do j = 0, nsp-1
         s(row, indexv(j+5))     = (km1_pl%Ys(j+1)-km1_mi%Ys(j+1))/varmod%Ys(j+1)
         s(row, indexv(j+5+nsp)) = (km1_pl%Ncol(j+1)-km1_mi%Ncol(j+1))/varmod%Ncol(j+1)
      end do
      s(row, ne+indexv(1)) = (k_pl%v-k_mi%v)/varmod%v
      s(row, ne+indexv(2)) = (k_pl%z-k_mi%z)/varmod%z
      s(row, ne+indexv(3)) = (k_pl%rho-k_mi%rho)/varmod%rho
      s(row, ne+indexv(4)) = (k_pl%T-k_mi%T)/varmod%T
      do j = 0, nsp-1
         s(row, ne+indexv(j+5))     = (k_pl%Ys(j+1)-k_mi%Ys(j+1))/varmod%Ys(j+1)
         s(row, ne+indexv(j+5+nsp)) = (k_pl%Ncol(j+1)-k_mi%Ncol(j+1))/varmod%Ncol(j+1)
      end do
      s(row, jsf) = eqn
      end subroutine fill_row

      !---------------------------------------------------------------!
      subroutine get_bc_evals(varmod, eqn_eval, eqn_pl, eqn_mi, k, y, ymod, eqnnum, ne)
      ! Numerical derivative of `eqnnum` wrt each variable at a single
      ! boundary cell (ymod(:,1) only). Mirrors difeq.c:get_bc_evals.
      integer, intent(in)  :: k, eqnnum, ne
      real*8,  intent(in)  :: y(ne, wae_m)
      real*8,  intent(inout) :: ymod(0:ne,0:1)
      type(wae_varlist), intent(out) :: varmod, eqn_pl, eqn_mi
      real*8, intent(out) :: eqn_eval
      integer :: j

      call zero_ymod(ymod, ne)
      eqn_eval = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)

      varmod%rho = y(3,k)/DERIVDIV
      call zero_ymod(ymod, ne); ymod(3,1) = varmod%rho/2.0d0
      eqn_pl%rho = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(3,1) = -varmod%rho/2.0d0
      eqn_mi%rho = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)

      varmod%v = y(1,k)/DERIVDIV
      call zero_ymod(ymod, ne); ymod(1,1) = varmod%v/2.0d0
      eqn_pl%v = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(1,1) = -varmod%v/2.0d0
      eqn_mi%v = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)

      varmod%T = y(4,k)/DERIVDIV
      call zero_ymod(ymod, ne); ymod(4,1) = varmod%T/2.0d0
      eqn_pl%T = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(4,1) = -varmod%T/2.0d0
      eqn_mi%T = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)

      do j = 0, nsp-1
         varmod%Ys(j+1) = y(j+5,k)/DERIVDIV
         call zero_ymod(ymod, ne); ymod(j+5,1) = varmod%Ys(j+1)/2.0d0
         eqn_pl%Ys(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)
         call zero_ymod(ymod, ne); ymod(j+5,1) = -varmod%Ys(j+1)/2.0d0
         eqn_mi%Ys(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)
      end do

      do j = 0, nsp-1
         varmod%Ncol(j+1) = y(j+5+nsp,k)/DERIVDIV
         call zero_ymod(ymod, ne); ymod(j+5+nsp,1) = varmod%Ncol(j+1)/2.0d0
         eqn_pl%Ncol(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)
         call zero_ymod(ymod, ne); ymod(j+5+nsp,1) = -varmod%Ncol(j+1)/2.0d0
         eqn_mi%Ncol(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)
      end do

      varmod%z = y(2,k)/DERIVDIV
      call zero_ymod(ymod, ne); ymod(2,1) = varmod%z/2.0d0
      eqn_pl%z = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(2,1) = -varmod%z/2.0d0
      eqn_mi%z = wae_eval_eqn(k, x, y, ymod, eqnnum, 1, wae_m, ne)
      end subroutine get_bc_evals

      !---------------------------------------------------------------!
      subroutine get_eqn_evals(varmod, eqn_eval, km1_pl, km1_mi, k_pl, k_mi, &
                               avg, k, y, ymod, eqnnum, species, ne)
      ! Numerical derivative of `eqnnum` wrt each variable at cells k-1
      ! (ymod(:,0)) and k (ymod(:,1)). Step sizes use avg.*/DERIVDIV.
      ! Mirrors difeq.c:get_eqn_evals.
      integer, intent(in) :: k, eqnnum, species, ne
      real*8,  intent(in) :: y(ne, wae_m)
      type(wae_i_eqnvars), intent(in) :: avg
      real*8, intent(inout) :: ymod(0:ne,0:1)
      type(wae_varlist), intent(out) :: varmod, km1_pl, km1_mi, k_pl, k_mi
      real*8, intent(out) :: eqn_eval
      integer :: j

      call zero_ymod(ymod, ne)
      eqn_eval = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)

      ! rho
      varmod%rho = avg%rho/DERIVDIV
      call zero_ymod(ymod, ne); ymod(3,0) = varmod%rho/2.0d0
      km1_pl%rho = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(3,0) = -varmod%rho/2.0d0
      km1_mi%rho = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(3,1) = varmod%rho/2.0d0
      k_pl%rho = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(3,1) = -varmod%rho/2.0d0
      k_mi%rho = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)

      ! v
      varmod%v = avg%v/DERIVDIV
      call zero_ymod(ymod, ne); ymod(1,0) = varmod%v/2.0d0
      km1_pl%v = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(1,0) = -varmod%v/2.0d0
      km1_mi%v = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(1,1) = varmod%v/2.0d0
      k_pl%v = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(1,1) = -varmod%v/2.0d0
      k_mi%v = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)

      ! T
      varmod%T = avg%T/DERIVDIV
      call zero_ymod(ymod, ne); ymod(4,0) = varmod%T/2.0d0
      km1_pl%T = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(4,0) = -varmod%T/2.0d0
      km1_mi%T = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(4,1) = varmod%T/2.0d0
      k_pl%T = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(4,1) = -varmod%T/2.0d0
      k_mi%T = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)

      ! Ys
      do j = 0, nsp-1
         varmod%Ys(j+1) = avg%Ys(j+1)/DERIVDIV
         call zero_ymod(ymod, ne); ymod(j+5,0) = varmod%Ys(j+1)/2.0d0
         km1_pl%Ys(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
         call zero_ymod(ymod, ne); ymod(j+5,0) = -varmod%Ys(j+1)/2.0d0
         km1_mi%Ys(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
         call zero_ymod(ymod, ne); ymod(j+5,1) = varmod%Ys(j+1)/2.0d0
         k_pl%Ys(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
         call zero_ymod(ymod, ne); ymod(j+5,1) = -varmod%Ys(j+1)/2.0d0
         k_mi%Ys(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      end do

      ! Ncol
      do j = 0, nsp-1
         varmod%Ncol(j+1) = avg%Ncol(j+1)/DERIVDIV
         call zero_ymod(ymod, ne); ymod(j+5+nsp,0) = varmod%Ncol(j+1)/2.0d0
         km1_pl%Ncol(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
         call zero_ymod(ymod, ne); ymod(j+5+nsp,0) = -varmod%Ncol(j+1)/2.0d0
         km1_mi%Ncol(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
         call zero_ymod(ymod, ne); ymod(j+5+nsp,1) = varmod%Ncol(j+1)/2.0d0
         k_pl%Ncol(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
         call zero_ymod(ymod, ne); ymod(j+5+nsp,1) = -varmod%Ncol(j+1)/2.0d0
         k_mi%Ncol(j+1) = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      end do

      ! z
      varmod%z = avg%z/DERIVDIV
      call zero_ymod(ymod, ne); ymod(2,0) = varmod%z/2.0d0
      km1_pl%z = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(2,0) = -varmod%z/2.0d0
      km1_mi%z = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(2,1) = varmod%z/2.0d0
      k_pl%z = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      call zero_ymod(ymod, ne); ymod(2,1) = -varmod%z/2.0d0
      k_mi%z = wae_eval_eqn(k, x, y, ymod, eqnnum, species, wae_m, ne)
      end subroutine get_eqn_evals

      end module wae_difeq
