      module wae_eqns
      ! Residual equations + eval_eqn dispatcher, ported from soe.c.
      ! ymod(0:ne,0:1): ymod(:,1) perturbs cell k, ymod(:,0) perturbs k-1,
      ! matching the C ymod[i][j] convention (row i = y row: 1=v,2=z,
      ! 3=rho,4=T,5..=Ys,5+nsp..=Ncol).
      use wae_config, only: nsp => wae_nspecies, wae_m
      use wae_types,  only: wae_i_eqnvars
      use wae_params, only: par => wae_par
      use wae_soe,    only: wae_get_mu, wae_get_gamma, wae_get_dvdr,      &
                            wae_get_drhodr, wae_get_dYsdr, wae_get_dTdr,  &
                            wae_get_dNcoldr, wae_get_spQ, wae_erf_norm
      implicit none

      ! eqn numbers (defs.h)
      integer, parameter :: VEQN=1, RHOEQN=2, IONEQN=3, NCOLEQN=4
      integer, parameter :: TEQN=5, SPVEQN=6, SPCRITEQN=7
      contains

      real*8 function wae_eval_eqn(k, x, y, ymod, eqnnum, species, m, ne)
      integer, intent(in) :: k, eqnnum, species, m, ne
      real*8,  intent(in) :: x(m), y(ne,m), ymod(0:ne,0:1)
      select case (eqnnum)
      case (VEQN);      wae_eval_eqn = v_eqn(k, x, y, ymod, m, ne)
      case (RHOEQN);    wae_eval_eqn = rho_eqn(k, x, y, ymod, m, ne)
      case (TEQN);      wae_eval_eqn = T_eqn(k, x, y, ymod, m, ne)
      case (SPVEQN);    wae_eval_eqn = spv_eqn(k, x, y, ymod, m, ne)
      case (SPCRITEQN); wae_eval_eqn = spcrit_eqn(k, x, y, ymod, m, ne)
      case (IONEQN);    wae_eval_eqn = ion_eqn(k, x, y, ymod, species, m, ne)
      case (NCOLEQN);   wae_eval_eqn = Ncol_eqn(k, x, y, ymod, species, m, ne)
      case default
         write(*,*) '(wae_eval_eqn) Unrecognized eqnnum:', eqnnum
         stop 701
      end select
      end function wae_eval_eqn

      !---------------------------------------------------------------!
      subroutine set_vars(k, x, y, kv, km1, avg, ymod, m, ne)
      integer, intent(in) :: k, m, ne
      real*8,  intent(in) :: x(m), y(ne,m), ymod(0:ne,0:1)
      type(wae_i_eqnvars), intent(out) :: kv, km1, avg
      integer :: j
      kv%q = x(k)
      kv%v = y(1,k)+ymod(1,1); kv%z = y(2,k)+ymod(2,1)
      kv%rho = y(3,k)+ymod(3,1); kv%T = y(4,k)+ymod(4,1)
      do j = 1, nsp
         kv%Ys(j)   = y(j+4,k)+ymod(j+4,1)
         kv%Ncol(j) = y(j+4+nsp,k)+ymod(j+4+nsp,1)
      end do
      km1%q = x(k-1)
      km1%v = y(1,k-1)+ymod(1,0); km1%z = y(2,k-1)+ymod(2,0)
      km1%rho = y(3,k-1)+ymod(3,0); km1%T = y(4,k-1)+ymod(4,0)
      do j = 1, nsp
         km1%Ys(j)   = y(j+4,k-1)+ymod(j+4,0)
         km1%Ncol(j) = y(j+4+nsp,k-1)+ymod(j+4+nsp,0)
      end do
      avg%q = 0.5d0*(kv%q+km1%q);   avg%rho = 0.5d0*(kv%rho+km1%rho)
      avg%v = 0.5d0*(kv%v+km1%v);   avg%z = 0.5d0*(kv%z+km1%z)
      avg%T = 0.5d0*(kv%T+km1%T)
      do j = 1, nsp
         avg%Ys(j)   = 0.5d0*(kv%Ys(j)+km1%Ys(j))
         avg%Ncol(j) = 0.5d0*(kv%Ncol(j)+km1%Ncol(j))
      end do
      end subroutine set_vars

      subroutine set_erf_norm_k2(kv)
      type(wae_i_eqnvars), intent(in) :: kv
      wae_erf_norm = 1.0d0 - erf((kv%v-par%erf_drop(1))/par%erf_drop(2))
      end subroutine set_erf_norm_k2

      !---------------------------------------------------------------!
      real*8 function spv_eqn(k, x, y, ymod, m, ne)
      integer, intent(in) :: k, m, ne
      real*8,  intent(in) :: x(m), y(ne,m), ymod(0:ne,0:1)
      type(wae_i_eqnvars) :: kv
      real*8 :: mu, gamma
      integer :: j
      kv%q = x(k)
      kv%v = y(1,k)+ymod(1,1); kv%z = y(2,k)+ymod(2,1)
      kv%rho = y(3,k)+ymod(3,1); kv%T = y(4,k)+ymod(4,1)
      do j = 1, nsp
         kv%Ys(j) = y(j+4,k)+ymod(j+4,1); kv%Ncol(j) = y(j+4+nsp,k)+ymod(j+4+nsp,1)
      end do
      call wae_get_mu(mu, kv); call wae_get_gamma(gamma)
      spv_eqn = kv%v - sqrt(kv%T*gamma/mu)*par%breezeparam
      end function spv_eqn

      real*8 function spcrit_eqn(k, x, y, ymod, m, ne)
      integer, intent(in) :: k, m, ne
      real*8,  intent(in) :: x(m), y(ne,m), ymod(0:ne,0:1)
      type(wae_i_eqnvars) :: kv
      real*8 :: T_term, grav_term, Q_term, spQ, r, mu, gamma, qinv, norm_a
      integer :: j
      kv%q = x(k)
      kv%v = y(1,k)+ymod(1,1); kv%z = y(2,k)+ymod(2,1)
      kv%rho = y(3,k)+ymod(3,1); kv%T = y(4,k)+ymod(4,1)
      do j = 1, nsp
         kv%Ys(j) = y(j+4,k)+ymod(j+4,1); kv%Ncol(j) = y(j+4+nsp,k)+ymod(j+4+nsp,1)
      end do
      r = par%Rmin + kv%q*kv%z
      grav_term = -1.0d0/r
      norm_a = par%semimajor/par%Rp
      qinv = par%Mstar/par%Mp
      grav_term = grav_term + (qinv*r*                                   &
         (-(norm_a - r + r/qinv)/norm_a**3 + 1.0d0/(norm_a-r)**2))*par%tidalforce
      grav_term = grav_term*par%Rp/par%H0
      call wae_get_mu(mu, kv); call wae_get_gamma(gamma)
      T_term = 2.0d0*gamma*kv%T/mu
      call wae_get_spQ(spQ, kv, dble(k))
      Q_term = -(gamma-1.0d0)*spQ*r/kv%v
      spcrit_eqn = T_term + grav_term + Q_term
      end function spcrit_eqn

      real*8 function ion_eqn(k, x, y, ymod, species, m, ne)
      integer, intent(in) :: k, species, m, ne
      real*8,  intent(in) :: x(m), y(ne,m), ymod(0:ne,0:1)
      type(wae_i_eqnvars) :: kv, km1, avg
      real*8 :: delta_Ys, delta_r, dYsdr_avg(nsp)
      call set_vars(k, x, y, kv, km1, avg, ymod, m, ne)
      if (k .eq. 2) call set_erf_norm_k2(kv)
      delta_r = avg%z*(kv%q-km1%q)
      call wae_get_dYsdr(dYsdr_avg, avg, dble(k))
      delta_Ys = kv%Ys(species) - km1%Ys(species)
      delta_r = avg%z*(kv%q-km1%q)
      ion_eqn = delta_Ys - dYsdr_avg(species)*delta_r
      end function ion_eqn

      real*8 function Ncol_eqn(k, x, y, ymod, species, m, ne)
      integer, intent(in) :: k, species, m, ne
      real*8,  intent(in) :: x(m), y(ne,m), ymod(0:ne,0:1)
      type(wae_i_eqnvars) :: kv, km1, avg
      real*8 :: delta_N, delta_r, dNdr_avg(nsp)
      call set_vars(k, x, y, kv, km1, avg, ymod, m, ne)
      if (k .eq. 2) call set_erf_norm_k2(kv)
      delta_r = avg%z*(kv%q-km1%q)
      call wae_get_dNcoldr(dNdr_avg, avg)
      delta_N = kv%Ncol(species) - km1%Ncol(species)
      Ncol_eqn = delta_N - dNdr_avg(species)*delta_r
      end function Ncol_eqn

      real*8 function v_eqn(k, x, y, ymod, m, ne)
      integer, intent(in) :: k, m, ne
      real*8,  intent(in) :: x(m), y(ne,m), ymod(0:ne,0:1)
      type(wae_i_eqnvars) :: kv, km1, avg
      real*8 :: delta_v, delta_r, dvdr_avg
      call set_vars(k, x, y, kv, km1, avg, ymod, m, ne)
      if (k .eq. 2) call set_erf_norm_k2(kv)
      delta_v = kv%v - km1%v
      delta_r = avg%z*(kv%q-km1%q)
      call wae_get_dvdr(dvdr_avg, avg)
      v_eqn = delta_v - dvdr_avg*delta_r
      end function v_eqn

      real*8 function rho_eqn(k, x, y, ymod, m, ne)
      integer, intent(in) :: k, m, ne
      real*8,  intent(in) :: x(m), y(ne,m), ymod(0:ne,0:1)
      type(wae_i_eqnvars) :: kv, km1, avg
      real*8 :: delta_rho, delta_r, drhodr_avg, dvdr_avg
      call set_vars(k, x, y, kv, km1, avg, ymod, m, ne)
      if (k .eq. 2) call set_erf_norm_k2(kv)
      delta_rho = kv%rho - km1%rho
      delta_r = avg%z*(kv%q-km1%q)
      call wae_get_dvdr(dvdr_avg, avg)
      call wae_get_drhodr(drhodr_avg, avg, dvdr_avg)
      rho_eqn = delta_rho - drhodr_avg*delta_r
      end function rho_eqn

      real*8 function T_eqn(k, x, y, ymod, m, ne)
      integer, intent(in) :: k, m, ne
      real*8,  intent(in) :: x(m), y(ne,m), ymod(0:ne,0:1)
      type(wae_i_eqnvars) :: kv, km1, avg
      real*8 :: delta_T, delta_r, dTdr_avg, dYsdr_avg(nsp)
      real*8 :: drhodr_avg, dvdr_avg
      call set_vars(k, x, y, kv, km1, avg, ymod, m, ne)
      if (k .eq. 2) call set_erf_norm_k2(kv)
      delta_T = kv%T - km1%T
      delta_r = avg%z*(kv%q-km1%q)
      call wae_get_dvdr(dvdr_avg, avg)
      call wae_get_drhodr(drhodr_avg, avg, dvdr_avg)
      call wae_get_dYsdr(dYsdr_avg, avg, dble(k))
      call wae_get_dTdr(dTdr_avg, avg, drhodr_avg, dYsdr_avg, dble(k))
      T_eqn = delta_T - dTdr_avg*delta_r
      end function T_eqn

      end module wae_eqns
