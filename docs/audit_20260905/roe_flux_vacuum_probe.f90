program roe_flux_vacuum_probe
  use global_parameters, only: gamma_ad,T0,use_plm,flux
  use Numerical_Fluxes, only: Num_flux,Phys_flux
  use Conversion, only: W_to_U_comp
  use caloric_eos, only: caloric_mixture_active
  use S_estimate_ROE, only: speed_estimate_ROE
  use, intrinsic :: ieee_arithmetic
  implicit none
  real*8 :: wl(3),wr(3),nf(3),fl(3),u(3),unew(3),pf,vs,cl,cr
  real*8 :: speed,tail_l,tail_r,dt_dx,p_exact
  integer :: k
  T0=1000.d0
  use_plm=.true.
  caloric_mixture_active=.false.
  flux='ROE'
  do k=1,2
    speed=1.d0
    if(k==2) speed=10.d0
    wl=[1.d0,-speed,1.d0]
    wr=[1.d0, speed,1.d0]
    call speed_estimate_ROE(wl,wr,vs,cl,cr,gamma_ad)
    call Num_flux(wl,wr,nf,pf,1,2)
    write(*,'(a,f6.1)') 'Symmetric separating speed: ',speed
    write(*,'(a,2l3)') 'Auxiliary sound speeds finite: ', &
      ieee_is_finite(cl),ieee_is_finite(cr)
    write(*,'(a,3es23.14)') 'Production ROE flux: ',nf
    if(k==1) then
      p_exact=(1.d0-(gamma_ad-1.d0)*speed/(2.d0*sqrt(gamma_ad))) &
        **(2.d0*gamma_ad/(gamma_ad-1.d0))
      write(*,'(a,es23.14)') 'Exact symmetric nonvacuum interface pressure: ',p_exact
    else
      tail_l=-speed+2.d0*sqrt(gamma_ad)/(gamma_ad-1.d0)
      tail_r= speed-2.d0*sqrt(gamma_ad)/(gamma_ad-1.d0)
      write(*,'(a,2es23.14)') 'Exact vacuum edges: ',tail_l,tail_r
      write(*,'(a)') 'Exact interface flux inside this vacuum: 0, 0, 0'
      call W_to_U_comp(wl,u,1)
      call Phys_flux(wl,fl,1)
      dt_dx=0.01d0
      unew=u-dt_dx*(nf-fl)
      write(*,'(a,es23.14)') 'Courant number: ',dt_dx*(speed+sqrt(gamma_ad))
      write(*,'(a,3es23.14)') 'One Euler update, rho,momentum,energy: ',unew
      write(*,'(a,es23.14)') 'Returned internal energy density: ', &
        unew(3)-0.5d0*unew(2)**2/unew(1)
    endif
  enddo
end program roe_flux_vacuum_probe
