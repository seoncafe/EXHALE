program audit_probe
  use global_parameters
  use Numerical_Fluxes
  use S_estimate_ROE
  use caloric_eos
  use h3p_cooling
  use mol_rates
  use molecular_infrared_cooling
  use h2_photo_channels
  use Cross_sections, only: sigma_H2
  use, intrinsic :: ieee_arithmetic
  implicit none
  real*8 :: wl(3),wr(3),nf(3),pf,cs1,cs2,vs,va,ca,cb
  real*8 :: tk,ur,cv,uth,eps,sg(4),threshold_chemical,scale
  real*8 :: temps(6), joins(3), rel, emit
  real*8 :: u_before(3),u_after(3)
  integer :: i
  T0=1000.d0
  use_plm=.true.
  caloric_mixture_active=.false.
  wl=(/1.d0,0.d0,1.d0/)
  call W_to_U_comp(wl,u_before,1)
  wl(3)=1.5d0
  call W_to_U_comp(wl,u_after,1)
  write(*,'(a,3es22.12)') 'Fixed-T H ionization projection: E_before, E_after, delta: ', &
    u_before(3),u_after(3),u_after(3)-u_before(3)
  wl=(/1.d0,-2.d0,0.6d0/)
  wr=(/2.d0,-2.d0,1.2d0/)
  call lax_friedrichs_flux(wl,wr,nf,pf,1,2)
  write(*,'(a,3es22.12)') 'LLF mass flux, inferred alpha, required alpha: ', &
    nf(1),(-6.d0-2.d0*nf(1)),3.d0
  wl=(/1.d0,-1.d0,1.d0/)
  wr=(/1.d0,1.d0,0.1d0/)
  call speed_estimate_ROE(wl,wr,vs,cs1,cs2,gamma_ad)
  scale=10.d0
  wl(1)=wl(1)*scale; wl(3)=wl(3)*scale
  wr(1)=wr(1)*scale; wr(3)=wr(3)*scale
  call speed_estimate_ROE(wl,wr,va,ca,cb,gamma_ad)
  write(*,'(a,4es22.12)') 'ROE rarefaction sound speeds, original / scaled: ',cs1,cs2,ca,cb
  wl=(/2.d0,5.d0,2.d0/)
  wr=(/8.d0,-5.d0,0.2d0/)
  call speed_estimate_ROE(wl,wr,vs,cs1,cs2,gamma_ad)
  wl(1)=wl(1)*scale; wl(3)=wl(3)*scale
  wr(1)=wr(1)*scale; wr(3)=wr(3)*scale
  call speed_estimate_ROE(wl,wr,va,ca,cb,gamma_ad)
  write(*,'(a,4es22.12)') 'ROE shock sound speeds, original / scaled: ',cs1,cs2,ca,cb
  wl=(/1.d0,-10.d0,1.d0/)
  wr=(/1.d0,10.d0,1.d0/)
  call speed_estimate_ROE(wl,wr,vs,cs1,cs2,gamma_ad)
  write(*,'(a,2l3)') 'ROE vacuum sound speeds finite: ',ieee_is_finite(cs1),ieee_is_finite(cs2)
  write(*,'(a,3es22.12)') 'H3+ cooling T=1000, nH3+=1, nH2=0,1,1e6: ', &
    h3p_cooling_rate(1000.d0,1.d0,0.d0),h3p_cooling_rate(1000.d0,1.d0,1.d0), &
    h3p_cooling_rate(1000.d0,1.d0,1.d6)
  joins=(/300.d0,800.d0,1800.d0/)
  do i=1,3
    tk=joins(i)
    rel=h3p_emission_lte(tk+1.d-7)/h3p_emission_lte(tk-1.d-7)-1.d0
    write(*,'(a,f9.2,es22.12)') 'H3+ LTE join T, relative jump: ',tk,rel
  enddo
  temps=(/300.d0,1000.d0,2000.d0,4000.d0,8000.d0,15000.d0/)
  eps=1.d-5
  do i=1,6
    tk=temps(i)
    call h2_rovibrational_energy_and_heat_capacity(tk,ur,cv)
    uth=tk*(log(q_rovib_H2(tk*exp(eps)))-log(q_rovib_H2(tk*exp(-eps))))/(2.d0*eps)
    write(*,'(a,f9.2,3es22.12)') 'H2 T, EOS u_rv/k, chemistry u_rv/k, relative difference: ', &
      tk,ur,uth,(uth-ur)/max(ur,1.d-99)
  enddo
  threshold_chemical=2.d0*13.598434599d0+h2_dissociation_energy_eV()
  write(*,'(a,3es22.12)') 'H2 double threshold: code, product energy, unassigned eV: ', &
    e_th_H2_dd,threshold_chemical,e_th_H2_dd-threshold_chemical
  h2_double_ionization_model='chung80'
  call h2_channel_cross_sections(80.d0,sigma_H2(80.d0),sg)
  write(*,'(a,es22.12)') 'H2 double event fraction at 80 eV: ',sg(ICH_D)/sum(sg)
  call molecular_infrared_init(700.d0)
  emit=h2_line_emission_lte(700.d0)
  write(*,'(a,es22.12)') 'H2 net/emission at Tgas=Trad=700,W=1: ', &
    h2_line_net_cooling_rate(700.d0,1.d0,1.d0)/emit
  emit=h2o_band_emission_lte(700.d0)
  write(*,'(a,es22.12)') 'H2O net/emission at Tgas=Trad=700,W=1: ', &
    h2o_band_net_cooling_rate(700.d0,1.d0,1.d0)/emit
  emit=co_band_emission_lte(700.d0)
  write(*,'(a,es22.12)') 'CO net/emission at Tgas=Trad=700,W=1: ', &
    co_band_net_cooling_rate(700.d0,1.d0,1.d0)/emit
end program audit_probe
