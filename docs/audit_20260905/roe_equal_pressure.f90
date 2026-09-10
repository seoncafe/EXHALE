program roe_equal_pressure
  use global_parameters, only: gamma_ad
  use S_estimate_ROE, only: speed_estimate_ROE
  use h3p_cooling, only: h3p_nonlte_factor,h3p_emission_lte
  use, intrinsic :: ieee_arithmetic
  implicit none
  real*8 :: wl(3), wr(3), ustar, clstar, crstar, sound, p_exact
  wl = [1.d0, -1.d0, 1.d0]
  wr = [1.d0,  1.d0, 1.d0]
  call speed_estimate_ROE(wl,wr,ustar,clstar,crstar,gamma_ad)
  sound = sqrt(gamma_ad)
  p_exact = (1.d0-(gamma_ad-1.d0)*(wr(2)-wl(2))/(4.d0*sound)) &
            **(2.d0*gamma_ad/(gamma_ad-1.d0))
  write(*,'(a,3es23.14)') 'gamma, exact p_star, exact c_star: ', &
     gamma_ad,p_exact,sound-(gamma_ad-1.d0)/2.d0
  write(*,'(a,3es23.14)') 'routine u_star, cL_star, cR_star: ',ustar,clstar,crstar
  write(*,'(a,2l3)') 'finite returned sound speeds: ', &
     ieee_is_finite(clstar),ieee_is_finite(crstar)
  write(*,'(a,2es23.14)') 'H3+ s at 5000 K, just below and at 1e14: ', &
     h3p_nonlte_factor(5000.d0,1.d14*(1.d0-1.d-8)), &
     h3p_nonlte_factor(5000.d0,1.d14)
  write(*,'(a,es23.14)') 'H3+ E_LTE(5000), W molecule^-1 sr^-1: ', &
     h3p_emission_lte(5000.d0)
end program roe_equal_pressure
