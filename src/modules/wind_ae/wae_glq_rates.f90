      module wae_glq_rates
      ! Multifrequency ionization + heating rates, ported verbatim from
      ! wind-ae glq_rates.c : calc_gql_rates() (Shull & van Steenberg 1985
      ! energy partition; Dere 2007 secondary-ionization coefficients).
      ! Faithful to the C arithmetic and evaluation order -- behavior
      ! changes (e.g. Stage-2 spectrum binning) belong in a later phase.
      !
      ! Rate layout matches ionization_rate[]: for 0-based species j,
      ! ion(j*(nsp+1)+1)        = primary ionization rate,
      ! ion(j*(nsp+1)+1+m+1)    = secondary rate induced in species m.
      use wae_config,   only: nsp => wae_nspecies, NCOL0 => wae_NCOL0,    &
                              xray_lo => wae_xray_lo
      use wae_spectrum, only: npts => wae_npts, hc => wae_hc_over_wl,     &
                              wPhi => wae_wPhi_wl, sigma_wl => wae_sigma_wl,&
                              ion_pot => wae_ion_pot, HX => wae_HX,        &
                              amass => wae_atomic_mass, Ftot => wae_Ftot
      use wae_rate_coeffs, only: Rtab => wae_R
      implicit none

      ! glq rate cache (mirrors the static last_N / ionization_rate /
      ! heating_rate of calc_gql_rates).
      logical :: wae_glq_have_cache = .false.
      real*8, allocatable :: wae_last_N(:), wae_ion_cache(:), wae_heat_cache(:)
      contains

      subroutine wae_glq_rates_eval(N, Ys, ion, heat)
      ! N(:)  column densities in NCOL0 units (windsoln Ncol columns)
      ! Ys(:) neutral fractions
      ! ion(:) length nsp*(nsp+1); heat(:) length nsp
      real*8, intent(in)  :: N(nsp), Ys(nsp)
      real*8, intent(out) :: ion(nsp*(nsp+1)), heat(nsp)

      integer :: i, j, m, jj
      real*8  :: tau, denom, sigma, f, E_0, eta_m
      real*8  :: background, frac_in_excite, frac_in_heat, frac_in_ion_tot
      real*8  :: frac_in_ion_m, primary_ion_rate_j, R_tot, n0_m
      real*8  :: n_tot, n_ion_tot, n_j, bg
      real*8  :: tau_arr(npts), f_denom(npts), Phi(npts)
      logical :: recalc

      ! --- cache, faithful to calc_gql_rates(): the prior rates are reused
      ! whenever the column densities N are unchanged. NOTE this ignores
      ! any change in Ys (a latent C quirk that matters when get_*_derivs
      ! perturbs Ys at fixed Ncol); preserved here for C-comparability. ---
      recalc = (.not. wae_glq_have_cache)
      if (wae_glq_have_cache) then
         do j = 1, nsp
            if (N(j) .ne. wae_last_N(j)) then
               recalc = .true.
               exit
            end if
         end do
      end if
      if (.not. recalc) then
         ion  = wae_ion_cache
         heat = wae_heat_cache
         return
      end if

      ion  = 0.0d0
      heat = 0.0d0

      ! Total and ionized number densities (per rho; the rho factor cancels)
      n_tot = 0.0d0
      n_ion_tot = 0.0d0
      do j = 1, nsp
         n_j = HX(j)/amass(j)
         n_tot = n_tot + n_j
         n_ion_tot = n_ion_tot + n_j*(1.0d0 - Ys(j))
      end do

      ! Optical depth, on-the-spot opacity denominator, attenuated flux
      do i = 1, npts
         tau = 0.0d0
         denom = 0.0d0
         do j = 1, nsp
            tau   = tau   + sigma_wl(i,j)*N(j)*NCOL0
            denom = denom + sigma_wl(i,j)*HX(j)*Ys(j)/amass(j)
         end do
         tau_arr(i) = tau
         f_denom(i) = denom
         Phi(i) = wPhi(i)*exp(-tau)
      end do

      ! Per-species ionization (primary + secondary) and heating
      do j = 1, nsp
         do i = 1, npts
            background = n_ion_tot/n_tot
            frac_in_excite = 0.4766d0 *                                   &
                 (1.0d0 - background**0.2735d0)**1.5221d0
            frac_in_heat = 0.9971d0 *                                     &
                 (1.0d0 - (1.0d0 - background**0.2663d0)**1.3163d0)
            frac_in_ion_tot = 1.0d0 - frac_in_heat - frac_in_excite
            sigma = sigma_wl(i,j)
            f = sigma*(HX(j)*Ys(j)/amass(j))/f_denom(i)
            E_0 = hc(i) - ion_pot(j)
            if (sigma .eq. 0.0d0) exit          ! C: break

            primary_ion_rate_j = Ftot*sigma*f*Phi(i)

            if (E_0 .gt. xray_lo) then
               bg = background
               if (bg .lt. 0.0d0) bg = 1.0d-10  ! (clamp kept; unused below)
               do m = 1, nsp
                  R_tot = 0.0d0
                  do jj = 1, nsp
                     R_tot = R_tot +                                      &
                        Rtab(nsp*(jj-1) + (m-1) + 1, i) *                 &
                        Ys(jj)*HX(jj)/amass(jj)
                  end do
                  n0_m = Ys(m)*HX(m)/amass(m)
                  frac_in_ion_m = Rtab(nsp*(j-1) + (m-1) + 1, i)*n0_m/R_tot
                  eta_m = E_0*(frac_in_ion_tot*frac_in_ion_m)/ion_pot(m)
                  ! secondary rate into species m from ionizing species j
                  ion((j-1)*(nsp+1) + 1 + m) =                            &
                     ion((j-1)*(nsp+1) + 1 + m) + eta_m*primary_ion_rate_j
               end do
            else
               frac_in_heat = 1.0d0
            end if

            ion((j-1)*(nsp+1) + 1) = ion((j-1)*(nsp+1) + 1)              &
                                     + primary_ion_rate_j
            heat(j) = heat(j) + frac_in_heat*E_0*sigma*f*Phi(i)*Ftot
         end do
      end do

      ! save state for the next cache check (as calc_gql_rates does)
      if (.not. allocated(wae_last_N)) then
         allocate(wae_last_N(nsp), wae_ion_cache(nsp*(nsp+1)),            &
                  wae_heat_cache(nsp))
      end if
      wae_last_N = N
      wae_ion_cache = ion
      wae_heat_cache = heat
      wae_glq_have_cache = .true.
      end subroutine wae_glq_rates_eval

      end module wae_glq_rates
