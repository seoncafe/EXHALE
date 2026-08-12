   module excited_hydrogen
   ! In-code excited hydrogen H(n=2).
   !
   ! Solves the 2s/2p rate-equilibrium of Christie, Arras & Li (2013, ApJ
   ! 772, 144; their Eqs. 12-13) -- collisional 1s<->2s/2p, 2s<->2p l-mixing,
   ! recombination cascade, two-photon decay, and Ly-alpha radiative pumping --
   ! with the Ly-alpha mean intensity J_lya either (a) estimated as in Huang
   ! et al. (2017, ApJ 851, 150) Eq. (6), J_lya ~ 0.1 F_LyC/Dnu_D, attenuated
   ! by 1/(1+tau_lya) below the Ly-alpha photosphere, or (b) read from an
   ! external Ly-alpha RT profile (selected by jlya_mode).
   ! The resulting H(n=2) population is then (i) photoionized by the stellar
   ! Balmer continuum (E > 3.4 eV), adding a proton source to the H ionization
   ! balance, and (ii) heated by the photoelectron excess energy (photoelectric
   ! heating) and by collisional de-excitation of the Ly-alpha-pumped n=2 atoms
   ! (incl_deexc_heat). This is the Fortran counterpart of the n=2 model in
   ! EXHALE_transit.py (which uses it for H-alpha/H-beta line opacity).
   !
   ! All feedback is gated by use_excited_H; when off, the module is inert and
   ! the global feedback arrays stay zero, so the build reproduces the no-excited-H result.
   !
   ! Coupling is DECOUPLED (lagged-explicit): excited_H_update is called once
   ! per timestep from EXHALE_main BEFORE the ionization/energy solve, fills the
   ! module-level feedback/diagnostic arrays in global_parameters from the
   ! current state, and those frozen arrays are read inside ioniz_eq -- never
   ! inside its Newton iteration. The single relaxation then converges hydro,
   ! ionization and the Balmer feedback together.

   use global_parameters
   use species_table, only: n_mion, mion_fsp
   use utils, only: calc_ne
   use lya_rt, only: jlya_escape_prob, jint_arr, jstar_arr,                  &
                     lya_line_center_optical_depth
   ! n=2 / Ly-alpha atomic data and collisional rate coefficients
   ! (Christie+2013 Table 2, Draine 2011); one definition, shared with lya_rt.
   use hydrogen_n2_rates

   implicit none

   ! ----- Balmer-continuum (n=2 photoionization) data ----- !
   real*8, parameter :: eV2Hz    = 2.417989242d14      ! Hz per eV
   real*8, parameter :: E2_thr   = 3.40d0              ! n=2 ionization edge [eV]
   real*8, parameter :: nu2_thr  = E2_thr*eV2Hz        ! [s^-1]
   real*8, parameter :: s2_thr   = 1.4d-17             ! sigma_2 at threshold [cm^2]

   ! ----- Ly-alpha pumping (parameterized J_lya) data ----- !
   real*8, parameter :: sigma_LyC = 6.3d-18            ! H photoion. xsec at LyC [cm^2]

   ! RT-supplied J_lya(r) profile (jlya_mode=1), interpolated onto the grid once.
   logical, save :: jlya_rt_loaded = .false.
   real*8, dimension(1-Ng:N+Ng) :: jlya_rt_grid = 0.0d0

   ! Saved context for the diagnostic dump (filled in excited_H_update).
   real*8, dimension(1-Ng:N+Ng) :: Tdiag = 0.0d0, nhidiag = 0.0d0, nediag = 0.0d0
   real*8, dimension(1-Ng:N+Ng) :: nhiidiag = 0.0d0   ! nHII for the recomb. sink
   real*8, dimension(1-Ng:N+Ng) :: taulya = 0.0d0      ! top-down Ly-alpha optical depth

   contains

   ! --------------------------------------------------------------- !

   subroutine excited_H_update(T_in, n_in, f_sp_in, v_in, rel_change)
   ! Fill the global excited-H feedback (gph_balmer_HI, heat_balmer) and
   ! diagnostic arrays from a converged equilibrium state. Inputs are the
   ! dimensionless ATES profiles (same convention as ioniz_eq): T_in*T0 = T[K],
   ! n_in*n0 = total number density [cm^-3], f_sp_in = ion fractions.
   ! rel_change returns the max relative change in heat_balmer vs the previous
   ! call (used by the outer iteration to test convergence).

   real*8, dimension(1-Ng:N+Ng), intent(in) :: T_in, n_in, v_in
   real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_in
   real*8, intent(out) :: rel_change

   integer :: j, im
   real*8, dimension(1-Ng:N+Ng) :: T_K, n_dim, nhi, nhii
   real*8, dimension(1-Ng:N+Ng) :: nhei, nheii, nheiii, ne
   real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
   real*8, dimension(1-Ng:N+Ng) :: heat_prev
   real*8 :: F_LyC, F_inc, xi, abs_frac, N_HI_tot, a_cm
   real*8 :: Dnu_D, n2s, n2p, n2tot, relc

   ! Remember the previous heating for the outer-iteration convergence test.
   heat_prev = heat_balmer

   ! Dimensionalize exactly as ioniz_eq does.
   T_K   = T_in*T0
   n_dim = n_in*n0
   nhi   = f_sp_in(:,1)*n_dim
   nhii  = f_sp_in(:,2)*n_dim
   if (thereis_He) then
      nhei   = f_sp_in(:,3)*n_dim
      nheii  = f_sp_in(:,4)*n_dim
      nheiii = f_sp_in(:,5)*n_dim
   else
      nhei = 0.0d0; nheii = 0.0d0; nheiii = 0.0d0
   endif
   do im = 1,n_mion
      nm(:,im) = f_sp_in(:,mion_fsp(im))*n_dim
   enddo
   ! Molecular ions are deliberately omitted as trace electron donors
   ! (negligible in the hot, atomic n=2 layer this module models).
   call calc_ne(nhii, nheii, nheiii, ne, nm)

   ! ----- Day-night / 2D dilution factor xi ----- !
   ! Same factor the code applies to ground-state EUV photoionization (the ATES
   ! "2D approximate method"): Rate/4 -> 1/4, Rate/2 -> 1/2, else full.
   if      (index(appx_mth,'Rate/4') .gt. 0) then
      xi = 0.25d0
   else if (index(appx_mth,'Rate/2') .gt. 0) then
      xi = 0.5d0
   else
      xi = 1.0d0
   endif

   ! ----- Scalar Balmer-continuum rates (depend only on T_star, R_star/a) ----- !
   ! Computed once per update; reduced by xi for the dayside hemisphere average,
   ! consistent with the ground-state EUV ionization treatment.
   gamma2_bal = xi*gamma_n2_balmer(T_star_eff, R_star/max(a_orb,1.0d-30))
   hpe2_bal   = xi*heat_n2_balmer(T_star_eff,  R_star/max(a_orb,1.0d-30))

   ! ----- Deposited Ly-continuum flux F_LyC (Huang+2017 Eq. 6 input) ----- !
   ! Single deposited flux from the total neutral-H column, matching
   ! EXHALE_transit.py:
   ! each absorbed LyC photon balanced by a recombination -> Ly-alpha photon.
   a_cm  = a_orb                                   ! a_orb already in cm
   F_inc = 10.0d0**LEUV/(4.0d0*pi*a_cm**2.0)       ! incident stellar LyC [erg cm^-2 s^-1]
   ! Vertical neutral-H column [cm^-2] by trapezoid over the physical grid.
   N_HI_tot = 0.0d0
   do j = 1, N
      N_HI_tot = N_HI_tot                                                  &
               + 0.5d0*(nhi(j) + nhi(j+1))*(r(j+1) - r(j))*R0
   enddo
   abs_frac = 1.0d0 - exp(-sigma_LyC*N_HI_tot)
   F_LyC    = xi*F_inc*abs_frac

   ! ----- Ly-alpha mean intensity J_lya(r) ----- !
   if (jlya_mode .eq. 1) then
      ! (b) Externally-computed Ly-alpha RT profile (already physical); loaded
      ! once and reused. No depth attenuation -- the RT carries it. tau_Lya is
      ! still evaluated, as a diagnostic only: it says where the imported field
      ! is optically thick, and nothing in this mode consumes it.
      if (.not. jlya_rt_loaded) call load_jlya_rt()
      call lya_line_center_optical_depth(T_K, nhi, taulya)
      do j = 1-Ng, N+Ng
         Jlya_arr(j) = jlya_rt_grid(j)
      enddo
   else if (jlya_mode .eq. 2) then
      ! (c) In-line escape-probability RT (Neufeld/Harrington wing escape),
      ! evaluated from the current state every timestep. Fills taulya too.
      call jlya_escape_prob(T_K, nhi, nhii, ne, v_in, Jlya_arr, taulya)
   else
      ! (a) Parameterized J_lya = 0.1 F_LyC/Dnu_D (Huang+2017 Eq. 6), attenuated
      ! by 1/(1+tau_lya) with tau_lya the top-down line-center Ly-alpha optical
      ! depth, so the pumping vanishes below the Ly-alpha photosphere. This is
      ! an approximate staging estimate; the accurate field is jlya_mode=1.
      call lya_line_center_optical_depth(T_K, nhi, taulya)
      do j = 1-Ng, N+Ng
         Dnu_D = nu_lya*sqrt(2.0d0*kb_erg*max(T_K(j),1.0d0)/mu)/c_light
         Jlya_arr(j) = 0.1d0*F_LyC/max(Dnu_D,1.0d-30)/(1.0d0 + taulya(j))
      enddo
   endif

   ! ----- Cell-by-cell n=2 populations + feedback ----- !
   do j = 1-Ng, N+Ng
      call n2_populations(T_K(j), max(nhi(j),0.0d0), max(nhii(j),0.0d0),    &
                          max(ne(j),0.0d0), Jlya_arr(j),                    &
                          gamma2_bal, gamma2_bal, n2s, n2p)
      n2tot = n2s + n2p

      ! Diagnostics
      n2s_arr(j)  = n2s
      n2p_arr(j)  = n2p
      Tdiag(j)    = T_K(j)
      nhidiag(j)  = nhi(j)
      nhiidiag(j) = nhii(j)
      nediag(j)   = ne(j)

      ! (i) Balmer photoionization of H(n=2): proton source [cm^-3 s^-1].
      Sproton_arr(j) = gamma2_bal*n2tot
      ! Folded into the H balance as an EFFECTIVE extra HI photoion. rate
      ! [s^-1] so that nhi*gph reproduces the volumetric source (ioniz_eq
      ! adds this to P_HI; see its injection point).
      gph_balmer_HI(j) = gamma2_bal*n2tot/max(nhi(j),1.0d-30)

      ! (ii) Photoelectric heating: photoelectron excess energy [erg cm^-3 s^-1].
      Hpe_arr(j) = n2tot*hpe2_bal
      ! (ii') Optional collisional de-excitation heating of the n=2 atoms.
      if (incl_deexc_heat) then
         Hdx_arr(j) = ne(j)*E21_erg                                        &
                    *( c2s1s_rate(T_K(j))*n2s + c2p1s_rate(T_K(j))*n2p )
      else
         Hdx_arr(j) = 0.0d0
      endif
      heat_balmer(j) = Hpe_arr(j) + Hdx_arr(j)
   enddo

   ! Outer-iteration convergence metric: max relative change in heating
   ! over the physical domain (ignore cells with negligible heating).
   rel_change = 0.0d0
   do j = 1, N
      if (heat_balmer(j) .gt. 1.0d-30) then
         relc = abs(heat_balmer(j) - heat_prev(j))/heat_balmer(j)
         if (relc .gt. rel_change) rel_change = relc
      endif
   enddo

   end subroutine excited_H_update

   ! --------------------------------------------------------------- !

   subroutine n2_populations(T, n1s, nHII_l, ne_l, Jlya,                     &
                             gam_ion_2s, gam_ion_2p, n2s, n2p)
   ! Christie+2013 Eqs. 12-13: solve the 2x2 2s/2p rate equilibrium for the
   ! H(n=2) populations [cm^-3]. T in K, densities in cm^-3, Jlya in cgs
   ! (erg s^-1 cm^-2 Hz^-1 sr^-1). gam_ion_2s/gam_ion_2p = photoionization
   ! rates [s^-1] out of 2s and 2p (the stellar Balmer continuum).
   !
   ! The rate coefficients come from hydrogen_n2_rates so that the statistical
   ! weights g1s/g2s/g2p are evaluated once, at module scope, and cannot be
   ! shadowed by a dummy argument of this routine (Fortran is case-insensitive,
   ! and the dummies used to be called G2s/G2p).

   real*8, intent(in)  :: T, n1s, nHII_l, ne_l, Jlya, gam_ion_2s, gam_ion_2p
   real*8, intent(out) :: n2s, n2p

   real*8 :: Tl, a2s, a2p
   real*8 :: C1s2s, C1s2p, C2s2p, C2s1s, C2p1s, C2p2s
   real*8 :: Ppump, Pstim, L2p, L2s, S2p, S2s, M12, M21, det

   Tl = max(T, 1.0d0)

   ! Level-resolved recombination (Draine 2011; Table 2 R2,R8,R9).
   a2s = alpha_2s_hydrogen(Tl)
   a2p = alpha_2p_hydrogen(Tl)

   ! Collisional excitation 1s->2s, 1s->2p, 2s<->2p l-mixing and the reverse
   ! 2s/2p->1s de-excitation (R3,R4,R5 and their detailed-balance partners).
   C1s2s = c1s2s_rate(Tl)
   C1s2p = c1s2p_rate(Tl)
   C2s2p = c2s2p_rate(Tl)
   C2s1s = c2s1s_rate(Tl)
   C2p1s = c2p1s_rate(Tl)
   C2p2s = c2p2s_rate(Tl)

   ! Ly-alpha pump (1s->2p) and stimulated emission (2p->1s).
   Ppump = B12_lya*Jlya
   Pstim = B21_lya*Jlya

   ! 2x2 system [[L2p,-M12],[-M21,L2s]] [n2p,n2s]^T = [S2p,S2s]^T.
   ! The cascade source is alpha_2l * ne * nHII: the recombining partner of an
   ! electron is a proton, not another electron. The two differ wherever the
   ! electrons come from helium and metals while hydrogen is still neutral,
   ! which is the case through the base.
   L2p = A_2p1s + Pstim + (C2p1s + C2p2s)*ne_l + gam_ion_2p
   L2s = (C2s1s + C2s2p)*ne_l + gam_ion_2s + A_2s1s
   S2p = (Ppump + C1s2p*ne_l)*n1s + a2p*ne_l*nHII_l
   S2s = (C1s2s*ne_l)*n1s + a2s*ne_l*nHII_l
   M12 = C2s2p*ne_l
   M21 = C2p2s*ne_l
   det = L2p*L2s - M12*M21
   if (abs(det) .le. 0.0d0) det = 1.0d0

   n2p = max((S2p*L2s + M12*S2s)/det, 0.0d0)
   n2s = max((L2p*S2s + M21*S2p)/det, 0.0d0)

   end subroutine n2_populations

   ! --------------------------------------------------------------- !

   real*8 function gamma_n2_balmer(Tstar, R_over_a)
   ! n=2 photoionization rate [s^-1] from a diluted stellar blackbody Balmer
   ! continuum (E > 3.4 eV). Hydrogenic sigma_2(nu) = s2_thr*(nu2/nu)^3, flux
   ! F_nu = pi B_nu(Tstar) (R_star/a)^2. Trapezoid over 3.4-13.6 eV (above
   ! 13.6 eV the BB flux is negligible and ground-state H absorbs it).

   real*8, intent(in) :: Tstar, R_over_a
   ! n_nu, not ng: Ng is the global ghost-cell count and would be shadowed.
   integer, parameter :: n_nu = 400
   integer :: i
   real*8 :: Eg, nu, dnu, Bnu, Fnu, sig2, integ, nu_a, nu_b, x

   gamma_n2_balmer = 0.0d0
   if (Tstar .le. 0.0d0) return

   nu_a  = E2_thr  *eV2Hz
   nu_b  = e_th_HI *eV2Hz
   dnu   = (nu_b - nu_a)/dble(n_nu-1)
   integ = 0.0d0
   do i = 1, n_nu
      nu   = nu_a + dble(i-1)*dnu
      x    = hp_erg*nu/(kb_erg*Tstar)
      Bnu  = (2.0d0*hp_erg*nu**3.0/c_light**2.0)/(exp(x) - 1.0d0)
      Fnu  = pi*Bnu*R_over_a**2.0
      sig2 = s2_thr*(nu2_thr/nu)**3.0
      ! Integrand F_nu/(h nu) * sigma_2  ; trapezoid weights (endpoints 1/2).
      Eg   = Fnu/(hp_erg*nu)*sig2
      if (i .eq. 1 .or. i .eq. n_nu) Eg = 0.5d0*Eg
      integ = integ + Eg
   enddo
   gamma_n2_balmer = integ*dnu

   end function gamma_n2_balmer

   ! --------------------------------------------------------------- !

   real*8 function heat_n2_balmer(Tstar, R_over_a)
   ! Photoelectric heating per H(n=2) atom [erg s^-1]: the rate integrand of
   ! gamma_n2_balmer weighted by the photoelectron excess energy (h nu - 3.4 eV).

   real*8, intent(in) :: Tstar, R_over_a
   ! n_nu, not ng: Ng is the global ghost-cell count and would be shadowed.
   integer, parameter :: n_nu = 400
   integer :: i
   real*8 :: nu, dnu, Bnu, Fnu, sig2, integ, nu_a, nu_b, x, w, term

   heat_n2_balmer = 0.0d0
   if (Tstar .le. 0.0d0) return

   nu_a  = E2_thr  *eV2Hz
   nu_b  = e_th_HI *eV2Hz
   dnu   = (nu_b - nu_a)/dble(n_nu-1)
   integ = 0.0d0
   do i = 1, n_nu
      nu   = nu_a + dble(i-1)*dnu
      x    = hp_erg*nu/(kb_erg*Tstar)
      Bnu  = (2.0d0*hp_erg*nu**3.0/c_light**2.0)/(exp(x) - 1.0d0)
      Fnu  = pi*Bnu*R_over_a**2.0
      sig2 = s2_thr*(nu2_thr/nu)**3.0
      w    = hp_erg*(nu - nu2_thr)                  ! photoelectron excess [erg]
      term = Fnu/(hp_erg*nu)*sig2*w
      if (i .eq. 1 .or. i .eq. n_nu) term = 0.5d0*term
      integ = integ + term
   enddo
   heat_n2_balmer = integ*dnu

   end function heat_n2_balmer

   ! --------------------------------------------------------------- !

   subroutine load_jlya_rt
   ! Read an externally-computed Ly-alpha mean-intensity profile J_lya(r) from
   ! jlya_rt_file (two columns: r/Rp, J_lya [cgs]; '#'/blank lines skipped) and
   ! linearly interpolate onto the simulation grid (clamped at the endpoints).
   ! Cached after the first call (jlya_rt_loaded).
   integer, parameter :: mx = 200000
   integer :: io, nrt, k, j
   real*8  :: rr, jj
   real*8, allocatable :: r_rt(:), j_rt(:)
   character(len=300) :: line

   allocate(r_rt(mx), j_rt(mx))
   nrt = 0
   open(unit=73, file=trim(jlya_rt_file), status='old', action='read', iostat=io)
   if (io .ne. 0) then
      write(*,*) '(excited_hydrogen) ERROR: cannot open Jlya RT file ',      &
                 trim(jlya_rt_file), ' -- using J_lya = 0.'
      jlya_rt_grid   = 0.0d0
      jlya_rt_loaded = .true.
      deallocate(r_rt, j_rt)
      return
   endif
   do
      read(73,'(a)',iostat=io) line
      if (io .ne. 0) exit
      line = adjustl(line)
      if (len_trim(line) .eq. 0)  cycle
      if (line(1:1) .eq. '#')     cycle
      read(line,*,iostat=io) rr, jj
      if (io .ne. 0) cycle
      nrt = nrt + 1
      r_rt(nrt) = rr
      j_rt(nrt) = jj
   enddo
   close(73)

   if (nrt .lt. 2) then
      write(*,*) '(excited_hydrogen) ERROR: Jlya RT file ',                  &
                 trim(jlya_rt_file), ' has < 2 points -- using J_lya = 0.'
      jlya_rt_grid   = 0.0d0
      jlya_rt_loaded = .true.
      deallocate(r_rt, j_rt)
      return
   endif

   do j = 1-Ng, N+Ng
      if (r(j) .le. r_rt(1)) then
         jlya_rt_grid(j) = j_rt(1)
      else if (r(j) .ge. r_rt(nrt)) then
         jlya_rt_grid(j) = j_rt(nrt)
      else
         do k = 1, nrt-1
            if (r(j) .ge. r_rt(k) .and. r(j) .le. r_rt(k+1)) then
               jlya_rt_grid(j) = j_rt(k) + (j_rt(k+1) - j_rt(k))            &
                               *(r(j) - r_rt(k))/(r_rt(k+1) - r_rt(k))
               exit
            endif
         enddo
      endif
   enddo
   jlya_rt_loaded = .true.
   deallocate(r_rt, j_rt)
   write(*,'(a,i0,a,a)') ' (excited_hydrogen) loaded ', nrt,                 &
                         ' Jlya RT points from ', trim(jlya_rt_file)
   end subroutine load_jlya_rt

   ! --------------------------------------------------------------- !

   subroutine write_excited_H
   ! Dump the cell-by-cell excited-H diagnostic to output/Excited_H.txt for
   ! validation against Huang et al. (2023) Figs. 11/27 (proton source) and
   ! Figs. 10/26 (heating budget). All quantities cgs; written for the last
   ! converged outer pass.
   integer :: j
   open(unit = 72, file = './output/Excited_H.txt')
   write(72,'(a)') '# Excited hydrogen H(n=2) diagnostic.'
   write(72,'(a,es12.5)') '# Stellar T_eff [K]       = ', T_star_eff
   write(72,'(a,es12.5)') '# Stellar R_star [Rsun]   = ', R_star/Rsun
   write(72,'(a,es12.5)') '# Gamma_2 (Balmer) [s-1]  = ', gamma2_bal
   write(72,'(a,es12.5)') '# heat per n2 atom [erg/s]= ', hpe2_bal
   if (jlya_mode .eq. 1) then
      write(72,'(a,a)')   '# Jlya mode = 1 (RT profile from ', trim(jlya_rt_file)//')'
   else if (jlya_mode .eq. 2) then
      write(72,'(a)')     '# Jlya mode = 2 (in-line escape-probability RT, lya_rt.f90)'
   else
      write(72,'(a)')     '# Jlya mode = 0 (parameterized 0.1 F_LyC/Dnu_D / (1+tau_Lya))'
   endif
   write(72,'(a)') '# col1 r/Rp  col2 T[K]  col3 nHI  col4 ne  col5 Jlya'   &
                // '  col6 n2s  col7 n2p  col8 Sproton[cm-3 s-1]'           &
                // '  col9 Hpe[erg cm-3 s-1]  col10 Hdx[erg cm-3 s-1]'      &
                // '  col11 tau_Lya'                                        &
                // '  col12 Sgrnd_photoion[cm-3 s-1]'                       &
                // '  col13 Scoll_ion[cm-3 s-1]  col14 Srecomb[cm-3 s-1]'   &
                // '  col15 Jint[cgs]  col16 Jstar[cgs]'
   do j = 1-Ng, N+Ng
      write(72,*) r(j), Tdiag(j), nhidiag(j), nediag(j),                    &
                  Jlya_arr(j), n2s_arr(j), n2p_arr(j),                      &
                  Sproton_arr(j), Hpe_arr(j), Hdx_arr(j), taulya(j),        &
                  gph_ground_HI(j)*nhidiag(j),                              &
                  cion_HI(j)*nediag(j)*nhidiag(j),                          &
                  arec_HII(j)*nediag(j)*nhiidiag(j),                        &
                  jint_arr(j), jstar_arr(j)
   enddo
   close(72)
   end subroutine write_excited_H

   ! End of module
   end module excited_hydrogen
