   module excited_hydrogen
   ! In-code excited hydrogen H(n=2).
   !
   ! Solves the 2s/2p rate-equilibrium of Christie, Arras & Li (2013, ApJ
   ! 772, 144; their Eqs. 12-13) -- collisional 1s<->2s/2p, 2s<->2p l-mixing,
   ! recombination cascade, two-photon decay, Ly-alpha radiative pumping, and
   ! in a molecular run the dissociative recombinations H2+ + e and HeH+ + e,
   ! which leave one of their hydrogen fragments in n = 2 rather than in the
   ! ground state (the term is formed by molecular_reaction_heat, which
   ! subtracts the same excitation energy from the chemical heat) --
   ! with the Ly-alpha mean intensity J_lya either (a) estimated as in Huang
   ! et al. (2017, ApJ 851, 150) Eq. (6), J_lya ~ 0.1 F_LyC/Dnu_D, attenuated
   ! by 1/(1+tau_lya) below the Ly-alpha photosphere, or (b) read from an
   ! external Ly-alpha RT profile (selected by jlya_mode).  Every stellar
   ! Ly-alpha field a cell is pumped by is the MEAN of that field over the
   ! cell's own line-centre optical depth (lya_rt.f90).
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
   use species_table, only: n_mion, mion_fsp, isp_H2p, isp_HeHp
   use utils, only: calc_ne, write_row_layout_header, one_minus_exp_over_x
   ! The two dissociative recombinations of the molecular network that leave
   ! one hydrogen atom in n = 2, and the excitation energy the heat ledger
   ! subtracts for them. Reading the source from the module that owns the
   ! subtraction is what keeps the level production and the energy that left
   ! the gas the same event counted once.
   use molecular_reaction_heat, only:                                        &
                     dissociative_recombination_n2_source
   use lya_rt, only: jlya_escape_prob, jint_arr, jstar_arr,                  &
                     lya_line_center_optical_depth,                         &
                     lya_photosphere_attenuation_cell_mean
   ! n=2 / Ly-alpha atomic data and collisional rate coefficients (see the
   ! header of hydrogen_n2_rates for the sources); one definition, shared
   ! with lya_rt and with the H I cooling of Cool_coeff.
   use hydrogen_n2_rates
   ! The stellar field the Balmer continuum is photoionized by, the n = 2
   ! ionization threshold that heads it, and the frequency of a
   ! one-electronvolt photon: the single definitions of all three, shared
   ! with the photon grid (J_inc.f90). stellar_flux_eV is the field of
   ! the RUN'S spectrum type, so the Balmer band is built from the same
   ! spectrum as every other band, and e_th_HI_n2 is the same
   ! threshold that floors that grid whenever this coupling is armed
   ! (sed_read's photon_grid_floor_eV).
   use J_incident, only: stellar_flux_eV, spectrum_covers_eV, eV2Hz,        &
                         e_th_HI_n2

   implicit none

   ! ----- Balmer-continuum (n=2 photoionization) data ----- !
   ! sigma_2 at the n = 2 edge [cm^2], one value for 2s and 2p, taken to
   ! fall as (e_th_HI_n2/E)^3 above it (the Kramers frequency dependence).
   ! 1.4e-17 is Seaton (1959, MNRAS 119, 81) at the n = 2 edge: his eq. (3)
   ! Kramers cross section 1.5814e-17 times his eq. (10) Gaunt factor
   ! g_II(2, 0) = 0.8715 gives 1.378e-17 (DERIVED). The exact hydrogenic
   ! values are 1.478e-17 (2s) and 1.355e-17 (2p), statistical mean
   ! 1.386e-17 (DERIVED by integrating the bound-free dipole matrix
   ! elements). Osterbrock & Ferland (2006) give no n = 2 value (their
   ! eq. 2.4 is the 1s cross section), and Christie, Arras & Li (2013),
   ! whose level balance this module solves, leave sigma_2s,2p unstated.
   real*8, parameter :: s2_thr   = 1.4d-17             ! sigma_2 at threshold [cm^2]

   ! ----- Ly-alpha pumping (parameterized J_lya) data ----- !
   real*8, parameter :: sigma_LyC = 6.3d-18            ! H photoion. xsec at LyC [cm^2]

   ! Whether excited_H_update has run at least once, so that the level
   ! populations, the Ly-alpha field and the Balmer-continuum rates the
   ! level balance is measured against exist. Before that there is no field
   ! and no population: a residual would be a statement about zeros.
   logical, save :: excited_H_field_ready = .false.

   ! RT-supplied J_lya(r) profile (jlya_mode=1), interpolated onto the grid once.
   logical, save :: jlya_rt_loaded = .false.
   real*8, dimension(:), allocatable :: jlya_rt_grid

   ! Saved context for the diagnostic dump (filled in excited_H_update).
   real*8, dimension(:), allocatable :: Tdiag, nhidiag, nediag
   real*8, dimension(:), allocatable :: nhiidiag   ! nHII for the recomb. sink
   real*8, dimension(:), allocatable :: taulya     ! top-down Ly-alpha optical depth

   contains

   ! --------------------------------------------------------------- !

   subroutine excited_H_allocate_arrays
   ! Allocate the grid-sized module arrays once the number of cells N is
   ! known; called from EXHALE_main right after input_read. The zeros are
   ! the initializers the declarations used to carry.

   allocate(jlya_rt_grid(1-Ng:N+Ng))
   allocate(Tdiag(1-Ng:N+Ng), nhidiag(1-Ng:N+Ng), nediag(1-Ng:N+Ng))
   allocate(nhiidiag(1-Ng:N+Ng), taulya(1-Ng:N+Ng))

   jlya_rt_grid = 0.0d0
   Tdiag        = 0.0d0
   nhidiag      = 0.0d0
   nediag       = 0.0d0
   nhiidiag     = 0.0d0
   taulya       = 0.0d0

   end subroutine excited_H_allocate_arrays

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
   ! Cell-by-cell line-centre Ly-alpha depths: taulya_d(j) is cell j's own
   ! depth and taulya_out(j) the depth at its star-ward face, the pair the
   ! cell mean of the field is taken over.
   real*8, dimension(1-Ng:N+Ng) :: taulya_d, taulya_out
   real*8 :: F_LyC, F_inc, xi, abs_frac, N_HI_tot, a_cm
   real*8 :: Dnu_D, n2s, n2p, n2tot, relc
   ! The two molecular ions whose dissociative recombination leaves an H
   ! atom in n = 2, and the production that follows from them.
   real*8, dimension(1-Ng:N+Ng) :: n_h2p_mol, n_hehp_mol
   real*8 :: chem_n2

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

   ! H2+ and HeH+, the two carriers whose dissociative recombination leaves
   ! one hydrogen atom in n = 2. They exist only in a molecular run; in an
   ! atomic one the columns are absent and the term is identically zero, so
   ! nothing an atomic configuration computes is touched.
   if (thereis_mol) then
      n_h2p_mol  = f_sp_in(:,isp_H2p)*n_dim
      n_hehp_mol = f_sp_in(:,isp_HeHp)*n_dim
   else
      n_h2p_mol  = 0.0d0
      n_hehp_mol = 0.0d0
   endif

   ! ----- Day-night / 2D dilution factor xi ----- !
   ! global_parameters' dayside_dilution(), the single definition the ATES
   ! "2D approximate method" has; the ground-state EUV grid and the four FUV
   ! bands read the same one.
   xi = dayside_dilution()

   ! ----- Scalar Balmer-continuum rates (depend only on the stellar spectrum) ----- !
   ! Computed once per update; reduced by xi for the dayside hemisphere average,
   ! consistent with the ground-state EUV ionization treatment.
   gamma2_bal = xi*gamma_n2_balmer()
   hpe2_bal   = xi*heat_n2_balmer()

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
   ! Fraction of the incident Lyman continuum the column absorbs,
   ! 1 - e^-tau with tau = sigma_LyC N_HI, formed as tau (1 - e^-tau)/tau
   ! with the quotient from one_minus_exp_over_x (utilities.f90): the
   ! closed form cancels for a thin column, where the fraction tends to
   ! tau itself.
   abs_frac = (sigma_LyC*N_HI_tot)*one_minus_exp_over_x(sigma_LyC*N_HI_tot)
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
      call jlya_escape_prob(T_K, nhi, nhii, nheii, nheiii, ne, v_in,     &
                            Jlya_arr, taulya)
   else
      ! (a) Parameterized J_lya = 0.1 F_LyC/Dnu_D (Huang+2017 Eq. 6), attenuated
      ! by 1/(1+tau_lya) with tau_lya the top-down line-centre Ly-alpha optical
      ! depth, so the pumping vanishes below the Ly-alpha photosphere. This is
      ! an approximate staging estimate; the accurate field is jlya_mode=1.
      !
      ! The pumping rate of a cell is the MEAN of that field over the cell's
      ! own line-centre depth, not its value at a face: taulya(j) is the depth
      ! at the cell's INNER face and taulya_out(j) the depth at its star-ward
      ! face, and lya_photosphere_attenuation_cell_mean holds the closed-form
      ! mean of 1/(1+tau) between them.
      call lya_line_center_optical_depth(T_K, nhi, taulya, taulya_d,        &
                                         taulya_out)
      do j = 1-Ng, N+Ng
         Dnu_D = nu_lya*sqrt(2.0d0*kb_erg*max(T_K(j),1.0d0)/mu)/c_light
         Jlya_arr(j) = 0.1d0*F_LyC/max(Dnu_D,1.0d-30)                       &
                     *lya_photosphere_attenuation_cell_mean(taulya_out(j),  &
                                                            taulya_d(j))
      enddo
   endif

   ! ----- Cell-by-cell n=2 populations + feedback ----- !
   do j = 1-Ng, N+Ng
      chem_n2 = dissociative_recombination_n2_source(T_K(j),               &
                       n_h2p_mol(j), n_hehp_mol(j), max(ne(j),0.0d0))
      call n2_populations(T_K(j), max(nhi(j),0.0d0), max(nhii(j),0.0d0),    &
                          max(nheii(j),0.0d0), max(nheiii(j),0.0d0),        &
                          max(ne(j),0.0d0), Jlya_arr(j),                    &
                          gamma2_bal, gamma2_bal, chem_n2, n2s, n2p)
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

   excited_H_field_ready = .true.

   end subroutine excited_H_update

   ! --------------------------------------------------------------- !

   subroutine n2_populations(T, n1s, nHII_l, nHeII_l, nHeIII_l, ne_l, Jlya, &
                             gam_ion_2s, gam_ion_2p, chem_n2, n2s, n2p)
   ! Christie+2013 Eqs. 12-13: solve the 2x2 2s/2p rate equilibrium for the
   ! H(n=2) populations [cm^-3]. T in K, densities in cm^-3, Jlya in cgs
   ! (erg s^-1 cm^-2 Hz^-1 sr^-1). gam_ion_2s/gam_ion_2p = photoionization
   ! rates [s^-1] out of 2s and 2p (the stellar Balmer continuum).
   !
   ! The rate coefficients come from hydrogen_n2_rates so that the statistical
   ! weights g1s/g2s/g2p are evaluated once, at module scope, and cannot be
   ! shadowed by a dummy argument of this routine (Fortran is case-insensitive,
   ! and the dummies used to be called G2s/G2p).
   !
   ! chem_n2 is the chemical H(n=2) production of the molecular network
   ! [cm^-3 s^-1]; see n2_rate_matrix for what it is and how it is split
   ! between the two levels.

   real*8, intent(in)  :: T, n1s, nHII_l, nHeII_l, nHeIII_l, ne_l, Jlya
   real*8, intent(in)  :: gam_ion_2s, gam_ion_2p
   real*8, intent(in)  :: chem_n2
   real*8, intent(out) :: n2s, n2p

   real*8 :: L2p, L2s, S2p, S2s, M12, M21, det

   call n2_rate_matrix(T, n1s, nHII_l, nHeII_l, nHeIII_l, ne_l, Jlya,     &
                       gam_ion_2s, gam_ion_2p, chem_n2,                   &
                       L2p, L2s, S2p, S2s, M12, M21)

   det = L2p*L2s - M12*M21
   if (abs(det) .le. 0.0d0) det = 1.0d0

   n2p = max((S2p*L2s + M12*S2s)/det, 0.0d0)
   n2s = max((L2p*S2s + M21*S2p)/det, 0.0d0)

   end subroutine n2_populations

   ! --------------------------------------------------------------- !

   subroutine n2_rate_matrix(T, n1s, nHII_l, nHeII_l, nHeIII_l, ne_l,     &
                             Jlya, gam_ion_2s, gam_ion_2p, chem_n2,       &
                             L2p, L2s, S2p, S2s, M12, M21)
   ! The 2x2 rate matrix and source vector of the 2s/2p statistical
   ! equilibrium (Christie+2013 Eqs. 12-13), the ONE definition of the
   ! coefficients: n2_populations solves the system with them and
   ! excited_hydrogen_level_residual measures the imbalance of a given pair
   ! of populations against them, so the two cannot drift apart.
   !
   !   L2p n2p - M12 n2s = S2p ,   L2s n2s - M21 n2p = S2s
   !
   ! L is the total destruction rate of a level [s^-1], M the l-mixing
   ! transfer from the other one [s^-1] and S its production [cm^-3 s^-1].
   !
   ! chem_n2 is the CHEMICAL production of H(n=2) [cm^-3 s^-1]: the two
   ! dissociative recombinations of the molecular network, H2+ + e and
   ! HeH+ + e, which leave one of their hydrogen fragments in n = 2 rather
   ! than in the ground state (Takagi 2002, Phys. Scr. T96, 52; Guberman
   ! 1994, Phys. Rev. A 49, R4277; sources and ranges at the site that
   ! forms the term, molecular_reaction_heat). It is zero in an atomic run,
   ! where the network does not exist.
   !
   ! WHICH OF THE TWO LEVELS RECEIVES IT is not resolved by the sources and
   ! is assigned by statistical weight, g2s : g2p = 1 : 3. Takagi states the
   ! product as "n = 2" without an l, and Guberman's HeH+ C state likewise
   ! dissociates to "an excited n = 2 H atom". Giusti-Suzor, Bardsley &
   ! Derkits (1983), Phys. Rev. A 28, 682, do name the H(1s) + H(2s) limit,
   ! but for H2+ in its lowest three vibrational states and electron
   ! energies below 0.5 eV, and this code's H2+ is vibrationally hot (its
   ! lifetime against R5, R8 and R9 is about 1e-3 s against a radiative
   ! vibrational relaxation of about 1 s), so that condition does not hold
   ! here. THE SPLIT IS NOT NEUTRAL and is recorded as open: 2p decays by
   ! Lyman-alpha at A_2p1s, 2s only by the two-photon continuum, so an
   ! all-2s assignment would route the same energy out of the gas by a
   ! different channel and leave a larger n = 2 population behind.

   real*8, intent(in)  :: T, n1s, nHII_l, nHeII_l, nHeIII_l, ne_l, Jlya
   real*8, intent(in)  :: gam_ion_2s, gam_ion_2p
   real*8, intent(in)  :: chem_n2
   real*8, intent(out) :: L2p, L2s, S2p, S2s, M12, M21

   real*8 :: Tl, a2s, a2p
   real*8 :: C1s2s, C1s2p, C2s1s, C2p1s, Mix_2s2p, Mix_2p2s
   real*8 :: Ppump, Pstim

   Tl = max(T, 1.0d0)

   ! Level-resolved recombination: the balance's case-B coefficient split
   ! into 2s and 2p (hydrogen_n2_rates).
   a2s = alpha_2s_hydrogen(Tl)
   a2p = alpha_2p_hydrogen(Tl)

   ! Collisional excitation 1s->2s, 1s->2p and the reverse 2s/2p->1s
   ! de-excitation (hydrogen_n2_rates: the H I rate set of Cool_coeff and
   ! its detailed-balance partners), and the 2s <-> 2p l-mixing rates
   ! [s^-1] by electrons (Seaton 1955) and by H+, He+ and He2+ (Pengelly
   ! & Seaton 1964), which at 1e4 K the protons dominate by a factor 8.3
   ! in ionized gas.
   C1s2s = c1s2s_rate(Tl)
   C1s2p = c1s2p_rate(Tl)
   C2s1s = c2s1s_rate(Tl)
   C2p1s = c2p1s_rate(Tl)
   Mix_2s2p = l_mixing_rate_2s2p(Tl, ne_l, nHII_l, nHeII_l, nHeIII_l)
   Mix_2p2s = l_mixing_rate_2p2s(Tl, ne_l, nHII_l, nHeII_l, nHeIII_l)

   ! Ly-alpha pump (1s->2p) and stimulated emission (2p->1s).
   Ppump = B12_lya*Jlya
   Pstim = B21_lya*Jlya

   ! 2x2 system [[L2p,-M12],[-M21,L2s]] [n2p,n2s]^T = [S2p,S2s]^T.
   ! The cascade source is alpha_2l * ne * nHII: the recombining partner of an
   ! electron is a proton, not another electron. The two differ wherever the
   ! electrons come from helium and metals while hydrogen is still neutral,
   ! which is the case through the base.
   L2p = A_2p1s + Pstim + C2p1s*ne_l + Mix_2p2s + gam_ion_2p
   L2s = C2s1s*ne_l + Mix_2s2p + gam_ion_2s + A_2s1s
   S2p = (Ppump + C1s2p*ne_l)*n1s + a2p*ne_l*nHII_l                      &
       + max(chem_n2, 0.0d0)*g2p/(g2s + g2p)
   S2s = (C1s2s*ne_l)*n1s + a2s*ne_l*nHII_l                              &
       + max(chem_n2, 0.0d0)*g2s/(g2s + g2p)
   M12 = Mix_2s2p
   M21 = Mix_2p2s

   end subroutine n2_rate_matrix

   ! --------------------------------------------------------------- !

   subroutine excited_hydrogen_level_residual(T_in, n_in, f_sp_in,        &
                                              res, scale, ok, why)
   ! THE H(n=2) LEVEL BALANCE, MEASURED ON A STATE.
   !
   ! The n=2 populations are not unknowns of any solve: n2_populations
   ! closes them one outer pass BEHIND the composition they are a closure
   ! for, and their products enter the ionization system as a lagged rate
   ! and a lagged heat. What that lag is worth has never been a number.
   ! This is the number: production minus destruction of each of the two
   ! levels, for the populations the module currently holds, at the field
   ! it currently holds, and at the composition of the state handed in.
   !
   !   r2p = L2p n2p - M12 n2s - S2p ,   r2s = L2s n2s - M21 n2p - S2s
   !
   ! [cm^-3 s^-1]. NOTHING IS UPDATED: n2s_arr, n2p_arr, Jlya_arr,
   ! gph_balmer_HI and heat_balmer are read and left as they are.
   !
   ! res(j) and scale(j) are the row of the two that is furthest out in
   ! units of its OWN terms, so that the level whose rates are the smaller
   ! of the pair is not hidden by the other. The scale is the sum of the
   ! magnitudes of that row's terms; the floor beneath it, 1e-30 cm^-3
   ! s^-1, is numerical only -- one n=2 atom in 3e21 s -- and exists so
   ! that a cell with no hydrogen at all divides by something.
   !
   ! THE FIELD IS THE ONE THE STATE CARRIES, in the sense the code makes
   ! available: Jlya_arr and the Balmer-continuum rates gamma2_bal /
   ! hpe2_bal were built by the last excited_H_update, which is the update
   ! this state's own lagged coupling used. Recomputing them here would
   ! measure a different closure from the one the run solved.

   real*8, dimension(1-Ng:N+Ng),           intent(in)  :: T_in, n_in
   real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp_in
   real*8, dimension(1:N),                 intent(out) :: res, scale
   logical,                                intent(out) :: ok
   character(len=*),                       intent(out) :: why

   integer :: j, im
   real*8, dimension(1-Ng:N+Ng) :: T_K, n_dim, nhi, nhii
   real*8, dimension(1-Ng:N+Ng) :: nhei, nheii, nheiii, ne
   real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
   real*8, dimension(1-Ng:N+Ng) :: n_h2p_mol, n_hehp_mol
   real*8 :: L2p, L2s, S2p, S2s, M12, M21
   real*8 :: r2p, r2s, s2pd, s2sd, q2p, q2s, chem_n2

   res = 0.0d0;  scale = 1.0d0;  ok = .false.;  why = ''
   if (.not. use_excited_H) then
      why = 'the run does not carry H(n=2)'
      return
   endif
   if (.not. excited_H_field_ready) then
      why = 'no excited-H update has built the field and the populations'
      return
   endif

   ! Densities exactly as excited_H_update dimensionalizes them, so the
   ! residual is evaluated on the same quantities the closure was.
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
   call calc_ne(nhii, nheii, nheiii, ne, nm)
   ! The chemical n = 2 source of the same state, so that the residual
   ! measures the balance the closure solved and not one without that term.
   if (thereis_mol) then
      n_h2p_mol  = f_sp_in(:,isp_H2p)*n_dim
      n_hehp_mol = f_sp_in(:,isp_HeHp)*n_dim
   else
      n_h2p_mol  = 0.0d0
      n_hehp_mol = 0.0d0
   endif

   do j = 1, N
      chem_n2 = dissociative_recombination_n2_source(T_K(j),              &
                      n_h2p_mol(j), n_hehp_mol(j), max(ne(j),0.0d0))
      call n2_rate_matrix(T_K(j), max(nhi(j),0.0d0), max(nhii(j),0.0d0),  &
                          max(nheii(j),0.0d0), max(nheiii(j),0.0d0),      &
                          max(ne(j),0.0d0), Jlya_arr(j),                  &
                          gamma2_bal, gamma2_bal, chem_n2,                &
                          L2p, L2s, S2p, S2s, M12, M21)
      r2p  = L2p*n2p_arr(j) - M12*n2s_arr(j) - S2p
      r2s  = L2s*n2s_arr(j) - M21*n2p_arr(j) - S2s
      s2pd = abs(L2p*n2p_arr(j)) + abs(M12*n2s_arr(j)) + abs(S2p)
      s2sd = abs(L2s*n2s_arr(j)) + abs(M21*n2p_arr(j)) + abs(S2s)
      q2p  = abs(r2p)/max(s2pd, 1.0d-30)
      q2s  = abs(r2s)/max(s2sd, 1.0d-30)
      if (q2p .ge. q2s) then
         res(j) = r2p;  scale(j) = max(s2pd, 1.0d-30)
      else
         res(j) = r2s;  scale(j) = max(s2sd, 1.0d-30)
      endif
   enddo
   ok = .true.

   end subroutine excited_hydrogen_level_residual

   ! --------------------------------------------------------------- !

   real*8 function gamma_n2_balmer()
   ! n=2 photoionization rate [s^-1] in the stellar Balmer continuum, the
   ! band from the n=2 edge e_th_HI_n2 = 3.400 eV (3647 A) to the H I edge
   ! 13.6 eV, above which the ground state absorbs the photons.
   !
   !   gamma_2 = INT_{e_th_HI_n2}^{13.6 eV} F_E(E)/(h nu) sigma_2(E) dE ,
   !
   ! the photon number flux per unit energy times the hydrogenic cross
   ! section sigma_2(E) = s2_thr*(e_th_HI_n2/E)^3, with the photon energy
   ! written as h nu = hp_erg (E eV2Hz) [erg] from the two exact CODATA
   ! constants, and dE in eV.
   !
   ! THE FIELD IS THE RUN'S OWN SPECTRUM: stellar_flux_eV is the power law for
   ! "Spectrum type: Power-law", the photospheric blackbody for "Planck" and
   ! the loaded table for "Load". No band is built from a type the input did
   ! not select, so this integral is no longer a blackbody in a power-law
   ! run.
   !
   ! Undiluted: excited_H_update applies dayside_dilution() to the result.

   real*8 :: gam2, hpe2

   call balmer_band_integrals(gam2, hpe2)
   gamma_n2_balmer = gam2

   end function gamma_n2_balmer

   ! --------------------------------------------------------------- !

   real*8 function heat_n2_balmer()
   ! Photoelectric heating per H(n=2) atom [erg s^-1]: the integrand of
   ! gamma_n2_balmer weighted by the photoelectron excess energy
   ! h(nu - nu_2) = (E - e_th_HI_n2), in the same band and from the same
   ! field, the run's own spectrum type,
   !
   !   heat_2 = INT F_E(E) sigma_2(E) (1 - e_th_HI_n2/E) dE ,
   !
   ! on the same nodes as the rate. Undiluted, as gamma_n2_balmer is.

   real*8 :: gam2, hpe2

   call balmer_band_integrals(gam2, hpe2)
   heat_n2_balmer = hpe2

   end function heat_n2_balmer

   ! --------------------------------------------------------------- !

   subroutine balmer_band_integrals(gam2, hpe2)
   ! The two band integrals of the Balmer continuum, over the same nodes:
   !
   !   gam2 = Kg INT F_E(E) E^-4 dE               [s^-1]
   !   hpe2 = Kh INT F_E(E) (E^-3 - E_2 E^-4) dE  [erg s^-1]
   !
   ! with E_2 = e_th_HI_n2, Kg = s2_thr E_2^3/(hp_erg eV2Hz) and
   ! Kh = s2_thr E_2^3, which is the pair above with sigma_2(E) =
   ! s2_thr (E_2/E)^3 taken out of the integrand.
   !
   ! THE QUADRATURE RESOLVES THE FIELD IT INTEGRATES. A rule whose nodes are
   ! fixed independently of the spectrum cannot: a measured stellar table
   ! carries emission lines, Ly-alpha among them at 10.20 eV inside this
   ! band, that are narrower than any fixed step, and refining the table
   ! narrows the interpolated line further while the step stays put, so the
   ! error grows instead of falling. The nodes are therefore taken from the
   ! field:
   !
   !   Load        the table's OWN rows inside the band, plus the two band
   !               edges. Between two rows stellar_flux_eV is log-log, i.e.
   !               F = F_k (E/E_k)^p exactly, so each integrand above is a
   !               sum of powers of E and the segment integral is closed
   !               form (power_law_moment) rather than a rule. Where a row
   !               is not positive stellar_flux_eV is linear in E instead,
   !               and that segment is integrated in closed form too. The
   !               result is the exact integral of the field the code reads,
   !               to round-off, at any line width.
   !   Power-law   J_inc is a single power law on each side of e_mid, so the
   !               band (with e_mid inserted if it falls inside it) is the
   !               same closed form, again exact.
   !   Planck      pi B_nu is smooth and has no structure below its peak, so
   !               a composite 5-point Gauss-Legendre rule on a geometric
   !               subdivision of the band integrates it to round-off.
   !
   ! The band is [e_th_HI_n2, e_th_HI] exactly, including the partial
   ! intervals at both ends where a table row does not fall on the edge.
   real*8, intent(out) :: gam2, hpe2

   ! Geometric subdivisions of the band, and the 5-point Gauss-Legendre
   ! nodes and weights on [-1,1] (Abramowitz & Stegun Table 25.4), for a
   ! field that is not a power law on any segment.
   integer, parameter :: n_sub = 64
   integer, parameter :: n_gl  = 5
   real*8, parameter :: x_gl(n_gl) = (/ -9.0617984593866400d-1,            &
                                        -5.3846931010568309d-1,            &
                                         0.0d0,                            &
                                         5.3846931010568309d-1,            &
                                         9.0617984593866400d-1 /)
   real*8, parameter :: w_gl(n_gl) = (/  2.3692688505618909d-1,            &
                                         4.7862867049936647d-1,            &
                                         5.6888888888888889d-1,            &
                                         4.7862867049936647d-1,            &
                                         2.3692688505618909d-1 /)

   integer :: k, i, n_node, n_brk
   real*8  :: Ea, Eb, Kg, Kh, a, b, sg, sh, e1, e2, Em, Fm, brk(3)
   real*8  :: q_ratio, xm, xh, xx, F_E, sig2

   gam2 = 0.0d0
   hpe2 = 0.0d0
   call stop_if_balmer_band_unstated

   Ea = e_th_HI_n2
   Eb = e_th_HI
   if (Eb .le. Ea) return
   Kg = s2_thr*e_th_HI_n2**3.0d0/(hp_erg*eV2Hz)
   Kh = s2_thr*e_th_HI_n2**3.0d0

   if (do_read_sed .and. allocated(e_sed_node) .and.                       &
       allocated(F_sed_node)) then

      n_node = size(e_sed_node)
      if (n_node .lt. 2) return
      do k = 1, n_node-1
         e1 = e_sed_node(k)
         e2 = e_sed_node(k+1)
         a  = e1
         ! Below the lowest row stellar_flux_eV continues the first segment,
         ! over the less than one row between the grid floor and that row;
         ! the same continuation is integrated here.
         if (k .eq. 1) a = min(e1, Ea)
         a = max(a, Ea)
         b = min(e2, Eb)
         if (b .le. a) cycle
         call sed_segment_balmer(a, b, e1, F_sed_node(k), e2,              &
                                 F_sed_node(k+1), Kg, Kh, sg, sh)
         gam2 = gam2 + sg
         hpe2 = hpe2 + sh
      enddo

   else if (is_PL_sed) then

      n_brk    = 2
      brk(1)   = Ea
      brk(2)   = Eb
      if (e_mid .gt. Ea .and. e_mid .lt. Eb) then
         brk(2) = e_mid
         brk(3) = Eb
         n_brk  = 3
      endif
      do i = 1, n_brk-1
         a  = brk(i)
         b  = brk(i+1)
         ! One power law of index PLind on this piece; its normalization is
         ! read from the field itself at the geometric midpoint, so the
         ! branch J_inc takes there is the branch integrated.
         Em = sqrt(a*b)
         Fm = stellar_flux_eV(Em)
         sg = power_law_moment(a, b, Em, Fm, PLind, -4.0d0)
         sh = power_law_moment(a, b, Em, Fm, PLind, -3.0d0)
         gam2 = gam2 + Kg*sg
         hpe2 = hpe2 + Kh*(sh - Ea*sg)
      enddo

   else

      ! Geometric subdivision: the field falls by orders of magnitude across
      ! the band and the nodes follow it.
      q_ratio = (Eb/Ea)**(1.0d0/dble(n_sub))
      a = Ea
      do i = 1, n_sub
         b = a*q_ratio
         if (i .eq. n_sub) b = Eb
         xm = 0.5d0*(b + a)
         xh = 0.5d0*(b - a)
         do k = 1, n_gl
            xx   = xm + xh*x_gl(k)
            F_E  = stellar_flux_eV(xx)
            sig2 = s2_thr*(e_th_HI_n2/xx)**3.0d0
            sg   = w_gl(k)*xh*F_E*sig2
            gam2 = gam2 + sg/(hp_erg*xx*eV2Hz)
            hpe2 = hpe2 + sg*(1.0d0 - e_th_HI_n2/xx)
         enddo
         a = b
      enddo

   endif

   end subroutine balmer_band_integrals

   ! --------------------------------------------------------------- !

   subroutine sed_segment_balmer(a, b, e1, f1, e2, f2, Kg, Kh, sg, sh)
   ! The two Balmer integrands of a loaded table, integrated in closed form
   ! over [a,b] inside the table segment [e1,e2], from the SAME
   ! reconstruction stellar_flux_eV uses there: F = f1 (E/e1)^p with
   ! p = ln(f2/f1)/ln(e2/e1) where both rows are positive, and F linear in E
   ! otherwise.
   real*8, intent(in)  :: a, b, e1, f1, e2, f2, Kg, Kh
   real*8, intent(out) :: sg, sh
   real*8 :: p, slope, c0, c1, m4, m3, m2

   sg = 0.0d0
   sh = 0.0d0

   if (f1 .gt. 0.0d0 .and. f2 .gt. 0.0d0 .and. e2 .gt. e1) then

      p  = log(f2/f1)/log(e2/e1)
      m4 = power_law_moment(a, b, e1, f1, p, -4.0d0)
      m3 = power_law_moment(a, b, e1, f1, p, -3.0d0)
      sg = Kg*m4
      sh = Kh*(m3 - e_th_HI_n2*m4)

   else if (e2 .gt. e1) then

      ! F = c0 + c1 E, the linear reconstruction of a segment carrying a
      ! non-positive row.
      slope = (f2 - f1)/(e2 - e1)
      c1    = slope
      c0    = f1 - e1*slope
      m4    = power_moment(a, b, -4.0d0)
      m3    = power_moment(a, b, -3.0d0)
      m2    = power_moment(a, b, -2.0d0)
      sg    = Kg*(c0*m4 + c1*m3)
      sh    = Kh*(c0*m3 + c1*m2 - e_th_HI_n2*(c0*m4 + c1*m3))

   else

      ! Two rows at the same energy: stellar_flux_eV returns f1 there.
      m4 = power_moment(a, b, -4.0d0)
      m3 = power_moment(a, b, -3.0d0)
      sg = Kg*f1*m4
      sh = Kh*f1*(m3 - e_th_HI_n2*m4)

   endif

   end subroutine sed_segment_balmer

   ! --------------------------------------------------------------- !

   real*8 function power_law_moment(a, b, e0, f0, p, q)
   ! INT_a^b f0 (E/e0)^p E^q dE, the exact integral of a power-law field
   ! against a power-law cross section.
   !
   !   = f0 (a/e0)^p a^(q+1) ln(b/a) [exp(w) - 1]/w ,   w = (p+q+1) ln(b/a),
   !
   ! which is the elementary primitive rewritten so that no large power of a
   ! or of e0 is formed, and so that the near-cancellation of b^(m+1) -
   ! a^(m+1) over a narrow interval is carried by [exp(w)-1]/w instead. The
   ! index m = p+q = -1 is the logarithmic case and is the limit w -> 0 of
   ! the same expression, not a separate branch.
   real*8, intent(in) :: a, b, e0, f0, p, q
   real*8 :: u, w, ratio

   power_law_moment = 0.0d0
   if (a .le. 0.0d0 .or. b .le. a .or. e0 .le. 0.0d0) return

   u = log(b/a)
   w = (p + q + 1.0d0)*u
   if (abs(w) .lt. 1.0d-6) then
      ! [exp(w)-1]/w to double precision without the cancellation.
      ratio = 1.0d0 + w*(0.5d0 + w*(1.0d0/6.0d0 + w/24.0d0))
   else
      ratio = (exp(w) - 1.0d0)/w
   endif

   power_law_moment = f0*(a/e0)**p*a**(q + 1.0d0)*u*ratio

   end function power_law_moment

   ! --------------------------------------------------------------- !

   real*8 function power_moment(a, b, q)
   ! INT_a^b E^q dE, the same expression with a unit field.
   real*8, intent(in) :: a, b, q

   power_moment = power_law_moment(a, b, a, 1.0d0, 0.0d0, q)

   end function power_moment

   ! --------------------------------------------------------------- !

   subroutine stop_if_balmer_band_unstated
   ! The Balmer continuum can only be integrated over a band the run's
   ! spectrum states a field on. A monochromatic run states nothing there,
   ! and filling the band from another type is forbidden;
   ! the run therefore stops, with its message in the same
   ! form as the SED coverage stop of sed_read.f90.
   !
   ! THE LAST LINE OF DEFENCE, not the first. e_th_HI_n2 floors the photon
   ! grid whenever this coupling is armed (sed_read's photon_grid_floor_eV),
   ! so a loaded table that stops above 3647 A is refused by read_sed, at
   ! startup and naming the file, before any of this is reached. What
   ! survives to here is a spectrum type that states no field at all over
   ! the band.

   if (spectrum_covers_eV(e_th_HI_n2) .and. spectrum_covers_eV(e_th_HI)) return

   write(*,*) '(excited_hydrogen.f90) ERROR: the spectrum of this run '//  &
              'states no field over the Balmer continuum that '//         &
              'photoionizes H(n=2).'
   if (allocated(sp_type))                                                &
      write(*,'(A,A)')  '    spectrum type        : ', trim(sp_type)
   write(*,'(A,ES12.5,A,ES12.5,A)')                                       &
      '    band needed          : ', e_th_HI_n2, ' to ', e_th_HI, ' eV'
   write(*,'(A,ES12.5,A,ES12.5,A)')                                       &
      '                         = ', hp_eV*c_light*1.0d8/e_th_HI,         &
      ' to ', hp_eV*c_light*1.0d8/e_th_HI_n2, ' A'
   if (allocated(e_v) .and. Nl .ge. 2)                                    &
      write(*,'(A,ES12.5,A,ES12.5,A)')                                    &
         '    band the run has     : ', e_v(1), ' to ', e_v(Nl), ' eV'
   write(*,*) '   absorber that needs it:'
   write(*,'(A,ES12.5,A,ES12.5,A)')                                       &
      '       H(n=2), the Balmer continuum : ', e_th_HI_n2, ' eV = ',         &
      hp_eV*c_light*1.0d8/e_th_HI_n2, ' A'
   if (do_read_sed .and. allocated(sed_file))                             &
      write(*,'(A,A)')  '    SED file             : ', trim(sed_file)
   write(*,*) '   remedies (there is no key to continue):'
   write(*,*) '      remove "Stellar Teff [K]:" and "Stellar radius'//    &
              ' [R_sun]:" from input.inp,'
   write(*,*) '      which is what arms the excited-hydrogen coupling'//  &
              ' (input_read.f90 use_excited_H); or'
   write(*,'(A,ES12.5,A)')                                                &
      '       supply a spectrum that reaches ',                           &
      hp_eV*c_light*1.0d8/e_th_HI_n2, ' A'
   error stop 1

   end subroutine stop_if_balmer_band_unstated

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
   call write_row_layout_header(72)
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
