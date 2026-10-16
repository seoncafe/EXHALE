   module hydrogen_n2_rates
   ! Atomic data and rate coefficients of the hydrogen n=2 manifold and of the
   ! Ly-alpha line, in one place.
   !
   ! Both the H(n=2) statistical-equilibrium solver (excited_hydrogen.f90) and
   ! the Ly-alpha escape-probability transfer (lya_rt.f90) need the same
   ! constants (A_2p1s, A_2s1s, the statistical weights, the line opacity
   ! constant) and the same collisional rate coefficients. They used to carry
   ! private copies of both. The copies drifted: the reverse-rate expressions
   ! written out inside n2_populations evaluated their statistical weights
   ! against shadowing dummy arguments and so disagreed with the identical
   ! expressions a few lines further down in the same file. Keeping one
   ! definition here removes that failure mode.
   !
   ! Collisional excitation and de-excitation 1s <-> 2s, 2p: the one rate
   ! set of Cool_coeff.f90 (CHIANTI v11 / Anderson et al. 2000), which the
   ! H I collisional-excitation cooling uses as well. 2s <-> 2p l-mixing:
   ! by electrons, Seaton (1955, Proc. Phys. Soc. A 68, 457), eq. (55)
   ! (c2s2p_rate); by protons, He+ and He2+, Pengelly & Seaton (1964)
   ! (Cool_coeff: l_mixing_2s2p_pengelly_seaton). Recombination: the case-B coefficient
   ! of the ionization balance, split between 2s and 2p with the case-B
   ! fraction of Pengelly (1964, Table I) (Cool_coeff:
   ! case_b_2s_fraction_hydrogenic). All quantities cgs.

   use global_parameters, only: hp_erg, c_light, erg2eV, m_p, mu,         &
                                m_He_atom, m_e, kb_erg, pi
   use Cooling_Coefficients, only: alpha_rec_HII_B, E_HI_n2_eV,           &
                                   excitation_rate_HI_1s2s,               &
                                   excitation_rate_HI_1s2p,               &
                                   deexcitation_rate_HI_2s1s,             &
                                   deexcitation_rate_HI_2p1s,             &
                                   case_b_2s_fraction_hydrogenic,         &
                                   l_mixing_2s2p_pengelly_seaton

   implicit none

   ! ----- Ly-alpha line and n=2 level data ----- !
   real*8, parameter :: lA_lya   = 1215.6701d-8        ! Ly-alpha wavelength [cm]
   real*8, parameter :: nu_lya   = c_light/lA_lya      ! Ly-alpha frequency [s^-1]
   ! A(2p->1s) [s^-1].  NIST ASD (Wiese & Fuhr 2009, J. Phys. Chem. Ref.
   ! Data 38, 565) give 6.2649e8 for Ly-alpha; it is the nonrelativistic
   ! electric-dipole value (2^8/3^8) alpha^4 c/a_0 = 6.2683e8 carrying the
   ! reduced-mass factor mu/m_e = 0.9994557, since A goes as omega^3 a_0^2
   ! and so as mu.  ONE definition: the Voigt parameter of lya_rt.f90 and
   ! every n = 2 rate below read this one.
   real*8, parameter :: A_2p1s   = 6.2649d8
   ! A(2s->1s) [s^-1], the two-photon decay: Drake (1986, Phys. Rev. A 34,
   ! 2871, eq. 27, READ from the published article),
   !   w = 8.22938 Z^6 Z_r^4 (1 - m/M) [1 + 3.9448 (aZ)^2 - 2.040 (aZ)^4]
   !       / [1 + 4.6019 (aZ)^2],
   ! the nonrelativistic hydrogenic rate (his Table I, 8.2293810 Z^6,
   ! "substantially" more accurate than Klarsfeld's) with the finite
   ! nuclear mass (1 - m/M) and the effective radiative charge
   ! Z_r = (Z - 1) m/(M + m) + 1, times his fit to the relativistic
   ! calculations ("an error of less than +-0.005% in the range
   ! [1 <=] Z <= 92"). For H (Z = 1, Z_r = 1, m/M = m_e/m_p) it gives
   ! 8.22461 s^-1 (DERIVED); the relativistic factor is 0.999965. The
   ! Nussbaumer & Schmutz (1984, A&A 138, 495, eq. 6) value 8.2249 s^-1,
   ! nonrelativistic, is 0.0035% above it; Osterbrock & Ferland (2006,
   ! sect. 2.2, p. 20; sect. 4.3, p. 84) print 8.23 s^-1, and the 8.26 s^-1
   ! of Christie et al. (2013, Table 2, R10), attributed to them, is not
   ! the value they print.
   real*8, parameter :: A_2s1s   = 8.22461d0
   real*8, parameter :: g1s = 2.0d0, g2s = 2.0d0, g2p = 6.0d0
   real*8, parameter :: E21_erg  = E_HI_n2_eV/erg2eV   ! 1s-2s/2p gap, 10.2 eV [erg]
   ! Einstein-B in the J_nu (mean-intensity) convention: B*J gives s^-1.
   real*8, parameter :: B21_lya = A_2p1s*c_light**2/(2.0d0*hp_erg*nu_lya**3)
   real*8, parameter :: B12_lya = (g2p/g1s)*B21_lya
   ! Line-opacity constant: sqrt(pi) e^2/(m_e c) f_lya [cm^2 Hz], so that the
   ! line-center cross section is C_lya/Dnu_D.
   real*8, parameter :: f_lya    = 0.4162d0            ! Ly-alpha oscillator strength
   real*8, parameter :: C_lya    = 1.49736d-2*f_lya    ! [cm^2 Hz]

   ! 2s -> 2p l-mixing by charged particles (Pengelly & Seaton 1964): the
   ! separations of the final fine-structure states [erg], 2s1/2 - 2p1/2
   ! (the Lamb shift) 0.035 cm^-1 and 2s1/2 - 2p3/2 0.331 cm^-1 (NIST ASD
   ! as distributed in CHIANTI v11 h_1.elvlc, 3 decimals: the rate depends
   ! on them logarithmically, 1% in dE is 0.2% in the rate), and the
   ! reduced masses of the collisions of H with H+, He+ and He2+ [g].
   real*8, parameter :: dE_2s2p12 = 0.035d0*hp_erg*c_light
   real*8, parameter :: dE_2s2p32 = 0.331d0*hp_erg*c_light
   real*8, parameter :: mu_H_p    = m_p*mu/(m_p + mu)
   real*8, parameter :: mu_H_e    = m_e*mu/(m_e + mu)
   real*8, parameter :: mu_H_Hep  = (m_He_atom - m_e)*mu                 &
                                    /(m_He_atom - m_e + mu)
   real*8, parameter :: mu_H_He2p = (m_He_atom - 2.0d0*m_e)*mu           &
                                    /(m_He_atom - 2.0d0*m_e + mu)

   contains

   ! --------------------------------------------------------------- !
   ! Collisional rate coefficients [cm^3 s^-1]. T in K.

   real*8 function c1s2s_rate(T)
   ! 1s->2s collisional excitation (Cool_coeff: excitation_rate_HI_1s2s).
   real*8, intent(in) :: T
   c1s2s_rate = excitation_rate_HI_1s2s(max(T,1.0d0))
   end function c1s2s_rate

   real*8 function c1s2p_rate(T)
   ! 1s->2p collisional excitation, both 2p levels.
   real*8, intent(in) :: T
   c1s2p_rate = excitation_rate_HI_1s2p(max(T,1.0d0))
   end function c1s2p_rate

   real*8 function c2s1s_rate(T)
   ! 2s->1s collisional de-excitation, the detailed-balance reverse of
   ! c1s2s_rate with the Boltzmann factor cancelled analytically (no 0*inf
   ! at very low T).
   real*8, intent(in) :: T
   c2s1s_rate = deexcitation_rate_HI_2s1s(max(T,1.0d0))
   end function c2s1s_rate

   real*8 function c2p1s_rate(T)
   ! 2p->1s collisional de-excitation (detailed balance from c1s2p_rate).
   real*8, intent(in) :: T
   c2p1s_rate = deexcitation_rate_HI_2p1s(max(T,1.0d0))
   end function c2p1s_rate

   real*8 function c2s2p_rate(T)
   ! 2s -> 2p l-mixing by electron impact [cm^3 s^-1]: Seaton (1955, Proc.
   ! Phys. Soc. A 68, 457, READ from the published article), the paper
   ! written for these cross sections ("Cross Sections for 2s-2p
   ! Transitions in H ... Produced by Electron and by Proton Impact").
   ! His approximation V (Bethe approximation with allowance for strong
   ! coupling, sect. 3.3), eq. (55),
   !   Omega(dE) = 72 M^2 {ln(2 M v^2/dE) - mu},  mu = 2.21 for electrons,
   ! in atomic units (M the reduced mass; line strength S = 27), with the
   ! two final states 2p1/2 and 2p3/2 weighted 1/3 and 2/3 at their own
   ! separations dE (eq. 18). As a cross section, Q = pi a0^2 Omega/k^2 =
   ! 72 pi a0^2 (Ry/E) [ln(4E/dE) - mu] for E = M v^2/2, whose Maxwell
   ! average is closed:
   !   q = 72 pi a0^2 Ry (8 M/pi)^(1/2)/(m_e (kT)^(1/2))
   !       sum_j w_j [ln(4kT/dE_j) - gamma - mu],
   ! gamma = 0.5772 (Euler). It gives 2.23e-5 + 3.55e-5 = 5.78e-5 at 1e4 K
   ! (DERIVED; Seaton's Table 1, approximation V, prints 0.22e-4 and
   ! 0.35e-4 at 1e4 K and 0.17e-4 and 0.27e-4 at 2e4 K; Osterbrock &
   ! Ferland 2006, Table 4.10, list the same four numbers). Seaton: "the
   ! error in the approximation V results for W should not exceed +-20%".
   ! VALIDITY. The large-distance cut-off is the energy separation of the
   ! final states (the Lamb shift for 2p1/2, the fine structure for
   ! 2p3/2, as in l_mixing_2s2p_pengelly_seaton). Debye screening, which
   ! Seaton does not include and which is not applied here, would cut
   ! first where the Debye radius falls below 1.12 hbar v/dE: n_e above
   ! 5e9 cm^-3 for the 2p1/2 part and 5e11 cm^-3 for the 2p3/2 part,
   ! independent of T (DERIVED); there the electron term is overestimated
   ! logarithmically, and it is 1/8 of the proton term (which does carry
   ! the Debye cut-off) at 1e4 K. Seaton evaluates it at 1e4 and 2e4 K and treats mu as independent of
   ! energy "within the energy range of interest"; below ~5e3 K it is an
   ! extrapolation of the same closed form.
   ! THE RATE QUOTED ELSEWHERE. Christie, Arras & Li (2013, Table 2, R5)
   ! give 6.21e-5 (ln(T/1.02) - 0.577)/T^1/2 (5.35e-6 at 1e4 K), the
   ! Maxwell average of the cross section of Janev, Reiter & Samm (2003,
   ! Berichte des Forschungszentrums Juelich 4105, sect. 2.1.1 B, eq. 8,
   ! from Chibisov 1969, READ): "sigma(2s -> 2p) = 8.617/E ln(1.14 x 10^4
   ! E) (x10^-15 cm^2)", E in eV. Its logarithm is Seaton's for 2p3/2
   ! (ln(4E e^-2.21/dE'') = ln(1.08e4 E)), but its constant 8.617e-15 eV
   ! cm^2 is 1/10 of 72 pi a0^2 Ry = 8.618e-14 eV cm^2, the first-order
   ! dipole constant of the same theory, so the printed exponent appears
   ! to be one power of ten low; that rate is not used.
   real*8, intent(in) :: T
   real*8, parameter :: a_bohr = 0.529177210903d-8        ! [cm], CODATA 2018
   real*8, parameter :: Ry_erg = 2.1798723611035d-11      ! [erg], CODATA 2018
   real*8, parameter :: mu_strong = 2.21d0, gamma_euler = 0.5772156649d0
   real*8 :: kT, pref
   kT   = kb_erg*max(T, 1.0d0)
   pref = 72.0d0*pi*a_bohr**2*Ry_erg*sqrt(8.0d0*mu_H_e/pi)/(m_e*sqrt(kT))
   c2s2p_rate = pref*(                                                     &
        max(log(4.0d0*kT/dE_2s2p12) - gamma_euler - mu_strong, 0.0d0)/3.0d0 &
      + 2.0d0*max(log(4.0d0*kT/dE_2s2p32) - gamma_euler - mu_strong, 0.0d0) &
        /3.0d0)
   end function c2s2p_rate

   real*8 function l_mixing_rate_2s2p(T, ne_l, nHII_l, nHeII_l, nHeIII_l)
   ! THE 2s -> 2p TRANSFER RATE [s^-1] of one H(2s) atom: electrons
   ! (c2s2p_rate) and the ions H+, He+ and He2+ (Pengelly & Seaton 1964;
   ! Cool_coeff: l_mixing_2s2p_pengelly_seaton, where the physics and the
   ! validity are written). At 1e4 K the proton coefficient is 4.8e-4
   ! cm^3 s^-1, 8.3 times the electron one (5.8e-5): the 2s level of
   ! hydrogen is mixed into 2p faster than it decays by two photons
   ! (8.2246 s^-1) wherever n_e = n_p exceeds 1.5e4 cm^-3.
   real*8, intent(in) :: T, ne_l, nHII_l, nHeII_l, nHeIII_l
   l_mixing_rate_2s2p = ne_l*c2s2p_rate(T)                                &
      + nHII_l*l_mixing_2s2p_pengelly_seaton(T, ne_l, 1.0d0, 1.0d0,       &
                          mu_H_p, dE_2s2p12, dE_2s2p32, A_2s1s)           &
      + nHeII_l*l_mixing_2s2p_pengelly_seaton(T, ne_l, 1.0d0, 1.0d0,      &
                          mu_H_Hep, dE_2s2p12, dE_2s2p32, A_2s1s)         &
      + nHeIII_l*l_mixing_2s2p_pengelly_seaton(T, ne_l, 1.0d0, 2.0d0,     &
                          mu_H_He2p, dE_2s2p12, dE_2s2p32, A_2s1s)
   end function l_mixing_rate_2s2p

   real*8 function l_mixing_rate_2p2s(T, ne_l, nHII_l, nHeII_l, nHeIII_l)
   ! 2p->2s [s^-1], the detailed-balance reverse: the separations are
   ! 1e-5 of kT, so only the statistical weights enter.
   real*8, intent(in) :: T, ne_l, nHII_l, nHeII_l, nHeIII_l
   l_mixing_rate_2p2s = (g2s/g2p)                                         &
        *l_mixing_rate_2s2p(T, ne_l, nHII_l, nHeII_l, nHeIII_l)
   end function l_mixing_rate_2p2s

   ! --------------------------------------------------------------- !
   ! Recombination into the n=2 levels [cm^3 s^-1]. Every case-B
   ! recombination passes through n = 2, so the 2s and 2p sources sum to the
   ! case-B coefficient the ionization balance removes protons with
   ! (Cool_coeff: alpha_rec_HII_B, whichever "Atomic rate set:" selects);
   ! the n = 2 source and the balance then count the same captures into
   ! n >= 2. The ground captures that the balance adds where their photons
   ! escape (alpha_rec_HII_net) do not pass through n = 2 and are not a
   ! source here. The share into 2s is the case-B fraction of Pengelly
   ! (1964, Table I) (Cool_coeff: case_b_2s_fraction_hydrogenic), within 1%
   ! of f(2s) in Draine (2011, Table 14.3, from Brown & Mathews 1970) at
   ! 4e3, 1e4 and 2e4 K.

   real*8 function alpha_B_hydrogen(T)
   real*8, intent(in) :: T
   alpha_B_hydrogen = alpha_rec_HII_B(max(T,1.0d0))
   end function alpha_B_hydrogen

   real*8 function alpha_2s_hydrogen(T)
   real*8, intent(in) :: T
   alpha_2s_hydrogen = case_b_2s_fraction_hydrogenic(max(T,1.0d0), 1.0d0) &
                     *alpha_B_hydrogen(T)
   end function alpha_2s_hydrogen

   real*8 function alpha_2p_hydrogen(T)
   real*8, intent(in) :: T
   alpha_2p_hydrogen = alpha_B_hydrogen(T) - alpha_2s_hydrogen(T)
   end function alpha_2p_hydrogen

   ! --------------------------------------------------------------- !

   real*8 function n2p_destruction_rate(T, ne_l, nHII_l, nHeII_l, nHeIII_l, &
                                        gam_ion_2s, gam_ion_2p)
   ! Rate [s^-1] at which a 2p atom is removed WITHOUT emitting a Ly-alpha
   ! photon into the line. Three channels:
   !
   !   ne q(2p->1s)                 collisional de-excitation (10.2 eV back to
   !                                the electron gas)
   !   Gamma_2p                     photoionization of n=2 by the stellar
   !                                Balmer continuum
   !   M(2p->2s) * P_2gamma         l-mixing into 2s (electrons and ions,
   !                                l_mixing_rate_2p2s) followed by two-photon
   !                                decay, with the 2s branching ratio
   !                                P_2gamma = A_2s1s/(A_2s1s + ne q(2s->1s)
   !                                                   + M(2s->2p) + Gamma_2s)
   !
   ! Together with the escape rate beta*A_2p1s these close the 2p budget:
   ! n2p = P/(beta A_2p1s + n2p_destruction_rate). Dropping the destruction
   ! term over-estimates the trapped population wherever beta is small enough
   ! that beta*A_2p1s falls to the destruction rate, i.e. at the base.
   real*8, intent(in) :: T, ne_l, nHII_l, nHeII_l, nHeIII_l
   real*8, intent(in) :: gam_ion_2s, gam_ion_2p
   real*8 :: L2s, P_2gamma

   L2s = A_2s1s + ne_l*c2s1s_rate(T)                                        &
       + l_mixing_rate_2s2p(T, ne_l, nHII_l, nHeII_l, nHeIII_l) + gam_ion_2s
   P_2gamma = A_2s1s/max(L2s, 1.0d-30)
   n2p_destruction_rate = ne_l*c2p1s_rate(T) + gam_ion_2p                   &
        + l_mixing_rate_2p2s(T, ne_l, nHII_l, nHeII_l, nHeIII_l)*P_2gamma

   end function n2p_destruction_rate

   ! End of module
   end module hydrogen_n2_rates
