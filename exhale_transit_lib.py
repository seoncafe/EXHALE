"""Reusable, testable core for the EXHALE transit post-processor.

This module holds the pure pieces of ``EXHALE_transit.py``: physical
constants, the line rest-wavelength / oscillator-strength / Einstein-A
metadata, the self-contained physics/utility functions, the chord geometry
those functions and any census of them share, and the reader that decides
what an ``_adv`` profile says about the validity of its own rows.  Apart
from numpy/scipy and the profile readers of ``examples/exhale_io.py`` it
references only its own constants, so it imports cleanly on its own and can
be unit-tested without a simulation.  The main script imports these names
and keeps the top-down control flow (loading a run, density prep, the
loops over the lines, convolution, plotting, saving).
"""

import os
import numpy as np
from scipy.special import wofz
from astropy.convolution import convolve

# The Jupiter radius is defined once, in examples/exhale_io.py (RJ_CM), and
# imported rather than written down again, so that this tool and the run it
# analyses mean one planet by `Planet radius [R_J]`.  `__file__` is absent
# under exec(), hence the guard.
import sys
try:
    _HERE = os.path.dirname(os.path.realpath(__file__))
except NameError:
    _HERE = os.getcwd()
_EXAMPLES_DIR = os.path.join(_HERE, 'examples')
if _EXAMPLES_DIR not in sys.path:
    sys.path.append(_EXAMPLES_DIR)
from exhale_io import RJ_CM, loadtxt_cells               # noqa: E402


# ----- ENVIRONMENT OVERRIDES ----- #
# This module was renamed TPM.py -> EXHALE_transit.py.  Its run-time overrides
# are now EXHALE_TRANSIT_<NAME>, with the old TPM_<NAME> names kept as
# backward-compatible fallbacks.  `_tenv('PATH', '.')` returns
# EXHALE_TRANSIT_PATH if set, else TPM_PATH, else the default.
def _tenv(name, default=None):
    return os.environ.get('EXHALE_TRANSIT_' + name,
                          os.environ.get('TPM_' + name, default))


def _tenv_set(name):
    return ('EXHALE_TRANSIT_' + name) in os.environ or ('TPM_' + name) in os.environ


# ----- CONSTANTS ----- #

# Physicals constants
kb  = 1.380649e-23		# Boltzmann constant [J/K]
G   = 6.67430e-11		# Gravitational constant [m3/kg/s2] (CODATA 2018)
mp  = 1.672623e-27		# Proton mass [Kg]
me  = 9.109384e-31		# Electron mass [Kg]
mD  = 3.344497e-27		# Deuterium mass [Kg]
c_light = 2.99792458e8	        # Speed of light [m/s]
AU  = 1.495978707e11		# Astronomical unit [m]
E0  = 8.854188e-12		# Vacuum permittivity [F/m]
h   = 6.626070e-34		# Planck constant [J*sec]
ht  = 1.054572e-34		# Reduced Planck constant (h slash) [J*sec]
e   = -1.602176e-19		# Electron charge [C]
mHe = 6.64648157e-27		# Helium mass [Kg]
RJ  = RJ_CM*1.0e-2			# Jupiter radius [m], from the one Python
					# definition in examples/exhale_io.py, which
					# is in cm; this module works in SI
MJ  = 1.898e27				# Jupiter mass
R_sun = 6.96000000e8 	                # Sun radius [m]
M_sun = 1.989e30			# Sun mass [kg]
gray  = '#a0a0a0'			# Color gray for line plot


# ----- LINE METADATA ----- #

# Wavelengths in air for HeI metastable transitions (from NIST) [m]
l_He3_1 = 10830.33977e-10
l_He3_2 = 10830.25010e-10
l_He3_3 = 10829.09114e-10
# Wavelengths in vacuum for HI and Deuterium Ly-alpha (from NIST) [m]
lA  = 1215.6701e-10
lD = 1215.3379e-10

# Wave frequencies for HeI metastable transitions [s^-1]
nu_He3_1 = c_light/l_He3_1
nu_He3_2 = c_light/l_He3_2
nu_He3_3 = c_light/l_He3_3
# Wave frequencies for HI and D transitions [s^-1]
nu_HI = c_light/lA
nu_D = c_light/lD

# Oscillation strenghts for He
f10830_34 = 2.9958e-1
f10830_25 = 1.7974e-1
f10829_09 = 5.9902e-2
# Oscillation strenghts for HI and Deuterium
f_la = 4.1641e-1
f_D  = 4.1630e-1

# Einstein coeff. for Helium
A12_HeTR = 1.0216e7
# Einstein coeff. for HI and D
A12_HI = 4.6986e8
A12_D = 4.6999e8

# Common constant for Faddeeva integrals
Fadd_const = np.sqrt(np.pi)*e**2.0/(4.0*np.pi*E0*me*c_light)

# ----- H-alpha line + n=2 level-population atomic data ----- #
# The level balance is that of Christie, Arras & Li (2013, ApJ 772, 144,
# eqs. 12-13); its rate coefficients are the solver's (n2_populations below).

# H-alpha (n=2 -> n=3), air wavelength [m]
l_Ha   = 6562.8e-10
nu_Ha  = c_light/l_Ha
f_Ha   = 0.6407                 # multiplet oscillator strength (statistical 2s:2p = 1:3)
A12_Ha = 4.4101e7               # Einstein A(3->2) [s^-1]; line is Doppler-dominated

# H-beta (n=2 -> n=4), air wavelength [m]; same lower level (n=2) as H-alpha
l_Hb   = 4861.35e-10
nu_Hb  = c_light/l_Hb
f_Hb   = 0.11938                # multiplet oscillator strength n=2->n=4
A12_Hb = 8.4193e6               # Einstein A(4->2) [s^-1]; Doppler-dominated

# Sub-level absorption oscillator strengths. The 2s and 2p sub-levels have
# DIFFERENT absorption cross sections; the multiplet f_Ha/f_Hb above are only
# their statistical-weight (2:6) averages. Because the metastable 2s and the
# Ly-alpha-pumped 2p depart strongly from the 1:3 statistical ratio, the Balmer
# optical depth must be summed per sub-level, tau ~ f_2s n_2s + f_2p n_2p, not
# f_multiplet (n_2s + n_2p). NIST/Wiese absorption oscillator strengths:
f_Ha_2s = 0.4349      # 2s -> 3p
f_Ha_2p = 0.70941     # 2p -> 3s (0.01361) + 2p -> 3d (0.69580)
f_Hb_2s = 0.1028      # 2s -> 4p
f_Hb_2p = 0.125886    # 2p -> 4s (0.002986) + 2p -> 4d (0.12290)

# Statistical weights
g1s, g2s, g2p = 2.0, 2.0, 6.0

# Ly-alpha (1s<->2p) atomic data for the radiative pumping
nu_Lya  = c_light/lA            # lA = 1215.67 A (vacuum), defined above
A_2p1s  = 6.2649e8              # A(2p->1s) [s^-1], NIST ASD (Wiese and Fuhr 2009); one value with hydrogen_n2_rates.f90
# A(2s->1s), the two-photon decay [s^-1]: Drake (1986, Phys. Rev. A 34,
# 2871, eq. 27) for H with the finite nuclear mass, the value of
# hydrogen_n2_rates.f90 (A_2s1s), where the derivation is written.
A_2s1s  = 8.22461

# Einstein-B coefficients in the mean-intensity (J_nu) convention, in cgs
# so that B * J_lya [erg s^-1 cm^-2 Hz^-1 sr^-1] gives a rate in s^-1:
#   B21 = A21 c^2 / (2 h nu^3),   B12 = (g2/g1) B21
c_cgs, h_cgs, kb_cgs = 2.99792458e10, 6.62607015e-27, 1.380649e-16
eV2Hz = 2.417989242e14         # Hz per eV
B21_lya = A_2p1s*c_cgs**2.0/(2.0*h_cgs*nu_Lya**3.0)   # ~2.85e9
B12_lya = (g2p/g1s)*B21_lya                           # ~8.55e9 (1s->2p pump)

# Doppler / optical-depth constants for the auto-window sizing
_amu = 1.66053907e-27
_kB  = 1.380649e-23
_ec2 = 0.026540045                # pi e^2 / (m_e c)  [cm^2 Hz] (sqrt-pi form below)
_ccm = 2.99792458e10


# ----- THE H(n=2) RATE COEFFICIENTS OF THE SOLVER ----- #
# n2_populations below rebuilds the 2s/2p populations with the rate
# coefficients the solver's own level balance uses
# (src/modules/radiation/excited_hydrogen.f90: n2_rate_matrix), transcribed
# from hydrogen_n2_rates.f90 and the Cooling_Coefficients functions it reads
# (Cool_coeff.f90); each function names the Fortran it reproduces, and the
# sources and validity are written there. All in cgs, T in K. The constants
# are those of global_parameters (parameters.f90).
_KB_EV   = 8.617333262e-05        # kb_eV
_KB_ERG  = 1.380649e-16           # kb_erg
_HP_ERG  = 6.62607015e-27         # hp_erg
_C_CGS   = 2.99792458e10          # c_light
_ME_G    = 9.1093837015e-28       # m_e
_MP_G    = 1.67262192369e-24      # m_p
_MH_G    = 1.67353284e-24         # mu, the hydrogen atom
_MHE_G   = 6.6464790722e-24       # m_He_atom
_ERG2EV  = 6.241509075e11         # erg2eV
_E_TH_HI = 13.598434599           # e_th_HI [eV]
_E_HI_N2 = 10.19881               # E_HI_n2_eV, 1s-2s [eV]
# coll_rate_prefactor, h^2/((2 pi m_e)^(3/2) k^(1/2)) [cm^3 s^-1 K^(1/2)]
_COLL_PREF = _HP_ERG**2/((2.0*np.pi*_ME_G)**1.5*np.sqrt(_KB_ERG))
# Separations of 2s1/2 from 2p1/2 and 2p3/2 [erg] and the reduced masses
# of H with e, H+, He+ and He2+ [g] (hydrogen_n2_rates.f90).
_DE_2S2P12 = 0.035*_HP_ERG*_C_CGS
_DE_2S2P32 = 0.331*_HP_ERG*_C_CGS
_MU_H_E    = _ME_G*_MH_G/(_ME_G + _MH_G)
_MU_H_P    = _MP_G*_MH_G/(_MP_G + _MH_G)
_MU_H_HEP  = (_MHE_G - _ME_G)*_MH_G/(_MHE_G - _ME_G + _MH_G)
_MU_H_HE2P = (_MHE_G - 2.0*_ME_G)*_MH_G/(_MHE_G - 2.0*_ME_G + _MH_G)
# The 12 positive nodes and weights of the 24-point Gauss-Legendre rule of
# ground_capture_milne_moments.
_GL24_X = np.array([6.40568928626056300e-02, 1.91118867473616311e-01,
                    3.15042679696163397e-01, 4.33793507626045127e-01,
                    5.45421471388839563e-01, 6.48093651936975546e-01,
                    7.40124191578554358e-01, 8.20001985973902947e-01,
                    8.86415527004400960e-01, 9.38274552002732798e-01,
                    9.74728555971309474e-01, 9.95187219997021311e-01])
_GL24_W = np.array([1.27938195346752215e-01, 1.25837456346828303e-01,
                    1.21670472927803419e-01, 1.15505668053725613e-01,
                    1.07444270115965607e-01, 9.76186521041140648e-02,
                    8.61901615319532882e-02, 7.33464814110804109e-02,
                    5.92985849154367417e-02, 4.42774388174195510e-02,
                    2.85313886289337432e-02, 1.23412297999870909e-02])


def hydrogen_ground_capture(T):
	# alpha_1 of H II -> H I [cm^3 s^-1], the Milne relation on the
	# hydrogenic ground-state cross section of the transfer (Cool_coeff:
	# ground_capture_milne, ion 1, and cross_sec.f90: sigma). The solver
	# interpolates a table of this quadrature (capture_table_value); the
	# two agree to the table's interpolation error.
	T   = np.maximum(np.asarray(T, dtype=float), 1.0)
	kT  = _KB_EV*T
	I   = _E_TH_HI
	tmx = np.log(1.0 + 50.0*kT/I)
	acc = np.zeros_like(T)
	for x, w in zip(_GL24_X, _GL24_W):
		for sgn in (-1.0, 1.0):
			E   = I*np.exp(0.5*tmx*(1.0 + sgn*x))
			eps = np.sqrt(E/I - 1.0)
			sig = 6.3*(I/E)**4*np.exp(4.0 - 4.0*np.arctan(eps)/eps)            \
			      /(1.0 - np.exp(-2.0*np.pi/eps))              # [1e-18 cm^2]
			acc = acc + w*E**3*sig*np.exp(-(E - I)/kT)
	return (2.0*np.sqrt(2.0/np.pi)/(_C_CGS**2*(_ME_G*kT/_ERG2EV)**1.5)
	        *0.5*tmx*acc*1.0e-18/_ERG2EV**3)


def alpha_B_hydrogen(T, rate_set='default'):
	# Case-B recombination coefficient of H II [cm^3 s^-1], the one the
	# ionization balance removes protons with (Cool_coeff: alpha_rec_HII_B):
	#   'default'      Badnell (2006, 2023 update) total minus the Milne
	#                  ground capture (alphaB_HII_new);
	#   'koskinen2022' Koskinen et al. (2022) Table 1 R1 ("Atomic rate
	#                  set: Koskinen2022");
	#   'legacy'       Hui & Gnedin (1997) ("Legacy_HHe_rates: True").
	T = np.maximum(np.asarray(T, dtype=float), 1.0)
	if rate_set == 'koskinen2022':
		return 4.0e-12*(300.0/T)**0.64
	if rate_set == 'legacy':
		xl = 2.0*157807.0/T
		return 2.753e-14*xl**1.5/(1.0 + (xl/2.740)**0.407)**2.242
	tt  = np.sqrt(T/2.965)
	rr  = 8.318e-11/(tt*(1.0 + tt)**(1.0 - 0.7472)
	                 *(1.0 + np.sqrt(T/7.001e5))**(1.0 + 0.7472))
	return rr - hydrogen_ground_capture(T)


def case_b_2s_fraction(T, Z=1.0):
	# Share of the case-B captures of a hydrogenic ion that end in 2s,
	# from Pengelly (1964, MNRAS 127, 145, Table I), interpolated linearly
	# in log t' of the logarithms of the two coefficients and held at the
	# table ends (Cool_coeff: case_b_2s_fraction_hydrogenic).
	a2s = np.array([30.6, 20.4, 13.3, 8.37, 5.07, 2.93, 1.61])
	a2p = np.array([97.5, 56.8, 32.1, 17.6, 9.27, 4.68, 2.28])
	T   = np.maximum(np.asarray(T, dtype=float), 1.0)
	pos = np.log(1.0e-4*T/(Z*Z))/np.log(2.0) + 4.0
	k   = np.clip(np.floor(pos), 1, 6).astype(int)
	f   = np.where(pos > 1.0, np.where(pos >= 7.0, 1.0, pos - k), 0.0)
	k   = np.where(pos >= 7.0, 6, np.where(pos > 1.0, k, 1))
	l2s = np.exp(np.log(a2s[k-1]) + f*(np.log(a2s[k]) - np.log(a2s[k-1])))
	l2p = np.exp(np.log(a2p[k-1]) + f*(np.log(a2p[k]) - np.log(a2p[k-1])))
	return l2s/(l2s + l2p)


def _upsilon_quartic_logT(c, T):
	# Cool_coeff: upsilon_quartic_logT, T held inside 1e3-1e5 K.
	x = np.log10(np.clip(T, 1.0e3, 1.0e5)/1.0e4)
	return 10.0**(c[0] + x*(c[1] + x*(c[2] + x*(c[3] + x*c[4]))))


# Cool_coeff: upsilon_HI_1s2s, upsilon_HI_1s2p (CHIANTI v11.0.2 h_1, the
# rate set of the H I cooling as well).
_UPS_1S2S = (-5.41913497e-01, 1.45680880e-01, -5.14786365e-02,
             -4.89245376e-02, 4.41143616e-02)
_UPS_1S2P = (-3.02352097e-01, 3.07054967e-01, 2.18544412e-01,
             2.79631531e-02, -3.12591467e-02)


def hydrogen_n2_collision_rates(T):
	# (C1s2s, C1s2p, C2s1s, C2p1s) [cm^3 s^-1]: Cool_coeff
	# excitation_rate_HI_1s2s/_1s2p and their detailed-balance reverses
	# deexcitation_rate_HI_2s1s/_2p1s with the Boltzmann factor cancelled.
	T  = np.maximum(np.asarray(T, dtype=float), 1.0)
	u2s = _upsilon_quartic_logT(_UPS_1S2S, T)
	u2p = _upsilon_quartic_logT(_UPS_1S2P, T)
	bz  = np.exp(-_E_HI_N2/(_KB_EV*T))
	pre = _COLL_PREF/np.sqrt(T)
	return (pre/2.0*u2s*bz, pre/2.0*u2p*bz, pre/2.0*u2s, pre/6.0*u2p)


def l_mixing_2s2p_electron(T):
	# 2s -> 2p by electron impact [cm^3 s^-1], Seaton (1955, Proc. Phys.
	# Soc. A 68, 457, eq. 55, approximation V), Maxwell-averaged in closed
	# form (hydrogen_n2_rates: c2s2p_rate). 5.78e-5 at 1e4 K.
	a_bohr = 0.529177210903e-8
	Ry_erg = 2.1798723611035e-11
	kT  = _KB_ERG*np.maximum(np.asarray(T, dtype=float), 1.0)
	pre = 72.0*np.pi*a_bohr**2*Ry_erg*np.sqrt(8.0*_MU_H_E/np.pi)/(_ME_G*np.sqrt(kT))
	g, m = 0.5772156649, 2.21
	return pre*(np.maximum(np.log(4.0*kT/_DE_2S2P12) - g - m, 0.0)/3.0
	            + 2.0*np.maximum(np.log(4.0*kT/_DE_2S2P32) - g - m, 0.0)/3.0)


def _lower_gamma_3(x):
	with np.errstate(over='ignore', invalid='ignore'):
		return np.where(x < 1.0e-2, x**3*(1.0/3.0 - x/4.0 + x*x/10.0),
		                2.0 - np.exp(-x)*(x*x + 2.0*x + 2.0))


def _lower_gamma_2(x):
	with np.errstate(over='ignore', invalid='ignore'):
		return np.where(x < 1.0e-2, x**2*(0.5 - x/3.0 + x*x/8.0),
		                1.0 - np.exp(-x)*(1.0 + x))


def l_mixing_2s2p_ion(T, n_e, z, Z_c, mu_g, dE1, dE2, A_2q):
	# 2s -> 2p by a charged heavy particle [cm^3 s^-1], Pengelly & Seaton
	# (1964, MNRAS 127, 165, eqs. 34-41) with their cut-off
	# min(1.12 hbar v/dE, 0.72 v tau, R_D), Maxwell-averaged in closed form
	# (Cool_coeff: l_mixing_2s2p_pengelly_seaton, where the physics and the
	# validity are written).
	from scipy.special import exp1
	e_esu  = 1.602176634e-19*_C_CGS/10.0
	hbar   = _HP_ERG/(2.0*np.pi)
	a_bohr = hbar**2/(_ME_G*e_esu**2)
	Tk   = np.maximum(np.asarray(T, dtype=float), 1.0)
	R_D  = np.sqrt(_KB_ERG*Tk/(4.0*np.pi*np.maximum(n_e, 1.0e-30)*e_esu**2))
	K_1  = np.sqrt(6.0*(Z_c/z)**2*12.0)*e_esu**2/hbar*a_bohr
	v_th = np.sqrt(2.0*_KB_ERG*Tk/mu_g)
	q = np.zeros(np.broadcast(Tk, R_D).shape)
	with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
		for dE, wt in ((dE1, 1.0/3.0), (dE2, 2.0/3.0)):
			a_c = min(1.12*hbar/dE, 0.72/A_2q)
			x_1 = K_1/(a_c*v_th**2)
			x_a = R_D**2/(a_c*v_th)**2
			x_D = (K_1/(R_D*v_th))**2
			B0  = 0.5 + np.log(R_D*v_th/K_1)
			A0  = 0.5 + np.log(a_c*v_th**2/K_1)
			near = (0.5*np.pi*a_c**2*v_th**3*_lower_gamma_3(x_1)
			        + np.pi*K_1**2/v_th*(A0*(np.exp(-x_1) - np.exp(-x_a))
			          + np.log(x_1)*np.exp(-x_1) - np.log(x_a)*np.exp(-x_a)
			          + exp1(x_1) - exp1(x_a))
			        + np.pi*K_1**2/v_th*(B0*np.exp(-x_a)
			          + 0.5*(np.log(x_a)*np.exp(-x_a) + exp1(x_a))))
			far  = (0.5*np.pi*a_c**2*v_th**3*_lower_gamma_3(x_a)
			        + 0.5*np.pi*R_D**2*v_th*(_lower_gamma_2(x_D) - _lower_gamma_2(x_a))
			        + np.pi*K_1**2/v_th*(B0*np.exp(-x_D)
			          + 0.5*(np.log(x_D)*np.exp(-x_D) + exp1(x_D))))
			q = q + wt*np.where(x_1 < x_a, near, far)
	return 2.0/np.sqrt(np.pi)*q


def l_mixing_rate_2s2p(T, ne, nHII, nHeII, nHeIII):
	# THE 2s -> 2p TRANSFER RATE [s^-1] of one H(2s) atom: electrons and
	# the ions H+, He+ and He2+ (hydrogen_n2_rates: l_mixing_rate_2s2p).
	# The reverse, 2p -> 2s, is (g2s/g2p) times it.
	ion = lambda Zc, mu: l_mixing_2s2p_ion(T, ne, 1.0, Zc, mu, _DE_2S2P12,
	                                       _DE_2S2P32, A_2s1s)
	return (ne*l_mixing_2s2p_electron(T) + nHII*ion(1.0, _MU_H_P)
	        + nHeII*ion(1.0, _MU_H_HEP) + nHeIII*ion(2.0, _MU_H_HE2P))


def dissociative_recombination_n2_source(T, nH2p, nHeHp, ne):
	# Chemical production of H(n=2) [cm^-3 s^-1] by H2+ + e (Koskinen et
	# al. 2022, R5) and HeH+ + e (R16), each leaving one H atom in n = 2
	# (molecular_reaction_heat.f90: dissociative_recombination_n2_source,
	# mol_rates.f90: rk_R5_H2p_dr, rk_R16_HeHp_dr); zero, as there, when
	# EXHALE_REACTION_HEAT_RECIPIENTS=0.
	if os.environ.get('EXHALE_REACTION_HEAT_RECIPIENTS', '').strip() == '0':
		return np.zeros_like(np.asarray(T, dtype=float))
	T = np.maximum(np.asarray(T, dtype=float), 1.0)
	return ((2.3e-8*(300.0/T)**0.4*np.maximum(nH2p, 0.0)
	         + 1.0e-8*(300.0/T)**0.6*np.maximum(nHeHp, 0.0))*np.maximum(ne, 0.0))


def n2_populations(T, n1s, nHII, nHeII, nHeIII, ne, Jlya, G2s=0.0, G2p=0.0,
                   nH2p=0.0, nHeHp=0.0, rate_set='default'):
	# The 2s/2p statistical equilibrium (Christie, Arras & Li 2013, eqs.
	# 12-13) as the solver writes it (excited_hydrogen.f90: n2_rate_matrix
	# and n2_populations), with the solver's rate coefficients (above):
	# A(2s) of Drake (1986), the CHIANTI 1s-2s/2p collision rates of the
	# H I cooling, the l-mixing by electrons (Seaton 1955) and by H+, He+
	# and He2+ (Pengelly & Seaton 1964), the balance's own case-B
	# coefficient split with the 2s share of Pengelly (1964), and the
	# chemical n = 2 source of a molecular run. Densities [cm^-3] (negative
	# entries are clamped to zero, as the solver does), T [K], Jlya [cgs],
	# G2s, G2p the Balmer-continuum photoionization rates [s^-1], rate_set
	# as for alpha_B_hydrogen. Returns (n2s, n2p, n2s + n2p).
	#
	# WHAT THIS DOES NOT SHARE WITH THE SOLVER is the field: J_lya and
	# G2s/G2p are the caller's (EXHALE_transit.py: a J_lya file or the
	# Huang et al. 2017 estimate, and a diluted blackbody Balmer continuum),
	# because the solver writes neither to the profiles this tool reads.
	T   = np.maximum(np.asarray(T, dtype=float), 1.0)
	n1s, nHII, nHeII, nHeIII, ne = (np.maximum(q, 0.0) for q in
	                                (n1s, nHII, nHeII, nHeIII, ne))
	aB  = alpha_B_hydrogen(T, rate_set)
	a2s = case_b_2s_fraction(T, 1.0)*aB
	a2p = aB - a2s
	C1s2s, C1s2p, C2s1s, C2p1s = hydrogen_n2_collision_rates(T)
	M12 = l_mixing_rate_2s2p(T, ne, nHII, nHeII, nHeIII)   # 2s -> 2p [s^-1]
	M21 = (g2s/g2p)*M12                                     # 2p -> 2s [s^-1]
	chem = dissociative_recombination_n2_source(T, nH2p, nHeHp, ne)
	# Ly-alpha radiative pump (1s->2p) and stimulated emission (2p->1s)
	Ppump = B12_lya*Jlya
	Pstim = B21_lya*Jlya
	# 2x2 linear system  [[L2p, -M12],[-M21, L2s]] [n2p,n2s]^T = [S2p,S2s]^T
	L2p = A_2p1s + Pstim + C2p1s*ne + M21 + G2p
	L2s = C2s1s*ne + M12 + G2s + A_2s1s
	# Cascade source: the electron recombines onto a proton, so the rate is
	# alpha_2l*ne*nHII. ne and nHII part company wherever helium and metals
	# supply the electrons while hydrogen is still neutral, i.e. at the base.
	S2p = (Ppump + C1s2p*ne)*n1s + a2p*ne*nHII + chem*g2p/(g2s + g2p)
	S2s = (C1s2s*ne)*n1s + a2s*ne*nHII + chem*g2s/(g2s + g2p)
	det = L2p*L2s - M12*M21
	det = np.where(np.abs(det) > 0.0, det, 1.0)   # guard (det>0 physically)
	n2p = np.maximum((S2p*L2s + M12*S2s)/det, 0.0)
	n2s = np.maximum((L2p*S2s + M21*S2p)/det, 0.0)
	return n2s, n2p, n2s + n2p


def gamma_n2_balmer(T_star, R_over_a):
	# n=2 photoionization rate [s^-1] from a diluted stellar blackbody
	# Balmer continuum, over the band the solver integrates, from the n=2
	# edge e_th_HI/4 = 3.3996 eV (3647 A) to the H I edge 13.5984 eV
	# (J_incident: e_th_HI_n2; excited_hydrogen: balmer_band_integrals).
	# Hydrogenic cross section sigma_2(nu) = sigma2_th*(nu2/nu)^3, one value
	# for 2s and 2p, the solver's s2_thr. R_over_a = R_star/a (dimensionless).
	# For a stellar beam, the rate is int F_nu/(h nu)*sigma dnu with
	# F_nu = pi B_nu(T_star)*(R_star/a)^2.
	#
	# THE FIELD DIFFERS FROM THE SOLVER'S: the solver integrates the run's
	# own spectrum type (the loaded table, the power law or the Planck
	# field), and neither that spectrum nor the rate is in the profiles this
	# tool reads, so a blackbody at T_star stands in for it here.
	#
	# sigma2_th = 1.4e-17 cm^2 is the n = 2 threshold value of the
	# Kramers cross section with the bound-free Gaunt factor of Seaton
	# (1959, MNRAS 119, 81, eqs. 3 and 10, READ): 2^6 alpha pi a_0^2 n /
	# (3 sqrt(3) Z^2) g_II(n, 0) = 1.5814e-17 x 0.8715 = 1.378e-17 cm^2
	# for n = 2 (DERIVED), rounded. The exact nonrelativistic hydrogenic
	# values at the edge are 1.478e-17 (2s) and 1.355e-17 (2p), whose
	# statistical average (weights 2 : 6) is 1.386e-17 (DERIVED by direct
	# integration of the bound-free dipole matrix elements, reproducing
	# 6.304e-18 for 1s). The (nu2/nu)^3 dependence with a constant
	# threshold value leaves out the rise of g_II above the edge.
	if T_star <= 0.0:
		return 0.0
	E2   = 13.598434599/4.0              # eV, n=2 ionization threshold (e_th_HI_n2)
	E1   = 13.598434599                  # eV, H I edge (e_th_HI); above it H(1s) absorbs
	nu2  = E2*eV2Hz
	s2th = 1.4e-17                       # cm^2, sigma at the n=2 threshold (s2_thr)
	Eg   = np.linspace(E2, E1, 400)      # eV
	nu   = Eg*eV2Hz
	Bnu  = (2.0*h_cgs*nu**3.0/c_cgs**2.0)/(np.exp(h_cgs*nu/(kb_cgs*T_star)) - 1.0)
	Fnu  = np.pi*Bnu*R_over_a**2.0       # flux at planet [erg s^-1 cm^-2 Hz^-1]
	sig2 = s2th*(nu2/nu)**3.0            # cm^2
	return np.trapz(Fnu/(h_cgs*nu)*sig2, nu)   # [s^-1]


# ----- USER DEFINED FUNCTIONS ----- #

# Read word_number-th word from string_in
def get_word(string_in,word_number):

	# Initialize counters and string
	count = 0
	c_string = ''
	c_word_counter = 0
	string = string_in.strip()

	# Keep reading through string
	while count >= 0 and count <= len(string):

		if count == len(string):
			out_word = c_string
			c_word_counter += 1
			# Return word
			if c_word_counter == word_number:
				return out_word

		# If blank space or at the end of the string
		if string[count:count+1] != ' ':

			c_string += string[count:count+1]
			count += 1
		else:

			if string[count-1:count] != ' '  or count == len(string):
				# Update counters and get word
				out_word = c_string
				c_word_counter += 1
				# Reset reading string
				c_string = ''
				# Return word
				# Update counter
				count += 1
				if c_word_counter == word_number:
					return out_word
			else:
				count += 1
				continue


# ----- input.inp label matching (mirrors input_read.f90) ----- #

def _is_sep(c):
	# Value separator after a label: ':', '?', whitespace, tab or '='
	# (same set as input_read.f90 is_sep).
	return c in (':', '?', ' ', '\t', '=')


def lbl_match(line, key):
	"""True iff `line` (after lstrip) begins with `key` followed by a value
	separator (':', '?', whitespace, '=') or the end of the line. This mirrors
	the anchored label match in input_read.f90 (lbl_match), so line order is
	irrelevant and a substring buried mid-line does not false-match."""
	t = line.rstrip('\n').lstrip()
	if not t.startswith(key):
		return False
	rest = t[len(key):]
	if rest == '':
		return True
	return _is_sep(rest[0])


def find_input_label(lines, key):
	"""Return the LAST line in `lines` matching `key` as an anchored label, or
	None. Last-occurrence-wins matches input_read.f90 find_lbl / the keyword
	loop, which keep the last match."""
	found = None
	for ln in lines:
		if lbl_match(ln, key):
			found = ln
	return found


def read_input_params(path):
	"""Read the header parameters EXHALE_transit needs from input.inp by LABEL
	(order-independent), matching the input_read.f90 core block. Rp/Mp/T0/
	a_orb/Mstar are taken from the labeled lines the Fortran uses, at the same
	word positions as the legacy positional reader; LEUV and the 2D approximate
	method stay content-matched. The stellar radius and effective temperature
	are read from the same "Stellar radius [R_sun]:" / "Stellar Teff [K]:"
	lines the Fortran core uses (word 4 in both cases), and are None when the
	line is absent -- they are optional in input.inp. Returns a dict:
	  Rp [m], Mp [kg], T0 [K], a_orb [m], Mstar [kg], LEUV, appx_mth,
	  R_star_Rsun [R_sun or None], T_star [K or None], h_rate_set ('default',
	  'legacy' or 'koskinen2022', the H II recombination set of the run)."""
	with open(path, 'r') as f:
		lines = f.readlines()

	def word(key, n):
		ln = find_input_label(lines, key)
		if ln is None:
			raise ValueError('input.inp: mandatory line "%s" not found in %s'
			                 % (key, path))
		return get_word(ln, n)

	def optional_word(key, n):
		ln = find_input_label(lines, key)
		if ln is None:
			return None
		try:
			return float(get_word(ln, n))
		except (TypeError, ValueError):
			return None

	Rp    = float(word('Planet radius', 4)) * RJ
	Mp    = float(word('Planet mass', 4)) * MJ
	T0    = float(word('Equilibrium temperature', 4))
	a_orb = float(word('Orbital distance', 4)) * AU
	Mstar = float(word('Parent star mass', 5)) * M_sun

	# Stellar radius / effective temperature (same labels and word positions as
	# input_read.f90). A non-positive entry means "not specified" there too.
	R_star_Rsun = optional_word('Stellar radius', 4)
	T_star      = optional_word('Stellar Teff', 4)
	if R_star_Rsun is not None and R_star_Rsun <= 0.0:
		R_star_Rsun = None
	if T_star is not None and T_star <= 0.0:
		T_star = None

	LEUV = None
	appx_mth = ''
	leuv_ln = find_input_label(lines, 'Log10 of EUV luminosity')
	if leuv_ln is not None:
		LEUV = float(leuv_ln.split(':')[-1])
	appx_ln = find_input_label(lines, '2D approximate method')
	if appx_ln is not None:
		appx_mth = appx_ln.split(':')[-1].strip()

	# The H II recombination coefficient of the run (input_read.f90, the
	# same labels and word positions): "Atomic rate set: Koskinen2022"
	# takes precedence over "Legacy_HHe_rates: True", as in Cool_coeff:
	# alpha_rec_HII_B. n2_populations splits that coefficient into 2s, 2p.
	h_rate_set = 'default'
	legacy_ln = find_input_label(lines, 'Legacy_HHe_rates')
	if legacy_ln is not None and get_word(legacy_ln, 2) in ('True', 'true'):
		h_rate_set = 'legacy'
	set_ln = find_input_label(lines, 'Atomic rate set')
	if set_ln is not None and get_word(set_ln, 4) in ('Koskinen2022',
	                                                  'koskinen2022'):
		h_rate_set = 'koskinen2022'

	params = dict(Rp=Rp, Mp=Mp, T0=T0, a_orb=a_orb, Mstar=Mstar,
	              LEUV=LEUV, appx_mth=appx_mth,
	              R_star_Rsun=R_star_Rsun, T_star=T_star,
	              h_rate_set=h_rate_set, resolved=False)

	# Prefer the wind solver's resolved configuration when present.
	# EXHALE_resolved.out is written by write_setup_report.f90
	# (write_resolved_config) AFTER the base.inp handoff overrides are
	# applied, so it is the radius/temperature/He ratio the wind actually
	# used -- input.inp alone is stale whenever base.inp overrides them.
	resolved_path = os.path.join(os.path.dirname(os.path.abspath(path)),
	                             'EXHALE_resolved.out')
	if os.path.isfile(resolved_path):
		rv = {}
		with open(resolved_path, 'r') as f:
			for ln in f:
				ln = ln.strip()
				if not ln or ln.startswith('#'):
					continue
				parts = ln.split()
				if len(parts) >= 2:
					rv[parts[0]] = parts[1]
		try:
			params['Rp']    = float(rv['planet_radius_RJ']) * RJ
			params['Mp']    = float(rv['planet_mass_MJ']) * MJ
			params['T0']    = float(rv['equilibrium_temperature_K'])
			params['a_orb'] = float(rv['orbital_distance_AU']) * AU
			params['Mstar'] = float(rv['star_mass_Msun']) * M_sun
			params['HeH']   = float(rv['HeH_number_ratio'])
			params['resolved'] = True
		except (KeyError, ValueError) as err:
			raise ValueError('EXHALE_resolved.out present but unreadable '
			                 '(%s); refusing to silently fall back to '
			                 'input.inp' % err)

	return params


def orbital_period_days(a_orb, Mstar, Mp):
	"""Keplerian orbital period [days] of the two-body system,
	P = 2 pi sqrt(a^3 / (G (M_star + M_p))).  a_orb [m], Mstar and Mp [kg].
	Used as the planet spin period under the tidal-locking assumption."""
	return 2.0*np.pi*np.sqrt(a_orb**3.0/(G*(Mstar + Mp)))/86400.0


def parameter_with_source(env_name, input_value, default,
                          input_label='input.inp',
                          default_label='built-in default'):
	"""Resolve one run parameter and report where it came from.

	Precedence: EXHALE_TRANSIT_<env_name> (explicit run-time override, with
	the legacy TPM_<env_name> spelling) > the value read from ./input.inp >
	the built-in default.  Returns (value, source_string)."""
	raw = _tenv(env_name)
	if raw is not None:
		src = 'EXHALE_TRANSIT_' + env_name
		if ('EXHALE_TRANSIT_' + env_name) not in os.environ:
			src = 'TPM_' + env_name
		return float(raw), 'env ' + src
	if input_value is not None:
		return float(input_value), input_label
	return float(default), default_label


def _line_halfwidth(lam0_A, m_atom_amu, fosc, A21, N_col_cm2, T, v):
    """Return the physical line half-width [m/s].

    ``T`` and ``v`` are the temperature and velocity profiles of the run
    (v in cm/s, as in the EXHALE output); passing them explicitly keeps the
    function pure so it can be tested without the module globals."""
    T_max = float(np.nanmax(T))
    v_max = float(np.nanmax(np.abs(v)))*1.0e-2                  # m/s
    v_th  = np.sqrt(2.0*_kB*T_max/(m_atom_amu*_amu))           # m/s
    v_rot = 2.0e4                                              # ~20 km/s
    dv_kin = v_max + 5.0*v_th + v_rot
    # damping-wing extent from the vertical line-center optical depth
    lam0_cm = lam0_A*1.0e-8
    nu0     = _ccm/lam0_cm
    dnuD    = nu0*(v_th*1.0e2)/_ccm                            # Doppler width [Hz]
    sig0    = np.sqrt(np.pi)*(4.803204e-10**2/(9.109384e-28*_ccm))*fosc/max(dnuD,1e-30)
    tau0    = sig0*max(N_col_cm2, 0.0)
    avoigt  = (A21/(4.0*np.pi))/max(dnuD, 1e-30)
    x_w     = np.sqrt(max(tau0*avoigt/np.sqrt(np.pi), 0.0))    # Doppler widths to tau=1
    dv_wing = v_max + x_w*v_th
    # A very optically thick line (Lya, tau0 ~ 1e8) is black far into its
    # 1/x^2 damping wings; capping the wing extent keeps the window generous
    # enough to show the profile turning over without an absurd velocity span.
    DV_CAP = 2.5e6                                             # 2500 km/s
    return 1.25*min(max(dv_kin, dv_wing), DV_CAP)


def _apply_window(lam0_A, dv_half, lmin0, lmax0, forced=False):
    if forced:
        return lmin0, lmax0
    dlam = lam0_A*dv_half/2.99792458e8
    return min(lmin0, lam0_A - dlam), max(lmax0, lam0_A + dlam)


# adequate, ODD sampling after widening (astropy convolution needs odd kernels)
def _odd(n):
    n = int(n); return n if n % 2 == 1 else n+1


# Optional turbulence broadening, opt-in through EXHALE_TRANSIT_TURB=1.
# Same definition p_winds.transit uses (after Lampon et al. 2020):
# v_turb = sqrt(5/6 kT/m), added in quadrature to the thermal width. Off by
# default, so every spectrum synthesized before this option is unchanged.
def _turb_factor():
	"""Multiplier on v_th when turbulence broadening is enabled.

	p_winds widens its Gaussian as sqrt(kT/m + v_turb^2) with
	v_turb^2 = 5/6 kT/m, i.e. kT/m -> (11/6) kT/m.  This module carries the
	Doppler b-parameter v_th = sqrt(2kT/m) instead, so the same physics is
	b -> sqrt(2*(11/6) kT/m) = sqrt(11/6) * v_th.
	"""
	if _tenv('TURB', '0') not in ('1', 'True', 'true'):
		return 1.0
	return (11.0/6.0)**0.5

# ----- CHORD GEOMETRY: THE ONE SHELL SELECTION ----- #
# The selection below is shared by every line-of-sight integral in this module
# and by the census of the material those integrals sample, in
# EXHALE_transit.py.  Two selections would let the diagnostic describe a
# different atmosphere from the one the spectrum was built from.
def chord_shell_indices(shell_r, b):
	"""Rows of a radial profile that the chord at impact parameter b crosses.

	The chord meets every shell of radius r >= b, at line-of-sight coordinate
	|x| = sqrt(r^2 - b^2); there is no upper radial limit.  A shell at r = 12
	crosses the ray at b = 2 at |x| = sqrt(140) whether or not the ray grid
	stops at b = 10, so a shell lying outside the stellar disk still absorbs on
	the rays that cross the disk.  Any radial census that imposes an upper limit
	therefore describes less material than the spectrum contains.

	shell_r  row radii.  The signed, mirrored (night + day) chord array of the
	         integrals is accepted unchanged: only |r| enters.
	b        impact parameter, in the unit of shell_r.
	"""
	return np.where(np.abs(shell_r) >= b)[0]


def refused_line_center_tau_share(lines, r_grid, Rp, data_r, data_v, data_T,
                                  data_refused):
	"""Share of each ray's line-center optical depth carried by flagged rows.

	The spectra integrate the Voigt absorption coefficient along the chord with
	the trapezoid rule, so at one wavelength the optical depth of the ray at
	impact parameter b is exactly tau(b) = sum_k w_k I_k, with I_k the integrand
	of row k and w_k that row's trapezoid weight.  The sum is additive over
	rows, so a subset of the rows carries the exact share
	sum_{k in subset} w_k I_k / tau(b).  Here the subset is the rows whose
	steady advective correction was refused (`adv_T_status != 0`), which carry
	the run's own temperature and its equilibrium composition instead.

	The result is a CONTRIBUTION, not an uncertainty on transit depth: the
	disk-averaged transmission is an average of exp(-tau) and is nonlinear in
	tau, so a saturated ray moves little when its optical depth changes, and a
	share near one does not mean the depth is wrong by a comparable amount.

	lines         [(label, lam_eval_A, components), ...] with
	              components = [(lam0_A, f_osc, A21, mass, n_lower), ...],
	              n_lower in m^-3 on the mirrored chord array, exactly as
	              resonance_spectrum receives it.  The optical depth is
	              evaluated at lam_eval_A with every component of the line
	              contributing there, since components of one multiplet overlap.
	data_refused  boolean array on the mirrored chord array, True where the row
	              is not the steady advective correction.
	returns       {label: (tau, share)}, both arrays over the ray grid r_grid.
	              share is zero where the ray carries no optical depth.
	"""
	out = {}
	for (label, lam_eval_A, components) in lines:
		nu_eval = c_light/(lam_eval_A*1e-10)
		tau     = np.zeros(len(r_grid))
		tau_ref = np.zeros(len(r_grid))
		for p in range(len(r_grid)):
			b     = r_grid[p]*Rp
			arg   = chord_shell_indices(data_r, b)
			if arg.size < 2:
				continue
			r_LOS = data_r[arg]
			x_LOS = np.sqrt(r_LOS**2.0 - b**2.0)*np.sign(r_LOS)
			dx    = np.abs(x_LOS[1:] - x_LOS[:-1])
			v_x   = x_LOS*data_v[arg]/r_LOS
			I     = np.zeros(arg.size)
			for (lam0_A, f_osc, A21, mass, n_lower) in components:
				nu0  = c_light/(lam0_A*1e-10)
				v_th = np.sqrt(2.0*kb*data_T[arg]/mass)*_turb_factor()
				Dnu  = nu0*v_th/c_light
				a_v  = A21/(4.0*np.pi*Dnu)
				X    = (nu_eval - nu0)/Dnu
				I   += n_lower[arg]*f_osc*Fadd_const/Dnu \
				       * wofz(X - v_x/v_th + 1j*a_v).real
			# The trapezoid weight of each row, so that the quadrature the
			# spectrum uses can be split over a subset of the rows:
			# sum_k w_k I_k == sum_k dx_k (I_k + I_{k+1})/2.
			w = np.zeros(arg.size)
			w[:-1] += dx/2.0
			w[1:]  += dx/2.0
			sel = np.asarray(data_refused, dtype=bool)[arg]
			tau[p]     = np.sum(w*I)
			tau_ref[p] = np.sum(w[sel]*I[sel])
		share = np.divide(tau_ref, tau, out=np.zeros_like(tau),
		                  where=(tau > 0.0))
		out[label] = (tau, share)
	return out


# ----- THE VALIDITY THE SPECTRUM INHERITS ----- #
# A transit spectrum is an integral over the rows of an `_adv` profile, so it
# is worth exactly what those rows are worth.  The `_adv` writer states, row by
# row, whether the steady advective correction was adopted, and the schema of
# that statement is versioned in the file's own header.  The reader below is
# the only place this tool decides what a profile says about itself, and the
# metadata writer below it is the only place that statement is passed on to a
# saved spectrum.

TRANSIT_SCHEMA = 1                # version of the comment block that
                                  # transit_metadata_block writes into every
                                  # saved tpm_ file

# Row validity, schema 2: two fields, one for temperature and one for
# composition, so a row may carry a corrected temperature on a retained
# composition.
ADV_STATUS_NAMES_2 = ('corrected', 'retained', 'failed', 'unsupported',
                      'not_evaluated')
# Row validity, schema 1: one field, whose values mix the two.  The physics of
# each value is the writer's own legend.
ADV_STATUS_NAMES_1 = ('corrected',
                      'T kept: local radiative balance sets it',
                      'T and composition kept: the mass flux is not stationary',
                      'composition kept: an ionization validity condition',
                      'T kept: the cell solve did not converge')


def header_comment_statements(path):
	"""Comment lines at the head of a profile file as (continuation, text).

	The writer states one thing per '# <key> ...' line and wraps the rest of
	that statement over further lines it indents by more than one space
	('#   ...'), which carry no key of their own.  A reader that takes each
	line separately gets a sentence cut in half, so the indentation is what
	says which lines belong to the line above.
	"""
	out = []
	with open(path) as fh:
		for row in fh:
			s = row.strip()
			if not s:
				continue
			if not s.startswith('#'):
				break
			body = s[1:]
			out.append((body.startswith('  '), body.strip()))
	return out


def first_data_row_ncol(path):
	"""Fields in the first data row of a profile file (0 if there is none)."""
	with open(path) as fh:
		for row in fh:
			s = row.strip()
			if s and not s.startswith('#'):
				return len(s.split())
	return 0


def read_adv_validity(path):
	"""What an `_adv` profile says about the validity of its own rows.

	Three cases, decided by the file's header and by nothing else:

	schema 2  the two integer columns `adv_T_status` and `adv_comp_status`,
	          each 0 corrected, 1 retained, 2 failed, 3 unsupported,
	          4 not_evaluated, and the `adv_mass_row` column, the measure both
	          were decided by.  A row is refused for the census when its
	          temperature field is not zero; its composition field is reported
	          beside it.  A corrected row is a CONDITIONAL correction accurate
	          to `adv_conditional_tol` of itself in the mass flux, and both the
	          fraction and the row's own measure come back so that the spectrum
	          carries the condition of the rows it stands on.
	schema 1  the single `adv_status` column, one field mixing temperature and
	          composition; a row is refused when it is not zero.
	schema 0  neither an `adv_schema` line nor a status column: the validity of
	          every row is UNKNOWN.  It is not 'corrected'.  A profile written
	          by a writer that did not state its refusals says nothing about
	          them, and a census of it can only say so.

	Columns are located by the `# columns` header line, so a file that gains or
	reorders columns is read by name.

	returns a dict; `refused` is None exactly when `schema` is 0.
	"""
	head         = header_comment_statements(path)
	cols         = []
	schema       = None
	certified    = None
	cond_tol     = None
	provenance   = []
	coupling     = None
	counts       = None
	legend_seen  = False
	# The writer wraps a statement of the block over several comment lines,
	# the continuations carrying no key of their own, so a statement read one
	# line at a time is a sentence cut in half.  `text` names the field the
	# continuations belong to.
	text = {'adv_stationarity_operator': None, 'adv_conditional_tol': None,
	        'adv_model_restrictions': None}
	key = None
	for cont, line in head:
		w = line.split()
		if not w:
			continue
		if cont:
			if key is not None:
				text[key] = text[key] + ' ' + line
			continue
		if w[0] == 'columns':
			cols = w[1:]
			key = None
		elif w[0] == 'adv_schema' and len(w) > 1:
			key = None
			try:
				schema = int(w[1])
			except ValueError:
				schema = None
		elif w[0] == 'adv_input_certified':
			certified = ' '.join(w[1:])
			key = None
		elif w[0] in text:
			key = w[0]
			text[key] = ' '.join(w[1:])
			if key == 'adv_conditional_tol':
				try:
					cond_tol = float(w[1])
				except (ValueError, IndexError):
					cond_tol = None
		elif w[0] == 'adv_status_counts':
			counts = ' '.join(w[1:])
			key = None
		elif w[0] in ('adv_status', 'adv_T_status', 'adv_comp_status',
		              'adv_mass_row', 'adv_product'):
			legend_seen = True
			key = None
		elif w[0].rstrip(':') == 'provenance':
			provenance.append(line)
			key = None
		elif w[0].rstrip(':') == 'coupling':
			coupling = line
			key = None
		else:
			key = None

	operator     = text['adv_stationarity_operator']
	cond_text    = text['adv_conditional_tol']
	restrictions = text['adv_model_restrictions']

	if schema is None:
		schema = 1 if ('adv_status' in cols or legend_seen) else 0

	T_status    = None
	comp_status = None
	mass_row    = None
	if schema >= 2:
		if 'adv_T_status' in cols:
			T_status = loadtxt_cells(path,
			                         usecols=(cols.index('adv_T_status'),))
		if 'adv_comp_status' in cols:
			comp_status = loadtxt_cells(path,
			                            usecols=(cols.index('adv_comp_status'),))
		if 'adv_mass_row' in cols:
			mass_row = loadtxt_cells(path,
			                         usecols=(cols.index('adv_mass_row'),))
	elif schema == 1:
		icol = cols.index('adv_status') if 'adv_status' in cols else -1
		if icol < 0 and first_data_row_ncol(path) >= 8:
			# The position the writer uses when the header does not name it:
			# the seven physical columns, then the status.
			icol = 7
		if icol >= 0:
			T_status = loadtxt_cells(path, usecols=(icol,))

	if T_status is None:
		# Nothing to read: the file's row validity is UNKNOWN whatever its
		# header claimed.
		schema = 0

	if certified is None:
		certified = 'unknown'
		if coupling is not None and 'certified=' in coupling:
			certified = ('unknown; the input carries no adv_input_certified '
			             'line, and its own "%s" is a certification of the run '
			             'state, not of the post-processed rows' % coupling)

	return {'path':                   path,
	        'schema':                 schema,
	        'columns':                cols,
	        'T_status':               T_status,
	        'comp_status':            comp_status,
	        'refused':                None if T_status is None
	                                  else (np.asarray(T_status) != 0),
	        'status_names':           ADV_STATUS_NAMES_2 if schema >= 2
	                                  else ADV_STATUS_NAMES_1,
	        'input_certified':        certified,
	        'stationarity_operator':  operator,
	        'mass_row':               mass_row,
	        'conditional_tol':        cond_tol,
	        'conditional_tol_text':   cond_text,
	        'model_restrictions':     restrictions,
	        'status_counts':          counts,
	        'provenance':             provenance}


def transit_tool_identity(paths):
	"""Identity of the code that produced a spectrum.

	The repository HEAD names the committed text; the modification time of each
	file names the text actually executed, which is not the same thing when the
	working tree carries changes.  Both are reported.
	"""
	import subprocess
	import time as _time
	out = []
	head = ''
	try:
		here = os.path.dirname(os.path.realpath(paths[0]))
		head = subprocess.run(['git', '--no-optional-locks', '-C', here,
		                       'rev-parse', '--short=12', 'HEAD'],
		                      capture_output=True, text=True,
		                      timeout=20).stdout.strip()
	except Exception:                                       # noqa: BLE001
		head = ''
	out.append('git=%s' % (head if head else 'unavailable'))
	for p in paths:
		try:
			out.append('%s mtime=%s'
			           % (os.path.basename(p),
			              _time.strftime('%Y-%m-%dT%H:%M:%S',
			                             _time.localtime(os.path.getmtime(p)))))
		except OSError:
			out.append('%s mtime=unavailable' % os.path.basename(p))
	return ' '.join(out)


def transit_environment_overrides():
	"""The EXHALE_TRANSIT_ / TPM_ overrides in effect, name and value."""
	return sorted((k, v) for k, v in os.environ.items()
	              if k.startswith('EXHALE_TRANSIT_') or k.startswith('TPM_'))


def transit_state_files(path, selection, require=True):
	"""The (hydro, ioniz) pair a transit spectrum is synthesized from.

	Two states of one run describe the same column and are NOT the same
	numbers, so which one a curve stands on is part of the curve:

	``solution``  ``output/Hydro_ioniz.txt`` and ``output/Ion_species.txt``,
	              the state the wind solver converged and the certification
	              judged.  Its composition is the one the solve carries: with
	              the ionization stages transported it is the transported
	              partition, and with them at local equilibrium it is that
	              equilibrium.
	``adv``       ``output/Hydro_ioniz_adv.txt`` and
	              ``output/Ion_species_adv.txt``, the advection-corrected
	              profile, an INDEPENDENT discretization of the same column
	              (second-order BDF2 marching at the bulk velocity, no eddy
	              term, no element drift).  It is what this tool has always
	              read, and it stays the default.

	The two are kept apart on purpose: reading one for the temperature and
	the other for the composition would synthesize a line from a state that
	solves neither set of equations.

	A stationary solve writes no `_adv` products: the advection
	correction is a post-process of a state, and the solve does not run
	it; the evaluation of a solved state (`Restart intent: stationary
	evaluate`) does run it and writes the pair.  With
	``require`` the pair a selection names must be on disk, and a missing
	one is refused by name; the other state is never read in its place,
	because a curve carrying the name of one state and the numbers of the
	other says something false about both.  ``require=False`` is for a
	caller that is only asking where the other state would be.

	Returns (hydro_path, ioniz_path, selection); an unknown selection raises
	ValueError and a missing required pair raises FileNotFoundError.
	"""
	sel = (selection or 'adv').strip().lower()
	if sel in ('adv', 'post_process', 'corrected'):
		pair = (os.path.join(path, 'output', 'Hydro_ioniz_adv.txt'),
		        os.path.join(path, 'output', 'Ion_species_adv.txt'), 'adv')
	elif sel in ('solution', 'solved', 'state'):
		pair = (os.path.join(path, 'output', 'Hydro_ioniz.txt'),
		        os.path.join(path, 'output', 'Ion_species.txt'), 'solution')
	else:
		raise ValueError('(EXHALE_transit) EXHALE_TRANSIT_STATE = %r is '
		                 'neither "adv" nor "solution"' % (selection,))
	if require:
		missing = [q for q in pair[:2] if not os.path.exists(q)]
		if missing:
			other = ('solution' if pair[2] == 'adv' else 'adv')
			raise FileNotFoundError(
			    '(EXHALE_transit) EXHALE_TRANSIT_STATE = %s names the %s '
			    'state of this run, and these files of it are not there: '
			    '%s. The advection-corrected profile is a post-process of a '
			    'marching state and a stationary solve writes none, so a '
			    'solved run has no "adv" state to read; nothing here reads '
			    'the %s state in its place. Set EXHALE_TRANSIT_STATE=%s to '
			    'synthesize the line from the state the solver certified, '
			    'and do not create `_adv` files to satisfy a file name.'
			    % (pair[2], pair[2], ', '.join(missing), other, other))
	return pair


def file_identity(path):
	"""What file a number came from: its path, its md5 and its size.

	A run directory is written again by the next run, so a path alone does
	not identify the numbers a spectrum was built from.  The digest does.
	"""
	import hashlib
	try:
		h = hashlib.md5()
		with open(path, 'rb') as fh:
			for chunk in iter(lambda: fh.read(1 << 20), b''):
				h.update(chunk)
		return '%s md5=%s bytes=%d' % (path, h.hexdigest(),
		                               os.path.getsize(path))
	except OSError:
		return '%s absent' % path


# The header keys by which a profile states which state it holds: where the
# numbers came from (`provenance`, `source`), which boundary the solve stood
# on (`boundary_model`, `boundary_reservoir`), which physics was active and
# whether the state was certified (`coupling`, and `derived_from`, the same
# fields as they stand in the advection-corrected products this tool reads),
# and which option was changed at a restart (`option_change`).
STATE_PROVENANCE_KEYS = ('provenance', 'source', 'boundary_model',
                         'boundary_reservoir', 'coupling', 'derived_from',
                         'option_change')


def state_provenance_statements(path):
	"""The statements a profile makes about the state it holds.

	Read from the file's own header and copied out unchanged, so a spectrum
	carries the identity of the state it was synthesized from whichever
	selection was made.  An unreadable file gives an empty list.
	"""
	out = []
	try:
		head = header_comment_statements(path)
	except OSError:
		return out
	for cont, line in head:
		if cont:
			continue
		w = line.split()
		if w and w[0].rstrip(':') in STATE_PROVENANCE_KEYS:
			out.append(line)
	return out


def transit_state_oi_levels(path, selection):
	"""The O I ground-term level file belonging to a state selection.

	The three-level statistical equilibrium behind the O I 1302 triplet is
	written for both states, and the level populations a line integrates
	must come from the same state as the density and the temperature it
	integrates them with.
	"""
	name = ('OI_levels_adv.txt' if selection == 'adv' else 'OI_levels.txt')
	return os.path.join(path, 'output', name)


def state_pair_difference(hydro_a, ioniz_a, hydro_b, ioniz_b):
	"""How far the two states of one run stand apart, measured.

	The largest relative difference, over the rows both files carry, of the
	temperature, the neutral hydrogen density and the He 2^3S density -- the
	three quantities a He I 10830 or a Balmer curve is built from.  Returns
	None when either file is missing or the two are not on one grid, so a
	caller can state that the comparison was not available rather than a zero.
	"""
	for q in (hydro_a, ioniz_a, hydro_b, ioniz_b):
		if not os.path.exists(q):
			return None
	try:
		Ta = loadtxt_cells(hydro_a, usecols=(4,))
		Tb = loadtxt_cells(hydro_b, usecols=(4,))
		ia = loadtxt_cells(ioniz_a, usecols=(1, 6), unpack=True)
		ib = loadtxt_cells(ioniz_b, usecols=(1, 6), unpack=True)
	except Exception:
		return None
	Ta = np.asarray(Ta); Tb = np.asarray(Tb)
	if Ta.size != Tb.size or np.asarray(ia[0]).size != np.asarray(ib[0]).size:
		return None

	def rel(x, y):
		x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
		d = np.abs(x - y)/np.maximum(np.maximum(np.abs(x), np.abs(y)), 1.0e-99)
		return float(np.nanmax(d)) if d.size else float('nan')

	return {'T': rel(Ta, Tb), 'HI': rel(ia[0], ib[0]),
	        'HeI_2_3S': rel(ia[1], ib[1]), 'rows': int(Ta.size)}


def transit_metadata_block(adv, tool_identity, overrides, census,
                           line_census=None, state=None):
	"""Comment lines that travel with a saved transit curve.

	Comments only: every line returned here is written behind a '#', so
	`np.loadtxt` and any reader that skips comment lines sees exactly the
	numerical columns it saw before.

	adv           the record of read_adv_validity for the profile the spectrum
	              was built from, or None if none was read.
	census        {'sampled', 'refused', 'above_cap', 'reasons', 'comp'} over
	              the rows the chords sample, or None when the input states no
	              row validity.
	line_census   (b [Rp], share, max share, b of the maximum) for THIS line:
	              the share of the line-center optical depth each ray takes
	              from refused rows.  None when there is nothing to report.
	state         {'selection', 'hydro', 'ioniz', 'difference', 'identity',
	              'provenance'}: WHICH of the run's two states this curve
	              stands on, named here so a curve always says so, how far
	              the other one stands from it, the md5 and size of each
	              file actually read, and the statements that state makes
	              about itself (its `provenance` lines, the boundary model
	              and reservoir it stood on, the physics that was active and
	              whether it was certified).  None when the caller states no
	              selection.
	"""
	L = ['transit_schema %d' % TRANSIT_SCHEMA,
	     'transit_product: a one-way transmission of the profile named below; '
	     'this block adds comments only, and the columns are the spectrum']
	if state is not None:
		L.append('state_selection %s' % state['selection'])
		L.append('state_files %s %s' % (state['hydro'], state['ioniz']))
		# The identity of the numbers, not of a name: a run directory is
		# overwritten by the next run, so the digest is what ties this
		# spectrum to the state it was synthesized from.
		for q in state.get('identity') or []:
			L.append('state_file_identity %s' % q)
		# Copied out of the state's own header, unchanged.
		for q in state.get('provenance') or []:
			L.append('state_%s' % q)
		d = state.get('difference')
		if d is None:
			L.append('state_difference: the other state of this run was not '
			         'read (absent, or not on this grid), so how far the two '
			         'stand apart is UNKNOWN, which is not zero')
		else:
			L.append('state_difference over %d rows, largest relative: '
			         'T %.3e, n(H I) %.3e, n(He 2^3S) %.3e -- the '
			         'advection-corrected profile is an INDEPENDENT '
			         'discretization of this column and these are its '
			         'distances from the state named above'
			         % (d['rows'], d['T'], d['HI'], d['HeI_2_3S']))
	if adv is None:
		L.append('adv_input: none read; row validity UNKNOWN')
	else:
		L.append('adv_input %s' % adv['path'])
		L.append('adv_schema %s'
		         % (adv['schema'] if adv['schema'] > 0 else
		            'absent (legacy file); row validity UNKNOWN'))
		L.append('adv_input_certified %s' % adv['input_certified'])
		if adv['stationarity_operator']:
			L.append('adv_stationarity_operator %s'
			         % adv['stationarity_operator'])
		# The fraction a corrected row of the input is accurate to, copied
		# out of the profile's own header: a spectrum built on corrected
		# rows is accurate to no better than that fraction of the mass flux
		# of those rows, and a reader of the curve alone cannot know it
		# otherwise.  The largest measure over the rows the chords sample
		# says how close to the fraction this profile came.
		if adv.get('conditional_tol_text'):
			L.append('adv_conditional_tol %s' % adv['conditional_tol_text'])
		if adv.get('mass_row') is not None:
			_m = np.asarray(adv['mass_row'])
			if _m.size:
				L.append('adv_mass_row max %.3e over the %d rows of the '
				         'profile; the measure of every row is the '
				         'adv_mass_row column of that file'
				         % (float(np.nanmax(_m)), _m.size))
		if adv['model_restrictions']:
			L.append('adv_model_restrictions %s' % adv['model_restrictions'])
		for p in adv['provenance']:
			L.append('input_%s' % p)
	L.append('transit_tool %s' % tool_identity)
	if census is None:
		L.append('census: the input states no row validity, so the rows the '
		         'chords sample are neither corrected nor refused here; the '
		         'count is UNKNOWN, which is not zero refused')
	else:
		L.append('census sampled_rows %d refused_rows %d above_disk_cap %d'
		         % (census['sampled'], census['refused'], census['above_cap']))
		if census.get('reasons'):
			L.append('census_refusal_reasons %s' % census['reasons'])
		if census.get('comp'):
			L.append('census_composition %s' % census['comp'])
	if line_census is not None:
		b, share, share_max, b_max = line_census
		L.append('census_share_b[Rp] ' + ' '.join('%.4f' % x for x in b))
		L.append('census_share       ' + ' '.join('%.4f' % x for x in share))
		L.append('census_share_max %.4f at b = %.4f Rp' % (share_max, b_max))
	L.append('census_note a share of line-center optical depth carried by '
	         'refused rows is a CONTRIBUTION diagnostic and not an '
	         'uncertainty on the transit depth: the disk-averaged '
	         'transmission averages exp(-tau) and is nonlinear in tau, so a '
	         'saturated ray moves little when its optical depth moves')
	if overrides:
		for k, v in overrides:
			L.append('override %s=%s' % (k, v))
	else:
		L.append('override none in effect: every window, resolving power and '
		         'geometry setting is the built-in default or input.inp')
	return L


def resonance_depth(lam0_A, f_osc, A21, mass, n_lower, instr_res,
                    Grid_Number, r_grid, Rp, data_r, data_v, data_T,
                    A_star, A_atm, A_planet,
                    half_A=4.0, nlam=401, band_A=4.0):
	"""Spherical transit depth for a single resonance line.

	The same integration as the He and Lya lines: spherical chord
	sqrt(r^2-b^2),
	LOS velocity v_los = v*x/r, Voigt profile via wofz, trapezoid optical
	depth, projected-area disk average, instrument convolution.

	lam0_A   line-center wavelength [Angstrom]
	n_lower  lower-level (ion ground) density on the data_r chord [m^-3]
	returns  (depth_linecenter_%, depth_4A_band_%, l_plot_A, transmission)."""
	lam0   = lam0_A*1e-10
	nu0    = c_light/lam0
	l_onde = np.linspace((lam0_A-half_A)*1e-10, (lam0_A+half_A)*1e-10, nlam)
	nu_l   = c_light/l_onde
	exp_tau = np.zeros((Grid_Number, nlam))
	for p in range(Grid_Number):
		r_temp = r_grid[p]*Rp
		arg    = chord_shell_indices(data_r, r_temp)
		r_LOS  = data_r[arg]
		x_LOS  = np.sqrt(r_LOS**2.0 - r_temp**2.0)*np.sign(r_LOS)
		dx     = np.abs(x_LOS[1:] - x_LOS[:-1])
		v_th   = np.sqrt(2.0*kb*data_T[arg]/mass)*_turb_factor()
		v_x    = x_LOS*data_v[arg]/r_LOS
		n_lo   = n_lower[arg]
		Dnu    = nu0*v_th/c_light
		a_v    = A21/(4.0*np.pi*Dnu)
		# Vectorized over wavelength: line-of-sight on axis 0, wavelength on axis 1.
		X     = (nu_l[None, :] - nu0)/Dnu[:, None]
		Voigt = f_osc*Fadd_const/Dnu[:, None] \
		        * wofz(X - (v_x/v_th)[:, None] + 1j*a_v[:, None]).real
		I     = n_lo[:, None]*Voigt
		exp_tau[p, :] = np.exp(-np.sum(dx[:, None]/2.0
		                               * (I[:-1, :] + I[1:, :]), axis=0))
	prob_tot = np.array([np.trapz(x=r_grid, y=2.0*exp_tau[:, li]*r_grid)
	                     *A_planet/(A_atm - A_planet) for li in range(nlam)])
	avg = ((A_star - A_atm) + (A_atm - A_planet)*prob_tot)/A_star
	avg = avg*A_star/(A_star - A_planet)
	FWHM = lam0/instr_res
	sig  = FWHM/(2.0*np.sqrt(2.0*np.log(2.0)))
	vg   = np.linspace(-(l_onde[-1]-l_onde[0])*0.5,
	                    (l_onde[-1]-l_onde[0])*0.5, nlam)
	conv = convolve(avg, np.exp(-0.5*(vg/sig)**2.0), boundary='extend')
	l_plot = l_onde*1e10
	band   = np.abs(l_plot - lam0_A) <= band_A*0.5
	return (1.0 - conv.min())*100.0, (1.0 - conv[band].mean())*100.0, l_plot, conv


def band_integrated_depth(l_plot_A, transmission, band_A, mask_A=()):
	"""Mean excess absorption [%] over an EXPLICIT wavelength band.

	l_plot_A     wavelength grid of the model curve [Angstrom]
	transmission T_lambda on that grid (1 = no absorption)
	band_A       (lam_lo, lam_hi) integration limits [Angstrom]
	mask_A       iterable of (lam_lo, lam_hi) sub-intervals to EXCLUDE,
	             for reproducing a measurement that had to throw pixels
	             away (interstellar cores, airglow emission).  Excluded
	             pixels are dropped from the mean, they are not set to
	             one.

	Returned as (1 - <T>)*100 over the surviving pixels.  This exists as a
	named function because a resonance line whose components are optically
	thick has a line-center depth that is a property of the line profile
	and of the instrument, not of the atmosphere: only a band-integrated
	quantity can be put next to a published band-integrated measurement.
	The band and the mask are arguments, never defaults, so that every
	comparison has to state which band it used."""
	sel = (l_plot_A >= band_A[0]) & (l_plot_A <= band_A[1])
	for lo, hi in mask_A:
		sel &= ~((l_plot_A >= lo) & (l_plot_A <= hi))
	if not np.any(sel):
		raise ValueError('band_integrated_depth: the band %s minus the mask '
		                 '%s contains no wavelength points'
		                 % (str(band_A), str(mask_A)))
	return (1.0 - transmission[sel].mean())*100.0


def resonance_spectrum(components, mass, instr_res, window_A, nlam,
                       Grid_Number, r_grid, Rp, data_r, data_v, data_T,
                       A_star, A_atm, A_planet, rotate_disk_average):
	"""Disk-averaged transmission spectrum of a multi-component resonance
	line: spherical chords, Voigt tau summed over the components, instrument
	convolution, and planet rotation via the exact projected-disk integral
	(rotate_disk_average), as for the He, Lya and Balmer lines.

	Each component carries its OWN lower-level density along the chord,
	   components = [(lam0_A, f_osc, A21, n_lower), ...]   n_lower in m^-3,
	because a component is defined by the level it absorbs out of.  For a
	doublet out of one ion ground state the same array is simply repeated;
	for the O I 1302/1304/1306 triplet the three arrays are the three
	ground-term fine-structure populations, which is the whole reason the
	density belongs to the component and not to the line."""
	l_onde = np.linspace(window_A[0]*1e-10, window_A[1]*1e-10, nlam)
	nu_l   = c_light/l_onde
	exp_tau = np.zeros((Grid_Number, nlam))
	for p in range(Grid_Number):
		r_temp = r_grid[p]*Rp
		arg    = chord_shell_indices(data_r, r_temp)
		r_LOS  = data_r[arg]
		x_LOS  = np.sqrt(r_LOS**2.0 - r_temp**2.0)*np.sign(r_LOS)
		dx     = np.abs(x_LOS[1:] - x_LOS[:-1])
		v_th   = np.sqrt(2.0*kb*data_T[arg]/mass)*_turb_factor()
		v_x    = x_LOS*data_v[arg]/r_LOS
		I = np.zeros((arg.size, nlam))
		for (lam0_A, f_osc, A21, n_lower) in components:
			n_lo = n_lower[arg]
			nu0 = c_light/(lam0_A*1e-10)
			Dnu = nu0*v_th/c_light
			a_v = A21/(4.0*np.pi*Dnu)
			X   = (nu_l[None, :] - nu0)/Dnu[:, None]
			I  += n_lo[:, None]*f_osc*Fadd_const/Dnu[:, None] 			      * wofz(X - (v_x/v_th)[:, None] + 1j*a_v[:, None]).real
		exp_tau[p, :] = np.exp(-np.sum(dx[:, None]/2.0
		                               * (I[:-1, :] + I[1:, :]), axis=0))
	prob_tot = np.array([np.trapz(x=r_grid, y=2.0*exp_tau[:, li]*r_grid)
	                     *A_planet/(A_atm - A_planet) for li in range(nlam)])
	avg = ((A_star - A_atm) + (A_atm - A_planet)*prob_tot)/A_star
	avg = avg*A_star/(A_star - A_planet)
	# instrument convolution
	lam_ref = np.mean([co[0] for co in components])*1e-10
	FWHM = lam_ref/instr_res
	sig  = FWHM/(2.0*np.sqrt(2.0*np.log(2.0)))
	vg   = np.linspace(-(l_onde[-1]-l_onde[0])*0.5,
	                    (l_onde[-1]-l_onde[0])*0.5, nlam)
	conv = convolve(avg, np.exp(-0.5*(vg/sig)**2.0), boundary='extend')
	# planet rotation: exact projected-disk integral, then instrument LSF
	avg_rot  = rotate_disk_average(exp_tau, l_onde)
	conv_rot = convolve(avg_rot, np.exp(-0.5*(vg/sig)**2.0), boundary='extend')
	# avg_rot is returned as well as conv_rot: it is the transmission with the
	# physical broadening (wind + rotation) but WITHOUT the instrument LSF,
	# which is what a band-integrated comparison needs.  Convolution with the
	# LSF conserves the integral, so it cancels from a ratio taken over a band
	# wide compared with the LSF; applying it before a weighted band average
	# with a narrow weight would smear absorption out of where the weight is.
	return dict(l_plot=l_onde*1e10, avg=avg, avg_rot=avg_rot,
	            conv=conv, conv_rot=conv_rot,
	            Tl=(1.0 - avg.min())*100.0,
	            Tl_conv=(1.0 - conv.min())*100.0,
	            Tl_conv_rot=(1.0 - conv_rot.min())*100.0)
