import numpy as np
import matplotlib.pyplot as plt
from scipy.special import wofz
from astropy.convolution import convolve
import os
import time

start = time.time()

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
G   = 6.67e-11			# Gravitational constant [m3/kg/s2]
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
RJ  = 6.9911e7				# Jupiter radius
MJ  = 1.898e27				# Jupiter mass
R_sun = 6.96000000e8 	                # Sun radius [m]
M_sun = 1.989e30			# Sun mass [kg]
gray  = '#a0a0a0'			# Color gray for line plot

# ----- NEEDED USER INPUTS ----- #

# Fill below here the necessary inputs to calculate the transmission spectra
#
# path       = path to the folder where ATES simulation has been performed
# Input_file = ATES' auto-generated input file
# Hydro_file = ATES' hydro output (with path)
# Ioniz_file = ATES' ionization output (with path)
# fig_name   = Name of the output figure (leave empty for not saving the figure)
# abs_file   = Name of the output file with absorption data(leave empty for not saving the file)

path = _tenv('PATH', '.')  # ATES' files destination folder (env override)
Input_file = path + '/input.inp'
Hydro_file = path + '/output/Hydro_ioniz_adv.txt'
Ioniz_file = path + '/output/Ion_species_adv.txt'
fig_name_hei = ''
fig_name_lya = ''
abs_file   = ''

# Data not in input_file
R_star    = float(_tenv('RSTAR_RSUN', '0.44'))*R_sun  # Stellar radius (env override)
# Instrument spectral resolving power R = lambda/Delta-lambda for the Gaussian
# line-spread convolution.  Each is env-overridable (EXHALE_TRANSIT_RES_* , with
# the TPM_RES_* fallback); defaults below match the instruments named inline.
Instr_res_HeTR = float(_tenv('RES_HETR', '8e4'))  # He I 10830: CARMENES 8e4 / GIANO-B 5e4
Instr_res_HI   = float(_tenv('RES_HI',   '5e4'))  # Ly-alpha 1215.67: HST-STIS ~ 1e4-1e5
Instr_res_Ha   = float(_tenv('RES_HA',   '1.15e5'))  # H-alpha 6562.8: HARPS/CARMENES-VIS ~ 1.1e5
Instr_res_Hb   = float(_tenv('RES_HB',   '1.15e5'))  # H-beta 4861.35
Instr_res_MgII = float(_tenv('RES_MGII', '3e4'))  # Mg II h&k 2796/2803 NUV: HST/STIS ~ 3e4
Instr_res_CaII = float(_tenv('RES_CAII', str(Instr_res_Ha)))  # Ca II H&K optical
Instr_res_NaI  = float(_tenv('RES_NAI',  str(Instr_res_Ha)))  # Na I D optical
# Planet rotation period [days]
rot_period = float(_tenv('ROTP', '4.88')) # [days] (env override)

# ----- H-alpha (n=2 -> n=3) inputs ----- #
# H-alpha absorption arises from the n=2 hydrogen population. Following
# Christie, Arras & Li (2013, ApJ 772, 144), the 2s/2p populations are
# set by the rate-equilibrium equations (their Eqs. 12-13) including
# Ly-alpha radiative pumping (1s<->2p). The Ly-alpha mean intensity
# J_lya(r) is NOT produced by ATES; it is obtained in one of two ways:
#
#  (1) If Jlya_file is set to an existing two-column text file, J_lya(r)
#      is read from it:
#         col 1 = r / R_p   (same radial coordinate as ATES output)
#         col 2 = J_lya     = Ly-alpha mean intensity J_nu at line
#                             center [erg s^-1 cm^-2 Hz^-1 sr^-1]
#      (linearly interpolated onto the ATES grid, clamped at the ends).
#
#  (2) Otherwise J_lya is estimated with the simple assumption of
#      Huang et al. (2017, ApJ 851, 150), Eq.(6) and related text:
#         J_lya(r) ~ 0.1 * F_LyC / Dnu_D(r)
#      where Dnu_D = nu_Lya*sqrt(2 kB T/m_p)/c is the local Ly-alpha
#      Doppler width and F_LyC is the DEPOSITED (absorbed) Lyman-continuum
#      flux (each LyC ionization balanced by a recombination -> Ly-alpha
#      photon; Huang+2017 page 11). F_LyC is computed from the actual
#      stellar input: the incident EUV-band (E>13.6 eV) flux at the
#      planet, 10^LEUV/(4*pi*a^2), times the fraction absorbed by the
#      atmosphere, (1 - exp(-tau_LyC)) with tau_LyC = sigma_LyC * N_HI
#      (the vertical neutral-H column). For the optically-thick atomic
#      layer the absorbed fraction ~ 1, i.e. F_LyC ~ incident stellar LyC.
Jlya_file  = ''            # path to J_lya(r) file; '' => Huang(2017) estimate
F_LyC_override = 0.0       # 0 => auto from stellar LyC; >0 => fixed [erg cm^-2 s^-1]
sigma_LyC      = 6.3e-18   # H photoionization cross section at LyC [cm^2]
# Day-night / 2D flux dilution (xi) applied to the incident stellar LyC,
# matching ATES's "2D approximate method": Rate/2 -> 0.5, Rate/4 -> 0.25,
# else 1.0. This is read from input.inp automatically; set xi_override>0
# to force a value (cf. the xi factor of Christie+2013 / Huang+2017).
xi_override    = 0.0       # 0 => auto from input.inp appx_mth
fig_name_ha = ''
# n=2 photoionization (Balmer continuum) rates Gamma_2s, Gamma_2p [s^-1].
# If T_star > 0 they are estimated from a diluted stellar blackbody
# Balmer continuum (E>3.4 eV) at the orbital distance; otherwise the
# manual values below are used (0 => neglected). T_star = stellar
# effective temperature [K] (e.g. ~6065 K for HD209458, ~5050 K for
# HD189733). R_star (above) and a_orb set the dilution.
T_star   = float(_tenv('TSTAR', '6065.0'))          # stellar effective temperature [K] (<=0 disables; env override)
Gamma_2s = 0.0             # used only if T_star <= 0
Gamma_2p = 0.0

# ------------------------------ #

# Number of discretization points
Grid_Number = 200

# Select range of considered wavelength and number of wavelengths
# (He window env-overridable: fast winds broaden the line past the default
# +/-1.5 A, so widen for high-velocity / extended-domain runs.)
lmin_HeTR = float(_tenv('HE_LMIN', '10828.2'))
lmax_HeTR = float(_tenv('HE_LMAX', '10831.2'))
number_lambda_HeTR = int(_tenv('HE_N', '201'))

lmin_HI = 1214.1
lmax_HI = 1217.0
number_lambda_HI = 201

lmin_Ha = 6561.0
lmax_Ha = 6564.6
number_lambda_Ha = 201

# H-beta (n=2 -> n=4) shares the same n=2 population as H-alpha
lmin_Hb = 4860.0
lmax_Hb = 4862.7
number_lambda_Hb = 201
fig_name_hb = ''

# ------------------------------ #

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
# (Christie, Arras & Li 2013, ApJ 772, 144; rate coeffs in their Table 2)

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
A_2p1s  = 6.3e8                 # A(2p->1s) [s^-1]   (Table 2, R11)
A_2s1s  = 8.26                  # A(2s->1s) two-photon [s^-1] (Table 2, R10)

# Einstein-B coefficients in the mean-intensity (J_nu) convention, in cgs
# so that B * J_lya [erg s^-1 cm^-2 Hz^-1 sr^-1] gives a rate in s^-1:
#   B21 = A21 c^2 / (2 h nu^3),   B12 = (g2/g1) B21
c_cgs, h_cgs, kb_cgs = 2.99792458e10, 6.62607015e-27, 1.380649e-16
eV2Hz = 2.417989242e14         # Hz per eV
B21_lya = A_2p1s*c_cgs**2.0/(2.0*h_cgs*nu_Lya**3.0)   # ~2.85e9
B12_lya = (g2p/g1s)*B21_lya                           # ~8.55e9 (1s->2p pump)

def n2_populations(T, n1s, ne, Jlya, G2s=0.0, G2p=0.0):
	# Solve the 2s/2p rate-equilibrium (Christie+2013 Eqs. 12-13) for the
	# n=2 populations [cm^-3]. All densities in cm^-3, T in K, Jlya in cgs.
	# Returns (n2s, n2p, n2tot). Forward collisional rates from Table 2;
	# reverse rates by detailed balance (g-weights; 2s,2p ~ degenerate).
	T  = np.maximum(T, 1.0)               # defensive floor (avoid 1/T, log(0))
	t4 = T/1.0e4
	# Case-B and level-resolved recombination (Draine 2011; Table 2 R2,R8,R9)
	aB  = 2.54e-13*t4**(-0.8163 - 0.0208*np.log(t4))
	a2s = (0.282 + 0.047*t4 - 0.006*t4**2.0)*aB
	a2p = aB - a2s
	# Collisional excitation 1s->2s, 1s->2p, and 2s<->2p l-mixing (R3,R4,R5)
	C1s2s = 1.21e-8*(1.0/t4)**0.455*np.exp(-118400.0/T)
	C1s2p = 1.71e-8*(1.0/t4)**0.077*np.exp(-118400.0/T)
	C2s2p = 6.21e-5*(np.log(T/1.02) - 0.57721)/np.sqrt(T)
	# Reverse (de-exciting) rates by detailed balance. Written in their
	# analytically-cancelled form: the Boltzmann exp(118400/T) factor
	# cancels the exp(-118400/T) inside C1s2s/C1s2p (super-elastic
	# collisions), so we avoid the numerical 0*inf at very low T.
	C2s1s = 1.21e-8*(1.0/t4)**0.455*(g1s/g2s)          # = C1s2s*exp(118400/T)
	C2p1s = 1.71e-8*(1.0/t4)**0.077*(g1s/g2p)          # = C1s2p*(g1s/g2p)*exp(..)
	C2p2s = C2s2p*(g2s/g2p)
	# Ly-alpha radiative pump (1s->2p) and stimulated emission (2p->1s)
	Ppump = B12_lya*Jlya
	Pstim = B21_lya*Jlya
	# 2x2 linear system  [[L2p, -M12],[-M21, L2s]] [n2p,n2s]^T = [S2p,S2s]^T
	L2p = A_2p1s + Pstim + (C2p1s + C2p2s)*ne + G2p
	L2s = (C2s1s + C2s2p)*ne + G2s + A_2s1s
	S2p = (Ppump + C1s2p*ne)*n1s + a2p*ne**2.0
	S2s = (C1s2s*ne)*n1s + a2s*ne**2.0
	M12 = C2s2p*ne
	M21 = C2p2s*ne
	det = L2p*L2s - M12*M21
	det = np.where(np.abs(det) > 0.0, det, 1.0)   # guard (det>0 physically)
	n2p = np.maximum((S2p*L2s + M12*S2s)/det, 0.0)
	n2s = np.maximum((L2p*S2s + M21*S2p)/det, 0.0)
	return n2s, n2p, n2s + n2p

def gamma_n2_balmer(T_star, R_over_a):
	# n=2 photoionization rate [s^-1] from a diluted stellar blackbody
	# Balmer continuum (E > 3.4 eV, lambda < 3647 A). Hydrogenic cross
	# section sigma_2(nu) = sigma2_th*(nu2/nu)^3. R_over_a = R_star/a
	# (dimensionless). For a stellar beam, the rate is int F_nu/(h nu)*sigma dnu
	# with F_nu = pi B_nu(T_star)*(R_star/a)^2.
	if T_star <= 0.0:
		return 0.0
	E2   = 3.40                          # eV, n=2 ionization threshold
	nu2  = E2*eV2Hz
	s2th = 1.4e-17                       # cm^2, sigma at the n=2 threshold
	Eg   = np.linspace(E2, 13.6, 400)    # eV (>13.6 eV: BB flux negligible + H1s absorbs)
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

# ------------------------- #

# Read useful parameters from the input file of ATES
LEUV = None      # log10 of EUV (Lyman-continuum-band) luminosity [erg/s]
appx_mth = ''    # ATES 2D flux approximation (sets the day-night xi factor)
with open(Input_file,'r') as f:

	data = f.readline()
	num = 1
	while data:
		data = f.readline()
		if num == 2:  Rp = float(get_word(data,4))*RJ
		if num == 3:  Mp = float(get_word(data,4))*MJ
		if num == 4:  T0 = float(get_word(data,4))
		if num == 5:  a_orb = float(get_word(data,4))*AU
		if num == 9:  Mstar = float(get_word(data,5))*M_sun
		# Content-based reads (robust to line shifts from spectrum options)
		if 'EUV luminosity'    in data:  LEUV = float(data.split(':')[-1])
		if 'approximate method' in data: appx_mth = data.split(':')[-1].strip()
		num += 1

f.close()

# Load profiles
r,rho,v,p,T,heat,cool = np.loadtxt(Hydro_file, unpack = True)
# Ion_species.txt: read only the first 7 columns (r + H/He). In EXHALE
# this file also carries trace-metal columns (C/N/O), so we slice rather
# than unpack all of them.
r,nhi,nhii,nhei,nheii,nheiii,nheiTR = \
    np.loadtxt(Ioniz_file, usecols = range(7), unpack = True)

# --------------------------------------------------------------------- #
# Auto-size the wavelength window for each line so the WHOLE line profile is
# captured (the line returns to the continuum inside the window) for any
# wind, instead of using fixed +/-few-A windows that clip fast/hot winds.
#
# The line half-width has two contributions, and we take whichever is
# larger:
#   (1) KINEMATIC:  dv_kin = |v_r|_max + K_th * v_thermal + v_rot
#       (the fastest line-of-sight gas plus its Doppler wing plus rotation);
#   (2) DAMPING WING: for an optically thick line the Lorentzian wings stay
#       black far past the kinematic edge.  The wing reaches optical depth
#       unity at x_w = sqrt(tau0 * a / sqrt(pi)) Doppler widths, where
#       tau0 = sigma0 * N_lower is the line-center optical depth of the
#       vertical column and a = A21 lambda / (4 pi c) / (v_th/c) is the
#       Voigt parameter.  dv_wing = x_w * v_thermal.
# The window is lambda0 * (1 +/- margin * max(dv_kin, dv_wing)), never
# narrower than the historical default, honoring TPM_HE_LMIN/LMAX overrides.
_amu = 1.66053907e-27
_kB  = 1.380649e-23
_ec2 = 0.026540045                # pi e^2 / (m_e c)  [cm^2 Hz] (sqrt-pi form below)
_ccm = 2.99792458e10

def _line_halfwidth(lam0_A, m_atom_amu, fosc, A21, N_col_cm2):
    """Return the physical line half-width [m/s]."""
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

# vertical columns of the lower-level absorbers [cm^-2] (rectangle rule; r in R_p)
_dr_cm = np.abs(np.gradient(r))*Rp*1.0e2       # cm (r in R_p, Rp in m -> *1e2)
_N_HI    = float(np.sum(nhi*_dr_cm))           # cm^-3 * cm -> cm^-2
_N_HeTR  = float(np.sum(nheiTR*_dr_cm))
_fHe = f10830_34 + f10830_25 + f10829_09
_he_forced = _tenv_set('HE_LMIN') or _tenv_set('HE_LMAX')
lmin_HeTR, lmax_HeTR = _apply_window(10830.34,
    _line_halfwidth(10830.34, 4.0, _fHe, 1.022e7, _N_HeTR), lmin_HeTR, lmax_HeTR, _he_forced)
lmin_HI, lmax_HI = _apply_window(1215.67,
    _line_halfwidth(1215.67, 1.0, 0.4162, 6.27e8, _N_HI), lmin_HI, lmax_HI)
# Balmer lines are optically thin in the n=2 population -> kinematic only
lmin_Ha, lmax_Ha = _apply_window(6562.80,
    _line_halfwidth(6562.80, 1.0, 0.6407, 4.41e7, 0.0), lmin_Ha, lmax_Ha)
lmin_Hb, lmax_Hb = _apply_window(4861.35,
    _line_halfwidth(4861.35, 1.0, 0.1193, 8.42e6, 0.0), lmin_Hb, lmax_Hb)

# adequate, ODD sampling after widening (astropy convolution needs odd kernels)
def _odd(n):
    n = int(n); return n if n % 2 == 1 else n+1
number_lambda_HeTR = _odd(max(number_lambda_HeTR, (lmax_HeTR-lmin_HeTR)/0.02))
number_lambda_HI   = _odd(max(number_lambda_HI,   (lmax_HI  -lmin_HI  )/0.02))
number_lambda_Ha   = _odd(max(number_lambda_Ha,   (lmax_Ha  -lmin_Ha  )/0.02))
number_lambda_Hb   = _odd(max(number_lambda_Hb,   (lmax_Hb  -lmin_Hb  )/0.02))
print('(TPM) auto windows [A]: He[%.1f,%.1f] Lya[%.1f,%.1f] Ha[%.1f,%.1f]'
      % (lmin_HeTR, lmax_HeTR, lmin_HI, lmax_HI, lmin_Ha, lmax_Ha))

# Save inverted profiles
r_I 	   = -np.flip(r)
T_I 	   =  np.flip(T)
v_I 	   = -np.flip(v)
nheiTR_I =  np.flip(nheiTR)
nH = nhi/(1 + 2.2e-5)
nH_I = np.flip(nH)
nD = nH*(2.2e-5)
nD_I = np.flip(nD)
Rib  = min(r[-1],R_star/Rp)

# ----- H-alpha: n=2 hydrogen population (Christie+2013) ----- #
# Electron density [cm^-3] (H + He contributions). n1s = neutral H.
ne_cm  = nhii + nheii + 2.0*nheiii
n1s_cm = nhi
do_Ha = True   # H-alpha is always computed (J_lya from file or estimate)
if do_Ha:
	if (len(Jlya_file) > 0) and os.path.exists(Jlya_file):
		# (1) Read J_lya(r) from file
		rj, Jlya_in = np.loadtxt(Jlya_file, usecols = (0,1), unpack = True)
		Jlya = np.interp(r, rj, Jlya_in)      # clamps to endpoints outside range
		jlya_src = 'file ' + Jlya_file
	else:
		if len(Jlya_file) > 0:
			print('(TPM) WARNING: Jlya_file "%s" not found; '
			      'using Huang(2017) estimate.' % Jlya_file)
		# (2) Simple estimate (Huang et al. 2017, Eq.(6) & related text):
		#     J_lya ~ 0.1 * F_LyC / Dnu_D, with Dnu_D the local Ly-alpha
		#     Doppler width (proton thermal speed) and F_LyC the deposited
		#     (absorbed) Lyman-continuum flux.
		# --- F_LyC = absorbed stellar LyC flux [erg cm^-2 s^-1] ---
		if F_LyC_override > 0.0:
			F_LyC = F_LyC_override
			print('(TPM)   F_LyC = %.3e erg/cm2/s (manual override)' % F_LyC)
		elif LEUV is not None:
			# Day-night / 2D dilution (xi), matching ATES's appx method
			if   xi_override > 0.0:        xi = xi_override
			elif 'Rate/4' in appx_mth:     xi = 0.25
			elif 'Rate/2' in appx_mth:     xi = 0.5
			else:                          xi = 1.0
			a_cm       = a_orb*1.0e2                        # orbital dist [cm]
			F_inc      = 10.0**LEUV/(4.0*np.pi*a_cm**2.0)   # incident stellar LyC flux
			r_cm       = r*Rp*1.0e2                         # radius [cm] (r in R_p)
			N_HI_tot   = np.trapz(nhi, r_cm)                # vertical HI column [cm^-2]
			abs_frac   = 1.0 - np.exp(-sigma_LyC*N_HI_tot)  # fraction of LyC absorbed
			F_LyC      = xi*F_inc*abs_frac
			print('(TPM)   F_LyC: xi=%.2f * incident %.3e * absorbed frac %.4f'
			      ' = %.3e erg/cm2/s' % (xi, F_inc, abs_frac, F_LyC))
		else:
			F_LyC = 1.0e4   # last-resort fallback if LEUV unavailable
			print('(TPM)   WARNING: LEUV not read; using F_LyC = %.2e' % F_LyC)
		Dnu_D_lya = nu_Lya*np.sqrt(2.0*kb*T/mp)/c_light  # Hz; F_LyC cgs -> J_lya cgs
		Jlya = 0.1*F_LyC/Dnu_D_lya
		jlya_src = 'Huang(2017) estimate (F_LyC = %.2e erg/cm2/s)' % F_LyC
	# n=2 photoionization (Balmer continuum) from a diluted stellar
	# blackbody, if a stellar T_eff is given (overrides manual Gamma_2*).
	if T_star > 0.0:
		Gamma_2s = gamma_n2_balmer(T_star, R_star/a_orb)
		Gamma_2p = Gamma_2s   # ~equal (Huang+2017 Table 2: 25.7 vs 21.5 s^-1)
		print('(TPM)   n=2 photoionization: Gamma_2s = Gamma_2p = %.2e s^-1'
		      ' (T_star = %g K)' % (Gamma_2s, T_star))
	n2s_cm, n2p_cm, n2_cm = n2_populations(T, n1s_cm, ne_cm, Jlya,
	                                       G2s = Gamma_2s, G2p = Gamma_2p)
	# Carry the 2s and 2p populations separately (different Balmer cross
	# sections); build the symmetric (inverted + normal) chord profiles in m^-3.
	data_n2  = np.concatenate((np.flip(n2_cm)*1.0e6,  n2_cm*1.0e6))
	data_n2s = np.concatenate((np.flip(n2s_cm)*1.0e6, n2s_cm*1.0e6))
	data_n2p = np.concatenate((np.flip(n2p_cm)*1.0e6, n2p_cm*1.0e6))
	print('(TPM) H-alpha: J_lya from %s; max n2 = %.3e cm^-3 (2s/2p max ratio %.2f)'
	      % (jlya_src, n2_cm.max(), (n2p_cm/np.maximum(n2s_cm,1e-99)).max()))

# Areas
A_star    = np.pi*R_star**2.0
A_planet  = np.pi*Rp**2.0
A_atm     = np.pi*(Rib*Rp)**2.0

# Instrument parameters
FWHM_HeTR 			= 1.0830e-6/Instr_res_HeTR # FWHM for convolution at 10830 Angstrom
sigma_gauss_HeTR = FWHM_HeTR/(2.0*np.sqrt(2.0*np.log(2.0)))
FWHM_HI 			= lA/Instr_res_HI # FWHM for convolution at 1215.67 Angstrom
sigma_gauss_HI = FWHM_HI/(2.0*np.sqrt(2.0*np.log(2.0)))
FWHM_Ha 			= l_Ha/Instr_res_Ha # FWHM for convolution at 6562.8 Angstrom
sigma_gauss_Ha = FWHM_Ha/(2.0*np.sqrt(2.0*np.log(2.0)))
FWHM_Hb 			= l_Hb/Instr_res_Hb # FWHM for convolution at 4861.4 Angstrom
sigma_gauss_Hb = FWHM_Hb/(2.0*np.sqrt(2.0*np.log(2.0)))

# ------------------------- #

# Construct grid 
r_grid = np.array([Rib**(j/(Grid_Number-1)) for j in range(Grid_Number)]) 

# Save wavelengths for figure title
lmin_lbl_HeTR = lmin_HeTR
lmax_lbl_HeTR = lmax_HeTR
lmin_lbl_HI = lmin_HI
lmax_lbl_HI = lmax_HI
lmin_lbl_Ha = lmin_Ha
lmax_lbl_Ha = lmax_Ha
lmin_lbl_Hb = lmin_Hb
lmax_lbl_Hb = lmax_Hb

# Convert to meter
lmin_HeTR = lmin_HeTR*1e-10
lmax_HeTR = lmax_HeTR*1e-10
lmin_HI = lmin_HI*1e-10
lmax_HI = lmax_HI*1e-10
lmin_Ha = lmin_Ha*1e-10
lmax_Ha = lmax_Ha*1e-10
lmin_Hb = lmin_Hb*1e-10
lmax_Hb = lmax_Hb*1e-10

# Generate vector of wavelength
l_onde_HeTR = np.linspace(lmin_HeTR, lmax_HeTR, number_lambda_HeTR)
l_plot_HeTR = l_onde_HeTR*1.0e10
l_onde_HI = np.linspace(lmin_HI, lmax_HI, number_lambda_HI)
l_plot_HI = l_onde_HI*1.0e10
l_onde_Ha = np.linspace(lmin_Ha, lmax_Ha, number_lambda_Ha)
l_plot_Ha = l_onde_Ha*1.0e10
l_onde_Hb = np.linspace(lmin_Hb, lmax_Hb, number_lambda_Hb)
l_plot_Hb = l_onde_Hb*1.0e10

# ------------------------- #

# Create double vectors 

# Radius
data_r = np.append(r_I*Rp, r*Rp)

# Temperature
data_T = np.append(T_I, T)

# Velocity
data_v = np.append(v_I, v)*1e-2

# HeI triplet density [m^3]
data_nheiTR = np.concatenate((nheiTR_I*1.0e6, nheiTR*1.0e6))

# HI and D density [m^3]
data_nHI = np.concatenate((nH_I*1.0e6, nH*1.0e6))
data_nD = np.concatenate((nD_I*1.0e6, nD*1.0e6))

# ------------------------- #

# Initialize temporary variables
exp_tau_HeTR = np.zeros((Grid_Number, number_lambda_HeTR))
prob_temp_HeTR = np.zeros((Grid_Number, number_lambda_HeTR))
prob_tot_HeTR = np.zeros(number_lambda_HeTR)
# --
exp_tau_HD = np.zeros((Grid_Number, number_lambda_HI))
prob_temp_HD = np.zeros((Grid_Number, number_lambda_HI))
prob_tot_HD = np.zeros(number_lambda_HI)
# --
if do_Ha:
	exp_tau_Ha = np.zeros((Grid_Number, number_lambda_Ha))
	prob_tot_Ha = np.zeros(number_lambda_Ha)
	exp_tau_Hb = np.zeros((Grid_Number, number_lambda_Hb))
	prob_tot_Hb = np.zeros(number_lambda_Hb)


# ----- LOOP ON GRID POINTS ----- #
for p in range(Grid_Number):

	# Center distance
	r_temp = r_grid[p]*Rp

	# My version of x_value + x_coordinate
	x_LOS = np.array([])	# Initialize as empty

	# Get index of current inner points
	data_arg = np.where( (abs(data_r) >= r_temp))[0]
	r_LOS    = data_r[data_arg]
	
	for rad in r_LOS: 
		
		# Get sign of radius
		Rad_sign = np.sign(rad)
		
		# Calculate distance
		dist     = np.sqrt(rad**2.0 - r_temp**2.0)*Rad_sign
			
		# Append to vector
		x_LOS = np.append(x_LOS,dist)
	
	# Get number of points on current LOS
	len_x = len(x_LOS)
	
	# Get x-grid spacing
	data_dx = np.abs(x_LOS[1:] - x_LOS[:-1])
	
	# ------------------------ #

	# Get values of temperature and velocity corresponding to the selected points 
	T_LOS = data_T[data_arg]
	v_th_HeTR = np.sqrt(2.0*kb*T_LOS/mHe)
	v_th_HI = np.sqrt(2.0*kb*T_LOS/mp)
	v_th_D = np.sqrt(2.0*kb*T_LOS/mD)
	n_HeTR = data_nheiTR[data_arg]
	n_HI = data_nHI[data_arg]
	n_D = data_nD[data_arg]
	v_LOS = data_v[data_arg]
	v_x = x_LOS*v_LOS/r_LOS

	# ------------------------ #

	# Doppler shifts parameters
	Dnu_He3_1 = (nu_He3_1*v_th_HeTR)/c_light 
	Dnu_He3_2 = (nu_He3_2*v_th_HeTR)/c_light 
	Dnu_He3_3 = (nu_He3_3*v_th_HeTR)/c_light

	Dnu_HI = (nu_HI*v_th_HI)/c_light
	Dnu_D = (nu_D*v_th_D)/c_light 

	# Arguments of Voigt Function 
	a_He3_1 = A12_HeTR/(4.0*np.pi*Dnu_He3_1)
	a_He3_2 = A12_HeTR/(4.0*np.pi*Dnu_He3_2)
	a_He3_3 = A12_HeTR/(4.0*np.pi*Dnu_He3_3)

	a_HI = A12_HI/(4.0*np.pi*Dnu_HI)
	a_D = A12_D/(4.0*np.pi*Dnu_D)
	
	# ------------------------ #
	
	# ----- Absorption integrals via Faddeeva method, vectorized over wavelength ----- #
	nu_HeTR_l = c_light/l_onde_HeTR
	vx_over_vth_HeTR = (v_x/v_th_HeTR)[:, None]

	# First line HeI3
	X_He3_1 = (nu_HeTR_l[None, :] - nu_He3_1)/Dnu_He3_1[:, None]
	Voigt_1 = f10830_34*Fadd_const/Dnu_He3_1[:, None] \
	          * wofz(X_He3_1 - vx_over_vth_HeTR + 1j*a_He3_1[:, None]).real

	# Second line HeI3
	X_He3_2 = (nu_HeTR_l[None, :] - nu_He3_2)/Dnu_He3_2[:, None]
	Voigt_2 = f10830_25*Fadd_const/Dnu_He3_2[:, None] \
	          * wofz(X_He3_2 - vx_over_vth_HeTR + 1j*a_He3_2[:, None]).real

	# Third line HeI3
	X_He3_3 = (nu_HeTR_l[None, :] - nu_He3_3)/Dnu_He3_3[:, None]
	Voigt_3 = f10829_09*Fadd_const/Dnu_He3_3[:, None] \
	          * wofz(X_He3_3 - vx_over_vth_HeTR + 1j*a_He3_3[:, None]).real

	# Integrand summed over the three components
	I_HeTR = n_HeTR[:, None]*(Voigt_1 + Voigt_2 + Voigt_3)

	# Optical depth via trapezoids along the line of sight
	exp_tau_HeTR[p, :] = np.exp(-np.sum(data_dx[:, None]/2.0
	                                    * (I_HeTR[:-1, :] + I_HeTR[1:, :]), axis=0))

# ----- End HeI metastable triplet ----- #


	# ----- Absorption integrals via Faddeeva method, vectorized over wavelength ----- #
	nu_HI_l = c_light/l_onde_HI

	# Hydrogen
	X_HI = (nu_HI_l[None, :] - nu_HI)/Dnu_HI[:, None]
	Voigt_HI = f_la*Fadd_const/Dnu_HI[:, None] \
	           * wofz(X_HI - (v_x/v_th_HI)[:, None] + 1j*a_HI[:, None]).real

	# Deuterium
	X_D = (nu_HI_l[None, :] - nu_D)/Dnu_D[:, None]
	Voigt_D = f_D*Fadd_const/Dnu_D[:, None] \
	          * wofz(X_D - (v_x/v_th_D)[:, None] + 1j*a_D[:, None]).real

	# Integrand summed over Hydrogen and Deuterium
	I_HD = n_HI[:, None]*Voigt_HI + n_D[:, None]*Voigt_D

	# Optical depth via trapezoids along the line of sight
	exp_tau_HD[p, :] = np.exp(-np.sum(data_dx[:, None]/2.0
	                                  * (I_HD[:-1, :] + I_HD[1:, :]), axis=0))

# ----- End Hydrogen and Deuterium ----- #

	# ----- H-alpha (n=2 -> n=3), sub-level-resolved absorption ----- #
	# 2s and 2p have different absorption oscillator strengths, so the optical
	# depth is f_2s*n_2s + f_2p*n_2p (not f_multiplet*n_2). The fine-structure
	# components share the same Doppler-dominated profile.
	if do_Ha:
		n_2s   = data_n2s[data_arg]
		n_2p   = data_n2p[data_arg]
		Dnu_Ha = (nu_Ha*v_th_HI)/c_light     # same thermal width as HI (proton)
		a_Ha   = A12_Ha/(4.0*np.pi*Dnu_Ha)

		# Voigt absorption profile via Faddeeva method, vectorized over wavelength
		nu_Ha_l = c_light/l_onde_Ha
		X_Ha = (nu_Ha_l[None, :] - nu_Ha)/Dnu_Ha[:, None]
		Voigt_Ha = Fadd_const/Dnu_Ha[:, None] \
		           * wofz(X_Ha - (v_x/v_th_HI)[:, None] + 1j*a_Ha[:, None]).real

		# Optical depth via trapezoids along the line of sight
		I_Ha = (f_Ha_2s*n_2s + f_Ha_2p*n_2p)[:, None]*Voigt_Ha
		exp_tau_Ha[p, :] = np.exp(-np.sum(data_dx[:, None]/2.0
		                                  * (I_Ha[:-1, :] + I_Ha[1:, :]), axis=0))

		# ----- H-beta (n=2 -> n=4), same 2s/2p populations ----- #
		Dnu_Hb = (nu_Hb*v_th_HI)/c_light
		a_Hb   = A12_Hb/(4.0*np.pi*Dnu_Hb)
		nu_Hb_l = c_light/l_onde_Hb
		X_Hb = (nu_Hb_l[None, :] - nu_Hb)/Dnu_Hb[:, None]
		Voigt_Hb = Fadd_const/Dnu_Hb[:, None] \
		           * wofz(X_Hb - (v_x/v_th_HI)[:, None] + 1j*a_Hb[:, None]).real
		I_Hb = (f_Hb_2s*n_2s + f_Hb_2p*n_2p)[:, None]*Voigt_Hb
		exp_tau_Hb[p, :] = np.exp(-np.sum(data_dx[:, None]/2.0
		                                  * (I_Hb[:-1, :] + I_Hb[1:, :]), axis=0))

# ----- End H-alpha / H-beta ----- #



# Integral over the planet's projected area metastable HeI triplet
for l in range(number_lambda_HeTR):
	prob_tot_HeTR[l] = np.trapz(x = r_grid, y = 2.0*exp_tau_HeTR[:,l]*r_grid)*A_planet/(A_atm - A_planet)

# Do geometric average with star area
avg_prob_HeTR = ((A_star - A_atm) + (A_atm - A_planet)*prob_tot_HeTR[:])/A_star

# Normalization to continuum
avg_prob_HeTR = avg_prob_HeTR[:]*A_star/(A_star - A_planet)


# Integral over the planet's projected area Hydrogen and Deuterium
for l in range(number_lambda_HI):
	prob_tot_HD[l] = np.trapz(x = r_grid, y = 2.0*exp_tau_HD[:,l]*r_grid)*A_planet/(A_atm - A_planet)
	               
# Do geometric average with star area
avg_prob_HD = ((A_star - A_atm) + (A_atm - A_planet)*prob_tot_HD[:])/A_star

# Normalization to continuum
avg_prob_HD = avg_prob_HD[:]*A_star/(A_star - A_planet)

# Physical transmission is in [0,1]; clip any non-finite entries (a saturated
# Lya damping wing can overflow exp(-tau) in a single grid cell) so the line
# profile and the auto-scaled plot axes stay well defined.
avg_prob_HD = np.clip(np.nan_to_num(avg_prob_HD, nan=1.0, posinf=1.0, neginf=0.0),
                      0.0, 1.0)

# Integral over the planet's projected area H-alpha and H-beta
if do_Ha:
	for l in range(number_lambda_Ha):
		prob_tot_Ha[l] = np.trapz(x = r_grid, y = 2.0*exp_tau_Ha[:,l]*r_grid)*A_planet/(A_atm - A_planet)
	avg_prob_Ha = ((A_star - A_atm) + (A_atm - A_planet)*prob_tot_Ha[:])/A_star
	avg_prob_Ha = avg_prob_Ha[:]*A_star/(A_star - A_planet)
	for l in range(number_lambda_Hb):
		prob_tot_Hb[l] = np.trapz(x = r_grid, y = 2.0*exp_tau_Hb[:,l]*r_grid)*A_planet/(A_atm - A_planet)
	avg_prob_Hb = ((A_star - A_atm) + (A_atm - A_planet)*prob_tot_Hb[:])/A_star
	avg_prob_Hb = avg_prob_Hb[:]*A_star/(A_star - A_planet)

# -------------------------------------------------------------------- #
# Planet-rotation broadening: exact projected-disk integral.
#
# For a tidally-/solid-body-rotating atmosphere the line-of-sight velocity of
# a projected disk patch is v = v_ang * x, where x = b*cos(phi) is the sky
# coordinate perpendicular to the spin axis (transit geometry), b the impact
# parameter, and phi the azimuth around the disk.  We Doppler-shift each
# impact-parameter transmission profile by its LOCAL rotational velocity,
# average over azimuth, and only then take the same projected-area disk
# average as the non-rotating profile.  This is the full rotational-broadening
# integral for an azimuthally-symmetric absorber -- it replaces the former
# approximation of convolving the disk-averaged profile with a single Gaussian
# kernel built from one effective-radius velocity v_ang*R_eff.
v_ang     = 2.0*np.pi/(60.0*60.0*24.0*rot_period)   # angular rate [rad/s]
N_phi_rot = int(_tenv('ROT_NPHI', '64'))            # azimuthal samples


def _disk_average(mat, nlam):
	"""Projected-area disk average of a transmission
	matrix mat[p, l] (p over r_grid, l over wavelength), identical to the
	non-rotating normalization used above."""
	prob = np.array([np.trapz(x=r_grid, y=2.0*mat[:, l]*r_grid)
	                 * A_planet/(A_atm - A_planet) for l in range(nlam)])
	avg = ((A_star - A_atm) + (A_atm - A_planet)*prob)/A_star
	return avg*A_star/(A_star - A_planet)


def _rotate_disk_average(exp_tau_mat, l_onde):
	"""Rotation-broadened, projected-area disk-averaged transmission.
	exp_tau_mat[p, l] is the line-of-sight transmission at impact parameter
	r_grid[p] (in R_p) and rest-frame wavelength l_onde[l] (m).  Each chord's
	profile is shifted by the solid-body rotational Doppler velocity
	v = v_ang*b*cos(phi) and averaged over azimuth phi before the disk
	average -- the exact rotational broadening for this azimuthally-symmetric
	absorber (no single-R_eff / Gaussian-kernel approximation)."""
	nlam = l_onde.size
	if v_ang <= 0.0:
		return _disk_average(exp_tau_mat, nlam)
	cphi = np.cos(np.linspace(0.0, 2.0*np.pi, N_phi_rot, endpoint=False))
	rot = np.empty_like(exp_tau_mat)
	for p in range(exp_tau_mat.shape[0]):
		b_m  = r_grid[p]*Rp                         # impact parameter [m]
		prof = exp_tau_mat[p, :]
		acc  = np.zeros(nlam)
		for cp in cphi:
			v = v_ang*b_m*cp                        # LOS rotational velocity
			# absorption seen at l_onde originates at rest wavelength
			# l_onde*(1 + v/c); continuum (transmission ~ 1) outside window.
			acc += np.interp(l_onde*(1.0 + v/c_light), l_onde, prof,
			                 left=prof[0], right=prof[-1])
		rot[p, :] = acc/N_phi_rot
	return _disk_average(rot, nlam)


# Create vector of wavelengths to convolve metastable HeI triplet profile
l_range_HeTR = l_onde_HeTR[-1] - l_onde_HeTR[0]
v_gauss_HeTR = np.linspace(-l_range_HeTR*0.5, l_range_HeTR*0.5, number_lambda_HeTR)

# Do convolution with gaussian nIR instrument Resolution
gaussian_HeTR = np.exp(-0.5*(v_gauss_HeTR[:]/sigma_gauss_HeTR)**2.0)   # Normalized at 1
convolved_avg_prob_HeTR = convolve(avg_prob_HeTR, gaussian_HeTR, boundary = 'extend')

# -----

# Create vector of wavelengths to convolve Hydrogen and Deuterium profile
l_range_HD = l_onde_HI[-1] - l_onde_HI[0]
v_gauss_HD = np.linspace(-l_range_HD*0.5, l_range_HD*0.5, number_lambda_HI)

# Do convolution with gaussian UV instrument Resolution
gaussian_HD = np.exp(-0.5*(v_gauss_HD[:]/sigma_gauss_HI)**2.0)   # Normalized at 1
convolved_avg_prob_HD = convolve(avg_prob_HD, gaussian_HD, boundary = 'extend')

# -----

# Create vector of wavelengths to convolve H-alpha / H-beta profiles
if do_Ha:
	l_range_Ha = l_onde_Ha[-1] - l_onde_Ha[0]
	v_gauss_Ha = np.linspace(-l_range_Ha*0.5, l_range_Ha*0.5, number_lambda_Ha)
	gaussian_Ha = np.exp(-0.5*(v_gauss_Ha[:]/sigma_gauss_Ha)**2.0)   # Normalized at 1
	convolved_avg_prob_Ha = convolve(avg_prob_Ha, gaussian_Ha, boundary = 'extend')
	l_range_Hb = l_onde_Hb[-1] - l_onde_Hb[0]
	v_gauss_Hb = np.linspace(-l_range_Hb*0.5, l_range_Hb*0.5, number_lambda_Hb)
	gaussian_Hb = np.exp(-0.5*(v_gauss_Hb[:]/sigma_gauss_Hb)**2.0)
	convolved_avg_prob_Hb = convolve(avg_prob_Hb, gaussian_Hb, boundary = 'extend')

# ------------------------- #

# Planet rotation (exact projected-disk integral), then instrument LSF.
# The intrinsic rotation broadening is applied at the disk-integration level
# via _rotate_disk_average; the instrument Gaussian is convolved afterwards.
transit_depth = (Rp/R_star)**2.0    # retained for the metal generic path below

# metastable HeI triplet
avg_prob_rot_HeTR = _rotate_disk_average(exp_tau_HeTR, l_onde_HeTR)
convolved_rot_prob_HeTR = convolve(avg_prob_rot_HeTR, gaussian_HeTR, boundary='extend')

# Hydrogen and Deuterium (Ly-alpha); guard the saturated-wing overflow as above
avg_prob_rot_HD = _rotate_disk_average(exp_tau_HD, l_onde_HI)
avg_prob_rot_HD = np.clip(np.nan_to_num(avg_prob_rot_HD, nan=1.0, posinf=1.0,
                                        neginf=0.0), 0.0, 1.0)
convolved_rot_prob_HD = convolve(avg_prob_rot_HD, gaussian_HD, boundary='extend')

# H-alpha / H-beta
if do_Ha:
	avg_prob_rot_Ha = _rotate_disk_average(exp_tau_Ha, l_onde_Ha)
	convolved_rot_prob_Ha = convolve(avg_prob_rot_Ha, gaussian_Ha, boundary='extend')
	avg_prob_rot_Hb = _rotate_disk_average(exp_tau_Hb, l_onde_Hb)
	convolved_rot_prob_Hb = convolve(avg_prob_rot_Hb, gaussian_Hb, boundary='extend')

# Transmission minima

Tl_HeI3 = (1.0 - avg_prob_HeTR.min())*100.0
Tl_HeI3_conv = (1.0 - convolved_avg_prob_HeTR.min())*100.0
Tl_HeI3_conv_rot = (1.0 - convolved_rot_prob_HeTR.min())*100.0

Tl_HI = (1.0 - avg_prob_HD.min())*100.0
Tl_HI_conv = (1.0 - convolved_avg_prob_HD.min())*100.0
Tl_HI_conv_rot = (1.0 - convolved_rot_prob_HD.min())*100.0

if do_Ha:
	Tl_Ha = (1.0 - avg_prob_Ha.min())*100.0
	Tl_Ha_conv = (1.0 - convolved_avg_prob_Ha.min())*100.0
	Tl_Ha_conv_rot = (1.0 - convolved_rot_prob_Ha.min())*100.0
	Tl_Hb = (1.0 - avg_prob_Hb.min())*100.0
	Tl_Hb_conv = (1.0 - convolved_avg_prob_Hb.min())*100.0
	Tl_Hb_conv_rot = (1.0 - convolved_rot_prob_Hb.min())*100.0

# ===================================================================== #
# Phase 5a: metal resonance lines (Mg II, Ca II, Na I D), spherical LOS
# --------------------------------------------------------------------- #
# Huang+2023's NUV/optical signal. These are resonance doublets whose
# lower level is the ION GROUND STATE; at ~1e4 K the excited fine-structure
# levels are Boltzmann-negligible, so the lower-level density is the ion
# density itself (n_lower ~= n_ion). The ion densities sit past the H/He
# block of Ion_species_adv.txt, in mion_fsp order (species_table.f90):
#   numpy col 17 = Mg II, 23 = Ca II, 25 = Na I  (0-indexed; r = col 0).
# This block is self-contained (it re-reads Ioniz_file and reuses the
# global chord/area grids) so the validated He/Lya/Ha/Hb output above is
# untouched.  Geometry here is spherical (Phase 5a gate vs Huang Case A);
# the triaxial Roche reconstruction is added in Phase 5b (roche_recon.py).
amu = 1.66053907e-27
mMg, mCa, mNa = 24.305*amu, 40.078*amu, 22.990*amu

# Metals-off runs write only the H/He columns (<=16): skip the metal
# resonance lines automatically (He/Lya/Ha/Hb above are unaffected).
_ncol_ion = np.loadtxt(Ioniz_file, max_rows=1).size
do_metals = (_ncol_ion >= 26)
if do_metals:
	nMgII_cm, nCaII_cm, nNaI_cm = np.loadtxt(Ioniz_file, usecols=(17, 23, 25),
	                                         unpack=True)              # cm^-3
else:
	_r0col = np.loadtxt(Ioniz_file, usecols=(0,))
	nMgII_cm = np.zeros_like(_r0col)   # metals-off: zero metal absorption,
	nCaII_cm = np.zeros_like(_r0col)   # so the metal lines are flat and the
	nNaI_cm  = np.zeros_like(_r0col)   # He/Lya/Ha output still saves.
	print('(TPM) metals-off run: skipping metal resonance lines.')
# symmetric (night + day) chord arrays, m^-3 (mirror data_nHI at line ~460)
data_nMgII = np.concatenate((np.flip(nMgII_cm), nMgII_cm))*1.0e6
data_nCaII = np.concatenate((np.flip(nCaII_cm), nCaII_cm))*1.0e6
data_nNaI  = np.concatenate((np.flip(nNaI_cm),  nNaI_cm ))*1.0e6


def resonance_depth(lam0_A, f_osc, A21, mass, n_lower, instr_res,
                    half_A=4.0, nlam=401, band_A=4.0):
	"""Spherical transit depth for a single resonance line.

	Mirrors the He/Lya pipeline exactly: spherical chord sqrt(r^2-b^2),
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
		arg    = np.where(np.abs(data_r) >= r_temp)[0]
		r_LOS  = data_r[arg]
		x_LOS  = np.sqrt(r_LOS**2.0 - r_temp**2.0)*np.sign(r_LOS)
		dx     = np.abs(x_LOS[1:] - x_LOS[:-1])
		v_th   = np.sqrt(2.0*kb*data_T[arg]/mass)
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


# Mg II is quoted in a 4 A bin; Ca II K / Na D2 line-center. Resolutions:
# Mg II NUV (HST/STIS ~ 3e4), Ca II / Na I optical (~ Instr_res_Ha).
metal_lines = [
	# label,        lam0_A,   f,      A21,     mass, n_lower,     R_instr
	('Mg II 2796',  2796.35,  0.608,  2.60e8,  mMg,  data_nMgII,  Instr_res_MgII),
	('Ca II K 3934', 3933.66, 0.6267, 1.47e8,  mCa,  data_nCaII,  Instr_res_CaII),
	('Na I D2 5890', 5889.95, 0.641,  6.16e7,  mNa,  data_nNaI,   Instr_res_NaI),
]
metal_depth = {}
print('')
print('(TPM) ===== Phase 5a: metal resonance lines (spherical) =====')
print('(TPM)   line           line-center [%]   4A-band [%]')
for lbl, lam_, f_, A_, m_, ncol_, R_ in metal_lines:
	d_lc, d_band, _, _ = resonance_depth(lam_, f_, A_, m_, ncol_, R_)
	metal_depth[lbl] = (d_lc, d_band)
	print('(TPM)   %-13s   %10.3f      %10.3f' % (lbl, d_lc, d_band))
print('')


# --------------------------------------------------------------------- #
# Full doublet transmission spectra for the metal resonance lines.
# Same spherical pipeline as resonance_depth, but with BOTH doublet
# components summed in one wavelength window, and with the instrument
# and planet-rotation convolutions applied exactly as for the H/He
# lines, so the metal lines are first-class TPM outputs (figures saved
# as PNG+PDF below). The Phase 5a table above keeps the validated
# single-component numbers. Skipped automatically for a metals-off run
# (all-zero ion columns). NIST atomic data (Kramida 2020).

fig_name_mgii = ''   # empty = compute but do not save
fig_name_caii = ''
fig_name_nai  = ''

METAL_DOUBLETS = [
	# key, label, components (lam0_A, f, A21), mass, chord density,
	#   instrument R, window [A], nlam
	('MgII', 'Mg II h&k',
	 [(2796.352, 0.608, 2.60e8), (2803.531, 0.303, 2.57e8)],
	 mMg, data_nMgII, Instr_res_MgII, (2790.0, 2810.0), 601, fig_name_mgii),
	('CaII', 'Ca II H&K',
	 [(3933.663, 0.6267, 1.47e8), (3968.469, 0.3116, 1.40e8)],
	 mCa, data_nCaII, Instr_res_CaII, (3927.0, 3975.0), 961, fig_name_caii),
	('NaI', 'Na I D',
	 [(5889.951, 0.641, 6.16e7), (5895.924, 0.320, 6.14e7)],
	 mNa, data_nNaI, Instr_res_NaI, (5884.0, 5902.0), 541, fig_name_nai),
]


def resonance_spectrum(components, mass, n_lower, instr_res, window_A, nlam):
	"""Disk-averaged transmission spectrum of a multi-component resonance
	line: spherical chords, Voigt tau summed over the components, instrument
	convolution, and planet rotation via the exact projected-disk integral
	(_rotate_disk_average), matching the He/Lya/Balmer pipeline."""
	l_onde = np.linspace(window_A[0]*1e-10, window_A[1]*1e-10, nlam)
	nu_l   = c_light/l_onde
	exp_tau = np.zeros((Grid_Number, nlam))
	for p in range(Grid_Number):
		r_temp = r_grid[p]*Rp
		arg    = np.where(np.abs(data_r) >= r_temp)[0]
		r_LOS  = data_r[arg]
		x_LOS  = np.sqrt(r_LOS**2.0 - r_temp**2.0)*np.sign(r_LOS)
		dx     = np.abs(x_LOS[1:] - x_LOS[:-1])
		v_th   = np.sqrt(2.0*kb*data_T[arg]/mass)
		v_x    = x_LOS*data_v[arg]/r_LOS
		n_lo   = n_lower[arg]
		I = np.zeros((arg.size, nlam))
		for (lam0_A, f_osc, A21) in components:
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
	avg_rot  = _rotate_disk_average(exp_tau, l_onde)
	conv_rot = convolve(avg_rot, np.exp(-0.5*(vg/sig)**2.0), boundary='extend')
	return dict(l_plot=l_onde*1e10, avg=avg, conv=conv, conv_rot=conv_rot,
	            Tl=(1.0 - avg.min())*100.0,
	            Tl_conv=(1.0 - conv.min())*100.0,
	            Tl_conv_rot=(1.0 - conv_rot.min())*100.0)


metal_spec = {}
for key_m, lbl_m, comps_m, m_m, ncol_m, R_m, win_m, nl_m, fnm_m 		in METAL_DOUBLETS:
	if ncol_m.max() <= 0.0:
		print('(TPM)   %s: ion column is zero (metals off?); skipped'
		      % lbl_m)
		continue
	metal_spec[key_m] = resonance_spectrum(comps_m, m_m, ncol_m, R_m,
	                                       win_m, nl_m)
	metal_spec[key_m].update(label=lbl_m, comps=comps_m, fig=fnm_m)


# ===================================================================== #
# Phase 5b: 3-D Roche-equipotential reconstruction + velocity broadening
# --------------------------------------------------------------------- #
# Set geometry = 'triaxial' to map the 1-D substellar profile onto the full
# 3-D Roche potential (roche_recon.py) and integrate the transit along the
# +x (substellar) line of sight at each sky point (y,z), with the
# LOS velocity v_sub(r_eff)*x/r - Omega*y (wind + tidally-locked rotation,
# Huang Eq. 16) over 20 angular sectors. Only meaningful for a Roche (Case
# B/D) run; for the spherical Case A run leave geometry = 'spherical'.
geometry = 'spherical'    # 'spherical' (5a) | 'triaxial' (5b); set per run

if geometry == 'triaxial':
	import roche_recon as rr
	from scipy.interpolate import interp1d

	q_rec  = Mstar/Mp                       # M*/M_p
	a_rec  = a_orb/Rp                       # orbital distance [R_p]
	Omega  = 2.0*np.pi/(rot_period*86400.0) # tidally-locked spin [rad/s]
	recon  = rr.ReconMap(r, q_rec, a_rec)   # r = substellar grid [R_p]
	r_lo, r_hi = recon.r_sub[0], recon.r_sub[-1]
	T_of_r = interp1d(r, T,        bounds_error=False, fill_value=(T[0], T[-1]))
	v_of_r = interp1d(r, v*1.0e-2, bounds_error=False, fill_value=0.0)   # m/s
	# substellar lower-level density interpolators [m^-3]
	nlo_of_r = {
		'Mg II 2796':  interp1d(r, nMgII_cm*1e6, bounds_error=False, fill_value=0.0),
		'Ca II K 3934': interp1d(r, nCaII_cm*1e6, bounds_error=False, fill_value=0.0),
		'Na I D2 5890': interp1d(r, nNaI_cm*1e6,  bounds_error=False, fill_value=0.0),
	}
	if do_Ha:
		nlo_of_r['Halpha 6563'] = interp1d(r, n2_cm*1e6, bounds_error=False, fill_value=0.0)
		nlo_of_r['Hbeta 4861']  = interp1d(r, n2_cm*1e6, bounds_error=False, fill_value=0.0)

	def triaxial_depth(lam0_A, f_osc, A21, mass, nlo, instr_res,
	                   half_A=4.0, nlam=161, band_A=4.0,
	                   n_sectors=20, n_b=80, n_x=250):
		"""Transit depth through the 3-D Roche reconstruction.  LOS along +x;
		state at (x,y,z) = substellar state at r_eff (roche_recon); LOS
		velocity v_sub(r_eff)*x/r - Omega*y.  Vectorized over wavelength."""
		lam0 = lam0_A*1e-10
		nu0  = c_light/lam0
		l_onde = np.linspace((lam0_A-half_A)*1e-10, (lam0_A+half_A)*1e-10, nlam)
		nu_l   = c_light/l_onde
		b_grid = np.array([r_hi**(j/(n_b-1)) for j in range(n_b)])   # 1 -> r_hi
		x_los  = np.linspace(-1.5*recon.L1, 1.5*recon.L1, n_x)        # [R_p]
		x_phys = x_los*Rp
		ang    = (np.arange(n_sectors) + 0.5)*2.0*np.pi/n_sectors
		radial = np.ones((n_sectors, n_b, nlam))
		for ia, aa in enumerate(ang):
			ca, sa = np.cos(aa), np.sin(aa)
			for ib, b in enumerate(b_grid):
				y, z = b*ca, b*sa
				reff = recon.r_eff(x_los, y, z)
				inside = np.isfinite(reff)
				if not inside.any():
					continue
				reffc = np.clip(np.where(inside, reff, r_hi), r_lo, r_hi)
				n_l = np.where(inside, nlo(reffc), 0.0)
				T_l = np.where(inside, T_of_r(reffc), 1.0)
				vw  = np.where(inside, v_of_r(reffc), 0.0)
				r3  = np.sqrt(x_los**2 + y*y + z*z)
				vth = np.sqrt(2.0*kb*T_l/mass)
				vlos = vw*x_los/r3 - Omega*(y*Rp)            # m/s
				Dnu = nu0*vth/c_light
				a_v = A21/(4.0*np.pi*Dnu)
				X   = (nu_l[None, :] - nu0)/Dnu[:, None]
				Vo  = f_osc*Fadd_const/Dnu[:, None] \
				      * wofz(X - (vlos/vth)[:, None] + 1j*a_v[:, None]).real
				tau = np.trapz(n_l[:, None]*Vo, x=x_phys, axis=0)
				radial[ia, ib, :] = np.exp(-np.abs(tau))
		# azimuthal mean of the radial (impact-parameter) integral / (r_hi^2-1)
		prob = np.array([np.mean([np.trapz(2.0*b_grid*radial[ia, :, l], b_grid)
		                          for ia in range(n_sectors)])
		                 for l in range(nlam)])/(r_hi**2 - 1.0)
		avg = ((A_star - A_atm) + (A_atm - A_planet)*prob)/A_star
		avg = avg*A_star/(A_star - A_planet)
		FWHM = lam0/instr_res
		sig  = FWHM/(2.0*np.sqrt(2.0*np.log(2.0)))
		vg   = np.linspace(-(l_onde[-1]-l_onde[0])*0.5,
		                    (l_onde[-1]-l_onde[0])*0.5, nlam)
		conv = convolve(avg, np.exp(-0.5*(vg/sig)**2.0), boundary='extend')
		l_plot = l_onde*1e10
		band = np.abs(l_plot - lam0_A) <= band_A*0.5
		return (1.0 - conv.min())*100.0, (1.0 - conv[band].mean())*100.0

	# triaxial radii prediction (at the base equipotential) + lobe relations
	Rx, Ry, Rz = recon.triaxial_radii(min(1.05, r_hi))
	print('(TPM) ===== Phase 5b: 3-D Roche reconstruction (triaxial) =====')
	print('(TPM)   L1 = %.3f R_p   R_RL = %.3f R_p   R_RL/L1 = %.3f (~2/3)   '
	      '4/9 = %.3f' % (recon.L1, recon.roche_lobe_radius(),
	                      recon.roche_lobe_radius()/recon.L1,
	                      (recon.roche_lobe_radius()/recon.L1)**2))
	print('(TPM)   triaxial radii @ r=%.2f: R_px=%.3f R_py=%.3f R_pz=%.3f  '
	      '(R_py*R_pz=%.3f)' % (min(1.05, r_hi), Rx, Ry, Rz, Ry*Rz))
	print('(TPM)   line           line-center [%]   4A-band [%]')
	tri_lines = [
		('Mg II 2796',  2796.35,  0.608,  2.60e8,  mMg, 3.0e4),
		('Ca II K 3934', 3933.66, 0.6267, 1.47e8,  mCa, Instr_res_Ha),
		('Na I D2 5890', 5889.95, 0.641,  6.16e7,  mNa, Instr_res_Ha),
	]
	if do_Ha:
		tri_lines += [('Halpha 6563', 6562.8, f_Ha, A12_Ha, mp, Instr_res_Ha),
		              ('Hbeta 4861',  4861.35, f_Hb, A12_Hb, mp, Instr_res_Hb)]
	for lbl, lam_, f_, A_, m_, R_ in tri_lines:
		d_lc, d_band = triaxial_depth(lam_, f_, A_, m_, nlo_of_r[lbl], R_)
		print('(TPM)   %-13s   %10.3f      %10.3f' % (lbl, d_lc, d_band))
	print('')

# ----- Optional: save model transmission curves for external overplotting --- #
# Triggered by env var TPM_SAVE_PREFIX (no effect when unset). Each file has
# columns: wavelength [A], T_lambda (theoretical), T_lambda (instr. conv.),
# T_lambda (planet-rot + instr. conv.).  Excess absorption [%] = (1 - T)*100.
_save_prefix = _tenv('SAVE_PREFIX', '')
if len(_save_prefix) > 0:
	np.savetxt(_save_prefix + 'tpm_He10830.txt',
	           np.c_[l_plot_HeTR, avg_prob_HeTR,
	                 convolved_avg_prob_HeTR, convolved_rot_prob_HeTR],
	           header='lambda[A]  T_theo  T_instr  T_rot+instr  (He I 10830)')
	if do_Ha:
		np.savetxt(_save_prefix + 'tpm_Halpha.txt',
		           np.c_[l_plot_Ha, avg_prob_Ha,
		                 convolved_avg_prob_Ha, convolved_rot_prob_Ha],
		           header='lambda[A]  T_theo  T_instr  T_rot+instr  (H-alpha 6562.8)')
		# Lyman-alpha (HI 1215.67) shares the n=2 / excited-H pipeline with
		# H-alpha, so it is written alongside it.
		np.savetxt(_save_prefix + 'tpm_Lya.txt',
		           np.c_[l_plot_HI, avg_prob_HD,
		                 convolved_avg_prob_HD, convolved_rot_prob_HD],
		           header='lambda[A]  T_theo  T_instr  T_rot+instr  (Ly-alpha 1215.67)')
	print('(TPM) saved model curves with prefix:', _save_prefix)

# ----- Setup of the figure ----- #

plt.figure(figsize = (8,7))

##### Figure metastable HeI triplet #####

# Line plots
plt.plot([l_He3_1*1.0e10,l_He3_1*1.0e10], [0.0, 1.1], '--', color = gray)
plt.plot([l_He3_2*1.0e10,l_He3_2*1.0e10], [0.0, 1.1], '--', color = gray)
plt.plot([l_He3_3*1.0e10,l_He3_3*1.0e10], [0.0, 1.1], '--', color = gray)

# Plot the curves of transmission
plt.plot(l_plot_HeTR, avg_prob_HeTR, '--', label = r'Theoretical T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HeI3, 2)))
plt.plot(l_plot_HeTR, convolved_avg_prob_HeTR, '-.', label = 'Instrument conv. T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HeI3_conv,2)))
plt.plot(l_plot_HeTR, convolved_rot_prob_HeTR, label = 'Planet rot. + Inst. conv T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HeI3_conv_rot,2)))

# Axis setup
plt.xlabel(r"Wavelength [$\AA{}$]", fontsize = 15)
plt.ylabel(r"T$_{\lambda}$", fontsize = 15)
plt.xlim([lmin_HeTR*1.0e10, lmax_HeTR*1.0e10])
plt.ylim([0.99*avg_prob_HeTR.min(),1.02*avg_prob_HeTR.max()])	
plt.legend(loc = 'best', labelspacing = 1)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

# Construct the title based on inputs
title = r'Avg. Transm. Prob. with $' + str(Grid_Number) + \
         '\\times' + str(Grid_Number) + '$ grid points' + \
         ' -- ' + str(number_lambda_HeTR) + ' pt in $\lambda$ -- $ ' + \
         str(lmin_lbl_HeTR) + ' < \lambda < ' + str(lmax_lbl_HeTR) + '~ \AA{}$'

plt.title(title)

# Print out max transmission probability
print('\n ----- Transmission probability at peak metastable HeI triplet ----- \n')
print(' - Theoretical: ', Tl_HeI3, '%')
print(' - Instrument convolution: ', Tl_HeI3_conv, '%')
print(' - Planet rot. + Inst. convolution: ', Tl_HeI3_conv_rot, '%')
print('\n')

# Save figure
if len(fig_name_hei) > 0 :
	plt.savefig(fig_name_hei)
	plt.savefig(fig_name_hei.rsplit('.',1)[0]+'.pdf')

plt.show(block=False)

##### Figure Hydrogen and Deuterium #####

plt.figure(figsize=(8,7))
# Line plots
plt.plot([lA*1.0e10,lA*1.0e10], [0.0, 1.2], '--', color = gray)
plt.plot([lD*1.0e10,lD*1.0e10], [0.0, 1.2], '--', color = gray)

# Plot the curves of transmission
plt.plot(l_plot_HI, avg_prob_HD, '--', label = 'Theoretical T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HI,2)))
plt.plot(l_plot_HI, convolved_avg_prob_HD, '-.', label = 'Instrument conv. T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HI_conv,2)))
plt.plot(l_plot_HI, convolved_rot_prob_HD, label = 'Planet rot. + Inst. conv. T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HI_conv_rot,2)))

dlm_ism = (1.0 - 3e4/c_light)*1.0e10
dlp_ism = (1.0 + 3e4/c_light)*1.0e10
plt.fill_between(np.arange(lA*dlm_ism,lA*dlp_ism,0.01),0, 1.2, \
                 color = gray, hatch = '//', alpha = 0.1,
                 label = r'ISM absorption $\pm 30$ km/s')

# Axis setup
plt.xlabel(r"Wavelength [$\AA{}$]", fontsize = 15)
plt.ylabel(r"T$_{\lambda}$", fontsize = 15)
plt.xlim([lmin_HI*1.0e10, lmax_HI*1.0e10])
plt.ylim([0.8*avg_prob_HD.min(),1.1*avg_prob_HD.max()])	
plt.legend(loc = 'best', labelspacing = 1)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

# Construct the title based on inputs
title_HD = r'Avg. Transm. Prob. with $' + str(Grid_Number) + \
           '\\times' + str(Grid_Number) + '$ grid points' +  \
           ' -- ' + str(number_lambda_HI) + ' pt in $\lambda$ -- $ ' + \
           str(lmin_lbl_HI) + ' < \lambda < ' + str(lmax_lbl_HI) + '~ \AA{}$'

plt.title(title_HD)

# Print out max transmission probability
print('\n ----- Transmission probability at peak HI Lya ----- \n')
print(' - Theoretical: ', Tl_HI, '%')
print(' - Instrument convolution: ', Tl_HI_conv, '%')
print(' - Planet rot. + Inst. convolution: ', Tl_HI_conv_rot, '%')
print('\n')

# Save figure (must happen while the Lya figure is still current --
# saving at the end of the script would capture the H-beta figure)
if len(fig_name_lya) > 0 :
	plt.savefig(fig_name_lya)
	plt.savefig(fig_name_lya.rsplit('.',1)[0]+'.pdf')


##### Figure H-alpha #####
if do_Ha:
	plt.figure(figsize=(8,7))
	# Line center marker
	plt.plot([l_Ha*1.0e10, l_Ha*1.0e10], [0.0, 1.2], '--', color = gray)

	# Plot the curves of transmission
	plt.plot(l_plot_Ha, avg_prob_Ha, '--', label = 'Theoretical T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_Ha,2)))
	plt.plot(l_plot_Ha, convolved_avg_prob_Ha, '-.', label = 'Instrument conv. T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_Ha_conv,2)))
	plt.plot(l_plot_Ha, convolved_rot_prob_Ha, label = 'Planet rot. + Inst. conv. T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_Ha_conv_rot,2)))

	# Axis setup
	plt.xlabel(r"Wavelength [$\AA{}$]", fontsize = 15)
	plt.ylabel(r"T$_{\lambda}$", fontsize = 15)
	plt.xlim([lmin_Ha*1.0e10, lmax_Ha*1.0e10])
	plt.ylim([0.95*avg_prob_Ha.min(), 1.02*avg_prob_Ha.max()])
	plt.legend(loc = 'best', labelspacing = 1)
	plt.xticks(fontsize = 14)
	plt.yticks(fontsize = 14)
	title_Ha = r'H$\alpha$ Avg. Transm. Prob. -- ' + str(number_lambda_Ha) + \
	           ' pt in $\lambda$ -- $ ' + str(lmin_lbl_Ha) + \
	           ' < \lambda < ' + str(lmax_lbl_Ha) + '~ \AA{}$'
	plt.title(title_Ha)

	print('\n ----- Transmission probability at peak HI H-alpha ----- \n')
	print(' - Theoretical: ', Tl_Ha, '%')
	print(' - Instrument convolution: ', Tl_Ha_conv, '%')
	print(' - Planet rot. + Inst. convolution: ', Tl_Ha_conv_rot, '%')
	print('\n')

	if len(fig_name_ha) > 0 :
		plt.savefig(fig_name_ha)
		plt.savefig(fig_name_ha.rsplit('.',1)[0]+'.pdf')

##### Figure H-beta #####
if do_Ha:
	plt.figure(figsize=(8,7))
	plt.plot([l_Hb*1.0e10, l_Hb*1.0e10], [0.0, 1.2], '--', color = gray)
	plt.plot(l_plot_Hb, avg_prob_Hb, '--', label = 'Theoretical T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_Hb,2)))
	plt.plot(l_plot_Hb, convolved_avg_prob_Hb, '-.', label = 'Instrument conv. T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_Hb_conv,2)))
	plt.plot(l_plot_Hb, convolved_rot_prob_Hb, label = 'Planet rot. + Inst. conv. T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_Hb_conv_rot,2)))
	plt.xlabel(r"Wavelength [$\AA{}$]", fontsize = 15)
	plt.ylabel(r"T$_{\lambda}$", fontsize = 15)
	plt.xlim([lmin_Hb*1.0e10, lmax_Hb*1.0e10])
	plt.ylim([0.95*avg_prob_Hb.min(), 1.02*avg_prob_Hb.max()])
	plt.legend(loc = 'best', labelspacing = 1)
	plt.xticks(fontsize = 14)
	plt.yticks(fontsize = 14)
	title_Hb = r'H$\beta$ Avg. Transm. Prob. -- ' + str(number_lambda_Hb) + \
	           ' pt in $\lambda$ -- $ ' + str(lmin_lbl_Hb) + \
	           ' < \lambda < ' + str(lmax_lbl_Hb) + '~ \AA{}$'
	plt.title(title_Hb)

	print('\n ----- Transmission probability at peak HI H-beta ----- \n')
	print(' - Theoretical: ', Tl_Hb, '%')
	print(' - Instrument convolution: ', Tl_Hb_conv, '%')
	print(' - Planet rot. + Inst. convolution: ', Tl_Hb_conv_rot, '%')
	print('\n')

	if len(fig_name_hb) > 0 :
		plt.savefig(fig_name_hb)
		plt.savefig(fig_name_hb.rsplit('.',1)[0]+'.pdf')


##### Figures: metal resonance doublets #####
for key_m, sp in metal_spec.items():
	plt.figure(figsize=(8, 7))
	for (lam0_A, f_osc, A21) in sp['comps']:
		plt.plot([lam0_A, lam0_A], [0.0, 1.2], '--', color=gray)
	plt.plot(sp['l_plot'], sp['avg'], '--',
	         label='Theoretical T$_{{\lambda}}$ = {} $\%$'.format(
	               round(sp['Tl'], 2)))
	plt.plot(sp['l_plot'], sp['conv'], '-.',
	         label='Instrument conv. T$_{{\lambda}}$ = {} $\%$'.format(
	               round(sp['Tl_conv'], 2)))
	plt.plot(sp['l_plot'], sp['conv_rot'],
	         label='Planet rot. + Inst. conv. T$_{{\lambda}}$ = {} $\%$'.format(
	               round(sp['Tl_conv_rot'], 2)))
	plt.xlabel(r"Wavelength [$\AA{}$]", fontsize=15)
	plt.ylabel(r"T$_{\lambda}$", fontsize=15)
	plt.xlim([sp['l_plot'][0], sp['l_plot'][-1]])
	plt.ylim([0.95*sp['avg'].min(), 1.02*sp['avg'].max()])
	plt.legend(loc='best', labelspacing=1)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.title(sp['label'].replace('&', r'\&') + r' Avg. Transm. Prob. -- '
	          + str(sp['l_plot'].size) + ' pt in $\lambda$')

	print('\n ----- Transmission probability at peak ' + sp['label']
	      + ' ----- \n')
	print(' - Theoretical: ', sp['Tl'], '%')
	print(' - Instrument convolution: ', sp['Tl_conv'], '%')
	print(' - Planet rot. + Inst. convolution: ', sp['Tl_conv_rot'], '%')
	print('\n')

	if len(sp['fig']) > 0:
		plt.savefig(sp['fig'])
		plt.savefig(sp['fig'].rsplit('.', 1)[0] + '.pdf')


print("--- Execution time: %s seconds ---" % (time.time() - start))

plt.show()




