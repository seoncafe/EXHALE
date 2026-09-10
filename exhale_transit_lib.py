"""Reusable, testable core for the EXHALE transit post-processor.

This module holds the pure pieces of ``EXHALE_transit.py``: physical
constants, the line rest-wavelength / oscillator-strength / Einstein-A
metadata, the self-contained physics/utility functions, the chord geometry
those functions and any census of them share, and the reader that decides
what an ``_adv`` profile says about the validity of its own rows.  Apart
from numpy/scipy and the profile readers of ``examples/exhale_io.py`` it
references only its own constants, so it imports cleanly on its own and can
be unit-tested without a simulation.  The main script imports these names
and keeps the top-down orchestration (loading a run, density prep, the
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
A_2p1s  = 6.2649e8              # A(2p->1s) [s^-1], NIST ASD (Wiese and Fuhr 2009); one value with hydrogen_n2_rates.f90
A_2s1s  = 8.26                  # A(2s->1s) two-photon [s^-1] (Table 2, R10)

# Einstein-B coefficients in the mean-intensity (J_nu) convention, in cgs
# so that B * J_lya [erg s^-1 cm^-2 Hz^-1 sr^-1] gives a rate in s^-1:
#   B21 = A21 c^2 / (2 h nu^3),   B12 = (g2/g1) B21
c_cgs, h_cgs, kb_cgs = 2.99792458e10, 6.62607015e-27, 1.380649e-16
eV2Hz = 2.417989242e14         # Hz per eV
B21_lya = A_2p1s*c_cgs**2.0/(2.0*h_cgs*nu_Lya**3.0)   # ~2.85e9
B12_lya = (g2p/g1s)*B21_lya                           # ~8.55e9 (1s->2p pump)

# Doppler / optical-depth helper constants for the auto-window sizing
_amu = 1.66053907e-27
_kB  = 1.380649e-23
_ec2 = 0.026540045                # pi e^2 / (m_e c)  [cm^2 Hz] (sqrt-pi form below)
_ccm = 2.99792458e10


def n2_populations(T, n1s, nHII, ne, Jlya, G2s=0.0, G2p=0.0):
	# Solve the 2s/2p rate-equilibrium (Christie+2013 Eqs. 12-13) for the
	# n=2 populations [cm^-3]. All densities in cm^-3, T in K, Jlya in cgs.
	# Returns (n2s, n2p, n2tot). Forward collisional rates from Table 2;
	# reverse rates by detailed balance (g-weights; 2s,2p ~ degenerate).
	# Mirrors src/modules/radiation/excited_hydrogen.f90::n2_populations.
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
	# Cascade source: the electron recombines onto a proton, so the rate is
	# alpha_2l*ne*nHII. ne and nHII part company wherever helium and metals
	# supply the electrons while hydrogen is still neutral, i.e. at the base.
	S2p = (Ppump + C1s2p*ne)*n1s + a2p*ne*nHII
	S2s = (C1s2s*ne)*n1s + a2s*ne*nHII
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
	  R_star_Rsun [R_sun or None], T_star [K or None]."""
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

	params = dict(Rp=Rp, Mp=Mp, T0=T0, a_orb=a_orb, Mstar=Mstar,
	              LEUV=LEUV, appx_mth=appx_mth,
	              R_star_Rsun=R_star_Rsun, T_star=T_star,
	              resolved=False)

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


def transit_metadata_block(adv, tool_identity, overrides, census,
                           line_census=None):
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
	"""
	L = ['transit_schema %d' % TRANSIT_SCHEMA,
	     'transit_product: a one-way transmission of the profile named below; '
	     'this block adds comments only, and the columns are the spectrum']
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
	(rotate_disk_average), matching the He/Lya/Balmer pipeline.

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
