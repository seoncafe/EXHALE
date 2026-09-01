"""Reusable, testable core for the EXHALE transit post-processor.

This module holds the pure pieces of ``EXHALE_transit.py``: physical
constants, the line rest-wavelength / oscillator-strength / Einstein-A
metadata, and the self-contained physics/utility functions.  It references
only its own constants plus numpy/scipy, so it imports cleanly on its own
and can be unit-tested without a simulation.  The main script imports these
names and keeps the top-down orchestration (file reading, density prep,
per-line loops, convolution, plotting, saving).
"""

import os
import numpy as np
from scipy.special import wofz
from astropy.convolution import convolve


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
RJ  = 6.9911e7				# Jupiter radius
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
A_2p1s  = 6.3e8                 # A(2p->1s) [s^-1]   (Table 2, R11)
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
		arg    = np.where(np.abs(data_r) >= r_temp)[0]
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
		arg    = np.where(np.abs(data_r) >= r_temp)[0]
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
