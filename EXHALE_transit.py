import numpy as np
import matplotlib.pyplot as plt
from scipy.special import wofz
from astropy.convolution import convolve
import os
import sys
import time

# WHERE THIS SCRIPT'S OWN MODULES COME FROM.  Three siblings are imported by
# name -- exhale_transit_lib just below, roche_recon for the triaxial geometry
# and he_line_metrics for the He 10830 fit -- and they live next to THIS file,
# not next to the run: a run directory holds `input.inp` and `output/`, and the
# script is started from there.  Started as `python3 <path>/EXHALE_transit.py`,
# CPython puts this script's own directory (symlinks resolved) first on
# sys.path and all three resolve, so that route was never broken.  Through
# `runpy.run_path`, IPython's `%run`, or `exec(open(...).read())` it does not:
# sys.path[0] is then the run directory, the siblings are not there, and the
# he_line_metrics failure surfaces only as a skipped He 10830 fit because that
# one import is wrapped in a try.  Naming this file's directory explicitly
# makes every one of those routes resolve the same three modules.  `__file__`
# does not exist under `exec()`, hence the guard.  A physical COPY of this
# script into a run directory is the one case nothing can repair: a copy has
# no siblings, and its own directory is all it can know.
try:
	_HERE = os.path.dirname(os.path.realpath(__file__))
except NameError:
	_HERE = os.getcwd()
if _HERE not in sys.path:
	sys.path.insert(0, _HERE)
# examples/ holds exhale_io, the ONE reader of EXHALE's output files: it is
# what knows which rows of a profile are solution cells and which are ghosts
# (loadtxt_cells / physical_cell_rows).  This script used to call np.loadtxt
# on those files itself and so carried the ghost rows into every column
# density, chord integral and depth; it now goes through that reader instead,
# so the rule lives in one place.
sys.path.insert(0, os.path.join(_HERE, 'examples'))
from exhale_io import loadtxt_cells

# Pure constants, line metadata, and physics/utility functions live in the
# importable library so they can be tested without a simulation.  The
# control flow (file reading, density prep, the loop over each line,
# convolution, plotting, saving) stays in this script.
from exhale_transit_lib import (
    _tenv, _tenv_set,
    kb, G, mp, me, mD, c_light, AU, E0, h, ht, e, mHe, RJ, MJ,
    R_sun, M_sun, gray,
    l_He3_1, l_He3_2, l_He3_3, lA, lD,
    nu_He3_1, nu_He3_2, nu_He3_3, nu_HI, nu_D,
    f10830_34, f10830_25, f10829_09, f_la, f_D,
    A12_HeTR, A12_HI, A12_D, Fadd_const,
    l_Ha, nu_Ha, f_Ha, A12_Ha, l_Hb, nu_Hb, f_Hb, A12_Hb,
    f_Ha_2s, f_Ha_2p, f_Hb_2s, f_Hb_2p,
    g1s, g2s, g2p, nu_Lya, A_2p1s, A_2s1s,
    c_cgs, h_cgs, kb_cgs, eV2Hz, B21_lya, B12_lya,
    _amu, _kB, _ec2, _ccm,
    n2_populations, gamma_n2_balmer, get_word, read_input_params,
    orbital_period_days, parameter_with_source,
    _line_halfwidth, _apply_window, _odd,
    resonance_depth, resonance_spectrum, band_integrated_depth, _turb_factor,
    chord_shell_indices, refused_line_center_tau_share,
    first_data_row_ncol,
    read_adv_validity, transit_metadata_block, transit_tool_identity,
    adv_derived_state_verdict, adv_derived_state_text,
    transit_environment_overrides, transit_state_files,
    state_pair_difference, file_identity, state_provenance_statements,
    transit_state_oi_levels,
)

start = time.time()

# ----- NEEDED USER INPUTS ----- #

# Fill below here the necessary inputs to calculate the transmission spectra
#
# path       = path to the folder where the EXHALE simulation has been performed
# Input_file = EXHALE's auto-generated input file
# Hydro_file = EXHALE's hydro output (with path)
# Ioniz_file = EXHALE's ionization output (with path)
#
# WHICH STATE OF THE RUN THE LINE IS SYNTHESIZED FROM, selectable through
# EXHALE_TRANSIT_STATE:
#
#   adv       (default) the advection-corrected profile, Hydro_ioniz_adv.txt
#             and Ion_species_adv.txt.  This is what the tool has always
#             read and what every published curve of this code stands on.
#   solution  the state the wind solver converged and the certification
#             judged, Hydro_ioniz.txt and Ion_species.txt.  With the
#             ionization stages transported (Ionization transport) the
#             composition of that state is the transported partition, and a
#             second, differently discretized answer for the same fractions
#             in the file the line reads would mix two states.
#
# What the two selections MEAN, since they are two answers for one column
# and not two formats of one answer:
#
#   the `_adv` profile is a POST-PROCESS of a marching state.  It takes the
#   marched temperature and composition and applies a steady advective
#   correction cell by cell, second order in the cell width (variable-step
#   BDF2 along the flow, each element stepped on the velocity of its own
#   nuclei, v plus its diffusive drift where He_diffusion is on; update log
#   stage 3 section 82), with no stage eddy term and no drift of a stage
#   against its element, and it refuses that correction where
#   an assumption of the steady equations fails in a cell.  It is an
#   independent discretization of the same column, and it is not what any
#   equation of the run was solved for.
#
#   the solved pair is the state the solver produced and the certification
#   judged: the state whose residuals were measured against a tolerance.
#   Where the ionization stages are transported, its composition is what the
#   stage transport itself produced, so applying the advective correction on
#   top of it would be a SECOND and different transport approximation of the
#   same stages, laid over the one that was solved.
#
# A stationary SOLVE therefore writes no `_adv` products, and none are to be
# manufactured for one by this tool; the evaluation of a solved state
# (`Restart intent: stationary evaluate`) runs the post-process on that
# state and writes them, which is how the LHS 1140 b catalog EW is made.
# A run without them is read with EXHALE_TRANSIT_STATE=solution, and the
# default selection refuses by name when the `_adv` pair is absent rather
# than reading the other state.
#
# The pair is a PAIR: the temperature of one state and the composition of
# the other solve neither set of equations, so both files come from the same
# selection and the selection is named in the header of every product, with
# the md5 of each file read and the state's own provenance and boundary
# model beside it.  Where the post-process exists it is kept either way, as
# the independent discretization the solution is measured against, and the
# measured distance between the two states travels into the same header.

path = _tenv('PATH', '.')  # EXHALE's files destination folder (env override)
Input_file = path + '/input.inp'
# A selection whose files are not there is refused with the reason, not with
# a traceback: the usual case is a stationary solve, which writes no `_adv`
# products (its evaluation, `Restart intent: stationary evaluate`, does).
try:
	Hydro_file, Ioniz_file, _state_selection = transit_state_files(
	    path, _tenv('STATE', 'adv'))
except (FileNotFoundError, ValueError) as _exc:
	sys.exit(str(_exc))

# ----- Output naming: one rule for every product of this script ----- #
# Every line computed here is identified by the same key -- He10830, Lya,
# Halpha, Hbeta, MgII, CaII, NaI -- and both products of a line are named from
# that key and land in the run directory:
#
#   <path>/<save prefix>tpm_<line>.txt   model curve   (always written)
#   <path>/<fig prefix><line>.png/.pdf   figure        (only if a prefix is set)
#
# Both prefixes are name prefixes resolved against the run directory, so an
# unattended run leaves its products next to the output they were built from
# and never depends on the working directory.  Give an absolute prefix to write
# somewhere else.  The curve prefix is empty by default -- the canonical
# tpm_<line>.txt -- and EXHALE_TRANSIT_SAVE_PREFIX only decorates that name.
# Figures are off by default because they are a presentation product, not a
# data product; set EXHALE_TRANSIT_FIG_PREFIX (e.g. "tpm_") to save them.
_save_prefix = os.path.join(path, _tenv('SAVE_PREFIX', ''))
_fig_prefix  = _tenv('FIG_PREFIX', '')


def _checked_transmission(line, prob, lam, clip=False, tol=1.0e-12):
	"""A disk-averaged transmission, checked before the spectrum uses it.

	lam is the wavelength array of the line in METRES, as the l_onde_*
	arrays of this script are; the messages quote it in Angstrom.

	exp(-tau) with tau >= 0 lies in [0, 1], and a large tau underflows it to 0
	rather than overflowing it, so a non-finite transmission is a defect of the
	optical depth or of its inputs (density, temperature, geometry, opacity)
	and not a saturated line. It stops the synthesis and names the first
	wavelength that carries one, instead of being replaced by a value the
	invalid input gives no basis for. A departure from [0, 1] larger than tol
	is a defect as well; with clip=True the departures within tol (the
	rounding of the disk average) are set to the bound and counted.
	"""
	prob = np.asarray(prob, dtype=float)
	bad = ~np.isfinite(prob)
	if bad.any():
		i = int(np.flatnonzero(bad)[0])
		raise ValueError('%s: %d non-finite transmission value(s); the first '
		                 'at wavelength index %d (%.6f A) is %r'
		                 % (line, int(bad.sum()), i, 1.0e10*float(lam[i]), prob[i]))
	far = (prob < -tol) | (prob > 1.0 + tol)
	if far.any():
		i = int(np.flatnonzero(far)[0])
		raise ValueError('%s: transmission %r outside [0, 1] at wavelength '
		                 'index %d (%.6f A)' % (line, prob[i], i, 1.0e10*float(lam[i])))
	if clip:
		n_out = int(np.count_nonzero((prob < 0.0) | (prob > 1.0)))
		if n_out > 0:
			print('   (%s) %d transmission value(s) within %.0e of [0, 1] set '
			      'to the bound' % (line, n_out, tol))
		prob = np.clip(prob, 0.0, 1.0)
	return prob


def _fig_name(line):
	"""Figure file name for a line key ('' = do not save this figure)."""
	return (os.path.join(path, _fig_prefix + line + '.png')
	        if len(_fig_prefix) > 0 else '')


fig_name_hei  = _fig_name('He10830')
fig_name_lya  = _fig_name('Lya')
fig_name_ha   = _fig_name('Halpha')
fig_name_hb   = _fig_name('Hbeta')
fig_name_mgii = _fig_name('MgII')
fig_name_caii = _fig_name('CaII')
fig_name_nai  = _fig_name('NaI')
fig_name_oi   = _fig_name('OI')

# ----- Parameters read from input.inp ----- #
# Every read is matched by LABEL (read_input_params -> input_read.f90
# semantics), so the header line order is irrelevant, matching the Fortran
# core block. Rp/Mp/T0/a_orb/Mstar keep their legacy word positions; the
# stellar radius / Teff use the same word positions as the Fortran; LEUV and
# the 2D approximate method stay content-matched (last occurrence wins).
_par = read_input_params(Input_file)
if _par['resolved']:
	print('(EXHALE_transit) Using the wind solver\'s resolved configuration '
	      '(EXHALE_resolved.out): Rp = %.5f R_J, T0 = %.1f K'
	      % (_par['Rp']/RJ, _par['T0']))
else:
	print('(EXHALE_transit) No EXHALE_resolved.out next to input.inp -- '
	      'using input.inp values (stale if base.inp overrides them).')
Rp    = _par['Rp']       # planet radius [m]
Mp    = _par['Mp']       # planet mass [kg]
T0    = _par['T0']       # equilibrium temperature [K]
a_orb = _par['a_orb']    # orbital distance [m]
Mstar = _par['Mstar']    # parent star mass [kg]
LEUV     = _par['LEUV']      # log10 of EUV (Lyman-continuum-band) luminosity [erg/s]
appx_mth = _par['appx_mth']  # EXHALE 2D flux approximation (sets the day-night xi factor)

# ----- Stellar and planet-rotation parameters ----- #
# The star sets the transit normalization (A_star, and the R_star/Rp cap on
# the absorbing annulus) and the diluted-blackbody Balmer continuum that
# photoionizes H(n=2), so it must describe the SAME system as the simulation.
# Each parameter is resolved as
#     EXHALE_TRANSIT_* environment override  >  ./input.inp  >  built-in default
# so an explicit override still wins, but an unattended run now inherits the
# star of the run directory instead of a fixed built-in one.
R_star_Rsun, _src_R_star = parameter_with_source('RSTAR_RSUN',
                                                 _par['R_star_Rsun'], 0.44)
R_star = R_star_Rsun*R_sun                     # stellar radius [m]
# Planet rotation period [days]. It enters as the solid-body spin of the
# atmosphere (rotational Doppler broadening, and the Roche-geometry spin).
# Close-in giants are tidally locked, so the default is the Keplerian orbital
# period built from the orbital distance and the masses in input.inp.
_P_orb = orbital_period_days(a_orb, Mstar, Mp)
rot_period, _src_rot_period = parameter_with_source(
    'ROTP', _P_orb, _P_orb,
    input_label='input.inp (tidally locked: P_orb from a, M_star + M_p)',
    default_label='input.inp (tidally locked: P_orb from a, M_star + M_p)')

# Instrument spectral resolving power R = lambda/Delta-lambda for the Gaussian
# line-spread convolution.  Each is env-overridable (EXHALE_TRANSIT_RES_* , with
# the TPM_RES_* names still read); defaults below match the instruments named inline.
Instr_res_HeTR = float(_tenv('RES_HETR', '8e4'))  # He I 10830: CARMENES 8e4 / GIANO-B 5e4
# NOTE: the 8e4 default is CARMENES. It is not universal: the LHS 1140 b
# He 10830 transit was taken with WINERED in HIRES-Y mode, R = 68,000
# (Cherubim et al. 2026, Supplement). Runs for that target set
# EXHALE_TRANSIT_RES_HETR=68000 by sourcing LHS1140b/winered_hires_y.sh.
# Leave the default alone --
# the stored results for the other planets were produced at 8e4.
Instr_res_HI   = float(_tenv('RES_HI',   '5e4'))  # Ly-alpha 1215.67: HST-STIS ~ 1e4-1e5
Instr_res_Ha   = float(_tenv('RES_HA',   '1.15e5'))  # H-alpha 6562.8: HARPS/CARMENES-VIS ~ 1.1e5
Instr_res_Hb   = float(_tenv('RES_HB',   '1.15e5'))  # H-beta 4861.35
Instr_res_MgII = float(_tenv('RES_MGII', '3e4'))  # Mg II h&k 2796/2803 NUV: HST/STIS ~ 3e4
Instr_res_CaII = float(_tenv('RES_CAII', str(Instr_res_Ha)))  # Ca II H&K optical
Instr_res_NaI  = float(_tenv('RES_NAI',  str(Instr_res_Ha)))  # Na I D optical
# O I 1302/1304/1306 FUV. The published HD 209458 b measurement (Vidal-Madjar
# et al. 2004, ApJ 604, L69) used HST/STIS G140L, for which that paper quotes
# only a resolution of "~2.5 A" (R ~ 520 at 1302 A) while Ballester &
# Ben-Jaffel (2015, ApJ 804, 116) give R ~ 1000 for the same grating. R = 1000
# is the default here; the band-integrated depth this is compared against is
# insensitive to the choice (the band is 10 A wide, four instrument elements).
Instr_res_OI   = float(_tenv('RES_OI',   '1.0e3'))

# ----- H-alpha (n=2 -> n=3) inputs ----- #
# H-alpha absorption arises from the n=2 hydrogen population. Following
# Christie, Arras & Li (2013, ApJ 772, 144), the 2s/2p populations are
# set by the rate-equilibrium equations (their Eqs. 12-13) including
# Ly-alpha radiative pumping (1s<->2p), with the rate coefficients of the
# solver's own n=2 model (exhale_transit_lib: n2_populations). The Ly-alpha
# mean intensity J_lya(r) is not in the profiles this tool reads (a run
# with the excited-hydrogen model writes its own field to
# output/Excited_H.txt, which is not read here); it is obtained in one of
# two ways:
#
#  (1) If Jlya_file is set to an existing two-column text file, J_lya(r)
#      is read from it:
#         col 1 = r / R_p   (same radial coordinate as EXHALE output)
#         col 2 = J_lya     = Ly-alpha mean intensity J_nu at line
#                             center [erg s^-1 cm^-2 Hz^-1 sr^-1]
#      (linearly interpolated onto the EXHALE grid, clamped at the ends).
#
#  (2) Otherwise J_lya is estimated with the order-of-magnitude estimate of
#      Huang et al. (2017, ApJ 851, 150), Section 2, in the text after
#      their Eq. (6) (page 3), which they state for the peak of the H
#      photoionization:
#         J_lya(r) ~ 0.1 * F_LyC / Dnu_D(r)
#      where Dnu_D = nu_Lya*sqrt(2 kB T/m_p)/c is the local Ly-alpha
#      Doppler width and F_LyC is the DEPOSITED (absorbed) Lyman-continuum
#      flux (each LyC ionization balanced by a recombination -> Ly-alpha
#      photon; Huang+2017 page 3). Applying it at every radius, and the
#      column-integrated F_LyC below, are this tool's choices, not the
#      paper's. F_LyC is computed from the actual
#      stellar input: the incident EUV-band (E>13.6 eV) flux at the
#      planet, 10^LEUV/(4*pi*a^2), times the fraction absorbed by the
#      atmosphere, (1 - exp(-tau_LyC)) with tau_LyC = sigma_LyC * N_HI
#      (the vertical neutral-H column). For the optically-thick atomic
#      layer the absorbed fraction ~ 1, i.e. F_LyC ~ incident stellar LyC.
Jlya_file  = ''            # path to J_lya(r) file; '' => Huang(2017) estimate
F_LyC_override = 0.0       # 0 => auto from stellar LyC; >0 => fixed [erg cm^-2 s^-1]
sigma_LyC      = 6.3e-18   # H photoionization cross section at LyC [cm^2]
# Day-night / 2D flux dilution (xi) applied to the incident stellar LyC,
# matching EXHALE's "2D approximate method": Rate/2 -> 0.5, Rate/4 -> 0.25,
# else 1.0. This is read from input.inp automatically; set xi_override>0
# to force a value (cf. the xi factor of Christie+2013 / Huang+2017).
xi_override    = 0.0       # 0 => auto from input.inp appx_mth
# n=2 photoionization (Balmer continuum) rates Gamma_2s, Gamma_2p [s^-1].
# If T_star > 0 they are estimated from a diluted stellar blackbody
# Balmer continuum (E>3.4 eV) at the orbital distance; otherwise the
# manual values below are used (0 => neglected). T_star = stellar
# effective temperature [K], resolved with the same
# env > input.inp > built-in precedence as R_star. R_star (above) and a_orb
# set the dilution.
T_star, _src_T_star = parameter_with_source('TSTAR', _par['T_star'], 6065.0)
Gamma_2s = 0.0             # used only if T_star <= 0
Gamma_2p = 0.0

# ----- Resolved-parameter report (cf. write_setup_report.f90) ----- #
print('(EXHALE_transit) resolved parameters  [value, source]')
print('  input file          : %s' % Input_file)
print('  planet radius       : %-12.4f R_J     input.inp' % (Rp/RJ))
print('  planet mass         : %-12.4f M_J     input.inp' % (Mp/MJ))
print('  equilibrium T       : %-12.1f K       input.inp' % T0)
print('  orbital distance    : %-12.5f AU      input.inp' % (a_orb/AU))
print('  parent star mass    : %-12.4f M_sun   input.inp' % (Mstar/M_sun))
print('  stellar radius      : %-12.4f R_sun   %s' % (R_star_Rsun, _src_R_star))
print('  stellar Teff        : %-12.1f K       %s' % (T_star, _src_T_star))
print('  planet spin period  : %-12.5f d       %s' % (rot_period, _src_rot_period))
print('  log10 L_EUV         : %-12s erg/s   input.inp'
      % ('%.3f' % LEUV if LEUV is not None else 'not read'))
print('  2D approx. method   : %-12s         input.inp' % (appx_mth or 'none'))
print('  model curves        : %stpm_<line>.txt   %s'
      % (_save_prefix,
         'EXHALE_TRANSIT_SAVE_PREFIX' if _tenv_set('SAVE_PREFIX')
         else 'run directory, canonical name (default)'))
print('  figures             : %s'
      % (os.path.join(path, _fig_prefix) + '<line>.png (+.pdf)   '
         'EXHALE_TRANSIT_FIG_PREFIX' if len(_fig_prefix) > 0
         else 'not saved (set EXHALE_TRANSIT_FIG_PREFIX)'))
print('  resolving powers    : He %.3g, Lya %.3g, Ha %.3g, Hb %.3g,'
      ' MgII %.3g, CaII %.3g, NaI %.3g   (EXHALE_TRANSIT_RES_* overrides)'
      % (Instr_res_HeTR, Instr_res_HI, Instr_res_Ha, Instr_res_Hb,
         Instr_res_MgII, Instr_res_CaII, Instr_res_NaI))
if _par['R_star_Rsun'] is None:
	print('  WARNING: "Stellar radius [R_sun]:" absent from %s -- the transit'
	      ' normalization uses a built-in radius that need not match this'
	      ' system.' % Input_file)
if _par['T_star'] is None:
	print('  WARNING: "Stellar Teff [K]:" absent from %s -- the n=2'
	      ' photoionization rate uses a built-in temperature that need not'
	      ' match this system.' % Input_file)
print('')

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

# ------------------------------ #

# ------------------------- #

# Load profiles.  The seven physical columns come first and are named by
# position; the row-validity columns after them are read by name from the
# file's own header.
r,rho,v,p,T,heat,cool = loadtxt_cells(Hydro_file, usecols = range(7),
                                      unpack = True)
# The validity of the rows the spectrum is built from.  The advection
# post-process REFUSES the steady correction where an assumption of the steady
# equations fails in a cell, and such a row carries the run's own temperature
# and its equilibrium composition, so a spectrum built from it is a spectrum
# of the uncorrected state there -- except for a row that was an unknown of a
# rejected column energy solve (the transport terms), whose temperature is the
# MARCHING profile of the post-process, a number no equation was solved for.
# read_adv_validity reads the profile's schema: schema 2 states temperature
# and composition validity separately (a row is refused here when its
# temperature field is nonzero or its composition field is failed, and a
# retained composition is reported beside it), schema 1 states one mixed
# field, and a profile whose header states neither leaves the validity of
# every row UNKNOWN, which is not the same as no row refused.  A corrected row
# is a CONDITIONAL correction, accurate to the fraction of itself in the mass
# flux that the profile's own adv_conditional_tol line states, and the measure
# of each row is its adv_mass_row column: both travel into the metadata of
# every saved curve, so the accuracy the spectrum inherits is on the file.  The census of
# the rows the chords sample, and the share of line-center optical depth those
# rows carry, is printed once the chord arrays and every line's lower-level
# density exist (search for CONTRIBUTION diagnostic).
_adv          = read_adv_validity(Hydro_file)
# How far the run's other state stands from the one this curve is built on,
# measured on the three quantities a He I 10830 or a Balmer curve integrates.
# The other state is only being located here, so it is not required to
# exist: a stationary solve has no `_adv` products and that is not an error.
_state_other  = transit_state_files(
    path, 'solution' if _state_selection == 'adv' else 'adv', require=False)
_state_record = {
    'selection': _state_selection,
    'hydro': Hydro_file, 'ioniz': Ioniz_file,
    # A spectrum is traceable to the numbers it was built from only through
    # the digest of the files: the path is reused by the next run.
    'identity': [file_identity(Hydro_file), file_identity(Ioniz_file)],
    # What the state says about itself, copied out of its own header.
    'provenance': (state_provenance_statements(Hydro_file)
                   + [q for q in state_provenance_statements(Ioniz_file)
                      if q not in state_provenance_statements(Hydro_file)]),
    'difference': state_pair_difference(Hydro_file, Ioniz_file,
                                        _state_other[0], _state_other[1])}
print('(TPM) state read: %s (%s, %s)'
      % (_state_selection, os.path.basename(Hydro_file),
         os.path.basename(Ioniz_file)))
for _q in _state_record['identity']:
    print('(TPM)   %s' % _q)
for _q in _state_record['provenance']:
    print('(TPM)   state %s' % _q)
if _state_record['difference'] is not None:
    _d = _state_record['difference']
    print('(TPM) distance to the %s state over %d rows, largest relative: '
          'T %.3e, n(H I) %.3e, n(He 2^3S) %.3e'
          % (_state_other[2], _d['rows'], _d['T'], _d['HI'],
             _d['HeI_2_3S']))
else:
    print('(TPM) the %s state of this run was not read; how far the two '
          'stand apart is UNKNOWN' % _state_other[2])
_adv_refused  = _adv['refused']
_adv_T_status = _adv['T_status']
_adv_comp     = _adv['comp_status']
# THE DERIVED PRODUCT AS A WHOLE.  The `_adv` profile carries a record of how
# the solves that produced it ended (the column energy solve of the transport
# terms, the chemistry, the outer iteration).  A profile whose record says
# rejected is not an advection-corrected state, and no spectrum is made from
# it as an ordinary result; EXHALE_TRANSIT_DIAGNOSTIC=1 makes it anyway, as
# diagnostic data, with the refusal written into the metadata of every saved
# curve.  A profile without the record (written before it existed) has an
# UNKNOWN derived state, which is not a pass either: it is refused the same
# way (exit 7) unless EXHALE_TRANSIT_DIAGNOSTIC=1, and then the metadata says
# the curve stands on a profile that states no derived state.
_adv_verdict = (adv_derived_state_verdict(_adv) if _state_selection == 'adv'
                else None)
if _adv_verdict is not None:
	print('(TPM) derived state of the adv profile: %s'
	      % adv_derived_state_text(_adv))
	if _adv_verdict == 'rejected':
		if _tenv('DIAGNOSTIC', '0').strip() == '1':
			print('(TPM) the adv product is REJECTED by its own record; '
			      'EXHALE_TRANSIT_DIAGNOSTIC=1, so the spectrum is made as '
			      'DIAGNOSTIC data and says so in its metadata.')
		else:
			print('(TPM) the adv product is REJECTED by its own record, so no '
			      'spectrum is made from it (EXHALE_TRANSIT_DIAGNOSTIC=1 makes '
			      'one as diagnostic data).')
			# 7, the status the evaluation itself exits with for the same
			# verdict (EXHALE_main, evaluate route).
			sys.exit(7)
	elif _adv_verdict == 'unknown':
		print('(TPM) the adv profile states no derived-state record (a file '
		      'written before the record existed): its derived state is '
		      'UNKNOWN, which is not a pass.')
		if _tenv('DIAGNOSTIC', '0').strip() == '1':
			print('(TPM) EXHALE_TRANSIT_DIAGNOSTIC=1, so the spectrum is made '
			      'as DIAGNOSTIC data and says so in its metadata.')
		else:
			print('(TPM) no spectrum is made from it '
			      '(EXHALE_TRANSIT_DIAGNOSTIC=1 makes one as diagnostic data).')
			sys.exit(7)
# Ion_species.txt: read only the first 7 columns (r + H/He). In EXHALE
# this file also carries trace-metal columns (C/N/O), so we slice rather
# than unpack all of them.
r,nhi,nhii,nhei,nheii,nheiii,nheiTR = \
    loadtxt_cells(Ioniz_file, usecols = range(7), unpack = True)

# Metal- (and molecular-) ion electron donors, so the free-electron density
# below is not metal-blind. Ion_species.txt carries the 27 trace-metal ion
# columns right after the 7 H/He columns (canonical species_table order), and,
# for a molecular run, 4 more molecular columns (H2 H2+ H3+ HeH+). We sum each
# ion's net charge (= electrons released); neutral stages contribute nothing.
# Metals-off files (only the 7 H/He columns) keep ne_metal = 0, so the legacy
# H/He-only electron count is reproduced.
_ncol_ion = first_data_row_ncol(Ioniz_file)
# Net ionic charge per metal column (C,O,N,Mg,Si,Ca,Fe: 0/1/2; Na,K,S: 0/1).
_metal_charge = [0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,
                 0, 1,     0, 1,     0, 1,     0, 1, 2]
ne_metal_cm = np.zeros_like(nhi)
# H2+ and HeH+, whose dissociative recombination leaves one H atom in n = 2
# (the chemical source of n2_populations); zero in an atomic run.
nh2p_cm  = np.zeros_like(nhi)
nhehp_cm = np.zeros_like(nhi)
if _ncol_ion >= 7 + len(_metal_charge):
    _nm = loadtxt_cells(Ioniz_file,
                        usecols = range(7, 7 + len(_metal_charge)), unpack = True)
    for _ic, _z in enumerate(_metal_charge):
        if _z > 0:
            ne_metal_cm = ne_metal_cm + _z*_nm[_ic]
    # Molecular ions (H2+, H3+, HeH+ each release one electron), if present.
    if _ncol_ion >= 7 + len(_metal_charge) + 4:
        _c0 = 7 + len(_metal_charge)
        _h2p, _h3p, _hehp = loadtxt_cells(
            Ioniz_file, usecols = (_c0 + 1, _c0 + 2, _c0 + 3), unpack = True)
        ne_metal_cm = ne_metal_cm + _h2p + _h3p + _hehp
        nh2p_cm, nhehp_cm = _h2p, _hehp

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
# (_line_halfwidth / _apply_window / _odd are imported; the run's T and v
# profiles are passed to _line_halfwidth explicitly.)

# vertical columns of the lower-level absorbers [cm^-2] (rectangle rule; r in R_p)
_dr_cm = np.abs(np.gradient(r))*Rp*1.0e2       # cm (r in R_p, Rp in m -> *1e2)
_N_HI    = float(np.sum(nhi*_dr_cm))           # cm^-3 * cm -> cm^-2
_N_HeTR  = float(np.sum(nheiTR*_dr_cm))
_fHe = f10830_34 + f10830_25 + f10829_09
_he_forced = _tenv_set('HE_LMIN') or _tenv_set('HE_LMAX')
lmin_HeTR, lmax_HeTR = _apply_window(10830.34,
    _line_halfwidth(10830.34, 4.0, _fHe, 1.022e7, _N_HeTR, T, v), lmin_HeTR, lmax_HeTR, _he_forced)
lmin_HI, lmax_HI = _apply_window(1215.67,
    _line_halfwidth(1215.67, 1.0, 0.4162, 6.27e8, _N_HI, T, v), lmin_HI, lmax_HI)
# Balmer lines are optically thin in the n=2 population -> kinematic only
lmin_Ha, lmax_Ha = _apply_window(6562.80,
    _line_halfwidth(6562.80, 1.0, 0.6407, 4.41e7, 0.0, T, v), lmin_Ha, lmax_Ha)
lmin_Hb, lmax_Hb = _apply_window(4861.35,
    _line_halfwidth(4861.35, 1.0, 0.1193, 8.42e6, 0.0, T, v), lmin_Hb, lmax_Hb)

# adequate, ODD sampling after widening (astropy convolution needs odd kernels)
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
# Electron density [cm^-3] (H + He + metal/molecular ion contributions).
# n1s = neutral H.
ne_cm  = nhii + nheii + 2.0*nheiii + ne_metal_cm
n1s_cm = nhi
do_Ha = True   # H-alpha is always computed (J_lya from file or estimate)
if do_Ha:
	if (len(Jlya_file) > 0) and os.path.exists(Jlya_file):
		# (1) Read J_lya(r) from file.  Plain np.loadtxt: this is an
		#     externally computed J_lya(r) table (the LaRT Monte Carlo
		#     coupling), not a profile write_output.f90 emitted, so it has
		#     no ghost rows; it is interpolated onto r anyway.
		rj, Jlya_in = np.loadtxt(Jlya_file, usecols = (0,1), unpack = True)
		Jlya = np.interp(r, rj, Jlya_in)      # clamps to endpoints outside range
		jlya_src = 'file ' + Jlya_file
	else:
		if len(Jlya_file) > 0:
			print('(TPM) WARNING: Jlya_file "%s" not found; '
			      'using Huang(2017) estimate.' % Jlya_file)
		# (2) Estimate of Huang et al. 2017, Section 2, text after Eq. (6):
		#     J_lya ~ 0.1 * F_LyC / Dnu_D, with Dnu_D the local Ly-alpha
		#     Doppler width (proton thermal speed) and F_LyC the deposited
		#     (absorbed) Lyman-continuum flux.
		# --- F_LyC = absorbed stellar LyC flux [erg cm^-2 s^-1] ---
		if F_LyC_override > 0.0:
			F_LyC = F_LyC_override
			print('(TPM)   F_LyC = %.3e erg/cm2/s (manual override)' % F_LyC)
		elif LEUV is not None:
			# Day-night / 2D dilution (xi), matching EXHALE's appx method
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
			F_LyC = 1.0e4   # used only if LEUV is unavailable
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
	# The solver's rate coefficients (exhale_transit_lib: n2_populations),
	# with the H II recombination set the run used; the field (J_lya and
	# the Balmer continuum) is this tool's, as stated above.
	n2s_cm, n2p_cm, n2_cm = n2_populations(T, n1s_cm, nhii, nheii, nheiii,
	                                       ne_cm, Jlya,
	                                       G2s = Gamma_2s, G2p = Gamma_2p,
	                                       nH2p = nh2p_cm, nHeHp = nhehp_cm,
	                                       rate_set = _par['h_rate_set'])
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
	# _turb_factor() is 1 unless EXHALE_TRANSIT_TURB=1, which adds the
	# Lampon et al. (2020) turbulence term the p-winds comparison uses.
	_tf = _turb_factor()
	v_th_HeTR = np.sqrt(2.0*kb*T_LOS/mHe)*_tf
	v_th_HI = np.sqrt(2.0*kb*T_LOS/mp)*_tf
	v_th_D = np.sqrt(2.0*kb*T_LOS/mD)*_tf
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

# The transmission is checked before use (_checked_transmission): a
# non-finite value stops the synthesis, and only the rounding of the disk
# average is clipped to [0, 1].
avg_prob_HD = _checked_transmission('Lya', avg_prob_HD, l_onde_HI, clip=True)

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
avg_prob_rot_HeTR = _checked_transmission('He 10830 (rotated)',
                                          avg_prob_rot_HeTR, l_onde_HeTR)
convolved_rot_prob_HeTR = convolve(avg_prob_rot_HeTR, gaussian_HeTR, boundary='extend')

# Hydrogen and Deuterium (Ly-alpha), checked as above
avg_prob_rot_HD = _rotate_disk_average(exp_tau_HD, l_onde_HI)
avg_prob_rot_HD = _checked_transmission('Lya (rotated)', avg_prob_rot_HD,
                                        l_onde_HI, clip=True)
convolved_rot_prob_HD = convolve(avg_prob_rot_HD, gaussian_HD, boundary='extend')

# H-alpha / H-beta
if do_Ha:
	avg_prob_rot_Ha = _rotate_disk_average(exp_tau_Ha, l_onde_Ha)
	avg_prob_rot_Ha = _checked_transmission('H-alpha (rotated)',
	                                        avg_prob_rot_Ha, l_onde_Ha)
	convolved_rot_prob_Ha = convolve(avg_prob_rot_Ha, gaussian_Ha, boundary='extend')
	avg_prob_rot_Hb = _rotate_disk_average(exp_tau_Hb, l_onde_Hb)
	avg_prob_rot_Hb = _checked_transmission('H-beta (rotated)',
	                                        avg_prob_rot_Hb, l_onde_Hb)
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
mO = 15.999*amu

# Metals-off runs write only the H/He columns (<=16): skip the metal
# resonance lines automatically (He/Lya/Ha/Hb above are unaffected).  The
# width comes from the first data row, so the comment header is not part of
# the count and no row selection is involved.
do_metals = (_ncol_ion >= 26)
if do_metals:
	nMgII_cm, nCaII_cm, nNaI_cm = loadtxt_cells(Ioniz_file,
	                                            usecols=(17, 23, 25),
	                                            unpack=True)          # cm^-3
else:
	_r0col = loadtxt_cells(Ioniz_file, usecols=(0,))
	nMgII_cm = np.zeros_like(_r0col)   # metals-off: zero metal absorption,
	nCaII_cm = np.zeros_like(_r0col)   # so the metal lines are flat and the
	nNaI_cm  = np.zeros_like(_r0col)   # He/Lya/Ha output still saves.
	print('(TPM) metals-off run: skipping metal resonance lines.')
# symmetric (night + day) chord arrays, m^-3 (mirror data_nHI at line ~460)
data_nMgII = np.concatenate((np.flip(nMgII_cm), nMgII_cm))*1.0e6
data_nCaII = np.concatenate((np.flip(nCaII_cm), nCaII_cm))*1.0e6
data_nNaI  = np.concatenate((np.flip(nNaI_cm),  nNaI_cm ))*1.0e6

# ----- O I 1302/1304/1306: the three ground-term fine-structure levels ----- #
# The O I resonance triplet does NOT come out of one lower level: 1302.168,
# 1304.858 and 1306.029 A absorb out of the 3P2, 3P1 and 3P0 levels of the
# 2p4 3P ground term respectively.  Applying the total O I density to all
# three components would count the same atoms three times over -- at the
# near-statistical populations these winds carry that is a factor 3 in the
# total column and a factor 1.8 / 3.1 / 9.3 in the three components.
# So the populations are read from the wind solver, which writes them out of
# the same three-level statistical equilibrium its [O I] 63/145/44um cooling
# is built on (the OI_levels file of the selected state, Cool_coeff.f90).
# A run whose output
# predates that file, or a metals-off run that never wrote it, simply skips
# the line the way the other metal lines are skipped.
# It follows the state selection: the level populations a line integrates
# and the density and temperature it integrates them with come from one
# state, never from two.
OI_file = transit_state_oi_levels(path, _state_selection)
do_OI = do_metals and os.path.exists(OI_file)
if do_OI:
	# cols 8,9,10 (0-indexed) = n(3P2), n(3P1), n(3P0) in cm^-3
	n3P2_cm, n3P1_cm, n3P0_cm = loadtxt_cells(OI_file, usecols=(8, 9, 10),
	                                          unpack=True)
	_nOI_file = loadtxt_cells(OI_file, usecols=(4,))
	if n3P2_cm.size != nMgII_cm.size:
		raise ValueError('(EXHALE_transit) %s has %d rows but %s has %d; '
		                 'they must be the same grid'
		                 % (OI_file, n3P2_cm.size, Ioniz_file, nMgII_cm.size))
	# Gate, re-checked here rather than trusted: the three level densities
	# must add up to the total O I of the same file.
	_lvl_sum = n3P2_cm + n3P1_cm + n3P0_cm
	_lvl_err = np.max(np.abs(_lvl_sum - _nOI_file)
	                  /np.maximum(_nOI_file, 1e-99))
	print('(TPM) O I ground-term levels: max |sum(levels)/n(O I) - 1| = %.2e'
	      % _lvl_err)
	if _lvl_err > 1.0e-10:
		raise ValueError('(EXHALE_transit) O I level densities do not sum to '
		                 'the total O I in %s (max rel. error %.3e)'
		                 % (OI_file, _lvl_err))
	# Column-regime knob. The O I triplet is a resonance line whose three
	# components span a factor 6.6 in gf and are optically thick in the
	# launch region, so the forward model has to be exercised in BOTH
	# regimes: EXHALE_TRANSIT_OI_NSCALE multiplies all three level
	# densities, leaving their ratios (and therefore the level physics)
	# alone. 1e-4 puts every component on the linear part of the curve of
	# growth, 1e4 saturates all three. Default 1.0 = the run's own densities.
	_oi_nscale = float(_tenv('OI_NSCALE', '1.0'))
	if _oi_nscale != 1.0:
		print('(TPM) O I level densities scaled by %.3e '
		      '(EXHALE_TRANSIT_OI_NSCALE): a curve-of-growth test, not a '
		      'physical model' % _oi_nscale)
	data_n3P2 = np.concatenate((np.flip(n3P2_cm), n3P2_cm))*1.0e6*_oi_nscale
	data_n3P1 = np.concatenate((np.flip(n3P1_cm), n3P1_cm))*1.0e6*_oi_nscale
	data_n3P0 = np.concatenate((np.flip(n3P0_cm), n3P0_cm))*1.0e6*_oi_nscale
else:
	_r0col_oi = loadtxt_cells(Ioniz_file, usecols=(0,))
	_zero_oi = np.zeros(2*_r0col_oi.size)
	data_n3P2 = _zero_oi
	data_n3P1 = _zero_oi
	data_n3P0 = _zero_oi
	if do_metals:
		print('(TPM) no %s (run predates the O I level output): skipping '
		      'the O I 1302 triplet.' % OI_file)


# resonance_depth is imported from exhale_transit_lib; the compute-specific
# grid/area arrays are passed to it explicitly at the call site below.

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
	d_lc, d_band, _, _ = resonance_depth(lam_, f_, A_, m_, ncol_, R_,
	                                     Grid_Number, r_grid, Rp, data_r,
	                                     data_v, data_T, A_star, A_atm, A_planet)
	metal_depth[lbl] = (d_lc, d_band)
	print('(TPM)   %-13s   %10.3f      %10.3f' % (lbl, d_lc, d_band))
print('')


# --------------------------------------------------------------------- #
# Full doublet transmission spectra for the metal resonance lines.
# The same spherical chord integration as resonance_depth, but with BOTH doublet
# components summed in one wavelength window, and with the instrument
# and planet-rotation convolutions applied exactly as for the H/He
# lines, so the metal lines are first-class TPM outputs, written under
# the same tpm_<line> naming. The Phase 5a table above keeps the validated
# single-component numbers. Skipped automatically for a metals-off run
# (all-zero ion columns). NIST atomic data (Kramida 2020).

METAL_DOUBLETS = [
	# key, label, components (lam0_A, f, A21, lower-level chord density
	#   [m^-3]), mass, instrument R, window [A], nlam, figure name.
	# The lower-level density belongs to the COMPONENT: for these doublets
	# both components absorb out of the ion ground state, so the same array
	# appears twice; the O I triplet below is where that stops being true.
	('MgII', 'Mg II h&k',
	 [(2796.352, 0.608, 2.60e8, data_nMgII),
	  (2803.531, 0.303, 2.57e8, data_nMgII)],
	 mMg, Instr_res_MgII, (2790.0, 2810.0), 601, fig_name_mgii),
	('CaII', 'Ca II H&K',
	 [(3933.663, 0.6267, 1.47e8, data_nCaII),
	  (3968.469, 0.3116, 1.40e8, data_nCaII)],
	 mCa, Instr_res_CaII, (3927.0, 3975.0), 961, fig_name_caii),
	('NaI', 'Na I D',
	 [(5889.951, 0.641, 6.16e7, data_nNaI),
	  (5895.924, 0.320, 6.14e7, data_nNaI)],
	 mNa, Instr_res_NaI, (5884.0, 5902.0), 541, fig_name_nai),
	# O I 2p4 3P -> 3s 3S resonance triplet.  One upper level (3S1, g = 3,
	# E = 76794.978 cm^-1), three lower levels: the 3P2 / 3P1 / 3P0 ground-term
	# fine-structure levels at 0 / 158.265 / 226.977 cm^-1 (g = 5, 3, 1).
	# Vacuum wavelengths, A_ul and f_lu from NIST ASD (accuracy A), retrieved
	# 2026-08-30.  Each component therefore carries ITS OWN lower level.
	# Window 1296-1312 A at 0.02 A brackets the observers' 1300-1310 A band
	# with room for the instrument convolution to be applied before the band
	# average is taken.
	('OI', 'O I 1302 triplet',
	 [(1302.168, 0.0520, 3.41e8, data_n3P2),
	  (1304.858, 0.0518, 2.03e8, data_n3P1),
	  (1306.029, 0.0519, 6.76e7, data_n3P0)],
	 mO, Instr_res_OI, (1296.0, 1312.0), 801, fig_name_oi),
]


# --------------------------------------------------------------------- #
# WHAT THE SPECTRUM IS BUILT FROM: the rows the chords sample, and how much
# line-center optical depth comes from rows whose advective correction was
# refused.
# --------------------------------------------------------------------- #
# A ray at impact parameter b takes every shell with r >= b
# (chord_shell_indices, the same selection every optical-depth integral above
# uses), so the rows the spectrum is built from are the ones the innermost ray
# reaches.  There is no upper radial limit: a shell above the stellar-disk cap
# Rib is outside the ray grid as an impact parameter yet still crosses the rays
# that do fall on the disk, and absorbs there.
#
# A row whose steady advective correction was refused carries the run's own
# temperature and the equilibrium composition at that temperature, so the
# transmission of the chords crossing it is built on the uncorrected state;
# a failed row of a rejected column energy solve carries the marching profile
# of the post-process instead, and only the diagnostic mode reads one.
# How much optical depth those rows carry is a CONTRIBUTION, printed as such:
# the disk-averaged transmission is an average of exp(-tau), nonlinear in tau,
# so a large share is not an error bar on the depth.
#
# This block sits here because it needs the chord arrays and the lower-level
# density of every line, metals and O I included.  The same census travels
# with every saved curve, in the comment metadata block of the file.
#
# What the metadata block of every saved curve reports: the census over the
# sampled rows, and, line by line, the share of the line-center optical depth
# that comes from refused rows.
_census_summary = None
_transit_census = {}
_b_show = np.unique(np.linspace(0, Grid_Number - 1, 8).astype(int))

_sampled   = chord_shell_indices(r, r_grid[0])
_n_sampled = int(_sampled.size)
_n_cap     = int(np.count_nonzero(r[_sampled] > Rib))
print('')
if _adv_refused is None:
	# Legacy profile: it states no row validity, so neither does the census.
	print('(TPM) row validity: the profile states no row validity (no '
	      'adv_schema and no status column), so the validity of the %d rows '
	      'the chords sample is UNKNOWN' % _n_sampled)
	print('      UNKNOWN is not zero refused: this profile does not say which '
	      'rows carry the steady advective correction.')
else:
	_refused_row = _adv_refused
	_n_refused   = int(np.count_nonzero(_refused_row[_sampled]))
	_names       = _adv['status_names']
	print('(TPM) row validity (adv_schema %d): %d of the %d rows the chords '
	      'sample are not the steady advective correction'
	      % (_adv['schema'], _n_refused, _n_sampled))
	print('      (rays 1 <= b <= %.4f Rp; a row is sampled when r >= b for '
	      'one of them, with no upper limit, so the %d rows above the cap '
	      'count too)' % (Rib, _n_cap))
	_reasons = ', '.join(
	    '%s: %d rows' % (_names[_k],
	                     int(np.count_nonzero(_adv_T_status[_sampled] == _k)))
	    for _k in range(1, len(_names))
	    if np.any(_adv_T_status[_sampled] == _k))
	if len(_reasons) > 0:
		print('      temperature: ' + _reasons)
	# What refuses each refused row, over both fields: a failed temperature
	# or composition (2), else the temperature field (retained, unsupported,
	# not_evaluated), so a numerical failure is not counted as a retention.
	_refused_by = ''
	if _adv.get('refused_reason') is not None:
		_rr = np.asarray(_adv['refused_reason'])[_sampled]
		_refused_by = ', '.join(
		    '%s: %d rows' % (_names[_k], int(np.count_nonzero(_rr == _k)))
		    for _k in range(1, len(_names)) if np.any(_rr == _k))
		if len(_refused_by) > 0:
			print('      refused by: ' + _refused_by)
	# The record of the derived product as a whole, next to the row counts.
	if _state_selection == 'adv':
		print('      derived state: ' + adv_derived_state_text(_adv))
	# Schema 2 states the composition validity of a row separately from its
	# temperature: a corrected temperature can sit on a retained composition,
	# and the lower-level densities of every line come from that composition.
	_comp = ''
	if _adv_comp is not None:
		_comp = ', '.join(
		    '%s: %d rows' % (_names[_k],
		                     int(np.count_nonzero(_adv_comp[_sampled] == _k)))
		    for _k in range(len(_names))
		    if np.any(_adv_comp[_sampled] == _k))
		print('      composition: ' + _comp)
	# What a corrected row of this profile is worth: the fraction it is
	# accurate to in the mass flux, and how close to that fraction the rows
	# the chords sample actually came.  A count of corrected rows without
	# their measure says nothing about the accuracy of the spectrum.
	if _adv.get('conditional_tol') is not None:
		_mrow = _adv.get('mass_row')
		if _mrow is not None:
			_mrow = np.asarray(_mrow)
			print('      a corrected row is accurate to %.1e of itself in '
			      'the mass flux; the sampled rows measure up to %.3e'
			      % (_adv['conditional_tol'],
			         float(np.nanmax(_mrow[_sampled]))))
		else:
			print('      a corrected row is accurate to %.1e of itself in '
			      'the mass flux' % _adv['conditional_tol'])
	_census_summary = {'sampled':    _n_sampled,
	                   'refused':    _n_refused,
	                   'above_cap':  _n_cap,
	                   'reasons':    _reasons,
	                   'refused_by': _refused_by,
	                   'comp':       _comp}

	# Same mirrored (night + day) chord array as every density above.
	_data_refused = np.concatenate((np.flip(_refused_row), _refused_row))
	_census_lines = [
		('He10830', 'He I 10830', l_He3_1*1e10,
		 [(l_He3_1*1e10, f10830_34, A12_HeTR, mHe, data_nheiTR),
		  (l_He3_2*1e10, f10830_25, A12_HeTR, mHe, data_nheiTR),
		  (l_He3_3*1e10, f10829_09, A12_HeTR, mHe, data_nheiTR)]),
		('Lya', 'Ly-alpha', lA*1e10,
		 [(lA*1e10, f_la, A12_HI, mp, data_nHI),
		  (lD*1e10, f_D,  A12_D,  mD, data_nD)]),
	]
	if do_Ha:
		# 2s and 2p absorb with different oscillator strengths out of the same
		# n = 2, exactly as the Balmer integrals above.
		_census_lines += [
			('Halpha', 'H-alpha', l_Ha*1e10,
			 [(l_Ha*1e10, f_Ha_2s, A12_Ha, mp, data_n2s),
			  (l_Ha*1e10, f_Ha_2p, A12_Ha, mp, data_n2p)]),
			('Hbeta', 'H-beta', l_Hb*1e10,
			 [(l_Hb*1e10, f_Hb_2s, A12_Hb, mp, data_n2s),
			  (l_Hb*1e10, f_Hb_2p, A12_Hb, mp, data_n2p)]),
		]
	_census_lines += [
		(_kc, _lbl_c, _comps_c[0][0],
		 [(_co[0], _co[1], _co[2], _m_c, _co[3]) for _co in _comps_c])
		for _kc, _lbl_c, _comps_c, _m_c, _Rc, _wc, _nc, _fc
		in METAL_DOUBLETS
		if max(_co[3].max() for _co in _comps_c) > 0.0
	]
	if _n_refused > 0:
		_census_share = refused_line_center_tau_share(
			[(_k_s, _lam_s, _comps_s)
			 for _k_s, _lbl_s, _lam_s, _comps_s in _census_lines],
			r_grid, Rp, data_r, data_v, data_T, _data_refused)
	else:
		# No refused row carries optical depth, so every share is zero and
		# the quadrature need not be repeated to say so.
		_census_share = {_k_s: (None, np.zeros(len(r_grid)))
		                 for _k_s, _l, _la, _c in _census_lines}
	for _k_s, _lbl_s, _lam_s, _comps_s in _census_lines:
		_sh_s = _census_share[_k_s][1]
		_imax = int(np.argmax(_sh_s))
		_transit_census[_k_s] = ([r_grid[_i] for _i in _b_show],
		                         [_sh_s[_i] for _i in _b_show],
		                         float(_sh_s[_imax]), float(r_grid[_imax]))
	if _n_refused > 0:
		print('      CONTRIBUTION diagnostic, not an uncertainty on the '
		      'depth: share of the line-center')
		print('      optical depth each ray takes from those rows '
		      '(transmission is nonlinear in tau).')
		print('        %-17s' % 'ray b [Rp]'
		      + ' '.join('%6.3f' % r_grid[_i] for _i in _b_show)
		      + '   maximum')
		for _k_s, _lbl_s, _lam_s, _comps_s in _census_lines:
			_b_s, _s_s, _max_s, _bmax_s = _transit_census[_k_s]
			print('        %-17s' % _lbl_s
			      + ' '.join('%6.3f' % _v for _v in _s_s)
			      + '   %.3f at b = %.4f' % (_max_s, _bmax_s))
print('')

# resonance_spectrum is imported from exhale_transit_lib; the compute-specific
# grid/area arrays and the rotation disk-average routine are passed to it at
# the call site below.

metal_spec = {}
for key_m, lbl_m, comps_m, m_m, R_m, win_m, nl_m, fnm_m 		in METAL_DOUBLETS:
	if max(co[3].max() for co in comps_m) <= 0.0:
		print('(TPM)   %s: lower-level column is zero (metals off?); skipped'
		      % lbl_m)
		continue
	metal_spec[key_m] = resonance_spectrum(comps_m, m_m, R_m,
	                                       win_m, nl_m,
	                                       Grid_Number, r_grid, Rp, data_r,
	                                       data_v, data_T, A_star, A_atm,
	                                       A_planet, _rotate_disk_average)
	metal_spec[key_m].update(label=lbl_m, comps=comps_m, fig=fnm_m)


# --------------------------------------------------------------------- #
# O I 1302 triplet: the band-integrated comparison with the published
# HD 209458 b measurement.
# --------------------------------------------------------------------- #
# The two published depths come from the SAME four HST/STIS G140L transits:
#
#   Vidal-Madjar et al. 2004, ApJ 604, L69, Table 1:
#       O I / O I* / O I**, 1300-1310 A, 12.8 (+4.5/-4.5) %
#   Ben-Jaffel & Hosseini 2010, ApJ 709, 1284, Table 2:
#       O I, 1299-1310 A, 10.5 +/- 4.4 %
#
# Four things have to be right for the comparison to mean anything, and all
# four are recorded in the metadata written beside the curve.
#
# 1. BAND, NOT LINE CENTER.  STIS G140L does not resolve the triplet -- "the
#    low resolution (~2.5 A) does not allow the stellar emission lines to be
#    resolved" (VM04) -- and both papers quote a depth over a ~10 A window.
#    The line-center depth of an optically thick resonance line is a property
#    of the line profile and the instrument, not of the atmosphere; it is
#    printed here but is never the comparison quantity.  Same trap as Mg II.
#
# 2. THE WEIGHT INSIDE THAT BAND IS THE STELLAR LINE PROFILE, NOT A FLAT
#    AVERAGE.  What was measured is a ratio of fluxes, and essentially all the
#    flux in 1300-1310 A of a G0 V star is in the three narrow chromospheric
#    O I emission lines: Ben-Jaffel & Hosseini measured their "average FWHM
#    ... ~0.2 A (or ~45 km/s)" and their relative peaks, "the lines' peaks in
#    the ratio 1:1.5:1.17, respectively, for the O i (1302.17 A, 1304.86 A,
#    and 1306.03 A) lines".  Averaging the model flatly over 10 A instead
#    dilutes the absorption by the ratio of the band width to the line widths,
#    a factor of order 20.  The stellar-profile weight is therefore the
#    default; EXHALE_TRANSIT_OI_WEIGHT=flat gives the flat average, and both
#    numbers are always reported so the difference is never hidden.
#
# 3. THE INTERSTELLAR MEDIUM is inside those peak ratios and must not be
#    applied twice.  Neither paper removed the ISM.  VM04 argued from it:
#    "The O i ground-level line is strongly absorbed by the interstellar
#    medium.  Therefore the ~13% absorption observed during the transit in the
#    full O i triplet must be due to the presence of O i* and O i** ...".
#    The 1:1.5:1.17 ratios are measured from the observed HD 209458 spectrum,
#    i.e. after the ISM has already eaten into 1302.17 A -- which is exactly
#    why the ground-level component is the WEAKEST of the three there.  So no
#    separate ISM screen is applied here.  The ISM never multiplies the
#    planet's transmission in any case: it absorbs in and out of transit alike
#    and cancels from the ratio; it only reweights which wavelengths the
#    measurement is sensitive to.
#
# 4. THE CONTINUUM.  Both papers quote (R_abs/R_*)^2 from an occultation-curve
#    fit, which includes the opaque planetary disk: their own continuum bands
#    give 2.0 (+0.5/-0.7) % (VM04, 1350-1700 A) and 1.96 +/- 0.42 % (BJ10,
#    1400-1700 A).  This script's T_lambda is normalized so that T = 1 is the
#    planet's disk alone, so (1 - T) is the excess over that continuum, and
#    each published value below has that paper's own continuum subtracted.
OI_BANDS = [
	# label, (lam_lo, lam_hi) [A], published excess-of-continuum depth [%],
	#   its uncertainty [%], source
	('VM04  1300-1310 A', (1300.0, 1310.0), 12.8 - 2.0, 4.5,
	 'Vidal-Madjar et al. 2004, ApJ 604, L69, Table 1 '
	 '(12.8 +/- 4.5 % total, minus their 2.0 % continuum)'),
	('BJ10  1299-1310 A', (1299.0, 1310.0), 10.5 - 1.96, 4.4,
	 'Ben-Jaffel & Hosseini 2010, ApJ 709, 1284, Table 2 '
	 '(10.5 +/- 4.4 % total, minus their 1.96 % continuum)'),
]

# Stellar O I emission-line profile used as the band weight (Ben-Jaffel &
# Hosseini 2010, section 3.2): three Gaussians at the triplet wavelengths with
# this FWHM [A] and these relative peak heights.
OI_STAR_FWHM_A  = float(_tenv('OI_STAR_FWHM_A', '0.2'))
OI_STAR_PEAKS   = [float(x) for x in
                   _tenv('OI_STAR_PEAKS', '1.0,1.5,1.17').split(',')]
OI_WEIGHT_MODE  = _tenv('OI_WEIGHT', 'star').lower()

oi_result = {}
if 'OI' in metal_spec:
	_sp = metal_spec['OI']
	_lam = _sp['l_plot']
	_lam0s = [co[0] for co in _sp['comps']]
	# The curve to integrate is the one WITHOUT the instrument LSF. A Gaussian
	# LSF conserves the integral, so it cancels from a ratio taken over a band
	# many LSF widths wide -- which is what both papers measured. Convolving
	# first and then weighting with the 0.2 A stellar profile would smear the
	# absorption out of the wavelengths that carry the weight and understate
	# the depth by the LSF-to-line width ratio. The wind and rotation
	# broadening IS kept (avg_rot), because that is the atmosphere.
	_Tcmp = _sp['avg_rot']
	_sig = OI_STAR_FWHM_A/(2.0*np.sqrt(2.0*np.log(2.0)))
	_wstar = np.zeros_like(_lam)
	for _l0, _pk in zip(_lam0s, OI_STAR_PEAKS):
		_wstar += _pk*np.exp(-0.5*((_lam - _l0)/_sig)**2.0)
	if OI_WEIGHT_MODE == 'flat':
		_wgt = np.ones_like(_lam)
		_wnote = ('weight = FLAT over the band (EXHALE_TRANSIT_OI_WEIGHT=flat). '
		          'This dilutes the model by the band-to-line width ratio and '
		          'is not what the published depths measure.')
	else:
		_wgt = _wstar
		_wnote = ('weight = stellar O I emission profile, three Gaussians of '
		          'FWHM %.3f A with peak ratios %s (Ben-Jaffel & Hosseini '
		          '2010, ApJ 709, 1284, section 3.2; those ratios already '
		          'carry the interstellar absorption of the 1302.17 A '
		          'component, so no separate ISM screen is applied)'
		          % (OI_STAR_FWHM_A,
		             ':'.join('%g' % p for p in OI_STAR_PEAKS)))
	print('')
	print('(TPM) ===== O I 1302 triplet vs. the published HD 209458 b depth =====')
	print('(TPM)   %s' % _wnote)
	print('(TPM)   band                 model [%]  (flat)   published [%]')
	for _lbl, _bnd, _obs, _err, _src in OI_BANDS:
		_sel = (_lam >= _bnd[0]) & (_lam <= _bnd[1])
		_d_flat = band_integrated_depth(_lam, _Tcmp, _bnd)
		_d_w = (1.0 - np.sum(_wgt[_sel]*_Tcmp[_sel])
		              /np.sum(_wgt[_sel]))*100.0
		oi_result[_lbl] = dict(band=_bnd, model_flat=_d_flat,
		                       model_weighted=_d_w, obs=_obs, obs_err=_err,
		                       source=_src)
		print('(TPM)   %-20s %10.4e (%10.4e)  %6.2f +/- %.2f   %s'
		      % (_lbl, _d_w, _d_flat, _obs, _err,
		         'consistent' if abs(_d_w - _obs) <= _err else 'DISCREPANT'))
	print('(TPM)   line-center depth (NOT the comparison quantity): %.4e %%'
	      % _sp['Tl_conv'])
	print('')
	_sp.update(bands=oi_result, weight_note=_wnote)


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
		# The impact parameters end where the annulus of the normalization
		# below ends, R_ib = min(r_N, R_star/R_p): the depth is
		# (A_atm - A_planet) times the mean transmission of that annulus,
		# and chords beyond the stellar radius do not cross the disk.
		b_top  = min(r_hi, Rib)
		b_grid = np.array([b_top**(j/(n_b-1)) for j in range(n_b)])  # 1 -> R_ib
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
				vth = np.sqrt(2.0*kb*T_l/mass)*_turb_factor()
				vlos = vw*x_los/r3 - Omega*(y*Rp)            # m/s
				Dnu = nu0*vth/c_light
				a_v = A21/(4.0*np.pi*Dnu)
				X   = (nu_l[None, :] - nu0)/Dnu[:, None]
				Vo  = f_osc*Fadd_const/Dnu[:, None] \
				      * wofz(X - (vlos/vth)[:, None] + 1j*a_v[:, None]).real
				tau = np.trapz(n_l[:, None]*Vo, x=x_phys, axis=0)
				radial[ia, ib, :] = np.exp(-np.abs(tau))
		# azimuthal mean of the radial (impact-parameter) integral / (R_ib^2-1)
		prob = np.array([np.mean([np.trapz(2.0*b_grid*radial[ia, :, l], b_grid)
		                          for ia in range(n_sectors)])
		                 for l in range(nlam)])/(b_top**2 - 1.0)
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

# ----- Save the model transmission curves ----- #
# One file per line, named <save prefix>tpm_<line>.txt (see the output-naming
# block at the top).  Columns: wavelength [A], T_lambda (theoretical),
# T_lambda (instr. conv.), T_lambda (planet-rot + instr. conv.).
# Excess absorption [%] = (1 - T)*100.  A metals-off run simply has no metal
# lines to write.
_curves = [('He10830', 'He I 10830', l_plot_HeTR, avg_prob_HeTR,
            convolved_avg_prob_HeTR, convolved_rot_prob_HeTR),
           ('Lya', 'Ly-alpha 1215.67', l_plot_HI, avg_prob_HD,
            convolved_avg_prob_HD, convolved_rot_prob_HD)]
if do_Ha:
	_curves += [('Halpha', 'H-alpha 6562.8', l_plot_Ha, avg_prob_Ha,
	             convolved_avg_prob_Ha, convolved_rot_prob_Ha),
	            ('Hbeta', 'H-beta 4861.35', l_plot_Hb, avg_prob_Hb,
	             convolved_avg_prob_Hb, convolved_rot_prob_Hb)]
_curves += [(key_m, sp['label'], sp['l_plot'], sp['avg'],
             sp['conv'], sp['conv_rot'])
            for key_m, sp in metal_spec.items()]

# Every saved file carries, as comments, what it was built from and how much
# of that was the uncorrected state, so a curve and its validity never travel
# separately.  Comments only: the four numerical columns are unchanged, and the
# column line stays immediately above the data.
try:
	_tool_dir   = os.path.dirname(os.path.realpath(__file__))
	_tool_paths = [os.path.realpath(__file__),
	               os.path.join(_tool_dir, 'exhale_transit_lib.py')]
except NameError:                       # __file__ is absent under exec()
	_tool_paths = ['EXHALE_transit.py', 'exhale_transit_lib.py']
_tool_identity = transit_tool_identity(_tool_paths)
_overrides     = transit_environment_overrides()

for _key, _lbl, _lam, _t0, _t1, _t2 in _curves:
	_meta = transit_metadata_block(_adv, _tool_identity, _overrides,
	                               _census_summary,
	                               _transit_census.get(_key),
	                               state=_state_record)
	np.savetxt(_save_prefix + 'tpm_%s.txt' % _key,
	           np.c_[_lam, _t0, _t1, _t2],
	           header='\n'.join(
	               _meta
	               + ['lambda[A]  T_theo  T_instr  T_rot+instr  (%s)' % _lbl]))
print('(TPM) saved model curves: %stpm_{%s}.txt'
      % (_save_prefix, ','.join(k for k, _, _, _, _, _ in _curves)))

# ----- O I 1302: the validation record, with its treatment stated ----- #
# Written next to the curve so that the number and the conditions under which
# it may be compared with the published measurement never travel separately.
if 'OI' in metal_spec and len(oi_result) > 0:
	with open(_save_prefix + 'tpm_OI_band_depths.txt', 'w') as _fh:
		for _ml in transit_metadata_block(_adv, _tool_identity, _overrides,
		                                  _census_summary,
		                                  _transit_census.get('OI'),
		                                  state=_state_record):
			_fh.write('# %s\n' % _ml)
		_fh.write('# O I 1302.168/1304.858/1306.029 A band-integrated transit '
		          'depths, HD 209458 b comparison\n')
		_fh.write('# Lower levels: the 3P2/3P1/3P0 ground-term fine-structure '
		          'populations from %s,\n' % os.path.basename(OI_file))
		_fh.write('#   i.e. the same three-level statistical equilibrium the '
		          '[O I] 63/145/44um cooling uses.\n')
		_fh.write('# Atomic data: NIST ASD (accuracy A), vacuum wavelengths.\n')
		_fh.write('# Instrument: Gaussian LSF at R = %.4g (HST/STIS G140L; '
		          'EXHALE_TRANSIT_RES_OI). The band-integrated depths below '
		          'are taken\n' % Instr_res_OI)
		_fh.write('#   from the curve WITHOUT that LSF: a Gaussian LSF '
		          'conserves the integral, so it cancels from\n')
		_fh.write('#   a ratio over a band many LSF widths wide. The wind and '
		          'rotation broadening is kept.\n')
		_fh.write('# Normalization: T = 1 is the opaque planetary disk, so '
		          '(1-T) is the EXCESS over the\n')
		_fh.write('#   continuum. The published depths are (R_abs/R_*)^2 and '
		          'INCLUDE the disk, so each\n')
		_fh.write('#   published value below has that paper\'s own continuum '
		          'band subtracted.\n')
		_fh.write('# Geocoronal O I 1302: subtracted by both papers in the '
		          'cross-dispersion direction, with\n')
		_fh.write('#   no spectral pixels masked, so no wavelength is excluded '
		          'from the band here either.\n')
		_fh.write('# Band weight: %s\n' % _sp['weight_note'])
		_fh.write('# Not included in the model: Ly-beta pumping of O I '
		          '1025.76 A (Bowen 1947, PASP 59, 196),\n')
		_fh.write('#   which feeds 3d 3D and cascades through 11287 and '
		          '8446 A onto the 3s 3S upper level of\n')
		_fh.write('#   this triplet and back onto all three 3P levels. Its '
		          'size is bracketed by the two\n')
		_fh.write('#   published regimes and by nothing closer: ~10% of direct '
		          '1304 excitation in the Earth\n')
		_fh.write('#   dayglow (Meier 1991, SSRv 58, 1, p. 99) but ~20x '
		          'collisional excitation in the solar\n')
		_fh.write('#   chromosphere (Skelton & Shine 1982, ApJ 259, 869). No '
		          'published calculation exists for\n')
		_fh.write('#   an escaping exoplanet atmosphere.\n')
		_fh.write('# columns: band_lo[A] band_hi[A] model_weighted[%] '
		          'model_flat[%] published[%] published_err[%] label\n')
		for _lbl, _r in oi_result.items():
			_fh.write('%9.3f %9.3f %14.6e %14.6e %10.3f %10.3f   %s | %s\n'
			          % (_r['band'][0], _r['band'][1], _r['model_weighted'],
			             _r['model_flat'], _r['obs'], _r['obs_err'],
			             _lbl, _r['source']))
	print('(TPM) saved O I band comparison: %stpm_OI_band_depths.txt'
	      % _save_prefix)

# ----- He 10830 line metrics (three-Gaussian, Cherubim et al. 2026) ----- #
# Astrophysical metrics of the modeled He triplet: blended-red depth, blue
# depth, red/blue amplitude ratio, FWHM of the blended feature, and the
# shared Doppler shift.  Fit on the planet-rotation and instrument-convolved
# curve (column 4 of tpm_He10830.txt, the one the equivalent-width and
# width-matching tools read; before 2026-10-06 the instrument-only column 3),
# in the AIR wavelength frame this script uses.  A fit failure is reported,
# not fatal.
try:
	from he_line_metrics import fit_metrics as _he_fit_metrics
	_he_excess = (convolved_rot_prob_HeTR.max() - convolved_rot_prob_HeTR) \
	    / convolved_rot_prob_HeTR.max() * 100.0
	_hm = _he_fit_metrics(l_plot_HeTR, _he_excess, frame='air')
	with open(_save_prefix + 'tpm_He10830_metrics.txt', 'w') as _fh:
		for _ml in transit_metadata_block(_adv, _tool_identity, _overrides,
		                                  _census_summary,
		                                  _transit_census.get('He10830'),
		                                  state=_state_record):
			_fh.write('# %s\n' % _ml)
		_fh.write('# He 10830 line metrics (three-Gaussian fit, air frame,\n'
		          '# planet-rotation and instrument-convolved curve; he_line_metrics.py)\n')
		for _k in ('red_depth', 'blue_depth', 'red_blue', 'fwhm_A',
		           'shift_A', 'sigma_A'):
			_fh.write('%-12s %14.6e\n' % (_k, _hm[_k]))
	print('(TPM) He 10830 metrics: red %.3f%%  blue %.3f%%  red/blue %.2f  '
	      'FWHM %.3f A  -> %stpm_He10830_metrics.txt'
	      % (_hm['red_depth'], _hm['blue_depth'], _hm['red_blue'],
	         _hm['fwhm_A'], _save_prefix))
except Exception as _err:
	print('(TPM) WARNING: He 10830 metric fit failed: %s' % _err)

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
	for co in sp['comps']:
		plt.plot([co[0], co[0]], [0.0, 1.2], '--', color=gray)
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




