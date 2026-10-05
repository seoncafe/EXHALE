"""The He I 10830 red-pair equivalent width of the LHS 1140 b comparison.

ONE MEASUREMENT OPERATOR FOR THE OBSERVATION AND FOR EVERY MODEL (decision
of 2026-09-23, md/PLAN_20260923_rev2.md section 4).  The statistic is the
trapezoidal integral of the excess absorption [%] over the observation's own
vacuum wavelength samples inside [EW_LO, EW_HI].  The observation is
integrated on those samples as it stands.  A model spectrum, already
convolved to the instrument resolution by EXHALE_transit.py, is interpolated
linearly onto the SAME samples and integrated with the SAME weights, so that
the two numbers differ only by the spectra and not by where the integral
starts and stops or how it is weighted.  Before this module each caller
integrated a model on its own finer grid inside the nominal window; on the
certified molecular K_zz = 1e9, He/H = 2.13 spectrum that gave 1.5603193380
against 1.5603912281 on the observation's samples (MEASURED 2026-09-23).

What the operator does NOT include: the detector-bin response (the CSV
gives point values, so a model is sampled, not bin-integrated), any
line-of-sight velocity (models carry none and the observation is in its
planet rest frame; every EW here is computed unshifted), and any covariance
of the observational error (the error below assumes independent samples).

The window and the vacuum/air ratio are the ones make_memo_figures.py has
used since the catalog began.  The observational error is the error of the
SAME linear estimator as the equivalent width: EW = sum_i w_i d_i with the
trapezoidal weights w_i of the samples, so sigma^2 = sum_i (w_i sigma_i)^2
for independent samples (until 2026-09-29 every sample was given the median
sample spacing as its weight, 0.0295 against 0.0289 %A; code audit of
2026-09-29, md/CODE_AUDIT_20260929.md F5).  Neither includes a covariance of
the observational errors or the uncertainty of the continuum fit, which the
released spectrum does not provide.
"""
import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OBS_CSV = os.path.join(HERE, 'Cherubim_2026',
                       'LHS1140b_He10833_Fig3B_spectrum.csv')

# Vacuum over air wavelength of the red pair (10832.057 A vacuum,
# 10829.09114 A air); EXHALE_transit.py writes air wavelengths.
AIR = 10832.057/10829.09114
# The integration window, vacuum wavelengths in the planet rest frame [A].
EW_LO, EW_HI = 10832.60, 10834.20


def _read_observation():
    rows = list(csv.DictReader(open(OBS_CSV)))

    def col(name):
        return np.array([float(r[name]) for r in rows])
    return (col('wavelength_vacuum_planet_rest_A'),
            col('absorption_depth_positive_percent'),
            col('uncertainty_1sigma_percent'))


OBS_VAC, OBS_DEPTH, OBS_SIGMA = _read_observation()
IN_WINDOW = (OBS_VAC >= EW_LO) & (OBS_VAC <= EW_HI)
# The abscissae and, through np.trapz, the weights of the operator.
SAMPLES = OBS_VAC[IN_WINDOW]


def trapezoidal_weights(x):
    """The weights w_i with np.trapz(y, x) = sum_i w_i y_i: half the sum of
    the two adjacent intervals inside, half the one interval at each end."""
    dx = np.diff(x)
    w = np.zeros_like(x)
    w[:-1] += 0.5*dx
    w[1:] += 0.5*dx
    return w


def observed_equivalent_width():
    """(EW, sigma) of the observation [%A].

    sigma = sqrt(sum (w_i sigma_i)^2) with the trapezoidal weights w_i of the
    estimator itself: independent samples, no covariance
    (md/PLAN_20260923_rev2.md section 4 states what this leaves out)."""
    ew = float(np.trapz(OBS_DEPTH[IN_WINDOW], SAMPLES))
    err = float(np.sqrt(np.sum((trapezoidal_weights(SAMPLES)
                                * OBS_SIGMA[IN_WINDOW])**2)))
    return ew, err


def model_equivalent_width(lam_vac, excess_percent):
    """EW [%A] of a model excess-absorption curve through the operator.

    lam_vac: vacuum wavelengths [A], planet rest frame; excess_percent: the
    excess absorption [%] at those wavelengths.  Returns None when the curve
    does not span the observation's samples in the window, since linear
    interpolation would then extrapolate a constant, and when the curve is
    not a curve: arrays that are not one-dimensional or not of one length,
    fewer than two samples, a non-finite value in either array (a NaN
    wavelength passes the ordering test below, since every comparison with
    it is false), or wavelengths not strictly ordered."""
    lam = np.asarray(lam_vac, dtype=float)
    exc = np.asarray(excess_percent, dtype=float)
    if lam.ndim != 1 or exc.ndim != 1 or lam.size != exc.size:
        return None
    if lam.size < 2:
        return None
    if not (np.all(np.isfinite(lam)) and np.all(np.isfinite(exc))):
        return None
    if lam[0] > lam[-1]:
        lam, exc = lam[::-1], exc[::-1]
    if np.any(np.diff(lam) <= 0.0):
        return None
    if lam[0] > SAMPLES[0] or lam[-1] < SAMPLES[-1]:
        return None
    return float(np.trapz(np.interp(SAMPLES, lam, exc), SAMPLES))


def transit_file_excess(path):
    """(vacuum wavelength [A], excess absorption [%]) from an EXHALE
    `tpm_He10830.txt`: column 1 is the air wavelength and column 3 the
    instrument-convolved transit depth, whose maximum over the file is the
    continuum."""
    s = np.loadtxt(path)
    if s.ndim != 2 or s.shape[1] < 3:
        return None
    if not np.all(np.isfinite(s[:, [0, 2]])):
        return None
    cont = s[:, 2].max()
    if not cont > 0.0:
        return None
    lam = s[:, 0]*AIR
    exc = (cont - s[:, 2])/cont*100.0
    return lam, exc


def transit_file_equivalent_width(path):
    """EW [%A] of an EXHALE `tpm_He10830.txt` through the operator, or None."""
    if not os.path.isfile(path):
        return None
    le = transit_file_excess(path)
    if le is None:
        return None
    return model_equivalent_width(*le)


# ----- THE ADDED VELOCITY BROADENING (width-matching kernel) -----
# What it stands for: a Gaussian line-of-sight velocity distribution that the
# 1-D steady model does not contain (the measured line is about three times
# wider than any solved profile).  It is NOT a computed physics term: it
# carries no energy or momentum budget and moves no absorber.  Its width is
# not a parameter of the model; it is MEASURED from the data, for each model,
# as the value at which the broadened model reaches the measured line width
# (docs/lhs1140b_exhale_vs_pwinds.tex, Sect. "The broadening the data demand,
# measured", sec:broadening; the definitions below were moved here unchanged
# from make_memo_figures.py, which imports them, so the memo and the EW chain
# apply one rule).
#
# The convention, step by step:
#   1. excess = (max T - T)/max T [%] from column 3 of tpm_He10830.txt
#      (T_instr), the synthetic spectrum convolved with the instrument
#      profile (R = 68,000) only, WITHOUT planet rotation (column 4,
#      T_rot+instr, carries it; for LHS 1140 b, v_rot ~ 0.03 km/s, the two
#      give EWs equal to 1e-6 %A, measured 2026-10-05 at He/H 0.1, s = 2.53), on its own uniform AIR wavelength grid (transit_file_excess
#      without the vacuum conversion);
#   2. that excess is convolved with a Gaussian of FWHM f in VELOCITY,
#      sigma_lambda = f/(2.35482 c) x 10830 A, by scipy's gaussian_filter1d
#      (reflecting boundary, kernel truncated at 4 sigma), planet rest frame,
#      no shift; the instrument kernel is already inside the curve and
#      Gaussians commute, so this is the same as broadening before it;
#   3. depth, width and red/blue are read with the three-Gaussian extractor
#      used on the data (he_line_metrics.fit_metrics, air frame); the width is
#      the FWHM of the blended red pair;
#   4. f is chosen so that this FWHM equals the measured 0.841 A
#      (FWHM_OBS_A, the same extractor on the released spectrum), by brentq
#      on f in [1, 60] km/s with xtol 1e-3;
#   5. the equivalent width is model_equivalent_width of the broadened curve,
#      i.e. interpolated onto the observation's vacuum samples in
#      [EW_LO, EW_HI] and integrated with their trapezoid weights.
# The reflecting boundary conserves the absorption summed over the grid; the
# window EW falls with f because the wings carry absorption out of the
# window, not because absorbing power is lost.
# Reference only: the width audit (md/lhs1140b_width_measurement_audit.md)
# put the like-for-like requirement at sigma = 9.33 km/s (FWHM 21.96 km/s);
# the memo found f = 22.22 km/s (sigma 9.44) on the 2026-09 solutions.
FWHM_OBS_A = 0.841                       # A, measured red-pair width
C_KMS = 2.99792458e5
FWHM_PER_SIGMA = 2.35482
NONTHERMAL_SIGMA_REFERENCE_KMS = 9.33    # width audit, for comparison only


def broaden_excess(lam, exc, fwhm_kms, lam0=10830.0):
    """Extra Gaussian velocity broadening, FWHM in km/s.

    lam0 converts the velocity width to a wavelength width; use the line
    position in the frame of `lam` (10830 A in air, 10833 A in vacuum).
    """
    from scipy.ndimage import gaussian_filter1d
    if fwhm_kms <= 0.0:
        return exc
    dl = np.median(np.diff(lam))
    return gaussian_filter1d(exc, fwhm_kms/C_KMS*lam0/FWHM_PER_SIGMA/dl)


def _fit_metrics():
    root = os.path.dirname(HERE)
    if root not in __import__('sys').path:
        __import__('sys').path.insert(0, root)
    from he_line_metrics import fit_metrics
    return fit_metrics


def broadened_line_metrics(lam, exc, fwhm_kms, frame='air', lam0=10830.0):
    """Three-Gaussian metrics and window EW of the broadened excess."""
    e = broaden_excess(lam, exc, fwhm_kms, lam0)
    d = _fit_metrics()(lam, e, frame=frame)
    lv = lam*AIR if frame == 'air' else lam
    v = model_equivalent_width(lv, e)
    d['ew'] = np.nan if v is None else v
    return d


def width_matching_fwhm_kms(lam, exc, frame='air', lam0=10830.0):
    """Added kernel FWHM [km/s] at which the line reaches the measured width."""
    from scipy.optimize import brentq
    return brentq(lambda f: broadened_line_metrics(lam, exc, f, frame,
                                                   lam0)['fwhm_A']
                  - FWHM_OBS_A, 1.0, 60.0, xtol=1e-3)


def transit_file_air_excess(path):
    """(air wavelength [A], excess [%]) of column 3 of a tpm_He10830.txt."""
    s = np.loadtxt(path)
    return s[:, 0], (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100.0


def width_matched_product(path, out_path, fwhm_kms=None):
    """Broaden a tpm_He10830.txt to the measured width (or by a given FWHM)
    and write the broadened excess with its record.  Returns the record."""
    lam, exc = transit_file_air_excess(path)
    criterion = 'given'
    if fwhm_kms is None:
        fwhm_kms = width_matching_fwhm_kms(lam, exc)
        criterion = ('width_match: three-Gaussian blended red-pair FWHM of '
                     'the broadened curve = %.3f A (measured)' % FWHM_OBS_A)
    m0 = broadened_line_metrics(lam, exc, 0.0)
    m = broadened_line_metrics(lam, exc, fwhm_kms)
    rec = dict(fwhm_kms=fwhm_kms, sigma_kms=fwhm_kms/FWHM_PER_SIGMA,
               criterion=criterion, ew0=m0['ew'], ew=m['ew'],
               red=m['red_depth'], blue=m['blue_depth'],
               red_blue=m['red_blue'], fwhm_A=m['fwhm_A'],
               red0=m0['red_depth'], fwhm0_A=m0['fwhm_A'])
    head = ['he10830_broadened: the excess of column 3 of %s convolved with '
            'an added Gaussian line-of-sight velocity distribution' % path,
            'standing for the velocity field of the outflow that the 1-D '
            'model does not contain (not a computed physics term)',
            'kernel FWHM_kms %.6f sigma_kms %.6f (reference: width audit '
            'sigma %.2f km/s)' % (rec['fwhm_kms'], rec['sigma_kms'],
                                   NONTHERMAL_SIGMA_REFERENCE_KMS),
            'criterion %s' % criterion,
            'convention he10830_equivalent_width.py (docs/'
            'lhs1140b_exhale_vs_pwinds.tex sec:broadening)',
            'unbroadened: EW %.6f red %.6f FWHM_A %.6f' % (rec['ew0'],
                                                         rec['red0'],
                                                         rec['fwhm0_A']),
            'broadened:   EW %.6f red %.6f blue %.6f red/blue %.4f FWHM_A '
            '%.6f' % (rec['ew'], rec['red'], rec['blue'], rec['red_blue'],
                      rec['fwhm_A']),
            'EW in percent-angstrom through the window operator; depths in %',
            'lambda_air[A]  excess[%]  excess_broadened[%]']
    np.savetxt(out_path, np.c_[lam, exc, broaden_excess(lam, exc, fwhm_kms)],
               header='\n'.join(head))
    return rec


if __name__ == '__main__':
    import sys
    if len(sys.argv) >= 4 and sys.argv[1] == '--width-matched':
        # --width-matched <tpm_He10830.txt> <output> [FWHM km/s]
        r = width_matched_product(sys.argv[2], sys.argv[3],
                                  float(sys.argv[4]) if len(sys.argv) > 4
                                  else None)
        print('fwhm_kms=%.4f sigma_kms=%.4f EW_kernel=%.4f EW_nokernel=%.4f '
              'red=%.4f blue=%.4f red_blue=%.3f fwhm_A=%.4f'
              % (r['fwhm_kms'], r['sigma_kms'], r['ew'], r['ew0'], r['red'],
                 r['blue'], r['red_blue'], r['fwhm_A']))
        sys.exit(0)
    ew, err = observed_equivalent_width()
    print('observation: EW = %.10f +/- %.10f %%A over %d samples, '
          '%.6f to %.6f A (vacuum)' % (ew, err, SAMPLES.size, SAMPLES[0],
                                        SAMPLES[-1]))
    for p in sys.argv[1:]:
        v = transit_file_equivalent_width(p)
        print('%s: %s' % (p, 'None' if v is None else '%.10f %%A' % v))
