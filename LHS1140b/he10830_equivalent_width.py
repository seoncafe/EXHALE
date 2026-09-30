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


if __name__ == '__main__':
    import sys
    ew, err = observed_equivalent_width()
    print('observation: EW = %.10f +/- %.10f %%A over %d samples, '
          '%.6f to %.6f A (vacuum)' % (ew, err, SAMPLES.size, SAMPLES[0],
                                        SAMPLES[-1]))
    for p in sys.argv[1:]:
        v = transit_file_equivalent_width(p)
        print('%s: %s' % (p, 'None' if v is None else '%.10f %%A' % v))
