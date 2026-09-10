#!/usr/bin/env python3
"""LHS 1140b: EXHALE against p-winds and against the 2024 measurement.

Same planet, same stellar spectrum, same composition grid -- the comparison
the tutorial-planet series could not make.  What still differs between the
two models is the physics, and that is the point:

  p-winds  isothermal Parker wind at an imposed (Mdot, T) taken from the
           paper's retrieval; ionization and the metastable population are
           solved on that fixed structure.
  EXHALE   the wind itself is solved -- mass, momentum and energy with
           photoionization heating and radiative cooling -- so Mdot and T(r)
           are outputs, not inputs.

Run from LHS1140b/ after the exhale/ scan has converged.
"""
import os, re, glob, sys
import numpy as np

# EXHALE's own constants [cm].  The Jupiter radius has one Python definition,
# RJ_CM of examples/exhale_io.py (the IAU 2015 nominal equatorial radius of
# parameters.f90); the directory is APPENDED so nothing in it can shadow a
# standard-library module.
_EXAMPLES_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'examples')
if _EXAMPLES_DIR not in sys.path:
    sys.path.append(_EXAMPLES_DIR)
from exhale_io import RJ_CM                                 # noqa: E402

RJ, RSUN = RJ_CM, 6.957e10
CASES = [('solar', 0.0833), ('heh1', 1.0), ('heh10', 10.0),
         ('heh100', 100.0), ('heh1000', 1000.0), ('heh10000', 1e4)]
OBS_CSV = 'Cherubim_2026/LHS1140b_He10833_Fig3B_spectrum.csv'

def observed():
    """The 2024 measurement, put through the same extractor as the models.

    The spectrum is the authors' own released file (Zenodo 10.5281/zenodo.
    20723095) reconstructed to their Fig. 3B, not a reading of the figure.
    Fitting it with he_line_metrics gives 1.254 % against the 1.24
    (+0.22/-0.23) % the paper quotes from its own MCMC, so one extractor can
    carry both the models and the measurement."""
    import csv, sys
    sys.path.insert(0, '..')
    from he_line_metrics import fit_metrics
    rows = list(csv.DictReader(open(OBS_CSV)))
    lam = np.array([float(r['wavelength_vacuum_planet_rest_A']) for r in rows])
    dep = np.array([float(r['absorption_depth_positive_percent']) for r in rows])
    m = fit_metrics(lam, dep, frame='vacuum')
    m['published_red'] = 1.24
    m['published_red_hi'], m['published_red_lo'] = 0.22, 0.23
    m['Reff_over_Rp'] = 1.52          # the paper's equivalent opaque radius
    return m

def rp_rstar():
    t = open('exhale/heh1000/input.inp').read()
    rp = float(re.search(r'^Planet radius.*$', t, re.M).group(0).split()[3])
    rs = float(re.search(r'^Stellar radius.*$', t, re.M).group(0).split()[3])
    return rp*RJ, rs*RSUN

def exhale_case(tag):
    d = f'exhale/{tag}'
    out = dict(tag=tag)
    lg = f'{d}/run.log'
    if os.path.isfile(lg):
        for ln in open(lg, errors='replace'):
            if 'steady-state Mdot' in ln:
                out['logMdot'] = float(ln.split('=')[-1].split()[0])
            if 'JFNK) done' in ln:
                out['jfnk'] = ln.strip()
    f = f'{d}/output/Hydro_ioniz.txt'
    if os.path.isfile(f):
        h = np.loadtxt(f)
        out.update(r=h[:, 0], T=h[:, 4], v=h[:, 2], n=h[:, 1])
        out['Tmax'] = h[:, 4].max()
        out['rTmax'] = h[np.argmax(h[:, 4]), 0]
    fi = f'{d}/output/Ion_species.txt'
    if os.path.isfile(fi):
        lab = [l for l in open(fi) if l.startswith('# columns')][0].split()[2:]
        k = {n: j for j, n in enumerate(lab)}
        i = np.loadtxt(fi)
        n3 = i[:, k['HeITR']]
        out['n3max'] = n3.max()
        out['r_n3max'] = i[np.argmax(n3), 0]
        H = i[:, k['HI']] + i[:, k['HII']]
        out['x_HII'] = np.interp(2.0, i[:, 0], i[:, k['HII']]/np.maximum(H, 1e-300))
    m = f'{d}/tpm_He10830_metrics.txt'
    if os.path.isfile(m):
        out['metrics'] = {ln.split()[0]: float(ln.split()[1])
                          for ln in open(m) if not ln.startswith('#')}
    return out

def main():
    Rp, Rs = rp_rstar()
    base = (Rp/Rs)**2*100
    print(f'LHS 1140b: R_p/R_star = {Rp/Rs:.5f},  bare transit depth = {base:.4f} %\n')

    pw = np.loadtxt('pwinds_oracle/scan_hhe.txt')
    pw_heh = {round(1.0/h): row for h, row in zip(pw[:, 0], pw)}

    hdr = (f"{'He/H':>8} {'logMdot':>8} {'T_max[K]':>9} {'n(2^3S)max':>11} "
           f"{'red[%]':>8} {'blue[%]':>8} {'r/b':>6} {'R_eff/Rp':>9} | "
           f"{'p-winds red[%]':>14}")
    print(hdr); print('-'*len(hdr))
    for tag, heh in CASES:
        c = exhale_case(tag)
        mt = c.get('metrics', {})
        red = mt.get('red_depth', np.nan)
        reff = np.sqrt(1.0 + red/100*(Rs/Rp)**2) if red == red else np.nan
        p = pw_heh.get(round(heh))
        pred = p[4] if p is not None else np.nan
        print(f"{heh:8.4g} {c.get('logMdot', float('nan')):8.2f} "
              f"{c.get('Tmax', float('nan')):9.0f} {c.get('n3max', float('nan')):11.3e} "
              f"{red:8.3f} {mt.get('blue_depth', float('nan')):8.3f} "
              f"{mt.get('red_blue', float('nan')):6.2f} {reff:9.3f} | {pred:14.3f}")
    print('-'*len(hdr))
    o = observed()
    print(f"{'observed':>8} {'':8} {'':9} {'':11} {o['red_depth']:8.3f} "
          f"{o['blue_depth']:8.3f} {o['red_blue']:6.2f} "
          f"{np.sqrt(1.0 + o['red_depth']/100*(Rs/Rp)**2):9.3f} |")
    print(f"\npaper's own fit to the same spectrum: {o['published_red']:.2f} "
          f"+{o['published_red_hi']:.2f}/-{o['published_red_lo']:.2f} %.")
    print('Depths are the safe comparison. The R_eff column is the excess-area')
    print('convention, delta = (R_eff^2 - R_p^2)/R_star^2, which is what the')
    print('plotted spectrum measures: its continuum sits at zero, so the')
    print("white-light transit is already divided out. The paper's stated")
    print(f"equivalent opaque radius of {o['Reff_over_Rp']:.2f} R_p is instead recovered by")
    print('delta = (R_eff/R_star)^2, and the two differ by 20 % in radius.')
    print(f"\np-winds holds Mdot = 2.03e8 g/s (log 8.31) and T = 6100 K fixed;"
          f"\nEXHALE solves for both, so a column-by-column difference in depth"
          f"\nis a difference in the wind, not only in the line formation.")

if __name__ == '__main__':
    main()
