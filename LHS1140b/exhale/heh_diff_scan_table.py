#!/usr/bin/env python3
"""Read a LHS 1140 b composition scan taken with diffusion on at one fixed
eddy coefficient, and print the decision table plus the EW crossing.

Same measurements and same definitions as kzz_scan_table.py -- whose readers
are imported rather than repeated -- with two differences that the scan axis
forces:

  * the reservoir ratio He/H and the eddy coefficient He_Kzz are both read
    from each case's input.inp instead of being assumed, since He/H is now
    the variable and the scan has been repeated at every decade of He_Kzz;
  * the He 10830 profile metrics are reported in full (red, blue, ratio,
    FWHM) from tpm_He10830_metrics.txt, and the composition at which the
    red-pair equivalent width crosses the measurement is solved for by
    Brent's method on a log-log interpolation of the scanned points.

usage: ./heh_diff_scan_table.py [0 | 1e6 | 1e8 | 1e9 | 1e10 | flat | <case> ...]
       (no argument = the adopted 1e9 scan)

The "flat" group is not a composition scan: it holds one composition solved
at every decade from He_Kzz = 0 to 1e6, which is how the plateau below 1e6
is measured.  The crossing solve is skipped for it, there being a single
He/H to interpolate on.
"""
import os
import sys

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq

import kzz_scan_table as K

HERE = os.path.dirname(os.path.abspath(__file__))
RLIST = [1.05, 2.0, 5.0, 10.0, 20.0]
EW_OBS, EW_SIG = 1.108, 0.030

SCANS = {
    '0': ['heh4p5_diff_kzz0',
          'heh4p6_diff_kzz0',
          'heh4p8_diff_kzz0',
          'heh4p9_diff_kzz0',
          'heh5_diff_kzz0'],
    '1e5': ['heh4p5_diff_kzz1e5',
            'heh4p6_diff_kzz1e5',
            'heh4p8_diff_kzz1e5',
            'heh4p9_diff_kzz1e5',
            'heh5_diff_kzz1e5'],
    '1e7': ['heh3p4_diff_kzz1e7',
            'heh3p6_diff_kzz1e7',
            'heh3p8_diff_kzz1e7',
            'heh3p9_diff_kzz1e7',
            'heh4_diff_kzz1e7',
            'heh4p2_diff_kzz1e7'],
    '1e6': ['heh2_diff_kzz1e6',
            'heh4p3_diff_kzz1e6',
            'heh4p4_diff_kzz1e6',
            'heh4p7_diff_kzz1e6',
            'heh5_diff_kzz1e6',
            'heh10_diff_kzz1e6'],
    'flat': ['heh5_diff_kzz0',
             'heh5_diff_kzz1e1',
             'heh5_diff_kzz1e2',
             'heh5_diff_kzz1e3',
             'heh5_diff_kzz1e4',
             'heh5_diff_kzz1e5',
             'heh5_diff_kzz1e6'],
    '1e8': ['heh0p55_diff_kzz1e8',
            'heh2p7_diff_kzz1e8',
            'heh3_diff_kzz1e8',
            'heh3p5_diff_kzz1e8',
            'heh4_diff_kzz1e8',
            'heh10_diff_kzz1e8'],
    '1e9': ['heh0p55_diff_kzz1e9',
            'heh1_diff_kzz1e9',
            'heh2_diff_kzz1e9',
            'heh2p13_diff_kzz1e9',
            'heh4_diff_kzz1e9'],
    '1e10': ['heh0p55_diff_kzz1e10',
             'heh1_diff_kzz1e10',
             'heh1p4_diff_kzz1e10',
             'heh1p5_diff_kzz1e10',
             'heh1p6_diff_kzz1e10'],
    '1e11': ['heh1_diff_kzz1e11',
             'heh1p2_diff_kzz1e11',
             'heh1p4_diff_kzz1e11',
             'heh1p5_diff_kzz1e11'],
}

# heh10_diff_kzz1e8 (He/H = 10) is listed in the table and excluded from the
# crossing fit.  The stored files of that case were taken while the advection
# correction was gated on the *hydrogen* Damkohler number alone, which kept
# 116 of 503 cells at the equilibrium He 2^3S against 27-32 in every other
# case here and left a step of an order of magnitude in n(2^3S); the gate now
# uses the slowest relaxation rate of the solved species vector
# (post_process_adv.f90, Update_EXHALE.md section 72) and re-post-processing
# the case removes the step, but its EW (5.69 %A against 23.45 %A here) is
# still far off the He/H = 2.7-4.0 sequence this fit interpolates.
EXCLUDE_FROM_CROSSING = {'heh10_diff_kzz1e8'}


def key_of(d, key, default=None):
    for line in open(os.path.join(d, 'input.inp')):
        if line.startswith(key + ':'):
            return float(line.split(':')[1])
    if default is not None:
        return default
    raise ValueError('no ' + key + ' in ' + d)


def metrics(d):
    """Every entry of the run's He 10830 metric file."""
    p = os.path.join(d, 'tpm_He10830_metrics.txt')
    out = {}
    if os.path.exists(p):
        for line in open(p):
            w = line.split()
            if len(w) >= 2 and not line.startswith('#'):
                try:
                    out[w[0]] = float(w[1])
                except ValueError:
                    pass
    return out


def main(cases):
    rows = []
    for tag in cases:
        d = os.path.join(HERE, tag)
        if not os.path.isdir(d):
            print('missing: ' + tag)
            continue
        heh = key_of(d, 'He/H number ratio')
        kzz = key_of(d, 'He_Kzz', 0.0)
        K.HEH = heh                       # element_ratio normalizes by it
        s = K.solver(d)
        r, x = K.element_ratio(d)
        rh, _ = K.homopause(d, kzz)
        dep, ew = K.he_line(d)
        rows.append(dict(tag=tag, heh=heh, kzz=kzz, s=s, mdot=K.mdot(d),
                         rh=rh, dep=dep, ew=ew, m=metrics(d),
                         ratio=[float(np.interp(rr, r, x)) for rr in RLIST]))

    kset = sorted({q['kzz'] for q in rows})
    print('\n## solver, homopause, Mdot, He 10830  (He_Kzz = '
          + ', '.join('%.0e' % k for k in kset) + ')\n')
    print('| case | He/H | info | ||R|| | outer passes | last drift '
          '| homopause r [R_p] | log10 Mdot | red depth [%] | red EW [% A] |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    for q in rows:
        hp = ('%.4f' % q['rh']) if isinstance(q['rh'], float) \
            and np.isfinite(q['rh']) else str(q['rh'])
        print('| %s | %g | %s | %s | %d | %s | %s | %.4f | %.4f | %.4f |'
              % (q['tag'], q['heh'], q['s'].get('info', '-'),
                 q['s'].get('resid', '-'), q['s'].get('npass', 0),
                 q['s'].get('drift', '-'), hp, q['mdot'], q['dep'], q['ew']))

    print('\n## (He/H)/HeH, elemental, from the _adv profiles\n')
    print('| case | ' + ' | '.join('%g' % rr for rr in RLIST) + ' |')
    print('|---' * (len(RLIST) + 1) + '|')
    for q in rows:
        print('| %s | ' % q['tag']
              + ' | '.join('%.4f' % v for v in q['ratio']) + ' |')

    print('\n## He 10830 profile metrics\n')
    keys = sorted({k for q in rows for k in q['m']})
    print('| case | ' + ' | '.join(keys) + ' |')
    print('|---' * (len(keys) + 1) + '|')
    for q in rows:
        print('| %s | ' % q['tag']
              + ' | '.join('%.4g' % q['m'][k] if k in q['m'] else '-'
                           for k in keys) + ' |')

    fit = [q for q in rows if q['tag'] not in EXCLUDE_FROM_CROSSING]
    for q in rows:
        if q not in fit:
            print('\nexcluded from the crossing fit: %s (see '
                  'EXCLUDE_FROM_CROSSING)' % q['tag'])
    h = np.array([q['heh'] for q in fit])
    e = np.array([q['ew'] for q in fit])
    if len(np.unique(h)) < 2:
        print('\nsingle composition in this group: no crossing to solve for')
        return
    o = np.argsort(h)
    h, e = h[o], e[o]
    f = interp1d(np.log10(h), np.log10(e), kind='linear')
    print('\n## EW crossing (log-log interpolation of the scanned points)\n')
    for lab, target in (('1 sigma low', EW_OBS - EW_SIG),
                        ('central', EW_OBS),
                        ('1 sigma high', EW_OBS + EW_SIG)):
        if not (e.min() <= target <= e.max()):
            print('%-13s %.3f %%A: outside the bracket EW = %.4f .. %.4f'
                  % (lab, target, e.min(), e.max()))
            continue
        j = int(np.searchsorted(e, target))
        root = brentq(lambda t: f(t) - np.log10(target),
                      np.log10(h[0]), np.log10(h[-1]))
        print('%-13s %.3f %%A -> He/H = %.4f   (bracketed by He/H = %g '
              '[EW %.4f] and %g [EW %.4f])'
              % (lab, target, 10.0**root, h[j - 1], e[j - 1], h[j], e[j]))


if __name__ == '__main__':
    argv = sys.argv[1:] or ['1e9']
    cases = []
    for a in argv:
        cases += SCANS[a] if a in SCANS else [a]
    main(cases)
