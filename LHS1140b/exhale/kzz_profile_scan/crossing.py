#!/usr/bin/env python3
"""The reservoir He/H at which the red-pair equivalent width crosses the
measured 1.108 +/- 0.030 %A, one row per K_zz.

Two independent solutions of the same crossing are printed so they can be
checked against each other:

  bracket   a straight log EW vs log(He/H) line through the two arms that
            straddle the target, solved by Brent's method;
  allpts    a least-squares line through log EW vs log(He/H) over every arm
            at that K_zz, solved the same way.

They agree only if the response is locally a power law over the arms used;
where they do not, the bracket is the number to quote and the disagreement is
the honest error bar on the interpolation.

Both the reservoir He/H (what goes into the chemistry) and the match-surface
He/H (what the profile actually hands the wind) are crossed, because at low
K_zz they need not be the same number.
"""
import glob
import os
import sys

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import measure as M

EW_OBS, EW_ERR = M.EW_OBS, M.EW_ERR


def arms(pattern):
    out = []
    for c in sorted(glob.glob(pattern)):
        if not os.path.isfile(os.path.join(c, 'closure_history.txt')):
            continue
        try:
            out.append(M.row(c))
        except Exception as exc:
            print('  (skipped %s: %s)' % (c, exc))
    return out


def _unique_sorted(x, y):
    """One point per reservoir value, ordered by equivalent width."""
    seen = {}
    for xi, yi in zip(np.asarray(x, float), np.asarray(y, float)):
        seen[round(float(xi), 6)] = float(yi)
    xs = np.array(sorted(seen))
    ys = np.array([seen[k] for k in xs])
    o = np.argsort(ys)
    return xs[o], ys[o]


def cross_bracket(x, y, target):
    """Log-log line through the two arms that straddle `target`."""
    x, y = _unique_sorted(x, y)
    if len(x) < 2 or target < y[0] or target > y[-1]:
        return float('nan'), None
    j = int(np.searchsorted(y, target))
    j = min(max(j, 1), len(y) - 1)
    lo, hi = j - 1, j
    lx, ly = np.log(x[[lo, hi]]), np.log(y[[lo, hi]])
    f = lambda t: (ly[0] + (ly[1] - ly[0])*(t - lx[0])/(lx[1] - lx[0])
                   - np.log(target))
    return float(np.exp(brentq(f, lx[0], lx[1]))), (x[lo], y[lo], x[hi], y[hi])


def cross_allpts(x, y, target, span=1.5):
    """Least-squares log-log line through the arms near the target.

    The response is a power law only locally -- over the whole scanned
    range of the reservoir ratio the local slope runs from 0.6 to 0.25 --
    so the arms more than a factor `span` away in equivalent width are left
    out; including them would make this a worse estimate, not a wider one.
    """
    x, y = _unique_sorted(x, y)
    m = np.abs(np.log(y/target)) < np.log(span)
    if m.sum() >= 2:
        x, y = x[m], y[m]
    if len(x) < 2:
        return float('nan')
    lx, ly = np.log(x), np.log(y)
    a, b = np.polyfit(lx, ly, 1)
    lo, hi = lx.min() - 2.0, lx.max() + 2.0
    f = lambda t: a*t + b - np.log(target)
    if f(lo)*f(hi) > 0:
        return float('nan')
    return float(np.exp(brentq(f, lo, hi)))


def report(tag, x, y, label):
    # ordered by equivalent width, one point per reservoir value, so that
    # y[0] and y[-1] are the ends of the ladder and not the ends of a
    # directory listing
    x, y = _unique_sorted(x, y)
    c, br = cross_bracket(x, y, EW_OBS)
    lo, _ = cross_bracket(x, y, EW_OBS - EW_ERR)
    hi, _ = cross_bracket(x, y, EW_OBS + EW_ERR)
    # The second solution is a check on the first, so it is only shown where
    # there is a first: outside the measured range of the equivalent width
    # the least-squares line is an extrapolation and is not reported.
    a = cross_allpts(x, y, EW_OBS) if br is not None else float('nan')
    s = ('%-9s %-9s  %8.3f  %8.3f - %-8.3f  %8.3f'
         % (tag, label, c, lo, hi, a))
    if br is not None:
        s += '   [%.3f EW %.4f | %.3f EW %.4f]' % br
    else:
        s += ('   not bracketed: the %d arms on this branch span EW %.4f'
              ' to %.4f' % (len(x), y.min(), y.max()))
        # When the nearest arm is within 10 % of the target, the local slope
        # of the two arms nearest it carries the crossing a short way past
        # the end of the ladder.  It is an extrapolation and is labelled one.
        near = x[0] if abs(np.log(y[0]/EW_OBS)) < abs(np.log(y[-1]/EW_OBS)) \
            else x[-1]
        ynear = y[0] if near == x[0] else y[-1]
        if len(x) >= 2 and abs(np.log(ynear/EW_OBS)) < np.log(1.10):
            i = 0 if near == x[0] else len(x) - 1
            j = 1 if i == 0 else len(x) - 2
            sl = np.log(y[i]/y[j])/np.log(x[i]/x[j])
            s += ('\n%-9s %-9s  extrapolated on the local slope %.3f from'
                  ' He/H = %.3f (EW %.4f): %.3f'
                  % ('', '', sl, near, ynear,
                     near*np.exp(np.log(EW_OBS/ynear)/sl)))
    print(s)
    return c


# The wind has two converged states over most of the scanned range: a cool
# one at log10 Mdot ~ 7.4-7.6 and a hot one at ~7.83-7.90, with the He 10830
# line four times stronger on the hot one.  A crossing may only be
# interpolated within one of them, so the arms are separated here.  The
# split is on log10 Mdot, which is bimodal with nothing between 7.60 and
# 7.82 in any arm measured.
BRANCH_MDOT = 7.70


def branch_of(d):
    return 'hot' if d['lgMdot'] >= BRANCH_MDOT else 'cool'


def main(patterns):
    want = None
    patterns = list(patterns)
    for k in ('cool', 'hot'):
        if '--' + k in patterns:
            patterns.remove('--' + k)
            want = k
    if not patterns:
        patterns = [os.path.join(HERE, 'kzz*')]
    print('crossing of the red-pair EW with %.3f +/- %.3f %%A%s'
          % (EW_OBS, EW_ERR,
             '' if not want else '   (%s branch only)' % want))
    print('%-9s %-9s  %8s  %-19s  %8s   %s'
          % ('K_zz', 'scale', 'crossing', '1 sigma range', 'allpts',
             'bracketing arms'))
    groups = {}
    for p in patterns:
        found = sorted(glob.glob(os.path.join(p, 'heh*')))
        if os.path.isfile(os.path.join(p, 'closure_history.txt')):
            found = [p]          # a single case directory named outright
        for c in found:
            if not os.path.isfile(os.path.join(c, 'closure_history.txt')):
                continue
            try:
                d = M.row(c)
            except Exception as exc:
                print('  (skipped %s: %s)' % (c, exc))
                continue
            if want and branch_of(d) != want:
                continue
            if not np.isfinite(d['ew']):
                continue          # no transit spectrum on this arm
            groups.setdefault(d['kzz'], []).append(d)
    for kz in sorted(groups, reverse=True):
        g = groups[kz]
        if len(g) < 2:
            print('%-9.1E %-9s  only %d arm(s); no bracket'
                  % (kz, '-', len(g)))
            continue
        report('%.2E' % kz, [d['resv'] for d in g], [d['ew'] for d in g],
               'reservoir')
        report('', [d['heh_match'] for d in g], [d['ew'] for d in g],
               'at match')


if __name__ == '__main__':
    main(sys.argv[1:])
