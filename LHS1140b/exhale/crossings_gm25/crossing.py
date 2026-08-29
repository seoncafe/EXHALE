#!/usr/bin/env python3
"""Solve one measured bracket for the He/H at which the red-pair equivalent
width equals the measured 1.108 +/- 0.030 %A.

Two independent solves are reported so that the answer does not depend on
one of them:

  chord      log-log straight line through the two arms that straddle the
             target -- an interpolation inside the bracket, nothing else;
  quadratic  log-log parabola through the three arms nearest the target,
             which uses the local curvature the chord ignores.

The 1-sigma interval is the same solve run against 1.108 -/+ 0.030, and is
reported only when both ends fall inside the measured arms.  Nothing here
extrapolates: a target outside the scanned range is reported as such.

usage: ./crossing.py <label> <He/H>=<EW> [<He/H>=<EW> ...]
"""
import sys
import numpy as np

EW_OBS, EW_ERR = 1.108, 0.030


def bracket(x, y, t):
    """Indices of the adjacent pair whose EW straddles the target."""
    for i in range(len(x)-1):
        if (y[i]-t)*(y[i+1]-t) <= 0.0:
            return i
    return None


def chord(x, y, t):
    i = bracket(x, y, t)
    if i is None:
        return float('nan'), None
    lx, ly, lt = np.log(x), np.log(y), np.log(t)
    f = (lt-ly[i])/(ly[i+1]-ly[i])
    return float(np.exp(lx[i] + f*(lx[i+1]-lx[i]))), (x[i], x[i+1])


def quadratic(x, y, t):
    """Log-log parabola through the three arms nearest the target, solved
    inside the bracket by bisection (the parabola is monotone there).
    Three arms are the minimum; with two, only the chord is defined."""
    i = bracket(x, y, t)
    if i is None or len(x) < 3:
        return float('nan')
    # the three points centred on the bracket
    j = min(max(i, 1), len(x)-2)
    idx = [j-1, j, j+1]
    lx, ly = np.log(np.asarray(x)[idx]), np.log(np.asarray(y)[idx])
    c = np.polyfit(lx, ly, 2)
    lo, hi = np.log(x[i]), np.log(x[i+1])
    lt = np.log(t)
    flo = np.polyval(c, lo)-lt
    for _ in range(200):
        mid = 0.5*(lo+hi)
        fm = np.polyval(c, mid)-lt
        if (flo < 0) == (fm < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return float(np.exp(0.5*(lo+hi)))


def solve(label, x, y):
    o = np.argsort(x)
    x, y = np.asarray(x)[o], np.asarray(y)[o]
    print('%s: %d arms, He/H %.4g-%.4g, EW %.4f-%.4f %%A'
          % (label, len(x), x[0], x[-1], y.min(), y.max()))
    mono = np.all(np.diff(y) > 0)
    if not mono:
        bad = [(x[i], y[i], x[i+1], y[i+1])
               for i in range(len(x)-1) if y[i+1] <= y[i]]
        print('   NOT monotonic in EW: '
              + '; '.join('%.4g->%.4g falls %.4f->%.4f' % b for b in bad))
    out = {}
    for tag, t in (('crossing', EW_OBS), ('-1sig', EW_OBS+EW_ERR),
                   ('+1sig', EW_OBS-EW_ERR)):
        c, br = chord(x, y, t)
        q = quadratic(x, y, t)
        out[tag] = (c, q)
        if br is None:
            print('   %-9s EW=%.3f  OUTSIDE the measured arms '
                  '(no bracket) -- not extrapolated' % (tag, t))
        elif q != q:
            print('   %-9s EW=%.3f  bracket [%.4g, %.4g]  chord %.4f  '
                  '(fewer than three arms: no quadratic)'
                  % (tag, t, br[0], br[1], c))
        else:
            print('   %-9s EW=%.3f  bracket [%.4g, %.4g]  chord %.4f  '
                  'quadratic %.4f  (differ %.2f%%)'
                  % (tag, t, br[0], br[1], c, q, 100*abs(q-c)/c))
    return out


if __name__ == '__main__':
    lab = sys.argv[1]
    xs, ys = [], []
    for a in sys.argv[2:]:
        k, v = a.split('=')
        xs.append(float(k)); ys.append(float(v))
    solve(lab, xs, ys)
