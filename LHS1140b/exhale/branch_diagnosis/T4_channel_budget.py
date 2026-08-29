#!/usr/bin/env python3
"""T4: which heating and cooling channels separate the two states?

Reads `Heating_breakdown.txt` and `Cooling_breakdown.txt` of the cool and the
hot solution and reports, at a set of radii and integrated over 4 pi r^2 dr in
bands, the share each channel holds of its own total, on both states.

The hypothesis under test is a positive feedback in which a cooling channel
carried by a NEUTRAL species collapses as the gas ionizes, so the hot state
loses the cooling that would have pulled it back.  The C I metal line term is
the candidate: on the photochemical base it carries most of the cooling near
1.5 R_p.  The measurement is reported whether or not it supports that.

Usage:  python3 T4_channel_budget.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCAN = os.path.join(HERE, '../kzz_profile_scan/kzz1e7')
COOL = os.path.join(SCAN, 'heh3p4/k02/output')
HOT = os.path.join(SCAN, 'heh3p399/k03/output')

METALS = ['CI', 'CII', 'CIII', 'OI', 'OII', 'OIII', 'NI', 'NII', 'NIII',
          'MgI', 'MgII', 'MgIII', 'SiI', 'SiII', 'SiIII', 'CaI', 'CaII',
          'CaIII', 'NaI', 'NaII', 'KI', 'KII', 'SI', 'SII',
          'FeI', 'FeII', 'FeIII']
COOL_CH = (['reco', 'coio', 'coex_HI[Lya]', 'coex_HeI', 'coex_HeII', 'brem',
            'H3p_IR'] + METALS)
HEAT_CH = ['HI', 'HeI', 'HeII', 'He23S', 'H2', 'metals', 'Hpe[excitedH]',
           'Hdx[Lya-deexc]', 'He_recomb', 'He23S_Penning',
           'He23S_H2_Penning', 'H2_LW']


def read_breakdown(path, channels):
    a = np.loadtxt(path)
    d = {'r': a[:, 0], 'T': a[:, 1], 'ne': a[:, 2], 'total': a[:, 3]}
    for k, nm in enumerate(channels):
        d[nm] = a[:, 4 + k]
    return d


def band_share(d, channels, lo, hi):
    """Share of the 4 pi r^2 dr integral of the total that each channel holds
    inside [lo, hi) R_p."""
    r = d['r']
    m = (r >= lo) & (r < hi)
    w = 4.0 * np.pi * r[m] ** 2 * np.gradient(r)[m]
    tot = np.sum(d['total'][m] * w)
    return tot, {c: np.sum(d[c][m] * w) / max(tot, 1e-300) for c in channels}


def main():
    out = []

    def p(s=''):
        print(s)
        out.append(s)

    for kind, fname, channels in (('COOLING', 'Cooling_breakdown.txt', COOL_CH),
                                  ('HEATING', 'Heating_breakdown.txt', HEAT_CH)):
        c = read_breakdown(os.path.join(COOL, fname), channels)
        h = read_breakdown(os.path.join(HOT, fname), channels)
        p('=' * 78)
        p('%s -- share of the volume-integrated total in each band [%%]' % kind)
        p('=' * 78)
        for lo, hi in ((1.0, 1.2), (1.2, 2.0), (1.0, 3.0), (3.0, 10.0),
                       (10.0, 40.0)):
            tc, sc = band_share(c, channels, lo, hi)
            th, sh = band_share(h, channels, lo, hi)
            p('')
            p('band %.1f - %.1f R_p    integrated total: cool %.4e  hot %.4e'
              ' erg/s   (hot/cool = %.3f)' % (lo, hi, tc, th, th / tc))
            p('  %-18s %10s %10s %10s' % ('channel', 'cool[%]', 'hot[%]',
                                          'hot/cool'))
            rows = sorted(channels, key=lambda k: -max(sc[k], sh[k]))
            for k in rows:
                if max(sc[k], sh[k]) < 5e-3:
                    continue
                ratio = (sh[k] * th) / max(sc[k] * tc, 1e-300)
                p('  %-18s %10.2f %10.2f %10.3f'
                  % (k, 100 * sc[k], 100 * sh[k], ratio))

        p('')
        p('%s at fixed radii: total and the two largest channels' % kind)
        for rq in (1.1, 1.3, 1.5, 2.0, 3.0):
            j = int(np.argmin(np.abs(c['r'] - rq)))
            jj = int(np.argmin(np.abs(h['r'] - rq)))
            top_c = sorted(channels, key=lambda k: -c[k][j])[:3]
            top_h = sorted(channels, key=lambda k: -h[k][jj])[:3]
            p('  r=%5.2f  cool T=%7.1f tot=%9.3e : %s' %
              (rq, c['T'][j], c['total'][j],
               ', '.join('%s %.0f%%' % (k, 100 * c[k][j] /
                                        max(c['total'][j], 1e-300))
                         for k in top_c)))
            p('          hot  T=%7.1f tot=%9.3e : %s' %
              (h['T'][jj], h['total'][jj],
               ', '.join('%s %.0f%%' % (k, 100 * h[k][jj] /
                                        max(h['total'][jj], 1e-300))
                         for k in top_h)))
        p('')

    # neutral fraction of the carbon that carries the C I cooling
    import sys
    sys.path.insert(0, HERE)
    from compare_branches import load
    lc = load(COOL)
    lh = load(HOT)
    p('=' * 78)
    p('carbon and oxygen ionization state (the neutral-line cooling carriers)')
    p('=' * 78)
    p('  %6s %12s %12s %12s %12s' % ('r[Rp]', 'f(CI) cool', 'f(CI) hot',
                                     'f(OI) cool', 'f(OI) hot'))
    for rq in (1.05, 1.1, 1.2, 1.3, 1.5, 2.0, 3.0, 5.0):
        j = int(np.argmin(np.abs(lc['r'] - rq)))
        def frac(d, el):
            tot = d['n_%sI' % el] + d['n_%sII' % el] + d['n_%sIII' % el]
            return d['n_%sI' % el][j] / max(tot[j], 1e-300)
        p('  %6.2f %12.4f %12.4f %12.4f %12.4f'
          % (lc['r'][j], frac(lc, 'C'), frac(lh, 'C'),
             frac(lc, 'O'), frac(lh, 'O')))

    with open(os.path.join(HERE, 'T4_channel_budget.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
