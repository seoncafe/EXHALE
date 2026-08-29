#!/usr/bin/env python3
"""Does the "hot branch" of the K_zz scan coincide with the metal dropout?

One row per arm of `kzz_profile_scan`, pairing the mass-loss rate and the
He 10830 equivalent width that `results.txt` reports with a single measured
property of the same solution: whether its C/N/O elemental density is exactly
zero over a band of cells.  Elemental carbon has no sink in these equations,
so a band of exactly zero carbon between two bands at the reservoir ratio is
not a state of the atmosphere.

Usage:  python3 branch_vs_dropout.py
"""
import os
import glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCAN = os.path.join(HERE, '../kzz_profile_scan')


def last_k(arm):
    ks = sorted(glob.glob(os.path.join(arm, 'k[0-9][0-9]')))
    return ks[-1] if ks else None


def dropout(outdir):
    a = np.loadtxt(os.path.join(outdir, 'Ion_species.txt'))
    r, C = a[:, 0], a[:, 7] + a[:, 8] + a[:, 9]
    z = np.where(C <= 0.0)[0]
    return (len(z), (r[z[0]], r[z[-1]]) if len(z) else None)


def measured_rows():
    """log10 Mdot and the red-pair EW as `results.txt` tabulated them."""
    rows = {}
    path = os.path.join(SCAN, 'results.txt')
    started = False
    for ln in open(path, errors='replace'):
        if ln.startswith('case '):
            started = True
            continue
        if not started:
            continue
        f = ln.split()
        if len(f) < 14 or not f[0].startswith(('kzz', '../')):
            if started and ln.startswith('='):
                break
            continue
        rows[f[0]] = (f[1], f[2], f[10], f[13])   # Kzz reserv lgMdot EW
    return rows


def main():
    tab = measured_rows()
    out = []
    for arm in sorted(glob.glob(os.path.join(SCAN, 'kzz*', '*'))):
        if not os.path.isdir(arm):
            continue
        k = last_k(arm)
        if k is None or not os.path.isfile(
                os.path.join(k, 'output', 'Ion_species.txt')):
            continue
        name = os.path.relpath(arm, SCAN)
        n, band = dropout(os.path.join(k, 'output'))
        row = tab.get(name, ('?', '?', '?', '?'))
        out.append((name, row, n, band))

    lines = ['one row per arm: the mass-loss rate and the He 10830 red-pair'
             ' equivalent width of the arm, against whether that solution has'
             ' a band of exactly zero C/N/O',
             '',
             '%-30s %9s %8s %9s %9s %6s  %s'
             % ('arm', 'K_zz', 'reserv', 'lgMdot', 'EW[%A]', 'ncell',
                'zero-C/N/O band [R_p]')]
    for name, row, n, band in sorted(out, key=lambda t: t[1][2]):
        lines.append('%-30s %9s %8s %9s %9s %6d  %s'
                     % (name, row[0], row[1], row[2], row[3], n,
                        '--' if band is None else '%.4f - %.4f' % band))

    md = [float(r[2]) for _, r, n, _ in out if n > 0 and r[2] not in ('?', 'nan')]
    mc = [float(r[2]) for _, r, n, _ in out if n == 0 and r[2] not in ('?', 'nan')]
    lines += ['',
              'log10 Mdot with a dropout   : n=%d, %.3f - %.3f (median %.3f)'
              % (len(md), min(md), max(md), float(np.median(md))),
              'log10 Mdot with none        : n=%d, %.3f - %.3f (median %.3f)'
              % (len(mc), min(mc), max(mc), float(np.median(mc))),
              '',
              'the two sets %s'
              % ('do not overlap' if min(md) > max(mc) else 'overlap')]
    txt = '\n'.join(lines)
    print(txt)
    open(os.path.join(HERE, 'branch_vs_dropout.txt'), 'w').write(txt + '\n')


if __name__ == '__main__':
    main()
