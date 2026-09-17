#!/usr/bin/env python3
"""Overlay a solved LHS 1140 b case on the archived solution of the same case.

    models/compare_archive.py <case dir> <archived dir> <out.pdf>

Six panels against radius: temperature, velocity, mass density, the hydrogen
ionized fraction, the He(2^3S) density, and -- where both directories carry a
`tpm_He10830.txt` -- the He I 10830 transmission profile against wavelength.
The archived curve is the 2026-08-30 solution as it stands (its code, its
grid, its transit reader); the new one is whatever the case directory holds
now.  The point of the figure is the comparison the user's rule asks for
before any number from a re-solved case is quoted: the whole profile, not one
or two representative values.

A table of the four measures the He I 10830 line is judged by -- the red-pair
equivalent width, the red-pair depth, the line width and log10 Mdot -- is
written to stdout for whichever of the two sides carries them.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# The measurement's own vacuum window and the air/vacuum factor of
# `LHS1140b/make_memo_figures.py`, so the equivalent width here is the number
# the memo and the runners quote.
AIR = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20


def species_columns(path):
    """The species names of an `Ion_species.txt`.

    Files written before the species header was separated carry the names on
    the `# columns` line instead, and both forms start with the radius.
    """
    names = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            if line.startswith('# species_columns'):
                return line.split()[3:]
            if line.startswith('# columns'):
                names = line.split()[2:]
    return names


def state(run_dir):
    """r, T, v, rho, x(H+) and n(He 2^3S) of a solved run."""
    hyd = os.path.join(run_dir, 'output', 'Hydro_ioniz.txt')
    ion = os.path.join(run_dir, 'output', 'Ion_species.txt')
    h = np.loadtxt(hyd)
    s = np.loadtxt(ion)
    names = species_columns(ion)
    idx = {n: i for i, n in enumerate(names)}
    phys = slice(2, h.shape[0] - 2)
    h, s = h[phys], s[phys]
    n_hi, n_hii = s[:, idx['HI']], s[:, idx['HII']]
    out = dict(r=h[:, 0], rho=h[:, 1], v=h[:, 2], p=h[:, 3], T=h[:, 4],
               x_hii=n_hii/np.maximum(n_hi + n_hii, 1e-300),
               n_he23s=s[:, idx['HeITR']])
    return out


def line_profile(run_dir):
    path = os.path.join(run_dir, 'tpm_He10830.txt')
    if not os.path.isfile(path):
        return None
    d = np.loadtxt(path)
    lam = d[:, 0]*AIR
    exc = (d[:, 2].max() - d[:, 2])/d[:, 2].max()*100.0
    return lam, exc


def line_measures(run_dir):
    """Equivalent width, red-pair depth, FWHM and log10 Mdot."""
    out = {}
    prof = line_profile(run_dir)
    if prof is not None:
        lam, exc = prof
        m = (lam >= EW_LO) & (lam <= EW_HI)
        out['EW'] = float(np.trapz(exc[m], lam[m]))
        out['FWHM'] = fwhm(lam[m], exc[m])
    # The depth and the width of the three-Gaussian fit beside the curve are
    # what the runners and status.py quote; where they exist they replace the
    # width measured off the sampled curve above, so one quantity does not
    # carry two definitions in two tables.
    metrics = os.path.join(run_dir, 'tpm_He10830_metrics.txt')
    if os.path.isfile(metrics):
        for line in open(metrics):
            word = line.split()
            if len(word) >= 2 and word[0] == 'red_depth':
                out['depth'] = float(word[1])
            elif len(word) >= 2 and word[0] == 'fwhm_A':
                out['FWHM'] = float(word[1])
    for name in ('pp.log', 'transit.log', 'run.log'):
        path = os.path.join(run_dir, name)
        if not os.path.isfile(path):
            continue
        for line in open(path, errors='replace'):
            if 'steady-state Mdot' in line:
                out['logMdot'] = float(line.split('=')[-1].split()[0])
    return out


def fwhm(lam, exc):
    """Full width of the strongest feature at half its own maximum."""
    if exc.size < 3 or exc.max() <= 0.0:
        return float('nan')
    half = 0.5*exc.max()
    above = np.flatnonzero(exc >= half)
    if above.size < 2:
        return float('nan')
    return float(lam[above[-1]] - lam[above[0]])


def main():
    if len(sys.argv) != 4:
        sys.stderr.write(__doc__)
        return 1
    new_dir, old_dir, out_pdf = sys.argv[1:4]
    new, old = state(new_dir), state(old_dir)

    plt.rcParams.update({'font.size': 8, 'axes.linewidth': 0.7,
                         'lines.linewidth': 1.1})
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.0))
    panels = [('T', r'$T$ [K]', 'log'),
              ('v', r'$v$ [cm s$^{-1}$]', 'symlog'),
              ('rho', r'$\rho$ [$m_{\rm H}$ cm$^{-3}$]', 'log'),
              ('x_hii', r'$n({\rm H^+})/n({\rm H})$', 'linear'),
              ('n_he23s', r'$n({\rm He}\,2^3{\rm S})$ [cm$^{-3}$]', 'log')]
    for ax, (key, label, scale) in zip(axes.ravel(), panels):
        ax.plot(old['r'], old[key], color='0.55', label='archived 2026-08-30')
        ax.plot(new['r'], new[key], color='crimson', label='this solution')
        ax.set_xscale('log')
        ax.set_yscale(scale)
        ax.set_xlabel(r'$r$ [$R_{\rm p}$]')
        ax.set_ylabel(label)
        ax.grid(alpha=0.25, lw=0.4)
    axes[0, 0].legend(frameon=False, fontsize=7)

    ax = axes[1, 2]
    drawn = False
    for run_dir, color, label in ((old_dir, '0.55', 'archived 2026-08-30'),
                                  (new_dir, 'crimson', 'this solution')):
        prof = line_profile(run_dir)
        if prof is None:
            continue
        lam, exc = prof
        ax.plot(lam, exc, color=color, label=label)
        drawn = True
    if drawn:
        ax.set_xlim(10830.0, 10836.0)
        ax.axvspan(EW_LO, EW_HI, color='0.85', zorder=0)
        ax.set_xlabel(r'$\lambda$ (air) [\AA]')
        ax.set_ylabel(r'excess absorption [\%]')
        ax.grid(alpha=0.25, lw=0.4)
    else:
        ax.text(0.5, 0.5, 'no He I 10830 profile on either side',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_axis_off()

    fig.tight_layout()
    fig.savefig(out_pdf)
    print('wrote %s' % out_pdf)

    rows = [('archived', line_measures(old_dir)),
            ('this solution', line_measures(new_dir))]
    print('%-16s %10s %10s %10s %10s'
          % ('', 'EW [%A]', 'depth [%]', 'FWHM [A]', 'log10 Mdot'))
    for name, m in rows:
        print('%-16s %10s %10s %10s %10s'
              % (name,
                 '%.5f' % m['EW'] if 'EW' in m else '-',
                 '%.5f' % m['depth'] if 'depth' in m else '-',
                 '%.4f' % m['FWHM'] if 'FWHM' in m else '-',
                 '%.2f' % m['logMdot'] if 'logMdot' in m else '-'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
