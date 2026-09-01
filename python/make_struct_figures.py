#!/usr/bin/env python3
"""make_struct_figures.py -- per-planet radius-structure figures for the EXHALE
method paper.

For each planet it writes a 2x3 multi-panel PDF ``fig_struct_<tag>.pdf`` to
``EXHALE/paper/`` showing, versus r/R_p:

    (a) temperature T(r)                                  [log-y]
    (b) radial velocity v(r) [km/s]                       [linear-y]
    (c) ion number densities: HI, HII, HeI, HeII, HeIII,
        He 2^3S, and the most abundant metal ions present [log-y]
    (d) cooling breakdown: each radiative channel from
        Cooling_breakdown.txt (recombination, collisional
        ionization, the collisional-excitation channel split
        into H I [Lya-dominated], He I and He II,
        bremsstrahlung, major metal-line coolants + "other
        metals"), and the total radiative cooling          [log-y]
    (e) heating breakdown: the photoionization heating of
        each absorber from Heating_breakdown.txt (H I, He I,
        He II, He 2^3S, metals), the excited-H photoelectric
        and Lya de-excitation heating, He-recombination and
        He 2^3S Penning heating, and the total heating      [log-y]
    (f) energy-budget totals: total heating, total radiative
        cooling, and the ADIABATIC (expansion) cooling
        computed offline from the hydro profile             [log-y]

Usage
-----
    python3 make_struct_figures.py                 # all four planet runs
    python3 make_struct_figures.py --rundir DIR --tag TAG --name "Name"

The run directories are the self-contained planet folders of the EXHALE tree
(HD209458b/, HD189733b/, WASP-52b/, WASP-121b/; see RUNDIR in paper_data.py);
pass --rundir for any other run.  Column maps and the adiabatic-cooling formula
live in paper_data.py.

matplotlib usetex is left ON (the machine default); all labels are ASCII or
LaTeX strings.
"""
import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import paper_data as pd

OUTDIR = os.path.join(HERE, '..', 'paper')

# LaTeX-safe legend labels.
ION_LABEL = {
    'HI': r'H\,{\sc i}', 'HII': r'H\,{\sc ii}',
    'HeI': r'He\,{\sc i}', 'HeII': r'He\,{\sc ii}', 'HeIII': r'He\,{\sc iii}',
    'HeITR': r'He\,$2^{3}\!S$',
    'CII': r'C\,{\sc ii}', 'CIII': r'C\,{\sc iii}',
    'OI': r'O\,{\sc i}', 'OII': r'O\,{\sc ii}', 'OIII': r'O\,{\sc iii}',
    'NII': r'N\,{\sc ii}', 'NIII': r'N\,{\sc iii}',
    'MgII': r'Mg\,{\sc ii}', 'MgIII': r'Mg\,{\sc iii}',
    'SiII': r'Si\,{\sc ii}', 'SiIII': r'Si\,{\sc iii}',
    'CaII': r'Ca\,{\sc ii}', 'CaIII': r'Ca\,{\sc iii}',
    'NaI': r'Na\,{\sc i}', 'NaII': r'Na\,{\sc ii}',
    'FeII': r'Fe\,{\sc ii}', 'FeIII': r'Fe\,{\sc iii}',
}

COOL_CH_LABEL = {
    'reco': 'recombination',
    'coio': 'collisional ionization',
    'coex_HI': r'coll.\ excitation H\,{\sc i} (Ly$\alpha$)',
    'coex_HeI': r'coll.\ excitation He\,{\sc i}',
    'coex_HeII': r'coll.\ excitation He\,{\sc ii}',
    'brem': 'bremsstrahlung',
    'H3p_IR': r'H$_3^+$ infrared',
    'H2_IR': r'H$_2$ infrared lines',
    'H2O_IR': r'H$_2$O infrared bands',
    'CO_IR': r'CO infrared bands',
}

# Ordered list of the fixed cooling channels drawn in panel (d).  H3p_IR is
# zero outside a molecular layer and is then dropped by the amplitude cut, as
# are the three molecular band channels without "Molecular IR bands".  Those
# three are NET rates and go negative where the band heats; panel (d) is
# logarithmic, so only the cooling part of such a channel is drawn.
COOL_CH_ORDER = ['reco', 'coio', 'coex_HI', 'coex_HeI', 'coex_HeII', 'brem',
                 'H3p_IR', 'H2_IR', 'H2O_IR', 'CO_IR']

HEAT_CH_LABEL = {
    'heat_HI': r'H\,{\sc i} photoion.',
    'heat_HeI': r'He\,{\sc i} photoion.',
    'heat_HeII': r'He\,{\sc ii} photoion.',
    'heat_He23S': r'He\,$2^{3}\!S$ photoion.',
    'heat_H2': r'H$_2$ photoion.',
    'heat_metals': 'metal photoion.',
    'heat_Hpe': r'H($n{=}2$) photoelectric',
    'heat_Hdx': r'Ly$\alpha$ de-excitation',
    'heat_He_recomb': r'He recomb.\ heating',
    'heat_He23S_Penning': r'He\,$2^{3}\!S$ Penning',
}

# Ordered list of the heating channels drawn in panel (e).
HEAT_CH_ORDER = ['heat_HI', 'heat_HeI', 'heat_HeII', 'heat_He23S',
                 'heat_H2', 'heat_metals', 'heat_Hpe', 'heat_Hdx',
                 'heat_He_recomb', 'heat_He23S_Penning']


def _ion_label(sp):
    return ION_LABEL.get(sp, sp)


def top_metal_ions(d, rmask, n=6):
    """Return the n metal ions with the largest peak density in the wind.

    The ranking is taken over r >= 1.1 R_p so that the extended-wind ions that
    actually carry the transit signal are selected, rather than the neutral
    metals that survive only in the thin high-density base layer.  The metal set
    is whatever the run actually carried (d.metal_ions), so a metals-off run
    simply yields an empty list.
    """
    wind = rmask & (d.r_ion >= 1.1)
    if not np.any(wind):
        wind = rmask
    peaks = []
    for sp in d.metal_ions:
        m = np.nanmax(d.ion[sp][wind])
        if np.isfinite(m) and m > 0:
            peaks.append((m, sp))
    peaks.sort(reverse=True)
    return [sp for _, sp in peaks[:n]]


def top_metal_coolants(d, n=4):
    """Return (majors, others) split of metal cooling channels by peak rate."""
    peaks = []
    for sp in d.metal_ions:
        v = d.cool_ch.get(sp)
        if v is None:
            continue
        m = np.nanmax(v)
        if np.isfinite(m) and m > 0:
            peaks.append((m, sp))
    peaks.sort(reverse=True)
    majors = [sp for _, sp in peaks[:n]]
    others = [sp for _, sp in peaks[n:]]
    return majors, others


def _peak(arr, mask):
    """Finite peak of arr over mask, or 0.0."""
    if arr is None:
        return 0.0
    v = np.asarray(arr)[mask]
    v = v[np.isfinite(v)]
    return float(np.nanmax(v)) if v.size and np.nanmax(v) > 0 else 0.0


def make_one(d, tag):
    # Show the structure out to the escape radius -- the escape model's physical
    # outer boundary, where the mass-loss rate is evaluated.  The grid extends a
    # little past it as padding; those cells are dropped from the plot.  NaN/Inf
    # ghost cells were already removed by the loader.
    rmax = min(d.r_esc, float(np.nanmax(d.r)))
    hmask = d.r <= rmax + 1e-9
    imask = d.r_ion <= rmax + 1e-9

    fig, axes = plt.subplots(2, 3, figsize=(13.6, 7.6))
    (axT, axV, axN), (axC, axH, axB) = axes

    # (a) Temperature.  Guard the y-scale against a single boundary-cell spike
    # (e.g. the WASP-121 b no-diffusion run) with a robust upper percentile.
    axT.semilogy(d.r[hmask], d.T[hmask], 'k-', lw=1.6)
    axT.set_ylabel(r'$T$ [K]')
    axT.set_title('(a) temperature')
    tlo = np.nanmin(d.T[hmask]) * 0.8
    thi = max(np.nanpercentile(d.T[hmask], 99.0) * 1.4,
              np.nanmedian(d.T[hmask]) * 2.0)
    axT.set_ylim(tlo, thi)

    # (b) Velocity [km/s]
    axV.plot(d.r[hmask], d.v[hmask] / 1e5, 'k-', lw=1.6)
    axV.set_ylabel(r'$v$ [km s$^{-1}$]')
    axV.set_title('(b) velocity')

    # (c) Ion densities
    base_species = [s for s in pd.HHE_SPECIES if s in d.ion]  # H/He (+He 2^3S)
    metals = top_metal_ions(d, imask, n=6)
    for sp in base_species:
        axN.semilogy(d.r_ion[imask], d.ion[sp][imask], lw=1.4,
                     label=_ion_label(sp))
    for sp in metals:
        axN.semilogy(d.r_ion[imask], d.ion[sp][imask], lw=1.1, ls='--',
                     label=_ion_label(sp))
    axN.set_ylabel(r'number density [cm$^{-3}$]')
    axN.set_title('(c) ion densities')
    # Keep a sensible dynamic range (8 decades below the peak).
    allmax = max(np.nanmax(d.ion[sp][imask]) for sp in base_species)
    axN.set_ylim(allmax * 1e-9, allmax * 3)
    axN.legend(fontsize=6.5, ncol=3, loc='lower center', framealpha=0.9)

    # Reference rate for the log y-range on the rate panels.
    ref = np.nanmax(d.heat[hmask])

    # (d) Cooling breakdown (collisional excitation split by absorber).
    if d.cool_ch:
        rc = d.r_cool
        cmask = rc <= rmax + 1e-9
        # Skip absorber channels that are negligible everywhere (e.g. coex_HeII
        # in a helium-poor wind) to keep the legend readable.
        for ch in COOL_CH_ORDER:
            if ch in d.cool_ch and _peak(d.cool_ch[ch], cmask) > ref * 1e-8:
                axC.semilogy(rc[cmask], d.cool_ch[ch][cmask], lw=1.2,
                             label=COOL_CH_LABEL[ch])
        majors, others = top_metal_coolants(d, n=4)
        for sp in majors:
            axC.semilogy(rc[cmask], d.cool_ch[sp][cmask], lw=1.0, ls='--',
                         label=pd.METAL_LABEL.get(sp, sp) + ' cooling')
        if others:
            osum = np.zeros_like(rc)
            for sp in others:
                osum = osum + np.nan_to_num(d.cool_ch[sp])
            axC.semilogy(rc[cmask], osum[cmask], lw=1.0, ls=':',
                         color='0.5', label='other metals')
        axC.semilogy(rc[cmask], d.cool_ch['cool_total'][cmask], 'k-', lw=1.8,
                     label='total cooling')
    axC.set_ylabel(r'rate [erg cm$^{-3}$ s$^{-1}$]')
    axC.set_title('(d) cooling breakdown')
    axC.set_ylim(ref * 1e-7, ref * 5)
    axC.legend(fontsize=6.0, ncol=2, loc='lower center', framealpha=0.9)

    # (e) Heating breakdown (photoionization per absorber + Balmer + He terms).
    if d.heat_ch:
        rh = d.r_heat
        emask = rh <= rmax + 1e-9
        for ch in HEAT_CH_ORDER:
            if ch in d.heat_ch and _peak(d.heat_ch[ch], emask) > ref * 1e-6:
                axH.semilogy(rh[emask], d.heat_ch[ch][emask], lw=1.2,
                             label=HEAT_CH_LABEL[ch])
        axH.semilogy(rh[emask], d.heat_ch['heat_total'][emask],
                     color='tab:red', lw=1.8, label='total heating')
    else:
        axH.semilogy(d.r[hmask], d.heat[hmask], color='tab:red', lw=1.8,
                     label='total heating')
    axH.set_ylabel(r'rate [erg cm$^{-3}$ s$^{-1}$]')
    axH.set_title('(e) heating breakdown')
    axH.set_ylim(ref * 1e-7, ref * 5)
    axH.legend(fontsize=6.0, ncol=2, loc='lower center', framealpha=0.9)

    # (f) Energy-budget totals: heating, radiative cooling, adiabatic cooling.
    axB.semilogy(d.r[hmask], d.heat[hmask], color='tab:red', lw=1.8,
                 label='total heating')
    if d.cool_ch:
        rc = d.r_cool
        cmask = rc <= rmax + 1e-9
        axB.semilogy(rc[cmask], d.cool_ch['cool_total'][cmask], color='tab:blue',
                     lw=1.8, label='total radiative cooling')
    axB.semilogy(d.r[hmask], d.adiab[hmask], color='tab:brown', lw=1.4,
                 ls='-.', label='adiabatic (expansion) cooling')
    axB.set_ylabel(r'rate [erg cm$^{-3}$ s$^{-1}$]')
    axB.set_title('(f) energy-budget totals')
    axB.set_ylim(ref * 1e-7, ref * 5)
    axB.legend(fontsize=6.5, loc='lower center', framealpha=0.9)

    for ax in (axT, axV, axN, axC, axH, axB):
        ax.set_xlabel(r'$r/R_{\rm p}$')
        ax.set_xlim(1.0, rmax)
        ax.grid(alpha=0.25)

    fig.suptitle(r'%s -- radial structure' % d.name.replace('_', r'\_'),
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = os.path.join(OUTDIR, 'fig_struct_%s.pdf' % tag)
    fig.savefig(out)
    plt.close(fig)
    print('wrote', out)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rundir', default=None,
                    help='single run dir (with input.inp + output/)')
    ap.add_argument('--tag', default=None, help='output tag for --rundir')
    ap.add_argument('--name', default=None, help='display name for --rundir')
    ap.add_argument('--base', default=None,
                    help='base dir holding the staged <tag>/ run dirs')
    args = ap.parse_args()

    if args.rundir:
        tag = args.tag or os.path.basename(os.path.normpath(args.rundir))
        d = pd.load_planet(args.rundir, name=args.name or tag)
        make_one(d, tag)
        return

    for name, tag, _ in pd.PLANETS:
        rundir = pd.planet_rundir(tag, base=args.base)
        d = pd.load_planet(rundir, name=name)
        make_one(d, tag)


if __name__ == '__main__':
    main()
