#!/usr/bin/env python3
"""make_transit_figures.py -- transmission-spectrum figures for the EXHALE
method paper.

For each diagnostic line it writes a 2x2 (one panel per planet) PDF to
``EXHALE/paper/`` comparing the EXHALE model transit spectrum against published
observations where they exist:

    fig_transit_He10830.pdf   He I 10830 triplet
    fig_transit_Halpha.pdf    H-alpha 6562.8
    fig_transit_Lya.pdf       Ly-alpha 1215.67  (velocity axis)

Common y-axis
-------------
Everything is shown as EXCESS ABSORPTION in percent, A = (1 - T) * 100, where T
is the in-/out-of-transit flux ratio.  The model curve is the rotation- plus
instrument-convolved transmission (column 4 of the tpm_*.txt files).  Each
observation is converted to the same axis from its own native convention (see
the OBS registry): dF/F and TS = F_in/F_out - 1 give A = -y*100; a normalized
in/out flux gives A = (1 - y)*100; a file already in "excess absorption [%]" is
used directly.  This is a fair common axis; it does NOT force a match.

Model curves
------------
The tpm_*.txt curves come from EXHALE_transit.py.  With --run this script
regenerates them by invoking that post-processor on each run directory, which
writes the canonical tpm_<line>.txt set there; without --run it loads whatever
tpm files are already present.

Usage
-----
    python3 make_transit_figures.py            # load tpm_*.txt, make figures
    python3 make_transit_figures.py --run      # regenerate tpm_*.txt first
    python3 make_transit_figures.py --run --tags hd189 wasp52 --no-figures
                                               # curves only, chosen planets

Run directories are the self-contained planet folders of the EXHALE tree
(RUNDIR in paper_data.py; relocate the set with --base).  The planet and
stellar parameters -- radius, T_eff, orbit, and the tidally locked spin period
derived from it -- are read by EXHALE_transit.py from each folder's input.inp,
which is the single source; this script passes none of them.  matplotlib usetex
is left ON; all labels are ASCII or LaTeX strings.
"""
import os
import sys
import glob
import argparse
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import paper_data as pd

EXHALE_DIR = pd.EXHALE_DIR
TRANSIT_PY = os.path.join(EXHALE_DIR, 'EXHALE_transit.py')
OUTDIR = os.path.join(HERE, '..', 'paper')
OBS_A = os.path.expanduser('~/Exoplanetary_Atmosphere')   # "A" root in the doc
OBS_E = EXHALE_DIR                                         # "E" root in the doc

C_KMS = 2.99792458e5
LYA0 = 1215.67   # Ly-alpha line center [Angstrom]

# Panel order and titles.  The run directory of each tag comes from
# paper_data.RUNDIR, and every physical parameter from that directory's
# input.inp -- nothing about the planets is duplicated here.
ORDER = ['hd209', 'hd189', 'wasp52', 'wasp121']
PLANET_NAME = {tag: name for name, tag, _ in pd.PLANETS}

# Observation registry.  Each entry:
#   file, xcol, ycol, errcol (or None), conv, xkind, frame, label
# conv maps the native y to excess absorption [%]:
#   'dFF'/'TS'  -> A = -y*100     (dF/F, or TS = F_in/F_out - 1)
#   'pct_abs'   -> A = y          (already excess absorption in %)
#   'norm_flux' -> A = (1-y)*100  (normalized in/out flux)
# xkind: 'wave' [Angstrom] or 'vel' [km/s].
# frame: wavelength reference of the tabulated x column.  The EXHALE model line
#   centers are AIR (He 10830.34, H-alpha 6562.80); Ly-alpha is quoted in vacuum
#   but plotted on a velocity axis so the reference is immaterial.  Data listed
#   in VACUUM are converted to air (vac_to_air) before plotting so obs and model
#   share one air wavelength scale; 'air' data are used as-is.  Verified per file
#   by checking the absorption feature sits at the air line center after the
#   conversion (He vac center 10833.2 -> air 10830.2; Salz He already at 10830.3;
#   Jensen H-alpha dip at 6562.8 = air).
OBS = {
    'He10830': {
        'hd189': (os.path.join(OBS_E, 'observational_data',
                               'salz2018_HD189733b_He.txt'),
                  0, 1, None, 'dFF', 'wave', 'air', 'Salz et al. 2018'),
        'wasp52': (os.path.join(OBS_A, 'WASP-52b',
                                'wasp52_He10830_trans_spec.dat'),
                   0, 1, 2, 'pct_abs', 'wave', 'vac', 'Kirk et al. 2022'),
        'wasp121': (os.path.join(OBS_A, 'WASP-121b',
                                 'WASP-121b_3D-Models_RT_code',
                                 'Transmission_spec_observation.dat'),
                    0, 1, 2, 'norm_flux', 'wave', 'vac', 'Czesla et al. 2024'),
    },
    'Halpha': {
        'hd209': (os.path.join(OBS_E, 'HD209458b',
                               'jensen2012_HD209458b_Ha.txt'),
                  0, 1, 2, 'dFF', 'wave', 'air', 'Jensen et al. 2012'),
        'hd189': (os.path.join(OBS_E, 'HD189733b',
                               'jensen2012_HD189733b_Ha.txt'),
                  0, 1, 2, 'dFF', 'wave', 'air', 'Jensen et al. 2012'),
        'wasp52': (os.path.join(OBS_A, 'WASP-52b',
                                'WASP52b_Halpha_transpec_0p10_AA_binned.txt'),
                   0, 1, 2, 'TS', 'wave', 'air', 'Chen et al. 2020'),
    },
    'Lya': {
        'hd209': (os.path.join(OBS_E, 'observational_data',
                               'vidalmadjar2003_HD209458b_Lya.txt'),
                  0, 1, None, 'dFF', 'vel', 'vac', 'Vidal-Madjar et al. 2003'),
    },
}

LINE_TITLE = {
    'He10830': r'He\,{\sc i} 10830 \AA',
    'Halpha': r'H$\alpha$ 6562.8 \AA',
    'Lya': r'Ly$\alpha$ 1215.67 \AA',
}


# --------------------------------------------------------------------------- #
def vac_to_air(lv):
    """Vacuum -> air wavelength [Angstrom] (Ciddor 1996 / Morton 2000 fit).

    Same relation used by the companion poster analysis.  The correction is
    about -2.9 A at the He I 10830 triplet, so vacuum-referenced He data must be
    converted before overplotting the air-referenced EXHALE model.
    """
    lv = np.asarray(lv, float)
    s = 1.0e4 / lv
    return lv / (1.0 + 8.34254e-5 + 0.02406147 / (130.0 - s**2)
                 + 0.00015998 / (38.9 - s**2))


def run_transit(base=None, tags=None):
    """Regenerate the canonical tpm_<line>.txt set in each planet folder via
    EXHALE_transit.py.

    Only the run directory is passed; the planet, star, and spin period come
    from that directory's input.inp, and the curves are written there under
    their canonical names, so no naming or parameter override is needed here.
    """
    for tag in (tags or ORDER):
        rundir = pd.planet_rundir(tag, base=base)
        env = dict(os.environ)
        env['MPLBACKEND'] = 'Agg'
        env['EXHALE_TRANSIT_PATH'] = rundir
        for _stale in ('EXHALE_TRANSIT_SAVE_PREFIX', 'TPM_SAVE_PREFIX',
                       'EXHALE_TRANSIT_RSTAR_RSUN', 'TPM_RSTAR_RSUN',
                       'EXHALE_TRANSIT_TSTAR', 'TPM_TSTAR',
                       'EXHALE_TRANSIT_ROTP', 'TPM_ROTP'):
            env.pop(_stale, None)   # input.inp is the single source
        print('running transit for', tag, '->', rundir)
        subprocess.run([sys.executable, TRANSIT_PY], env=env,
                       cwd=EXHALE_DIR, check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def load_tpm(rundir, line):
    """Load a tpm model curve; returns (lambda[A], A_model[%]) or None.

    The canonical name is 'tpm_<line>.txt', which EXHALE_transit.py now always
    writes into the run directory.  A doubled 'tpm_tpm_<line>.txt' from an
    earlier generation, when the save prefix had to be given by hand, is still
    accepted so an older run directory keeps loading.
    """
    for pat in ('tpm_%s.txt' % line, 'tpm_tpm_%s.txt' % line,
                '*tpm_%s.txt' % line):
        hits = sorted(glob.glob(os.path.join(rundir, pat)))
        if hits:
            m = np.loadtxt(hits[0])
            lam = m[:, 0]
            A = (1.0 - m[:, 3]) * 100.0   # col 4 = T_rot+instr
            return lam, A
    return None


def load_obs(entry):
    """Load an observation entry -> (x, A[%], Aerr[%] or None, label, xkind).

    Wavelength columns tagged 'vac' are converted to air (vac_to_air) so obs and
    the air-referenced EXHALE model share one air wavelength scale; 'air' data
    (and velocity axes) are left untouched.
    """
    path, xcol, ycol, ecol, conv, xkind, frame, label = entry
    if not os.path.exists(path):
        return None
    m = np.loadtxt(path, comments='#')
    x = m[:, xcol]
    if xkind == 'wave' and frame == 'vac':
        x = vac_to_air(x)
    y = m[:, ycol]
    e = m[:, ecol] if ecol is not None else None
    if conv in ('dFF', 'TS'):
        A = -y * 100.0
        Ae = e * 100.0 if e is not None else None
    elif conv == 'pct_abs':
        A = y
        Ae = e if e is not None else None
    elif conv == 'norm_flux':
        A = (1.0 - y) * 100.0
        Ae = e * 100.0 if e is not None else None
    else:
        raise ValueError('unknown conv %s' % conv)
    return x, A, Ae, label, xkind


def _to_axis_x(lam, line):
    """Model wavelength [A] -> panel x (velocity for Lya, else wavelength)."""
    if line == 'Lya':
        return (lam - LYA0) / LYA0 * C_KMS
    return lam


def make_line_figure(line, base=None):
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.2))
    axes = axes.ravel()
    xlabel = (r'velocity [km s$^{-1}$]' if line == 'Lya'
              else r'wavelength [\AA] (air)')
    for ax, tag in zip(axes, ORDER):
        rundir = pd.planet_rundir(tag, base=base)
        tpm = load_tpm(rundir, line)
        if tpm is not None:
            lam, A = tpm
            ax.plot(_to_axis_x(lam, line), A, 'k-', lw=1.6,
                    label='EXHALE (rot.+instr.)', zorder=3)
        obs_entry = OBS.get(line, {}).get(tag)
        if obs_entry is not None:
            obs = load_obs(obs_entry)
            if obs is not None:
                x, Ao, Ae, label, xkind = obs
                if xkind == 'wave' and line == 'Lya':
                    x = (x - LYA0) / LYA0 * C_KMS
                ax.errorbar(x, Ao, yerr=Ae, fmt='o', ms=2.6, lw=0.7,
                            color='tab:red', alpha=0.7, capsize=0,
                            label=label, zorder=2)
        ax.axhline(0.0, color='0.6', lw=0.6)
        ax.set_title(PLANET_NAME[tag])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'excess absorption [\%]')
        ax.legend(fontsize=7.5, loc='best')
        ax.grid(alpha=0.25)
        # Focus the x-window on the model line core.
        if tpm is not None:
            lam, A = tpm
            xm = _to_axis_x(lam, line)
            ipk = np.nanargmax(A)
            if line == 'Lya':
                ax.set_xlim(-500, 500)
            else:
                half = 3.0 if line == 'He10830' else 6.0
                ax.set_xlim(xm[ipk] - half, xm[ipk] + half)

    fig.suptitle(r'Transit transmission: %s' % LINE_TITLE[line], fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = os.path.join(OUTDIR, 'fig_transit_%s.pdf' % line)
    fig.savefig(out)
    plt.close(fig)
    print('wrote', out)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true',
                    help='regenerate tpm_*.txt via EXHALE_transit.py first')
    ap.add_argument('--base', default=None,
                    help='base dir holding the planet run dirs '
                         '(default: the EXHALE tree)')
    ap.add_argument('--tags', nargs='+', default=None, choices=ORDER,
                    help='restrict --run to these planets '
                         '(e.g. skip a folder with a run in progress)')
    ap.add_argument('--no-figures', action='store_true',
                    help='stop after --run; do not touch the paper figures '
                         '(use while some planets are still re-converging, so '
                         'a panel is never drawn from a stale curve)')
    args = ap.parse_args()
    if args.run:
        run_transit(base=args.base, tags=args.tags)
    if args.no_figures:
        return
    for line in ('He10830', 'Halpha', 'Lya'):
        make_line_figure(line, base=args.base)


if __name__ == '__main__':
    main()
