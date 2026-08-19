"""lya_insitu_emissivity.py -- in-situ Ly-alpha volume emissivity for EXHALE runs.

Builds the spherically-symmetric Ly-alpha volume emissivity eps(r) produced
*inside* an escaping planetary atmosphere (recombination + collisional
excitation), for use as the internal source of a LaRT diffuse-emissivity RT
run.  For each EXHALE benchmark run it writes:

  <run_dir>/insitu_emiss.txt : headerless 2-column  r[R_p]  eps_tot(r)
                               (LaRT reads this with read(unit,*); it multiplies
                                by r^2 internally, so only the SHAPE matters).
  <run_dir>/insitu_lum.txt   : L_rec / L_col / L_tot summary (erg/s, photons/s).
  <run_dir>/insitu_emiss.png : diagnostic 2-panel figure.

Run for all four default benchmarks in one invocation to also emit the combined
paper figure ATES/paper/figs/fig_lya_insitu.pdf.

Physics (matches EXHALE_transit.py::n2_populations exactly):
  t4  = T/1e4  (T clamped to >= 1 K)
  aB  = 2.54e-13 * t4**(-0.8163 - 0.0208*ln(t4))          case-B recomb
  a2s = (0.282 + 0.047*t4 - 0.006*t4**2) * aB
  a2p = aB - a2s                                            recombs landing in 2p
  C1s2p = 1.71e-8 * (1/t4)**0.077 * exp(-118400/T)          collisional 1s->2p
  n_e = n_HII + n_HeII + 2 n_HeIII
  eps_rec = a2p * n_e * n_HII * E_Lya
  eps_col = C1s2p * n_e * n_HI  * E_Lya
  eps_tot = eps_rec + eps_col          [erg cm^-3 s^-1]

Comparison to Yan et al. (2022, ApJ 936, 177), their Eq. 8:
  S_nu = (alpha_B * n_H+ * n_e + C^(e)_{1s->2} * n_1s * n_e) / (4 pi),
  with C^(e)_{1s->2} = C_{1s->2s} + C_{1s->2p}.
  i.e. Yan use the FULL case-B rate alpha_B for the recombination term and BOTH
  collisional channels (1s->2s AND 1s->2p).  We instead keep only the fraction
  of recombinations landing in 2p (a2p = aB - a2s) and only the 1s->2p
  collisional channel, because only 2p decays yield a Ly-alpha photon (2s decays
  via the two-photon continuum).  This is consistent with
  EXHALE_transit.py::n2_populations.  It also differs (a) from the in-code
  src/modules/radiation/lya_rt.f90, which uses the full alphaB, and (b) from the
  classic 0.68*alphaB rule of thumb.  Only the SHAPE of eps(r) matters for the
  LaRT sampling (normalization is external); the rec/col split and integrated
  luminosities below are reported physically.
"""

import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# make exhale_io importable regardless of cwd
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import exhale_io as aio

# --- physical constants (cgs) ---
RJ = 6.9911e9              # Jupiter radius [cm]
E_LYA_ERG = 1.634e-11      # 10.2 eV in erg

DEFAULT_RUNS = ['hd189', 'hd209', 'wasp121', 'wasp52']


def emissivity(run):
    """Return dict with r[R_p] and eps_rec/eps_col/eps_tot [erg cm^-3 s^-1].

    Matches EXHALE_transit.py::n2_populations for the atomic coefficients.
    """
    r = run.r
    T = np.maximum(run.T, 1.0)               # defensive floor (avoid 1/T, ln 0)
    t4 = T / 1.0e4

    aB = 2.54e-13 * t4 ** (-0.8163 - 0.0208 * np.log(t4))
    a2s = (0.282 + 0.047 * t4 - 0.006 * t4 ** 2.0) * aB
    a2p = aB - a2s
    C1s2p = 1.71e-8 * (1.0 / t4) ** 0.077 * np.exp(-118400.0 / T)

    n_HI = run.ion['HI']
    n_HII = run.ion['HII']
    n_e = run.ion['HII'] + run.ion['HeII'] + 2.0 * run.ion['HeIII']

    eps_rec = a2p * n_e * n_HII * E_LYA_ERG
    eps_col = C1s2p * n_e * n_HI * E_LYA_ERG
    eps_tot = eps_rec + eps_col
    return dict(r=r, eps_rec=eps_rec, eps_col=eps_col, eps_tot=eps_tot)


def luminosity(run, em):
    """Integrate L_X = int 4 pi (r R_p)^2 eps_X dr over the run radius [erg/s].

    Returns (dict of L_X [erg/s], dict of cumulative L_tot(<r) [erg/s]).
    """
    Rp_cm = run.inp['Rp_RJ'] * RJ
    r_cm = run.r * Rp_cm
    shell = 4.0 * np.pi * r_cm ** 2.0
    lum = {}
    for key in ('rec', 'col', 'tot'):
        integrand = shell * em['eps_' + key]
        lum['L_' + key] = float(np.trapz(integrand, r_cm))
    cum_tot = np.concatenate(([0.0], np.cumsum(
        0.5 * (shell[1:] * em['eps_tot'][1:] + shell[:-1] * em['eps_tot'][:-1])
        * np.diff(r_cm))))
    return lum, cum_tot


def write_emiss_file(path, r, eps_tot):
    """Write the headerless 2-column LaRT emissivity file."""
    data = np.column_stack([r, eps_tot])
    # headerless: LaRT reader does read(unit,*) and would choke on a comment line
    np.savetxt(path, data, fmt='%.6e')


def write_lum_file(path, name, lum):
    L_rec, L_col, L_tot = lum['L_rec'], lum['L_col'], lum['L_tot']
    frac_rec = L_rec / L_tot if L_tot > 0 else float('nan')
    frac_col = L_col / L_tot if L_tot > 0 else float('nan')
    lines = [
        '# in-situ Ly-alpha luminosity summary for %s' % name,
        '# X    L_X[erg/s]      photons/s(=L_X/E_Lya, E_Lya=%.4e erg)' % E_LYA_ERG,
        'rec    %.6e   %.6e' % (L_rec, L_rec / E_LYA_ERG),
        'col    %.6e   %.6e' % (L_col, L_col / E_LYA_ERG),
        'tot    %.6e   %.6e' % (L_tot, L_tot / E_LYA_ERG),
        '# rec:col fraction of L_tot = %.4f : %.4f' % (frac_rec, frac_col),
    ]
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def plot_single(path, name, em, cum_tot):
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = 'serif'
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9, 3.4))
    r = em['r']
    ax0.plot(r, em['eps_rec'], label='eps_rec (recomb)', color='C0')
    ax0.plot(r, em['eps_col'], label='eps_col (coll)', color='C1')
    ax0.plot(r, em['eps_tot'], label='eps_tot', color='k', lw=1.6)
    ax0.set_yscale('log')
    ax0.set_xlabel('r [R_p]')
    ax0.set_ylabel('eps [erg/cm^3/s]')
    ax0.set_title('%s in-situ Ly-alpha emissivity' % name)
    ax0.legend(fontsize=8)

    ax1.plot(r, cum_tot, color='k')
    ax1.set_xlabel('r [R_p]')
    ax1.set_ylabel('L(<r) [erg/s]')
    ax1.set_title('cumulative Ly-alpha luminosity')
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def plot_combined(path, results):
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = 'serif'
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(8, 3.2))
    for i, (name, em, cum_tot, lum) in enumerate(results):
        c = 'C%d' % i
        ax0.plot(em['r'], em['eps_tot'], color=c, label=name)
        L_tot = lum['L_tot']
        norm = cum_tot / L_tot if L_tot > 0 else cum_tot
        ax1.plot(em['r'], norm, color=c, label=name)
    ax0.set_yscale('log')
    ax0.set_xlabel('r [R_p]')
    ax0.set_ylabel('eps_tot [erg/cm^3/s]')
    ax0.legend(fontsize=8)
    ax1.set_xlabel('r [R_p]')
    ax1.set_ylabel('L(<r) / L_tot')
    ax1.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def process(run_dir):
    """Process one benchmark run directory; returns (name, em, cum_tot, lum)."""
    name = os.path.basename(os.path.normpath(run_dir))
    outdir = os.path.join(run_dir, 'output')
    inputfile = os.path.join(run_dir, 'input.inp')
    run = aio.load_run(outdir, inputfile, adv=True)
    em = emissivity(run)
    lum, cum_tot = luminosity(run, em)

    write_emiss_file(os.path.join(run_dir, 'insitu_emiss.txt'), em['r'], em['eps_tot'])
    write_lum_file(os.path.join(run_dir, 'insitu_lum.txt'), name, lum)
    plot_single(os.path.join(run_dir, 'insitu_emiss.png'), name, em, cum_tot)

    frac_rec = lum['L_rec'] / lum['L_tot'] if lum['L_tot'] > 0 else float('nan')
    frac_col = lum['L_col'] / lum['L_tot'] if lum['L_tot'] > 0 else float('nan')
    print('%-10s  L_rec=%.4e  L_col=%.4e  L_tot=%.4e erg/s | '
          'phot_tot=%.4e /s | rec:col=%.3f:%.3f'
          % (name, lum['L_rec'], lum['L_col'], lum['L_tot'],
             lum['L_tot'] / E_LYA_ERG, frac_rec, frac_col))
    return name, em, cum_tot, lum


def _resolve(arg):
    """Resolve a CLI arg (either a benchmark short name or a run dir path)."""
    if os.path.isdir(os.path.join(arg, 'output')):
        return arg
    cand = os.path.join(_HERE, '..', 'benchmarks', arg)
    return cand


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    repo = os.path.normpath(os.path.join(_HERE, '..'))
    if argv:
        run_dirs = [_resolve(a) for a in argv]
        all_four = False
    else:
        run_dirs = [os.path.join(repo, 'benchmarks', p) for p in DEFAULT_RUNS]
        all_four = True

    results = []
    for d in run_dirs:
        results.append(process(d))

    if all_four:
        # paper figs live inside the EXHALE tree: EXHALE/paper/figs
        figdir = os.path.join(repo, 'paper', 'figs')
        os.makedirs(figdir, exist_ok=True)
        pdf = os.path.join(figdir, 'fig_lya_insitu.pdf')
        plot_combined(pdf, results)
        print('combined figure -> %s' % pdf)
    return 0


if __name__ == '__main__':
    sys.exit(main())
