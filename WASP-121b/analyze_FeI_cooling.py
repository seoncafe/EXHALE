"""Fe I cooling validation for WASP-121b (Huang et al. 2023, Fig. 5 / Fig. 10).

Two panels:
  (left)  Fe I line-cooling COEFFICIENT Lambda_FeI(T) per (n_e n_FeI)
          [erg cm^3 s^-1] - the tabulated curve ported into Cool_coeff.f90.
          This is the quantity Huang shows in Fig. 5 (shape/magnitude check).
  (right) Per-channel volumetric cooling rate vs radius from the converged
          WASP-121b run (output/Cooling_breakdown.txt), Fe I highlighted
          among the dominant coolants (Huang Fig. 10 style).

Reads the cooling table straight from Cool_coeff.f90 so the left panel is
exactly what the Fortran interpolates (no offline re-derivation).
"""
import os
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
COOLF = os.path.join(HERE, 'output', 'Cooling_breakdown.txt')
COOL_F90 = os.path.normpath(os.path.join(
    HERE, '..', 'src', 'modules', 'radiation', 'Cool_coeff.f90'))


def read_f90_array(path, name):
    """Parse a `real*8, parameter :: name(NCOOLT) = [ ... ]` block."""
    txt = open(path).read()
    m = re.search(name + r'\s*\(NCOOLT\)\s*=\s*\[(.*?)\]', txt, re.S)
    if not m:
        raise RuntimeError('array %s not found in %s' % (name, path))
    vals = re.findall(r'[-+]?\d+\.\d+(?:[dDeE][-+]?\d+)?', m.group(1))
    return np.array([float(v.replace('d', 'e').replace('D', 'e'))
                     for v in vals])


def metal_ion_order(header_lines):
    """Recover the metal-ion column order from the file's comment header."""
    for ln in header_lines:
        if 'metal-ion columns' in ln:
            toks = ln.split(':', 1)[1].split()
            return toks
    return None


def main():
    # ---- left panel: Lambda_FeI(T) coefficient from the ported table ----
    logT = read_f90_array(COOL_F90, 'cool_logT')
    logL_FeI = read_f90_array(COOL_F90, 'cool_logL_FeI')
    logL_FeII = read_f90_array(COOL_F90, 'cool_logL_FeII')
    T = 10.0**logT

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.5, 5.0))

    axL.plot(T, 10.0**logL_FeI, 'b-', lw=2, label=r'Fe\,I (this work)')
    axL.plot(T, 10.0**logL_FeII, 'r--', lw=1.5, alpha=0.8,
             label=r'Fe\,II (coronal, ref.)')
    axL.set_xscale('log')
    axL.set_yscale('log')
    axL.set_xlim(2e3, 3e4)
    axL.set_ylim(1e-22, 3e-18)
    axL.set_xlabel(r'$T\ \mathrm{[K]}$')
    axL.set_ylabel(r'$\Lambda\ /\ (n_e\,n_{\rm ion})\ '
                   r'\mathrm{[erg\,cm^3\,s^{-1}]}$')
    axL.set_title(r'Fe\,I line-cooling coefficient (Huang Fig.~5)')
    axL.legend(loc='lower right', frameon=False)
    axL.grid(True, which='both', alpha=0.25)

    # ---- right panel: per-channel cooling vs radius ----
    if not os.path.exists(COOLF):
        axR.text(0.5, 0.5, 'Cooling_breakdown.txt not found\n(run not finished)',
                 ha='center', va='center', transform=axR.transAxes)
    else:
        header = [ln.rstrip('\n') for ln in open(COOLF) if ln.startswith('#')]
        order = metal_ion_order(header)
        d = np.loadtxt(COOLF)
        r = d[:, 0]
        cool_tot = d[:, 3]
        reco, coio, coex, brem = d[:, 4], d[:, 5], d[:, 6], d[:, 7]
        nmet = d.shape[1] - 8
        met = {order[i]: d[:, 8 + i] for i in range(nmet)} if order \
            else {str(i): d[:, 8 + i] for i in range(nmet)}

        sel = (r >= 1.0) & (r <= 1.5)
        axR.plot(r[sel], cool_tot[sel], 'k-', lw=2.5, label='total')
        axR.plot(r[sel], (reco + coio + coex + brem)[sel], color='0.5',
                 lw=1.2, ls=':', label=r'H/He (rec+ion+exc+ff)')

        # dominant metal coolants by peak contribution in the window
        peaks = {k: np.nanmax(v[sel]) for k, v in met.items()}
        top = sorted(peaks, key=peaks.get, reverse=True)[:6]
        if 'FeI' not in top:
            top = top[:5] + ['FeI']
        styles = {'FeI': dict(color='b', lw=2.5, ls='-', zorder=5)}
        for k in top:
            st = styles.get(k, dict(lw=1.3, ls='-'))
            axR.plot(r[sel], met[k][sel], label=k, **st)

        axR.set_yscale('log')
        axR.set_xlim(1.0, 1.5)
        ymax = np.nanmax(cool_tot[sel])
        axR.set_ylim(ymax * 1e-5, ymax * 3)
        axR.set_xlabel(r'$r/R_p$')
        axR.set_ylabel(r'cooling rate $\mathrm{[erg\,cm^{-3}\,s^{-1}]}$')
        axR.set_title(r'Per-channel cooling vs radius (Huang Fig.~10)')
        axR.legend(loc='upper right', frameon=False, ncol=2, fontsize=9)
        axR.grid(True, which='both', alpha=0.25)

        # report Fe I share at the base
        ib = np.argmax(sel)  # first in-window index
        print('# Fe I cooling share vs radius (innermost in-window cells):')
        print('# r/Rp     T[K]      cool_tot     FeI         FeI/tot   '
              'FeII/tot')
        for j in np.where(sel)[0][:8]:
            tot = cool_tot[j]
            fei = met.get('FeI', np.zeros_like(r))[j]
            feii = met.get('FeII', np.zeros_like(r))[j]
            print('%7.4f  %8.1f  %.3e  %.3e  %7.4f  %7.4f' %
                  (r[j], d[j, 1], tot, fei, fei / tot, feii / tot))

    fig.tight_layout()
    out = os.path.join(HERE, 'FeI_cooling_validation.png')
    fig.savefig(out, dpi=130)
    print('wrote', out)


if __name__ == '__main__':
    main()
