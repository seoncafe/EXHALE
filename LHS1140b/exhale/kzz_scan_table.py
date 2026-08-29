#!/usr/bin/env python3
"""Read the LHS 1140 b He_Kzz scan and print the decision table.

For each case directory the script reports, all measured from the files the
run left behind and nothing quoted from a log line that a later pass could
have superseded:

  * the solver outcome (info, ||R||, number of composition outer passes and
    the drift each ended on),
  * the elemental ratio (He/H)/HeH on a fixed radius list,
  * the homopause radius, where the molecular coefficient D_eff of the
    diffusion operator equals the eddy coefficient K_zz of the case,
  * log10 Mdot from the post-processing pass,
  * the He 10830 red-pair depth and equivalent width.
"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RLIST = [1.02, 1.05, 1.10, 1.20, 2.0, 5.0, 10.0, 20.0]
HEH = 0.55
AIR_TO_VAC = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20


def _cols(path):
    """Column names from the '# columns ...' schema header."""
    with open(path) as fh:
        for line in fh:
            if line.startswith('# columns'):
                return line.split()[2:]
    raise ValueError('no column header in ' + path)


def element_ratio(d):
    """(He/H)/HeH against radius, from the elemental sums of Ion_species."""
    p = os.path.join(d, 'output', 'Ion_species_adv.txt')
    if not os.path.exists(p):
        p = os.path.join(d, 'output', 'Ion_species.txt')
    names = _cols(p)
    a = np.loadtxt(p)
    idx = {n: i for i, n in enumerate(names)}
    r = a[:, idx['r[Rp]']]
    # He I already contains the 2^3S metastable (a level of He I, not a
    # separate species -- species_table.f90); adding HeITR here would count
    # the metastable helium nuclei twice.
    nHe = sum(a[:, idx[s]] for s in ('HeI', 'HeII', 'HeIII') if s in idx)
    nH = sum(a[:, idx[s]] for s in ('HI', 'HII') if s in idx)
    return r, (nHe / np.maximum(nH, 1e-99)) / HEH


def _face_flux_file(d):
    """Path and D_eff column of the element face-flux table of a run.

    The file is ./output/element_flux_profile.txt since 2026-08-27 and was
    ./diffusion_faceflux.txt in the run root before; runs made under either
    binary are read. The D_eff column is located by name from the
    "# columns:" schema line, so an added column cannot silently shift it
    (F_H was added between F_He and Mdot_face).
    """
    for rel in ('output/element_flux_profile.txt', 'diffusion_faceflux.txt'):
        p = os.path.join(d, rel)
        if os.path.isfile(p):
            break
    else:
        return None, None
    icol = 8
    with open(p) as f:
        for line in f:
            if not line.startswith('#'):
                break
            if 'columns:' in line:
                toks = line.split('columns:', 1)[1].split()
                if 'D_eff[cm2/s]' in toks:
                    icol = toks.index('D_eff[cm2/s]')
                break
    return p, icol


def homopause(d, kzz):
    """Radius where the molecular coefficient equals K_zz.

    D_eff is the stage-resolved binary coefficient the operator used; the
    eddy term is added to it, so the crossing is the homopause of the run.
    """
    p, icol = _face_flux_file(d)
    if kzz <= 0 or p is None:
        return None, None
    a = np.loadtxt(p, usecols=(1, icol))
    r, D = a[:, 0], a[:, 1]
    s = D - kzz
    for j in range(len(r) - 1):
        if s[j] * s[j + 1] < 0.0:
            f = s[j] / (s[j] - s[j + 1])
            return r[j] + f * (r[j + 1] - r[j]), (D[0], D[-1])
    # No crossing: D_eff rises with radius, so the homopause is above the top
    # of the domain when the eddy term still dominates at the last face, and
    # below the base when the molecular one already dominates at the first.
    return (float('inf') if s[0] < 0 else float('-inf')), (D[0], D[-1])


def solver(d):
    """info / ||R|| / outer-pass record, from the last solver log present."""
    logs = sorted([f for f in os.listdir(d)
                   if re.fullmatch(r'run(_cont)?\d*\.log', f)],
                  key=lambda f: os.path.getmtime(os.path.join(d, f)))
    if not logs:
        return {}
    txt = open(os.path.join(d, logs[-1])).read()
    info = re.findall(r'done info=(\d+)\s+\|\|R\|\|=\s*([0-9.E+-]+)', txt)
    drift = re.findall(r'outer pass\s+(\d+):.*?drift =\s*([0-9.E+-]+)', txt)
    return dict(log=logs[-1],
                info=info[-1][0] if info else '-',
                resid=info[-1][1] if info else '-',
                npass=len(drift),
                drift=drift[-1][1] if drift else '-')


def mdot(d):
    for f in ('pp.log', 'run.log'):
        p = os.path.join(d, f)
        if not os.path.exists(p):
            continue
        m = re.findall(r'steady-state Mdot\s*=\s*([0-9.+-]+)', open(p).read())
        if m:
            return float(m[-1])
    return float('nan')


def he_line(d):
    """Red-pair depth [%] and equivalent width [% A] of the He 10830 line."""
    p = os.path.join(d, 'tpm_He10830_metrics.txt')
    depth = float('nan')
    if os.path.exists(p):
        for line in open(p):
            if line.startswith('red_depth'):
                depth = float(line.split()[1])
    ew = float('nan')
    q = os.path.join(d, 'tpm_He10830.txt')
    if os.path.exists(q):
        # Same red-pair equivalent width as ../make_memo_figures.py, so that
        # the scan rows are directly comparable to the He/H scan and to the
        # measured 1.108 +/- 0.030 %A: the instrument-convolved curve
        # (column index 2) normalized to its own continuum, in percent, on
        # the vacuum window 10832.60-10834.20 A.
        a = np.loadtxt(q)
        lam = a[:, 0]*AIR_TO_VAC
        dep = (a[:, 2].max() - a[:, 2])/a[:, 2].max()*100.0
        m = (lam >= EW_LO) & (lam <= EW_HI)
        ew = np.trapz(dep[m], lam[m])
    return depth, ew


def main(cases):
    ratio_rows, main_rows = [], []
    for tag, kzz in cases:
        d = os.path.join(HERE, tag)
        if not os.path.isdir(d):
            print('missing: ' + tag)
            continue
        s = solver(d)
        r, x = element_ratio(d)
        rh, Dends = homopause(d, kzz)
        dep, ew = he_line(d)
        main_rows.append((tag, kzz, s, mdot(d), rh, dep, ew, Dends))
        ratio_rows.append((tag, [float(np.interp(rr, r, x)) for rr in RLIST]))

    print('\n## solver, homopause, Mdot, He 10830\n')
    print('| case | K_zz [cm2/s] | info | ||R|| | outer passes | last drift '
          '| homopause r [R_p] | log10 Mdot | He 10830 red depth [%] '
          '| red EW [% A] |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    for tag, kzz, s, md, rh, dep, ew, De in main_rows:
        if rh is None:
            hp = 'n/a (K_zz = 0)'
        elif rh == float('-inf'):
            hp = '< base (D_eff > K_zz everywhere)'
        elif rh == float('inf'):
            hp = '> 30 (D_eff < K_zz everywhere)'
        else:
            hp = '%.4f' % rh
        print('| %s | %.0e | %s | %s | %d | %s | %s | %.4f | %.4e | %.4e |'
              % (tag, kzz, s.get('info', '-'), s.get('resid', '-'),
                 s.get('npass', 0), s.get('drift', '-'), hp, md, dep, ew))

    print('\n## (He/H)/HeH\n')
    print('| case | ' + ' | '.join('%g' % rr for rr in RLIST) + ' |')
    print('|---' * (len(RLIST) + 1) + '|')
    for tag, vals in ratio_rows:
        print('| %s | ' % tag + ' | '.join('%.4f' % v for v in vals) + ' |')

    print('\n## D_eff at the base and at the top face [cm2/s]\n')
    for tag, kzz, s, md, rh, dep, ew, De in main_rows:
        if De:
            print('%-24s D_eff(first face) = %.4e   D_eff(last face) = %.4e'
                  % (tag, De[0], De[1]))


if __name__ == '__main__':
    CASES = [('heh0p55_diff_ctrl', 0.0),
             ('heh0p55_diff_kzz1e6', 1.0e6),
             ('heh0p55_diff_kzz1e7', 1.0e7),
             ('heh0p55_diff_kzz1e8', 1.0e8),
             ('heh0p55_diff_kzz1e9', 1.0e9),
             ('heh0p55_diff_kzz1e10', 1.0e10),
             ('heh0p55_diff_kzz1e11', 1.0e11)]
    main(CASES if len(sys.argv) == 1
         else [c for c in CASES if c[0] in sys.argv[1:]])
