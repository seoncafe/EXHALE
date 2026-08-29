#!/usr/bin/env python3
"""One row per arm of the K_zz x reservoir-He/H scan under the profile
boundary condition.

For every `kzz*/heh*/` case this reads

  closure_history.txt   the converged iterate k, its measured elemental
                        fluxes, the window spreads, the He/H the profile
                        hands over at the microbar match, log10 Mdot,
                        and the EXHALE info flag;
  k<NN>/run.log         the JFNK residual norm the wind stopped at;
  k<NN>/EXHALE_resolved.out
                        the cold trap and tropopause the climate solve gave;
  k<NN>/tpm_He10830*.txt
                        the He I 10830 line, with the red-pair equivalent
                        width taken on the same vacuum window
                        `make_memo_figures.py` and
                        `flux_closure/closure_heh_table.py` use;
  k<NN>/output/Hydro_ioniz.txt, Ion_species_adv.txt, element_flux_profile.txt
                        the wind: where T drops, how much of the metastable
                        column sits beyond 10 R_p, and the He/H the wind
                        itself carries at 2 and 10 R_p.

Usage:  python3 measure.py [kzz1e8/heh11p1 ...]   (default: every case found)
"""
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'flux_closure'))

AIR = 10832.057/10829.09114          # vacuum -> air for the red-pair window
EW_LO, EW_HI = 10832.60, 10834.20    # vacuum [A], the memo's red-pair window
EW_OBS, EW_ERR = 1.108, 0.030        # Cherubim et al. (2026), [% A]


def red_ew(path):
    if not os.path.isfile(path):
        return float('nan')
    s = np.loadtxt(path)
    lam = s[:, 0]*AIR
    dep = (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100.0
    m = (lam >= EW_LO) & (lam <= EW_HI)
    return float(np.trapz(dep[m], lam[m]))


def metrics(path):
    out = {}
    if not os.path.isfile(path):
        return out
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        k, v = line.split()
        out[k] = float(v)
    return out


def resid_norm(kdir):
    r = float('nan')
    p = os.path.join(kdir, 'run.log')
    if os.path.isfile(p):
        for line in open(p, errors='replace'):
            m = re.search(r'\(JFNK\) done info=\s*(\S+)\s+\|\|R\|\|=\s*(\S+)',
                          line)
            if m:
                r = float(m.group(2))
    return r


def resolved_notes(kdir):
    p = os.path.join(kdir, 'EXHALE_resolved.out')
    txt = ''
    for line in open(p, errors='replace'):
        if line.startswith('lower_profile_notes'):
            txt = line
    def grab(pat, default=float('nan')):
        m = re.search(pat, txt)
        return float(m.group(1)) if m else default
    return dict(
        OH=grab(r'O/H gas ([0-9.eE+-]+)'),
        CH=grab(r'C/H gas ([0-9.eE+-]+)'),
        NH=grab(r'N/H gas ([0-9.eE+-]+)'),
        p_trop=grab(r'tropopause ([0-9.eE+-]+) bar'),
        T_trop=grab(r'tropopause [0-9.eE+-]+ bar at ([0-9.eE+-]+) K'),
        T_deep=grab(r'deep boundary 20 bar at ([0-9.eE+-]+) K'))


def kzz_on_grid(kdir):
    """The range of K_zz the wind actually saw, as init echoed it."""
    lo = hi = float('nan')
    p = os.path.join(kdir, 'run.log')
    for line in open(p, errors='replace'):
        m = re.search(r'K_zz\(r\) on the grid:\s+(\S+) to\s+(\S+)', line)
        if m:
            lo, hi = float(m.group(1)), float(m.group(2))
    return lo, hi


def wind(kdir):
    """r_drop, the metastable column beyond 10 R_p, and the wind He/H."""
    o = os.path.join(kdir, 'output')
    h = np.loadtxt(os.path.join(o, 'Hydro_ioniz.txt'))
    r, T = h[:, 0], h[:, 4]
    T12 = float(np.interp(12.0, r, T))
    r_drop = float('nan')
    below = np.where((r > 2.0) & (T < 0.5*T12))[0]
    if below.size:
        r_drop = float(r[below[0]])
    a = np.loadtxt(os.path.join(o, 'Ion_species_adv.txt'))
    rr, n23 = a[:, 0], a[:, 6]
    Rp_cm = 1.0                                   # column ratio: units cancel
    tot = np.trapz(n23, rr*Rp_cm)
    out = np.trapz(n23[rr >= 10.0], rr[rr >= 10.0]*Rp_cm)
    f_out = float(out/tot) if tot > 0 else float('nan')
    e = np.loadtxt(os.path.join(o, 'element_flux_profile.txt'),
                   usecols=(1, 2, 10))
    X = e[:, 1]
    heh = (X/4.002602)/np.maximum(1.0 - X, 1e-300)*1.00794
    # the molecular coefficient the eddy term competes with, at the first
    # face (the base) and where the metastable line is made
    d_base = float(e[0, 2])
    d_2 = float(np.interp(2.0, e[:, 0], e[:, 2]))
    return (r_drop, f_out, float(np.interp(2.0, e[:, 0], heh)),
            float(np.interp(10.0, e[:, 0], heh)), float(T.max()),
            d_base, d_2)


def last_row(case):
    rows = [l.split() for l in open(os.path.join(case, 'closure_history.txt'))
            if l.strip() and not l.startswith('#')]
    return rows[-1]


def case_kzz_heh(case):
    cfg = open(os.path.join(case, 'closure.json')).read()
    kz = re.search(r'"--kzz-const",\s*"([^"]+)"', cfg).group(1)
    hh = re.search(r'"--abundances",\s*"He=([^"]+)"', cfg).group(1)
    return float(kz), float(hh)


def row(case):
    r = last_row(case)
    k = int(r[0])
    kdir = os.path.join(case, 'k%02d' % k)
    kz, resv = case_kzz_heh(case)
    mt = metrics(os.path.join(kdir, 'tpm_He10830_metrics.txt'))
    ew = red_ew(os.path.join(kdir, 'tpm_He10830.txt'))
    nt = resolved_notes(kdir)
    klo, khi = kzz_on_grid(kdir)
    r_drop, f_out, heh2, heh10, Tmax, d_base, d_2 = wind(kdir)
    return dict(case=os.path.relpath(case, HERE), kzz=kz, resv=resv, k=k,
                info=r[14], eps=max(float(r[5]), float(r[6])),
                spread_H=float(r[9]), spread_He=float(r[10]),
                heh_match=float(r[11]), lgMdot=float(r[12]),
                resid=resid_norm(kdir), kzz_lo=klo, kzz_hi=khi,
                red=mt.get('red_depth', float('nan')),
                fwhm=mt.get('fwhm_A', float('nan')), ew=ew, ratio=ew/EW_OBS,
                r_drop=r_drop, f_out=f_out, heh2=heh2, heh10=heh10, Tmax=Tmax,
                d_base=d_base, d_2=d_2, kd_base=kz/d_base, **nt)


HDR = ('%-18s %8s %8s %2s %4s %9s %8s %8s %9s %8s %8s %7s %7s %8s %8s %7s'
       ' %8s %8s %8s %9s %9s %9s %6s %6s %9s %8s'
       % ('case', 'K_zz', 'reserv', 'k', 'info', '||R||', 'sprd_H', 'sprd_He',
          'HeH@mtch', 'mtch/res', 'lgMdot', 'red[%]', 'FWHM', 'EW[%A]',
          'EW/obs', 'r_drop', 'f_23S>10', 'HeH@2Rp', 'HeH@10Rp', 'O/H_ct',
          'C/H_ct', 'N/H_ct', 'p_trp', 'T_trp', 'D_base', 'Kzz/Dbse'))
FMT = ('%-18s %8.1E %8.4f %2d %4s %9.3E %8.4f %8.4f %9.5f %8.5f %8.3f %7.3f'
       ' %7.4f %8.4f %8.4f %7.3f %8.4f %8.4f %8.4f %9.4E %9.4E %9.4E %6.3f'
       ' %6.1f %9.3E %8.4f')


def main(cases):
    if not cases:
        cases = sorted(glob.glob(os.path.join(HERE, 'kzz*', 'heh*')))
    print(HDR)
    for c in cases:
        c = c if os.path.isabs(c) else os.path.join(HERE, c)
        if not os.path.isfile(os.path.join(c, 'closure_history.txt')):
            continue
        try:
            d = row(c)
        except Exception as exc:                       # report, do not hide
            print('%-18s  FAILED: %s' % (os.path.relpath(c, HERE), exc))
            continue
        print(FMT % (d['case'], d['kzz'], d['resv'], d['k'], d['info'],
                     d['resid'], d['spread_H'], d['spread_He'],
                     d['heh_match'], d['heh_match']/d['resv'], d['lgMdot'],
                     d['red'], d['fwhm'], d['ew'], d['ratio'], d['r_drop'],
                     d['f_out'], d['heh2'], d['heh10'], d['OH'], d['CH'],
                     d['NH'], d['p_trop'], d['T_trop'], d['d_base'],
                     d['kd_base']))


if __name__ == '__main__':
    main(sys.argv[1:])
