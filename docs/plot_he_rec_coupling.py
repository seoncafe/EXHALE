#!/usr/bin/env python3
"""Comparison figures for the He recombination radiation -> H I ionization
coupling (input key He_rec_coupling, flag use_he_rec_coupling; Draine 2011
y/z on-the-spot treatment).

Run protocol (2026-07-17, figures docs/figures/he_rec_coupling_{atomic,heitr}.pdf):
  - HD209458b/input.inp + metals.inp, converged HD209458b/output profiles
    loaded as IC (Load IC? True; copy output/*.txt to *_IC.txt first);
  - four runs: {atomic, He23S} x {off, on He_rec_coupling}, each relaxed
    with EXHALE_MAXSTEPS=15000 from the same IC (differential comparison;
    the IC predates the 2026-07-17 rate update, so both branches share the
    same partial re-relaxation and only the on/off difference is meaningful);
  - point `base` below at the directory holding the four run dirs
    (A_off, A_on, T_off, T_on), then run from docs/figures/.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

base = './he_rec_coupling_runs'   # directory holding A_off, A_on, T_off, T_on

plt.rcParams.update({'text.usetex': True, 'font.size': 11, 'axes.labelsize': 12,
                     'legend.fontsize': 9.5, 'figure.dpi': 110})
C_OFF = '#000000'; C_ON = '#D55E00'; C_B = '#0072B2'   # Okabe-Ito

def load(tag):
    ion = np.loadtxt(f'{base}/{tag}/output/Ion_species.txt')
    hyd = np.loadtxt(f'{base}/{tag}/output/Hydro_ioniz.txt')
    return ion, hyd

def hefrac(ion):
    nHe = ion[:, 3] + ion[:, 4] + ion[:, 5] + ion[:, 6]
    return ion[:, 3]/nHe, ion[:, 4]/nHe

for mode, tags, fig_name in [('atomic', ('A_off', 'A_on'), 'he_rec_coupling_atomic'),
                             ('HeITR', ('T_off', 'T_on'), 'he_rec_coupling_heitr')]:
    (i0, h0), (i1, h1) = load(tags[0]), load(tags[1])
    r = i0[:, 0]; m = (r >= 1.0) & (r <= 8.0)
    f0I, f0II = hefrac(i0); f1I, f1II = hefrac(i1)
    fig, ax = plt.subplots(2, 2, figsize=(9.2, 7.0), constrained_layout=True)
    a = ax[0, 0]
    a.plot(r[m], f0I[m],  C_OFF, ls='-',  lw=1.6, label=r'He\,I, off')
    a.plot(r[m], f1I[m],  C_ON,  ls='-',  lw=1.6, label=r'He\,I, on')
    a.plot(r[m], f0II[m], C_OFF, ls='--', lw=1.6, label=r'He\,II, off')
    a.plot(r[m], f1II[m], C_ON,  ls='--', lw=1.6, label=r'He\,II, on')
    a.set_yscale('log'); a.set_ylim(1e-3, 1.3)
    a.set_xlabel(r'$r/R_p$'); a.set_ylabel(r'$n/n_{\rm He}$')
    a.legend(ncol=2, frameon=False); a.set_title(r'He ionization fractions')
    a = ax[0, 1]
    a.plot(r[m], (f1I/f0I)[m],   C_ON, lw=1.6, label=r'He\,I')
    a.plot(r[m], (f1II/f0II)[m], C_B,  lw=1.6, ls='--', label=r'He\,II')
    if mode == 'HeITR':
        a.plot(r[m], (i1[:, 6]/np.maximum(i0[:, 6], 1e-30))[m], C_OFF, lw=1.6,
               ls=':', label=r'He\,I\,$2^3S$')
    a.axhline(1.0, color='0.6', lw=0.8)
    a.set_xlabel(r'$r/R_p$'); a.set_ylabel(r'on\,/\,off')
    a.legend(frameon=False); a.set_title(r'Coupling on\,/\,off ratio')
    a = ax[1, 0]
    a.plot(r[m], h0[m, 4], C_OFF, lw=1.6, label='off')
    a.plot(r[m], h1[m, 4], C_ON,  lw=1.6, label='on')
    a.set_xlabel(r'$r/R_p$'); a.set_ylabel(r'$T$ [K]')
    a.legend(frameon=False); a.set_title(r'Temperature')
    a = ax[1, 1]
    if mode == 'HeITR':
        a.plot(r[m], i0[m, 6], C_OFF, lw=1.6, label='off')
        a.plot(r[m], i1[m, 6], C_ON,  lw=1.6, label='on')
        a.set_yscale('log'); a.set_xlabel(r'$r/R_p$')
        a.set_ylabel(r'$n(\mathrm{He\,I}\,2^3S)$ [cm$^{-3}$]')
        a.legend(frameon=False); a.set_title(r'Metastable $2^3S$ density')
    else:
        a.plot(r[m], h0[m, 5], C_OFF, lw=1.6, label='off')
        a.plot(r[m], h1[m, 5], C_ON,  lw=1.6, label='on')
        a.set_yscale('log'); a.set_xlabel(r'$r/R_p$')
        a.set_ylabel(r'heating [erg\,cm$^{-3}$\,s$^{-1}$]')
        a.legend(frameon=False); a.set_title(r'Photoheating rate')
    fig.suptitle((r'HD\,209458\,b, %s mode: He recombination'
                  % ('atomic' if mode == 'atomic' else r'He\,I\,$2^3S$'))
                 + r' radiation $\rightarrow$ H\,I ionization coupling', fontsize=12)
    fig.savefig(f'{fig_name}.pdf')
    print('saved', fig_name)
