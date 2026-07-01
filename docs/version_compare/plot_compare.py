#!/usr/bin/env python3
"""Compare EXHALE_v1.0 (no diffusion) vs current EXHALE (He/H + metal diffusion)
for HD 209458b and WASP-121b: vector-PDF figures for the doc."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = False
plt.rcParams['font.size'] = 10
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

here = '/nfs/mocafe/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/docs/version_compare/'

# planet: (v1 tag, v2 tag, He/H input, rmax, list of (name,cols) metals present)
metcol = dict(C=[7, 8, 9], N=[13, 14, 15], O=[10, 11, 12], Mg=[16, 17, 18],
              Ca=[22, 23, 24], Na=[25, 26], Fe=[31, 32, 33])
planets = {
    'HD 209458b': dict(v1='v1_nodiff', v2='v2_diff', HeH=0.0833333, rmax=5.0,
                       metals=['C', 'N', 'O'], tag='hd209'),
    'WASP-121b':  dict(v1='wasp_v1_nodiff', v2='wasp_v2_diff', HeH=0.0833333,
                       rmax=1.6, metals=['O', 'Mg', 'Fe'], tag='wasp',
                       lxlim=(10827.5, 10832.5)),
}

def load(tag):
    I = np.loadtxt(here + tag + '/output/Ion_species.txt')
    H = np.loadtxt(here + tag + '/output/Hydro_ioniz.txt')
    L = np.loadtxt(here + tag + '/tpm_tpm_He10830.txt')
    nH = I[:, 1] + I[:, 2]
    return dict(r=I[:, 0], I=I, nH=nH,
                heh=(I[:, 3] + I[:, 4] + I[:, 5]) / nH, heTR=I[:, 6],
                xHII=I[:, 2] / nH, T=H[:, 4], v=H[:, 2] / 1e5,
                lam=L[:, 0], Tobs=L[:, 3])

for pname, p in planets.items():
    d1, d2 = load(p['v1']), load(p['v2'])
    He0 = p['HeH']

    # ---- profiles figure ----
    fig, ax = plt.subplots(2, 3, figsize=(11, 6.4))
    for d, lab, c, ls in [(d1, 'v1.0 (no diffusion)', 'tab:blue', '--'),
                          (d2, 'current (diffusion)', 'tab:red', '-')]:
        ax[0, 0].plot(d['r'], d['heh'] / He0, c, ls=ls, label=lab)
        ax[0, 1].plot(d['r'], d['heTR'], c, ls=ls, label=lab)
        ax[0, 2].plot(d['r'], d['T'], c, ls=ls, label=lab)
        ax[1, 0].plot(d['r'], d['v'], c, ls=ls, label=lab)
        ax[1, 1].plot(d['r'], d['xHII'], c, ls=ls, label=lab)
    ax[0, 0].axhline(1.0, color='gray', lw=0.7)
    ax[0, 0].set(ylabel='(He/H)/(He/H)$_{base}$', title='He/H element ratio', ylim=(0, 1.15))
    ax[0, 1].set(ylabel='n(He 2$^3$S) [cm$^{-3}$]', title='metastable He 2$^3$S')
    ax[0, 1].set_yscale('log')
    ax[0, 2].set(ylabel='T [K]', title='temperature')
    ax[1, 0].set(ylabel='v [km s$^{-1}$]', title='outflow velocity')
    ax[1, 1].set(ylabel='n(H II)/n(H)', title='H ionization fraction')

    # metals panel (v2): element/H normalized to base, mass ordering
    r = d2['r']; rb = np.argmin(np.abs(r - 1.0))
    cols = {'He': 'tab:red', 'C': 'tab:green', 'N': 'tab:orange',
            'O': 'tab:purple', 'Mg': 'tab:brown', 'Fe': 'k'}
    y = d2['heh']; ax[1, 2].plot(r, y / y[rb], cols['He'], label='He (4)')
    masslab = dict(C='C (12)', N='N (14)', O='O (16)', Mg='Mg (24)', Fe='Fe (56)')
    for name in p['metals']:
        e = d2['I'][:, metcol[name]].sum(1) / d2['nH']
        base = e[rb] if e[rb] > 0 else 1e-30
        ax[1, 2].plot(r, np.clip(e / base, 0, 1.2), cols.get(name, 'gray'),
                      label=masslab.get(name, name))
    ax[1, 2].set(ylabel='element/H (norm. to base)', ylim=(0, 1.15),
                 title='element diffusion (v2)')
    ax[1, 2].legend(fontsize=7, ncol=2)

    for a in ax.flat:
        a.set_xlabel('r [R$_p$]'); a.set_xlim(1, p['rmax'])
    ax[0, 0].legend(fontsize=8)
    fig.suptitle('%s: EXHALE v1.0 (no diffusion) vs current (He/H + metal diffusion)'
                 % pname, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(here + 'fig_%s_profiles.pdf' % p['tag'])

    # ---- 10830 line figure ----
    fig2, ax2 = plt.subplots(1, 1, figsize=(6, 4.2))
    for d, lab, c, ls in [(d1, 'v1.0 (no diffusion)', 'tab:blue', '--'),
                          (d2, 'current (diffusion)', 'tab:red', '-')]:
        dep = (1 - d['Tobs']) * 100
        ax2.plot(d['lam'], dep, c, ls=ls, label=lab + '  (peak %.1f%%)' % dep.max())
    ax2.axvline(10830.3, color='gray', lw=0.7, ls=':')
    ax2.set(xlabel=r'wavelength [$\AA$]', ylabel='absorption 1 - T [%]',
            title='%s: He I 10830 transit line' % pname,
            xlim=p.get('lxlim', (10827, 10834)))
    ax2.legend(fontsize=9)
    fig2.tight_layout()
    fig2.savefig(here + 'fig_%s_He10830.pdf' % p['tag'])
    print('saved fig_%s_profiles.pdf, fig_%s_He10830.pdf' % (p['tag'], p['tag']))
