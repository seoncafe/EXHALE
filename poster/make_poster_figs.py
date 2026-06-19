#!/usr/bin/env python3
"""Generate the WASP-52b poster figures as vector PDFs (negative/transit-depth
convention).  Run from the poster/ directory:  python make_poster_figs.py
Outputs to poster/figs/*.pdf."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'text.usetex': False, 'font.size': 13,
                            'font.family': 'serif',
                            'font.serif': ['Liberation Serif', 'DejaVu Serif'],
                            'mathtext.fontset': 'stix',
                            'axes.titlesize': 14, 'axes.labelsize': 13,
                            'legend.fontsize': 10.5, 'lines.linewidth': 2.2})
import matplotlib.pyplot as plt

ROOT = os.path.abspath('..')                       # EXHALE/
sys.path.insert(0, os.path.join(ROOT, 'examples'))
import exhale_io as aio
import tpm_halpha_lart2d as h2d
import exhale_to_lart as e2l

W = os.path.join(ROOT, 'WASP-52b')
OBSDIR = '/home/kiseon/Exoplanetary_Atmosphere/WASP-52b'
FIG = 'figs'
os.makedirs(FIG, exist_ok=True)

he_obs = np.loadtxt(os.path.join(OBSDIR, 'wasp52_He10830_trans_spec.dat'))   # lam_vac, exc%, sig
ha_obs = np.loadtxt(os.path.join(OBSDIR, 'WASP52b_Halpha_transpec_0p10_AA_binned.txt'))


def vac_to_air(lv):
    s = 1.0e4 / lv
    return lv / (1.0 + 8.34254e-5 + 0.02406147/(130.0 - s**2) + 0.00015998/(38.9 - s**2))


def he_depth(folder):
    """(lam[A], transit depth[%] negative) from tpm_He10830.txt (col3 = rot+instr T)."""
    d = np.loadtxt(os.path.join(W, folder, 'tpm_He10830.txt'))
    return d[:, 0], -(1.0 - d[:, 3]) * 100.0


# ---------- 1. atmospheric structure (T, v, n) ----------
def fig_profiles():
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    cases = [('fxuv0p25_he98_L1', 'He 98/2', 'C1', '-'),
             ('fxuv0p25_solar_L1', 'solar 92/8', 'C0', '--')]
    for fol, lab, c, s in cases:
        d = np.loadtxt(os.path.join(W, fol, 'output', 'Hydro_ioniz_adv.txt'))
        r, n, v, T = d[:, 0], d[:, 1], d[:, 2], d[:, 4]
        kw = dict(color=c, ls=s, label=lab)
        ax[0].plot(r, T, **kw); ax[1].plot(r, v/1e5, **kw); ax[2].semilogy(r, n, **kw)
    for a in ax:
        a.set_xlabel(r'$R/R_{\rm pl}$'); a.grid(alpha=0.3); a.set_xlim(1, 2.53)
    ax[0].set_ylabel('T [K]'); ax[1].set_ylabel('v [km/s]'); ax[2].set_ylabel(r'n [cm$^{-3}$]')
    ax[0].set_title('Temperature'); ax[1].set_title('Velocity'); ax[2].set_title('Total density')
    ax[0].legend()
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'profiles.pdf')); plt.close(fig)


# ---------- 2. He 10830: extended (10 Rp) over-predicts vs L1 reproduces ----------
def fig_he_compare():
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    lair = vac_to_air(he_obs[:, 0])
    ax.errorbar(lair, -he_obs[:, 1], yerr=he_obs[:, 2], fmt='o', ms=4, color='0.35',
                ecolor='0.7', capsize=2, label='obs (Kirk+2022)', zorder=1)
    for fol, lab, c in [('fxuv0p25_he98', r'10 $R_p$ (extended): over-predicts', 'C3'),
                        ('fxuv0p25_he98_L1', r'$L_1$ (Roche, physical): matches', 'C0')]:
        lam, dep = he_depth(fol)
        ax.plot(lam, dep, color=c, label=lab, zorder=3)
    ax.axhline(0, color='0.6', lw=0.8); ax.set_xlim(10828, 10832)
    ax.set_xlabel(r'wavelength [$\AA$] (air)')
    ax.set_ylabel('transit depth [%] (neg. = absorption)')
    ax.set_title(r'He I $\lambda$10830 (H/He = 98/2, 0.5 $F_0$)')
    ax.legend(loc='lower left')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'he_compare.pdf')); plt.close(fig)


# ---------- 3. H/He scan ----------
def fig_hescan():
    pts = [(0.086957, 'fxuv0p25_solar_L1'), (0.050000, 'hescan_heh005_L1'),
           (0.020408, 'fxuv0p25_he98_L1'), (0.010101, 'hescan_heh0010101_L1'),
           (0.005025, 'hescan_heh0005025_L1')]
    heh = np.array([h for h, _ in pts])
    dep = np.array([-(1 - np.loadtxt(os.path.join(W, c, 'tpm_He10830.txt'))[:, 3].min())*100 for _, c in pts])
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(heh, dep, 'o-', color='C1', ms=8, label=r'EXHALE ($L_1$, 0.5 $F_0$)')
    ax.axhline(-3.44, color='0.4', ls='--', label='observed He 10830')
    ax.axhspan(-3.84, -3.04, color='0.8', alpha=0.5, zorder=0)
    ax.set_xscale('log'); ax.set_xlabel('He/H number ratio')
    ax.set_ylabel('He 10830 transit depth [%]')
    ax.set_title('He 10830 vs composition')
    for h, lab in [(0.086957, 'solar 92/8'), (0.020408, '98/2')]:
        i = np.argmin(abs(heh-h)); ax.annotate(lab, (heh[i], dep[i]),
                                               textcoords='offset points', xytext=(6, -14), fontsize=10)
    ax.grid(alpha=0.3, which='both'); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'hescan.pdf')); plt.close(fig)


# ---------- 4. incident Ly-alpha profile ----------
def fig_lya():
    c = 299792.458; lam0 = 1215.6701
    ls, ps = e2l.double_gaussian_lya(74, 49, nwav=601, width_is_fwhm=False)
    lf, pf = e2l.double_gaussian_lya(74, 49, nwav=601, width_is_fwhm=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.plot((ls-lam0)/lam0*c, ps, 'C0-', label=r'width = $\sigma$ = 49 km/s')
    ax.plot((lf-lam0)/lam0*c, pf, 'C1--', label='width = FWHM = 49 km/s')
    for vp in (-74, 74): ax.axvline(vp, color='0.7', ls=':', lw=1)
    ax.set_xlim(-250, 250); ax.set_xlabel(r'velocity [km/s] from Ly$\alpha$')
    ax.set_ylabel('normalized flux'); ax.set_title(r'Incident stellar Ly$\alpha$ (double-Gaussian)')
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'lya_profile.pdf')); plt.close(fig)


# ---------- 5. LaRT-coupled H-alpha (3 cases) ----------
def fig_halpha():
    cases = [(r'10 $R_p$, $\sigma$=49', 'fxuv0p25_he98', 'LaRT_lya/WASP-52b_lya.h5', 'C0'),
             (r'$L_1$ (Roche), $\sigma$=49', 'fxuv0p25_he98_L1', 'LaRT_lya_L1/WASP-52b_L1_lya.h5', 'C1'),
             (r'$L_1$ (Roche), FWHM=49', 'fxuv0p25_he98_L1', 'LaRT_lya_L1_fwhm/WASP-52b_L1_fwhm_lya.h5', 'C2')]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.errorbar(ha_obs[:, 0], ha_obs[:, 1]*100, yerr=ha_obs[:, 2]*100, fmt='o', ms=4,
                color='0.35', ecolor='0.7', capsize=2, label='obs (Chen+2020)', zorder=1)
    for lab, rd, h5, c in cases:
        lam, Tl, info = h2d.compute_halpha(os.path.join(W, rd), os.path.join(W, h5))
        ax.plot(lam, -(1-Tl)*100, color=c, label='%s : %.2f%%' % (lab, -info['depth_pct']), zorder=3)
    ax.axhline(0, color='0.6', lw=0.8); ax.set_xlim(6561, 6564.6)
    ax.set_xlabel(r'wavelength [$\AA$] (air)')
    ax.set_ylabel('transit depth [%] (neg. = absorption)')
    ax.set_title(r'H$\alpha$: EXHALE + LaRT 2D (cylindrical Ly$\alpha$)')
    ax.legend(loc='lower left')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'halpha.pdf')); plt.close(fig)


# ---------- 6. H-alpha matched by scaling F_Lya ----------
def fig_halpha_matched():
    rd = os.path.join(W, 'fxuv0p25_he98_L1')
    h5s = os.path.join(W, 'LaRT_lya_L1/WASP-52b_L1_lya.h5')
    h5f = os.path.join(W, 'LaRT_lya_L1_fwhm/WASP-52b_L1_fwhm_lya.h5')
    l1, T1, i1 = h2d.compute_halpha(rd, h5s, lya_scale=1.0)
    lsm, Tsm, ism = h2d.compute_halpha(rd, h5s, lya_scale=0.35)
    lfm, Tfm, ifm = h2d.compute_halpha(rd, h5f, lya_scale=0.25)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.errorbar(ha_obs[:, 0], ha_obs[:, 1]*100, yerr=ha_obs[:, 2]*100, fmt='o', ms=4,
                color='0.35', ecolor='0.7', capsize=2, label='obs (Chen+2020)', zorder=1)
    ax.plot(l1, -(1-T1)*100, '--', color='0.6', label=r'$\sigma$=49, full $F_{\rm Ly\alpha}$ : %.2f%%' % -i1['depth_pct'])
    ax.plot(lsm, -(1-Tsm)*100, 'C3-', label=r'$\sigma$=49, $\times$0.35 : %.2f%%' % -ism['depth_pct'])
    ax.plot(lfm, -(1-Tfm)*100, 'C0-', label=r'FWHM=49, $\times$0.25 : %.2f%%' % -ifm['depth_pct'])
    ax.axhline(0, color='0.6', lw=0.8); ax.set_xlim(6561, 6564.6)
    ax.set_xlabel(r'wavelength [$\AA$] (air)')
    ax.set_ylabel('transit depth [%] (neg. = absorption)')
    ax.set_title(r'H$\alpha$: Ly$\alpha$ flux scaled to match the data')
    ax.legend(loc='lower left')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'halpha_matched.pdf')); plt.close(fig)


# ---------- 7. Level-diagram schematics: He 2^3S and H-alpha Lya pumping ----------
def fig_schematics():
    """Schematic (NOT to energy scale) Grotrian-style diagrams."""
    fig, (axH, axA) = plt.subplots(1, 2, figsize=(13, 5.4))
    YC, Y2, Y3 = 5.0, 2.8, 3.9       # continuum, lower-excited, upper-excited (schematic)

    def lev(ax, x0, x1, y, color='k', lw=3.5):
        ax.hlines(y, x0, x1, color=color, lw=lw)

    # ===== He I 10830 : singlet (left) | triplet (right) =====
    axH.set_title(r'He I $\lambda$10830 : He $2^3S$ metastable', fontsize=15, weight='bold')
    axH.hlines(YC, 0.3, 6.7, color='0.55', ls='--', lw=2)
    axH.text(6.7, YC + 0.22, r'He$^+$ + e$^-$  ($24.6$ eV)', ha='right', fontsize=11.5, color='0.4')
    lev(axH, 0.5, 2.1, 0.0)                                          # 1^1S ground (singlet)
    axH.text(1.3, -0.48, r'$1^1S$ (ground)', ha='center', fontsize=12.5)
    lev(axH, 4.0, 5.6, Y2, color='C0')                              # 2^3S metastable
    lev(axH, 4.0, 5.6, Y3, color='C0')                             # 2^3P
    axH.text(5.7, Y2, r'$2^3S$', ha='left', va='center', fontsize=13, color='C0')
    axH.text(5.7, Y3, r'$2^3P$', ha='left', va='center', fontsize=13, color='C0')
    axH.text(4.8, Y2 - 0.5, r'metastable ($\tau\!\sim\!2.2$ hr)', ha='center', fontsize=10, color='C0')
    axH.text(1.3, 1.7, 'singlet\n(para)', ha='center', fontsize=11, color='0.45')
    axH.text(4.8, 1.0, 'triplet\n(ortho)', ha='center', fontsize=11, color='0.45')
    axH.annotate('', xy=(4.2, Y2), xytext=(3.3, YC),               # recomb + cascade
                 arrowprops=dict(arrowstyle='-|>', color='0.5', lw=2.2, connectionstyle='arc3,rad=0.3'))
    axH.text(3.05, 4.05, 'recomb.\n+ cascade', ha='right', fontsize=10.5, color='0.4')
    axH.annotate('', xy=(4.55, Y3), xytext=(4.55, Y2),            # He 10830 (2^3S -> 2^3P)
                 arrowprops=dict(arrowstyle='-|>', color='C3', lw=3.2))
    axH.text(4.42, (Y2 + Y3) / 2, r'He $\lambda$10830', ha='right', va='center',
             fontsize=12, color='C3', weight='bold')
    axH.annotate('', xy=(5.25, YC), xytext=(5.25, Y2),           # FUV photoion of 2^3S (sink)
                 arrowprops=dict(arrowstyle='-|>', color='C1', lw=2.6))
    axH.text(5.42, 4.55, 'FUV photoion.\n' r'$\lambda\!<\!2600\,\AA$',
             ha='left', va='center', fontsize=10, color='C1')
    axH.set_xlim(-0.3, 8.3); axH.set_ylim(-1.0, 6.1); axH.axis('off')

    # ===== H-alpha : Lya pumping of n=2 =====
    axA.set_title(r'H$\alpha$ : Ly$\alpha$ pumping of H ($n{=}2$)', fontsize=15, weight='bold')
    axA.hlines(YC, 0.3, 6.7, color='0.55', ls='--', lw=2)
    axA.text(6.7, YC + 0.22, r'H$^+$ + e$^-$  ($13.6$ eV)', ha='right', fontsize=11.5, color='0.4')
    lev(axA, 1.4, 3.4, 0.0)                                          # n=1
    axA.text(1.3, 0.0, r'$n{=}1$ ($1s$)', ha='right', va='center', fontsize=13)
    lev(axA, 1.4, 3.4, Y2, color='C0')                             # n=2
    axA.text(3.5, Y2, r'$n{=}2$ ($2s,2p$)', ha='left', va='center', fontsize=13, color='C0')
    lev(axA, 1.4, 3.4, Y3)                                          # n=3
    axA.text(3.5, Y3, r'$n{=}3$', ha='left', va='center', fontsize=13)
    axA.annotate('', xy=(2.0, Y2), xytext=(2.0, 0.0),            # Lya pump 1s <-> 2p
                 arrowprops=dict(arrowstyle='<|-|>', color='C0', lw=3.2))
    axA.text(1.85, Y2 / 2, r'Ly$\alpha$ $1216\,\AA$' '\n' r'pump $B_{1s2p}\bar{J}$',
             ha='right', va='center', fontsize=11.5, color='C0', weight='bold')
    axA.annotate('', xy=(2.8, Y3), xytext=(2.8, Y2),            # H-alpha n=2 -> n=3
                 arrowprops=dict(arrowstyle='-|>', color='C3', lw=3.2))
    axA.text(2.95, (Y2 + Y3) / 2, r'H$\alpha$ $6563\,\AA$', ha='left', va='center',
             fontsize=12, color='C3', weight='bold')
    axA.annotate('', xy=(3.2, Y2), xytext=(4.2, YC),           # recomb cascade -> n=2
                 arrowprops=dict(arrowstyle='-|>', color='0.5', lw=2.2, connectionstyle='arc3,rad=-0.3'))
    axA.text(4.4, 4.05, 'recomb.\ncascade', ha='left', fontsize=10.5, color='0.4')
    axA.set_xlim(-0.3, 8.3); axA.set_ylim(-1.0, 6.1); axA.axis('off')

    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'schematics.pdf')); plt.close(fig)


if __name__ == '__main__':
    fig_profiles();        print('profiles.pdf')
    fig_he_compare();      print('he_compare.pdf')
    fig_hescan();          print('hescan.pdf')
    fig_lya();             print('lya_profile.pdf')
    fig_halpha();          print('halpha.pdf')
    fig_halpha_matched();  print('halpha_matched.pdf')
    fig_schematics();      print('schematics.pdf')
    print('all poster figures written to figs/*.pdf')
