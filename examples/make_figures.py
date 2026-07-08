#!/usr/bin/env python3
"""Generate the figures used in the EXHALE user manual (and the example
notebook). Loads converged runs with exhale_io and writes PNGs to docs/figures/.

Run from the examples/ directory:  python3 make_figures.py
Each figure is wrapped in try/except so a missing run is skipped, not fatal.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import exhale_io as aio

ROOT = os.path.join(os.path.dirname(__file__), '..')
W = os.path.join(ROOT, 'WASP-121b')
FIG = os.path.join(ROOT, 'docs', 'figures')
os.makedirs(FIG, exist_ok=True)


def save(fig, name):
    # PNG (notebooks/READMEs) + vector PDF (the LaTeX docs include the PDF)
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=130, bbox_inches='tight')
    fig.savefig(os.path.splitext(p)[0] + '.pdf', bbox_inches='tight')
    plt.close(fig)
    print('  wrote', os.path.relpath(p, ROOT), '(+.pdf)')




# --- Fig 1: tutorial overview (the "what a run produces" figure) ----------- #
def load_tut(subdir):
    """Load a tutorial run, preferring the _adv files but falling back to the
    equilibrium snapshot if the slow spherical run has not fully converged."""
    d = os.path.join(ROOT, 'examples', subdir)
    adv = os.path.exists(os.path.join(d, 'output', 'Hydro_ioniz_adv.txt'))
    return aio.load_run(os.path.join(d, 'output'),
                        os.path.join(d, 'input.inp'), adv=adv)


try:
    t = load_tut('tutorial')
    fig, ax = plt.subplots(2, 2, figsize=(9, 7))
    ax[0, 0].semilogy(t.r, t.n); ax[0, 0].set_ylabel(r'$n$ [cm$^{-3}$]')
    ax[0, 1].plot(t.r, t.v_kms); ax[0, 1].set_ylabel(r'$v$ [km s$^{-1}$]')
    ax[1, 0].plot(t.r, t.T); ax[1, 0].set_ylabel(r'$T$ [K]')
    xh = t.x_ion(['HI', 'HII'])
    ax[1, 1].plot(t.r, xh['HI'], label='H I')
    ax[1, 1].plot(t.r, xh['HII'], label='H II')
    ax[1, 1].set_ylabel('H ionization fraction'); ax[1, 1].legend()
    for a in ax.flat:
        a.set_xlabel(r'$r$ [$R_p$]'); a.grid(alpha=0.3)
    fig.suptitle('Tutorial run: generic hot Jupiter (spherical, H/He + metals)')
    save(fig, 'tutorial_overview.png')
except Exception as e:
    print('  [skip] tutorial_overview:', e)


# --- Fig 2: metals ON vs OFF (clean controlled tutorial pair) -------------- #
try:
    on = load_tut('tutorial')
    off = load_tut('tutorial_nometals')
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(off.r, off.T, '--', label='metals off')
    ax[0].plot(on.r, on.T, '-', label='metals on')
    ax[0].set_ylabel(r'$T$ [K]'); ax[0].legend()
    ax[0].set_title('Metal-line cooling lowers the peak temperature')
    xo = off.x_ion(['HI', 'HII']); xn = on.x_ion(['HI', 'HII'])
    ax[1].plot(off.r, xo['HII'], '--', label='H II (off)')
    ax[1].plot(on.r, xn['HII'], '-', label='H II (on)')
    ax[1].set_ylabel('H II fraction'); ax[1].legend()
    ax[1].set_title('while the H ionization is barely affected')
    for a in ax:
        a.set_xlabel(r'$r$ [$R_p$]'); a.grid(alpha=0.3)
    save(fig, 'metals_on_off.png')
except Exception as e:
    print('  [skip] metals_on_off:', e)


# --- Fig 3: cooling breakdown (dominant coolants vs radius) ----------------- #
try:
    c = aio.load_cooling(os.path.join(W, 'output', 'Cooling_breakdown.txt'))
    fig, axx = plt.subplots(figsize=(8, 5))
    axx.semilogy(c['r'], c['cool_total'], 'k-', lw=2, label='total')
    show = ['coll_exc', 'rec', 'FeII', 'MgII', 'CII', 'OII', 'CaII']
    for nm in show:
        if nm in c['chan']:
            y = np.clip(c['chan'][nm], 1e-30, None)
            axx.semilogy(c['r'], y, label=nm)
    axx.set_xlabel(r'$r$ [$R_p$]')
    axx.set_ylabel(r'cooling [erg cm$^{-3}$ s$^{-1}$]')
    axx.set_ylim(c['cool_total'].max() * 1e-4, c['cool_total'].max() * 2)
    axx.set_title('Radiative cooling by channel (WASP-121b, metals on)')
    axx.legend(ncol=2, fontsize=8); axx.grid(alpha=0.3)
    save(fig, 'cooling_breakdown.png')
except Exception as e:
    print('  [skip] cooling_breakdown:', e)


# --- Fig 4: spherical (Case A) vs Roche (Case B) --------------------------- #
try:
    A = aio.load_run(os.path.join(W, 'output_caseA2058'),
                     os.path.join(W, 'input.inp.caseA2058'))
    B = aio.load_run(os.path.join(W, 'output'), os.path.join(W, 'input.inp'))
    # The _adv (advection-corrected) profiles are plotted directly. The
    # metal-cooled post-processor used to land on a spurious hot root + sawtooth
    # at the breathing Roche base; the option-(c) guard in post_process_adv.f90
    # now falls back to the converged eq solution wherever the flow is not a
    # clean outflow (v<=0), so the base is smooth (see the manual caveat).
    mA = aio.mdot_Mp_per_Gyr(A); mB = aio.mdot_Mp_per_Gyr(B)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    ax[0].plot(A.r, A.T, label='spherical (A)')
    ax[0].plot(B.r, B.T, label='Roche (B)')
    ax[0].set_ylabel(r'$T$ [K]'); ax[0].legend()
    ax[1].plot(A.r, A.v_kms); ax[1].plot(B.r, B.v_kms)
    ax[1].set_ylabel(r'$v$ [km s$^{-1}$]')
    ax[2].semilogy(A.r, A.n); ax[2].semilogy(B.r, B.n)
    ax[2].set_ylabel(r'$n$ [cm$^{-3}$]')
    for a in ax:
        a.set_xlabel(r'$r$ [$R_p$]'); a.set_xlim(1, 2.0); a.grid(alpha=0.3)
    fig.suptitle(r'Spherical vs.\ Roche (RLOF): '
                 r'$\dot M_A=%.3f$, $\dot M_B=%.3f$ $M_p$/Gyr' % (mA, mB))
    save(fig, 'spherical_vs_roche.png')
except Exception as e:
    print('  [skip] spherical_vs_roche:', e)


# --- Fig 5: transmission summary (R_eff/R_star, model vs Huang) ------------ #
try:
    # Case A spherical depths from EXHALE_transit.py (line-center % and 4 A-band %),
    # converted to R_eff/R_star = sqrt((Rp/R*)^2 + h). Rp=2.058 RJ, R*=1.458 Rsun.
    Rp = 2.058 * aio.RJ
    Rstar = 1.458 * 6.957e10
    td = (Rp / Rstar) ** 2
    reff = lambda h: np.sqrt(td + h / 100.0)
    lines = ['Mg II (4A)', 'Ca II K', 'Na D2']
    # Newton solution of Case A (output_caseA2058, 2026-06-11 TPM run,
    # CHIANTI C/N/O + FS-saturation default):
    # MgII 4A-band; CaII, Na line-center [% absorption]
    model = [reff(2.674), reff(7.121), reff(1.012)]
    huangA = [0.182, 0.199, 0.152]
    huangD = [0.302, 0.278, 0.147]
    x = np.arange(len(lines)); w = 0.27
    fig, axx = plt.subplots(figsize=(7, 4.2))
    axx.bar(x - w, model, w, label='this work (Case A)')
    axx.bar(x, huangA, w, label='Huang Case A')
    axx.bar(x + w, huangD, w, label='Huang Case D')
    axx.set_xticks(x); axx.set_xticklabels(lines)
    axx.set_ylabel(r'$R_{\rm eff}/R_\star$')
    axx.set_title('Metal-line transit radii (spherical Case A)')
    axx.legend(); axx.grid(alpha=0.3, axis='y')
    save(fig, 'transmission_metals.png')
except Exception as e:
    print('  [skip] transmission_metals:', e)

# --- Fig 6: Newton-solver convergence history (WASP-121b Case B) ----------- #
try:
    import re
    log = os.path.join(W, 'ATES_caseB_newton.out')
    steps, rmarch = [], []
    iters, rjfnk = [], []
    with open(log) as f:
        for ln in f:
            m = re.search(r'\[resid\] step\s+(\d+).*max=\s*([0-9.E+-]+)', ln)
            if m:
                steps.append(int(m.group(1))); rmarch.append(float(m.group(2)))
            m = re.search(r'\(JFNK\) (?:it|start)\s*(\d*)\s+\|\|R\|\|=\s*([0-9.E+-]+)', ln)
            if m:
                iters.append(int(m.group(1)) if m.group(1) else 0)
                rjfnk.append(float(m.group(2)))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].semilogy(steps, rmarch, '-')
    ax[0].axhline(5e-2, color='r', ls='--', lw=1,
                  label=r'$R_{\rm switch}=5\times10^{-2}$')
    ax[0].set_xlabel('time-marching step')
    ax[0].set_ylabel(r'$\max\,\|R\|_\infty$ (mass, mom., energy)')
    ax[0].set_title('Stage 1: two-stage marching warm-up')
    ax[0].legend()
    ax[1].semilogy(iters, rjfnk, 'o-')
    ax[1].set_xlabel('JFNK outer iteration')
    ax[1].set_ylabel(r'$\|R\|_\infty$')
    ax[1].set_title('Stage 2: JFNK finish')
    for a in ax:
        a.grid(alpha=0.3)
    fig.suptitle(r'\texttt{Solver: Newton} pipeline from a cold IC (WASP-121b, Roche)')
    save(fig, 'newton_convergence.png')
except Exception as e:
    print('  [skip] newton_convergence:', e)


# --- Fig 7: TPM transmission spectra montage (WASP-121b, He23S on) ---------- #
# Source: the regression/tpm_wasp probe (Roche Newton solution with
# Include He23S? True), so the He I 10830 panel shows real absorption;
# the Case A planet-folder run has He23S off and an empty He line.
try:
    import shutil
    TPMW = os.path.join(ROOT, 'regression', 'tpm_wasp')
    names = ['HeI_10830.png', 'Lya.png', 'Halpha.png', 'Hbeta.png']
    imgs = [plt.imread(os.path.join(TPMW, n)) for n in names]
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    for a, im in zip(ax.flat, imgs):
        a.imshow(im)
        a.axis('off')
    fig.subplots_adjust(wspace=0.02, hspace=0.02)
    # raster montage (README/notebook use); the manual instead includes the
    # four vector PDFs copied below in a 2x2 LaTeX block
    p = os.path.join(FIG, 'tpm_spectra.png')
    fig.savefig(p, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print('  wrote', os.path.relpath(p, ROOT))
    # vector PDFs for the LaTeX docs: the four H/He spectra (manual Fig 7
    # block) plus the three metal resonance doublets (transmission_spectrum)
    for n in names + ['MgII_hk.png', 'CaII_HK.png', 'NaI_D.png']:
        src = os.path.join(TPMW, n[:-4] + '.pdf')
        if not os.path.exists(src):
            print('  [skip] missing', os.path.relpath(src, ROOT))
            continue
        dst = os.path.join(FIG, 'tpm_' + n[:-4] + '.pdf')
        shutil.copy(src, dst)
        print('  copied', os.path.relpath(dst, ROOT))
except Exception as e:
    print('  [skip] tpm_spectra:', e)

print('done.')
