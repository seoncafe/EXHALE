#!/usr/bin/env python3
"""Figures for docs/lower_atmosphere_coupling.tex.

fig_lc_planets.pdf   Tier-1 analytic column: q_H2(1 ubar; T) + 4-planet base radii
fig_mol_structure.pdf Tier-2 molecular gates: H2 fraction / T / molecular ions
fig_vulcan_vs_eq.pdf Tier-3: VULCAN photochemical H2/H vs chemical equilibrium
fig_w52_he10830.pdf  WASP-52b He 10830: old vs corrected base radius (run separately
                     once /tmp/w52o finishes; see add_w52 below)

Data: data_g1m (HD209 molecular, 3000 steps), data_g1a (HD209 atomic),
data_g2 (hot-Uranus molecular), ../../..../VULCAN/output/HD189-photo.vul.
"""
import sys, os, math, pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = False
plt.rcParams['font.size'] = 9

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(EX, 'src', 'utils'))
from run_lower import q_h2_equilibrium, integrate_column, RJ, MJ, KB, MH, G


def integrate_atomic(Mp, R1bar, Teq, fhe):
    """Fully atomic mu column (the bracket)."""
    mu = (1.0 + 4.0 * fhe) / (1.0 + fhe)
    x, dx, r = math.log(1.0), (math.log(1e-6) - math.log(1.0)) / 2000, R1bar
    for _ in range(2000):
        def f(r_):
            return -KB * Teq / (mu * MH * (G * Mp / r_**2))
        k1 = f(r); k2 = f(r + 0.5 * dx * k1)
        k3 = f(r + 0.5 * dx * k2); k4 = f(r + dx * k3)
        r += dx * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        x += dx
    return r


# ------------------------------------------------------------------ #
# Fig 1: Tier-1 analytic column
# ------------------------------------------------------------------ #
PLANETS = [  # name, Mp[MJ], Teq, HeH, R1bar, inputR0
    ('HD 209458 b', 0.720, 1450., 0.083333, 1.36,  1.401),
    ('HD 189733 b', 1.237, 1183., 0.083333, 1.138, 1.193),
    ('WASP-121 b',  1.1824, 2358., 0.0851,  1.753, 2.2075),
    ('WASP-52 b',   0.46,  1304., 0.020408, 1.27,  1.270),
]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.4))
T = np.linspace(600, 3200, 400)
ax1.plot(T, [q_h2_equilibrium(1e-6, t) for t in T], 'k-', lw=1.5)
OFF = {'HD 209458 b': (8, -14), 'HD 189733 b': (-58, -6),
       'WASP-121 b': (8, 6), 'WASP-52 b': (8, 4)}
for name, mp, teq, heh, r1, r0in in PLANETS:
    q = q_h2_equilibrium(1e-6, teq)
    ax1.plot(teq, q, 'o', ms=6)
    ax1.annotate(name, (teq, q), textcoords='offset points',
                 xytext=OFF[name], fontsize=8)
ax1.set_xlabel('T [K]')
ax1.set_ylabel(r'$q_{\rm H_2}$ (chem. equilibrium, 1 $\mu$bar)')
ax1.set_title('(a) equilibrium H$_2$ fraction at the base')
ax1.axhline(0.3, color='0.6', ls=':', lw=0.8)
ax1.text(2650, 0.33, 'molecular-base\nwarning level', fontsize=7, color='0.4')

ys = np.arange(len(PLANETS))[::-1]
for y, (name, mp, teq, heh, r1, r0in) in zip(ys, PLANETS):
    req, tb, q2, qh, qhe = integrate_column(mp * MJ, r1 * RJ, teq, heh, 'iso', 100.)
    rat = integrate_atomic(mp * MJ, r1 * RJ, teq, heh)
    ax2.plot([req / RJ, rat / RJ], [y, y], '-', lw=6, color='C0', alpha=0.4,
             solid_capstyle='butt')
    ax2.plot(r0in, y, 'k*', ms=11, zorder=5)
    ax2.plot(r1, y, 'v', ms=6, color='0.4')
    xtxt = min(req / RJ, r1, r0in) - 0.04
    ax2.text(xtxt, y, name, va='center', ha='right', fontsize=8)
ax2.plot([], [], '-', lw=6, color='C0', alpha=0.4,
         label='derived $r_0$(1 $\\mu$bar): eq.--atomic bracket')
ax2.plot([], [], 'k*', ms=10, label='input "Planet radius"')
ax2.plot([], [], 'v', color='0.4', label='$R_{\\rm 1bar}$ (transit)')
ax2.set_yticks([])
ax2.set_xlabel(r'radius [$R_{\rm J}$]')
ax2.set_title('(b) base-radius consistency check')
ax2.set_xlim(0.72, 2.5)
ax2.legend(fontsize=7, loc='upper right')
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'fig_lc_planets.pdf'))
plt.close(fig)
print('fig_lc_planets.pdf')

# ------------------------------------------------------------------ #
# Fig 2: Tier-2 molecular structure (gates 1 and 2)
# ------------------------------------------------------------------ #
def mol_profiles(d):
    I = np.loadtxt(os.path.join(HERE, d, 'Ion_species.txt'))
    H = np.loadtxt(os.path.join(HERE, d, 'Hydro_ioniz.txt'))
    r = I[:, 0]
    h2, h2p, h3p, hehp = I[:, -4], I[:, -3], I[:, -2], I[:, -1]
    nh = I[:, 1] + I[:, 2] + 2 * h2 + 2 * h2p + 3 * h3p + hehp
    return r, 2 * h2 / nh, h2, h2p, h3p, hehp, H[:, 4]

r_m, x2_m, h2m, h2pm, h3pm, hehpm, T_m = mol_profiles('data_g1m')
Ia = np.loadtxt(os.path.join(HERE, 'data_g1a', 'Ion_species.txt'))
Ha = np.loadtxt(os.path.join(HERE, 'data_g1a', 'Hydro_ioniz.txt'))
r_u, x2_u, h2u, h2pu, h3pu, hehpu, T_u = mol_profiles('data_g2')

fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.2))
ax[0].semilogy(r_m, np.maximum(x2_m, 1e-22), 'C0-', label='HD 209458 b')
ax[0].semilogy(r_u, np.maximum(x2_u, 1e-22), 'C3-', label='hot Uranus (Teq=1140 K)')
for rr, c in [(1.019, 'C0'), (1.166, 'C3')]:
    ax[0].axvline(rr, color=c, ls=':', lw=0.8)
ax[0].set_xlim(1, 3); ax[0].set_ylim(1e-14, 2)
ax[0].set_xlabel(r'$r/R_{\rm p}$')
ax[0].set_ylabel(r'$x_{\rm H_2}$ (H nuclei in H$_2$)')
ax[0].set_title(r'(a) H$_2\to$H dissociation fronts')
ax[0].legend(fontsize=7)

ax[1].plot(Ha[:, 0], Ha[:, 4], 'k--', label='atomic run')
ax[1].plot(r_m, T_m, 'C0-', label='molecular run')
ax[1].set_xlim(1, 4.2)
ax[1].set_xlabel(r'$r/R_{\rm p}$'); ax[1].set_ylabel('T [K]')
ax[1].set_title('(b) HD 209458 b: T (3000-step relaxation snapshot)')
ax[1].legend(fontsize=7)

for y, lab, c in [(h2m, r'H$_2$', 'C0'), (h2pm, r'H$_2^+$', 'C1'),
                  (h3pm, r'H$_3^+$', 'C2'), (hehpm, r'HeH$^+$', 'C4')]:
    ax[2].semilogy(r_m, np.maximum(y, 1e-22), color=c, label=lab)
ax[2].set_xlim(1, 1.35); ax[2].set_ylim(1e-8, 1e14)
ax[2].set_xlabel(r'$r/R_{\rm p}$'); ax[2].set_ylabel(r'n [cm$^{-3}$]')
ax[2].set_title('(c) HD 209458 b: molecular species')
ax[2].legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'fig_mol_structure.pdf'))
plt.close(fig)
print('fig_mol_structure.pdf')

# ------------------------------------------------------------------ #
# Fig 3: VULCAN photochemistry vs chemical equilibrium (HD 189733 b)
# ------------------------------------------------------------------ #
VUL = os.path.join(EX, '..', '..', 'VULCAN', 'output', 'HD189-photo.vul')
VUL = os.path.abspath(VUL)
if not os.path.exists(VUL):
    VUL = '/nfs/mocafe/kiseon/RT_Codes/Exoplanetary_Atmospheres/VULCAN/output/HD189-photo.vul'
d = pickle.load(open(VUL, 'rb'))
species = list(d['variable']['species'])
ymix = d['variable']['ymix']
p = d['atm']['pco'] / 1e6
Tv = d['atm']['Tco']

fig, ax = plt.subplots(figsize=(5.4, 4.0))
for sp, c, ls in [('H2', 'C0', '-'), ('H', 'C3', '-'), ('He', '0.5', '-'),
                  ('H2O', 'C2', '--'), ('CO', 'C1', '--'), ('CH4', 'C4', '--')]:
    if sp in species:
        ax.loglog(np.maximum(ymix[:, species.index(sp)], 1e-14), p,
                  color=c, ls=ls, label=sp)
# chemical-equilibrium H2 along the SAME T(p)
qe = np.array([q_h2_equilibrium(pp, tt) for pp, tt in zip(p, Tv)])
ax.loglog(np.maximum(qe, 1e-14), p, 'C0:', lw=2,
          label=r'H$_2$ (chem. equilibrium)')
ax.axhline(1e-6, color='k', lw=0.8, ls=':')
ax.text(2e-14, 6e-7, r'1 $\mu$bar handoff', fontsize=8)
ax.set_ylim(1e2, 1e-8); ax.set_xlim(1e-14, 3)
ax.set_xlabel('volume mixing ratio')
ax.set_ylabel('p [bar]')
ax.set_title('HD 189733 b: VULCAN photochemistry vs equilibrium')
ax.legend(fontsize=7, loc='lower left', ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'fig_vulcan_vs_eq.pdf'))
plt.close(fig)
print('fig_vulcan_vs_eq.pdf')

# ------------------------------------------------------------------ #
# Fig 4: WASP-52b He 10830 old vs corrected base radius (if data ready)
# ------------------------------------------------------------------ #
old_tpm = os.path.join(HERE, 'data_w52', 'tpm_old_R1p270.txt')
new_tpm = os.path.join(HERE, 'data_w52', 'tpm_new_R1p437.txt')
if os.path.exists(old_tpm) and os.path.exists(new_tpm):
    Lo = np.loadtxt(old_tpm); Ln = np.loadtxt(new_tpm)
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.plot(Lo[:, 0], (1 - Lo[:, 3]) * 100, 'k--',
            label=r'base $R_0=1.270\,R_{\rm J}$ (transit-radius shortcut)')
    ax.plot(Ln[:, 0], (1 - Ln[:, 3]) * 100, 'C0-',
            label=r'base $R_0=1.437\,R_{\rm J}$ (Tier-1 column)')
    ax.set_xlabel('wavelength [A]')
    ax.set_ylabel('absorption [%]')
    ax.set_title('WASP-52 b He 10830: base-radius correction')
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'fig_w52_he10830.pdf'))
    plt.close(fig)
    print('fig_w52_he10830.pdf')
else:
    print('fig_w52_he10830.pdf SKIPPED (old-radius rerun not finished)')
