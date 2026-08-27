#!/usr/bin/env python3
"""Figures for docs/lhs1140b_exhale_vs_pwinds.tex.

Run from LHS1140b/. Writes PDF (vector) into ../docs/figures/.
"""
import csv, os, re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq

sys.path.insert(0, '..')
plt.rcParams.update({'text.usetex': True, 'font.family': 'serif',
                     'font.size': 9, 'axes.labelsize': 10})

OUT = '../docs/figures'
os.makedirs(OUT, exist_ok=True)
AIR = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20

# Neither our p-winds oracle nor EXHALE_transit.py applies a line-of-sight
# bulk velocity, so the model curves sit at the planetary rest frame while
# the measured line is redshifted. Shift both by the value C26 fit,
# v_wind = 2.26 km/s, for the overlay. Equivalent width is unaffected.
C_KMS = 2.99792458e5
V_WIND_KMS = 2.26
dlam_air = V_WIND_KMS/C_KMS*10830.0

# ---------- observation ----------
rows = list(csv.DictReader(open(
    'Cherubim_2026/LHS1140b_He10833_Fig3B_spectrum.csv')))
g = lambda n: np.array([float(r[n]) for r in rows])
o_air  = g('wavelength_air_planet_rest_A')
o_vac  = g('wavelength_vacuum_planet_rest_A')
o_flux = g('fig4_normalized_flux')
o_fsig = g('fig4_normalized_flux_uncertainty_1sigma')
o_dep  = g('absorption_depth_positive_percent')
o_sig  = g('uncertainty_1sigma_percent')

mo = (o_vac >= EW_LO) & (o_vac <= EW_HI)
EW_obs = np.trapz(o_dep[mo], o_vac[mo])
EW_err = np.sqrt(np.sum((o_sig[mo]*np.median(np.diff(o_vac)))**2))

def exhale_curve(tag, sub=''):
    f = os.path.join('exhale', tag, sub, 'tpm_He10830.txt')
    if not os.path.isfile(f):
        return None
    s = np.loadtxt(f)
    return s[:, 0], (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100   # air, %

def red_ew(tag, sub=''):
    c = exhale_curve(tag, sub)
    if c is None:
        return np.nan
    lam = c[0]*AIR
    m = (lam >= EW_LO) & (lam <= EW_HI)
    return np.trapz(c[1][m], lam[m])

# ============ Figure 1: Fig-4 style, EXHALE + p-winds vs observation ========
fig, axs = plt.subplots(1, 2, figsize=(7.1, 3.0), sharey=True)
for a in axs:
    a.errorbar(o_air, o_flux, yerr=o_fsig, fmt='ko', ms=2.0, lw=0.6,
               capsize=0, zorder=3, label=r'2024 (GP-corrected)')
    a.axhline(1.0 - 0.006, color='0.6', ls='--', lw=0.8,
              label=r'2025 detection limit')
    a.axhline(1.0, color='0.9', lw=0.6, zorder=0)
    a.set_xlim(10827, 10831.7)
    a.set_xlabel(r'air wavelength [\AA]')
    a.ticklabel_format(axis='x', useOffset=False, style='plain')
    a.grid(alpha=0.18)
axs[0].set_ylabel(r'normalized flux')
axs[0].set_ylim(0.982, 1.006)

# the authors' own best fit (MCMC medians), computed with their own settings
for fn, lab, col, ls in (
        ('pwinds_oracle/tspec_authors_turb.txt',
         r"C26 best fit, authors' settings", 'mediumpurple', '-'),
        ('pwinds_oracle/tspec_authors_noturb.txt',
         r'same, turbulence off', 'mediumpurple', ':')):
    t = np.loadtxt(fn)
    axs[0].plot(t[:, 0]/AIR + dlam_air, 1.0 - t[:, 3]/100, color=col,
                lw=1.4 if ls == '-' else 0.9, ls=ls,
                label=r'p-winds, ' + lab)
axs[0].legend(fontsize=6, loc='lower left', framealpha=0.9)
axs[0].set_title(r'(a) p-winds at the C26 best fit (MCMC medians)', fontsize=8)
axs[0].text(0.03, 0.06, r'models shifted by $v_{\rm wind}=+2.26$\,km\,s$^{-1}$',
            transform=axs[0].transAxes, fontsize=5.5, color='0.35')

for tag, col, lab in (('solar', 'C1', r'solar, $\mathrm{He/H}=0.083$'),
                      ('heh0p55', 'C2', r'$\mathrm{He/H}=0.55$ (EW match)'),
                      ('heh1', 'C0', r'$\mathrm{He/H}=1$')):
    c = exhale_curve(tag)
    if c is not None:
        axs[1].plot(c[0] + dlam_air, 1.0 - c[1]/100, color=col, lw=0.9,
                    ls=':', alpha=0.8)
    ct = exhale_curve(tag, 'tpm_turb')
    if ct is not None:
        axs[1].plot(ct[0] + dlam_air, 1.0 - ct[1]/100, color=col, lw=1.4,
                    label=lab)
axs[1].plot([], [], color='0.4', lw=0.9, ls=':', label=r'(dotted: no turbulence)')
axs[1].legend(fontsize=6, loc='lower left', framealpha=0.9)
axs[1].set_title(r'(b) EXHALE, wind solved; solid: turbulence on', fontsize=8)
plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_fig4style.pdf')
plt.close()

# ============ Figure 2: equivalent width vs He/H ===========================
SCAN = [('solar', 0.0833), ('heh0p25', 0.25), ('heh0p5', 0.5),
        ('heh0p55', 0.55), ('heh0p6', 0.6), ('heh0p7', 0.7),
        ('heh1', 1.0), ('heh10', 10.0), ('heh100', 100.0),
        ('heh1000', 1000.0), ('heh10000', 1e4)]
hh   = np.array([s[1] for s in SCAN])
ew_n = np.array([red_ew(s[0]) for s in SCAN])
ew_t = np.array([red_ew(s[0], 'tpm_turb') for s in SCAN])

fig, ax = plt.subplots(figsize=(3.5, 2.8))
cross = {}
for y, lab, col in ((ew_n, r'no turbulence', 'C3'),
                    (ew_t, r'turbulence on', 'C0')):
    ax.loglog(hh, y, 'o-', color=col, ms=3, lw=1.1, label=lab)
    sel = hh <= 1.0
    f = lambda t: np.interp(np.log10(t), np.log10(hh[sel]), y[sel])
    cross[lab] = brentq(lambda t: f(t) - EW_obs, hh[0], 1.0)
    ax.axvline(cross[lab], color=col, ls=':', lw=0.8)
ax.axhspan(EW_obs - EW_err, EW_obs + EW_err, color='0.85', zorder=0,
           label=r'observed')
ax.set_xlabel(r'$\mathrm{He/H}$ number ratio')
ax.set_ylabel(r'red-pair EW [\%\,\AA]')
ax.grid(alpha=0.22, which='both')
ax.legend(fontsize=6, loc='lower right')
plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_ew_vs_heh.pdf')
plt.close()

# ============ Figure 3: the structural difference ==========================
lab_i = [l for l in open('exhale/heh1000/output/Ion_species.txt')
         if l.startswith('# columns')][0].split()[2:]
k = {n: j for j, n in enumerate(lab_i)}
ion = np.loadtxt('exhale/heh1000/output/Ion_species.txt')
hyd = np.loadtxt('exhale/heh1000/output/Hydro_ioniz.txt')
pw  = np.loadtxt('pwinds_oracle/profile_matched_gj1132.txt')
rp_, vp_, rhop_, fHIIp_, fHe3p_, n3p_ = pw.T

fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.4))
ax[0].semilogy(ion[:, 0], np.maximum(ion[:, k['HeITR']], 1e-8), 'C3', lw=1.3,
               label=r'EXHALE')
ax[0].semilogy(rp_, np.maximum(n3p_, 1e-8), 'C0', lw=1.3, label=r'p-winds')
ax[0].set_ylabel(r'$n(2\,^3S)$ [cm$^{-3}$]'); ax[0].set_ylim(1e-6, 1e4)
ax[1].plot(hyd[:, 0], hyd[:, 4], 'C3', lw=1.3, label=r'EXHALE')
ax[1].axhline(6100, color='C0', lw=1.3, label=r'p-winds (imposed)')
ax[1].set_ylabel(r'$T$ [K]')
ax[2].plot(hyd[:, 0], hyd[:, 2]/1e5, 'C3', lw=1.3, label=r'EXHALE')
ax[2].plot(rp_, vp_, 'C0', lw=1.3, label=r'p-winds')
ax[2].set_ylabel(r'$v$ [km\,s$^{-1}$]')
for a in ax:
    a.set_xlim(1, 20); a.set_xlabel(r'$r$ [$R_{\rm p}$]')
    a.axvline(13.6, color='0.7', ls=':', lw=0.8)
    a.grid(alpha=0.2); a.legend(fontsize=6)
plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_structure.pdf')
plt.close()

print(f'EW_obs = {EW_obs:.4f} +/- {EW_err:.4f}')
for lab, x in cross.items():
    print(f'  crossing ({lab}): He/H = {x:.3f}  H:He = {1/x:.2f}  '
          f'{x/0.0833:.1f}x solar')
print('wrote 3 PDFs to', OUT)

# ============ Figure 4: the best model of each code, Fig-4 format ==========
fig, ax = plt.subplots(figsize=(4.6, 3.2))
ax.errorbar(o_air, o_flux, yerr=o_fsig, fmt='ko', ms=2.2, lw=0.6, capsize=0,
            zorder=3, label=r'LHS\,1140\,b, 2024 (GP-corrected)')
ax.axhline(1.0 - 0.006, color='0.6', ls='--', lw=0.8,
           label=r'2025 detection limit')
ax.axhline(1.0, color='0.9', lw=0.6, zorder=0)

t = np.loadtxt('pwinds_oracle/tspec_authors_turb.txt')
ax.plot(t[:, 0]/AIR + dlam_air, 1.0 - t[:, 3]/100, color='mediumpurple',
        lw=1.5, label=r'p-winds, C26 best fit ($\mathrm{H\!:\!He}=1.0\times10^{-3}$)')

c = exhale_curve('heh0p55', 'tpm_turb')
ax.plot(c[0] + dlam_air, 1.0 - c[1]/100, color='C2', lw=1.5,
        label=r'EXHALE, $\mathrm{He/H}=0.55$ (EW-matched)')

ax.set_xlim(10827, 10831.7); ax.set_ylim(0.982, 1.006)
ax.ticklabel_format(axis='x', useOffset=False, style='plain')
ax.set_xlabel(r'air wavelength [\AA]')
ax.set_ylabel(r'normalized flux')
ax.grid(alpha=0.18)
ax.legend(fontsize=5.8, loc='lower left', framealpha=0.9)
plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_bestfit.pdf')
plt.close()
print('wrote lhs1140b_bestfit.pdf')

# ====== Figure 5: structure of the EW-matched EXHALE solution ==============
hy = np.loadtxt('exhale/heh0p55/output/Hydro_ioniz.txt')
lab_i = [l for l in open('exhale/heh0p55/output/Ion_species.txt')
         if l.startswith('# columns')][0].split()[2:]
ki = {n: j for j, n in enumerate(lab_i)}
io = np.loadtxt('exhale/heh0p55/output/Ion_species.txt')
r, n, v, T = hy[:, 0], hy[:, 1], hy[:, 2]/1e5, hy[:, 4]

fig, axg = plt.subplots(2, 3, figsize=(7.1, 3.9))
ax = axg.ravel()
ax[0].plot(r, T, 'C2', lw=1.3);            ax[0].set_ylabel(r'$T$ [K]')
ax[1].semilogy(r, n, 'C2', lw=1.3)
ax[1].set_ylabel(r'$\rho/m_{\rm H}$ [cm$^{-3}$]')
ax[2].plot(r, v, 'C2', lw=1.3)
ax[2].set_ylabel(r'$v$ [km\,s$^{-1}$]')
ax[3].semilogy(r, np.maximum(io[:, ki['HeITR']], 1e-4), 'C3', lw=1.3)
ax[3].set_ylabel(r'$n(2\,^3S)$ [cm$^{-3}$]')
ax[3].set_ylim(1.0, 100.0)
H = io[:, ki['HI']] + io[:, ki['HII']]
He = io[:, ki['HeI']] + io[:, ki['HeII']] + io[:, ki['HeIII']]
ax[4].plot(r, io[:, ki['HII']]/np.maximum(H, 1e-300), 'C2', lw=1.3, label=r'H')
ax[4].plot(r, (io[:, ki['HeII']] + io[:, ki['HeIII']])/np.maximum(He, 1e-300),
           'C2', lw=1.3, ls='--', label=r'He')
ax[4].set_ylabel(r'ionization fraction'); ax[4].legend(fontsize=6)
ax[5].axis('off')
for a in ax[:5]:
    a.set_xlim(1, 20); a.set_xlabel(r'$r$ [$R_{\rm p}$]')
    a.axvline(13.6, color='0.7', ls=':', lw=0.8)
    a.grid(alpha=0.2); a.tick_params(labelsize=7)
    a.yaxis.label.set_size(8); a.xaxis.label.set_size(8)
plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_bestfit_structure.pdf')
plt.close()
k2 = int(np.argmax(io[:, ki['HeITR']]))
print('EW-matched solution (He/H=0.55): T_max %.0f K at r=%.2f; '
      'n(2^3S) peak %.3g at r=%.2f; v(20Rp) %.2f km/s'
      % (T.max(), r[int(np.argmax(T))], io[k2, ki['HeITR']], r[k2],
         np.interp(20, r, v)))
print('wrote lhs1140b_bestfit_structure.pdf')

# ====== Figure 6: how much extra broadening the measurement demands ========
# Convolve the EW-matched EXHALE profile with a Gaussian of FWHM f km/s and
# ask which f reproduces the measured line width.  Depth, width and the
# red/blue ratio are then read off with the same three-Gaussian extractor
# used on the data, and the equivalent width over the same vacuum window.
from scipy.ndimage import gaussian_filter1d
from he_line_metrics import fit_metrics

FWHM_OBS = 0.841                      # A, measured red-pair width
RED_OBS, RED_ERR_LO, RED_ERR_HI = 1.254, 0.23, 0.22   # %, C26 MCMC

def broaden(lam, exc, fwhm_kms, lam0=10830.0):
    """Extra Gaussian velocity broadening, FWHM in km/s.

    lam0 converts the velocity width to a wavelength width; use the line
    position in the frame of `lam` (10830 A in air, 10833 A in vacuum).
    """
    if fwhm_kms <= 0.0:
        return exc
    dl = np.median(np.diff(lam))
    return gaussian_filter1d(exc, fwhm_kms/C_KMS*lam0/2.35482/dl)

def broadened_metrics(lam, exc, fwhm_kms, frame='air', lam0=10830.0):
    e = broaden(lam, exc, fwhm_kms, lam0)
    d = fit_metrics(lam, e, frame=frame)
    lv = lam*AIR if frame == 'air' else lam
    m = (lv >= EW_LO) & (lv <= EW_HI)
    d['ew'] = np.trapz(e[m], lv[m])
    return d

def matched_kernel_curve(lam, exc, frame='air', lam0=10830.0):
    """Added kernel FWHM [km/s] at which the line reaches the measured width."""
    return brentq(lambda f: broadened_metrics(lam, exc, f, frame, lam0)['fwhm_A']
                  - FWHM_OBS, 1.0, 60.0, xtol=1e-3)

def matched_kernel(tag, sub=''):
    c = exhale_curve(tag, sub)
    if c is None:
        return np.nan
    return matched_kernel_curve(c[0], c[1])

lam0, exc0 = exhale_curve('heh0p55')
f_match = matched_kernel('heh0p55')
print('matched kernel (He/H=0.55, no turbulence): %.2f km/s '
      '(sigma %.2f km/s)' % (f_match, f_match/2.35482))
for tag, sub, lab in (('heh0p55', 'tpm_turb', 'He/H=0.55, turbulence on'),
                      ('solar', '', 'solar, no turbulence')):
    print('  same crossing, %-26s %.2f km/s' % (lab, matched_kernel(tag, sub)))
print('  metrics at matched kernel: ' + ', '.join(
    '%s=%.3f' % (k, broadened_metrics(lam0, exc0, f_match)[k])
    for k in ('red_depth', 'blue_depth', 'red_blue', 'fwhm_A', 'ew')))

kern = np.linspace(0.0, 30.0, 61)
mm = [broadened_metrics(lam0, exc0, f) for f in kern]
red_k = np.array([d['red_depth'] for d in mm])
fw_k = np.array([d['fwhm_A'] for d in mm])

fig, axs = plt.subplots(1, 2, figsize=(7.1, 3.0))

a = axs[0]
a.errorbar(o_air, o_flux, yerr=o_fsig, fmt='ko', ms=2.0, lw=0.6, capsize=0,
           zorder=3, label=r'LHS\,1140\,b, 2024 (GP-corrected)')
a.axhline(1.0 - 0.006, color='0.6', ls='--', lw=0.8,
          label=r'2025 detection limit')
a.axhline(1.0, color='0.9', lw=0.6, zorder=0)
a.plot(lam0 + dlam_air, 1.0 - exc0/100, color='C2', lw=1.0, ls=':',
       label=r'EXHALE, $\mathrm{He/H}=0.55$, as solved')
a.plot(lam0 + dlam_air, 1.0 - broaden(lam0, exc0, f_match)/100, color='C2',
       lw=1.5, label=r'same, $+' + '%.1f' % f_match
       + r'$\,km\,s$^{-1}$ FWHM Gaussian')
a.set_xlim(10827, 10831.7); a.set_ylim(0.982, 1.006)
a.ticklabel_format(axis='x', useOffset=False, style='plain')
a.set_xlabel(r'air wavelength [\AA]'); a.set_ylabel(r'normalized flux')
a.grid(alpha=0.18); a.legend(fontsize=5.8, loc='lower left', framealpha=0.9)
a.set_title(r'(a) one added broadening, depth and width together', fontsize=8)

a = axs[1]
a.plot(kern, red_k, 'C3', lw=1.4, label=r'red depth (left)')
a.axhspan(RED_OBS - RED_ERR_LO, RED_OBS + RED_ERR_HI, color='C3', alpha=0.13,
          zorder=0)
a.axhline(RED_OBS, color='C3', lw=0.8, ls='--')
a.plot([f_match], [np.interp(f_match, kern, red_k)], 'o', color='C3', ms=4.5)
a.set_ylabel(r'red-pair depth [\%]', color='C3')
a.tick_params(axis='y', labelcolor='C3')
a.set_xlabel(r'added Gaussian FWHM [km\,s$^{-1}$]')
a.set_xlim(0, 30); a.set_ylim(0, 4.6)
a.grid(alpha=0.18)

b = a.twinx()
b.plot(kern, fw_k, 'C0', lw=1.4, label=r'line FWHM (right)')
b.axhline(FWHM_OBS, color='C0', lw=0.8, ls='--')
b.plot([f_match], [FWHM_OBS], 'o', color='C0', ms=4.5)
b.set_ylabel(r'red-pair FWHM [\AA]', color='C0')
b.tick_params(axis='y', labelcolor='C0')
b.set_ylim(0, 1.15)
a.axvline(f_match, color='0.35', ls=':', lw=1.0)
a.text(f_match + 0.6, 0.30, r'$' + '%.1f' % f_match + r'$\,km\,s$^{-1}$',
       fontsize=6.5, color='0.35')
h1, l1 = a.get_legend_handles_labels()
h2, l2 = b.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, fontsize=6, loc='upper right', framealpha=0.9)
a.set_title(r'(b) both cross the measurement at the same kernel', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_broadened.pdf')
plt.close()
print('wrote lhs1140b_broadened.pdf')

# ====== Figure 6b: the same demand, with diffusion and the adopted K_zz ====
# The EW crossing moves from He/H = 0.55 to 2.09 once binary H/He element
# diffusion is on at He_Kzz = 1e9 (kzz_decision.md section 6).  Repeat the
# kernel measurement on the EW-matched diffusive solution and compare the
# matched kernel with the diffusion-off 22.3 km/s.
DIFF_TAG = 'heh2p13_diff_kzz1e9'
lamd, excd = exhale_curve(DIFF_TAG)
f_match_diff = matched_kernel_curve(lamd, excd)
print('matched kernel (%s, He/H=2.13, He_Kzz=1e9): %.2f km/s (sigma %.2f), '
      'diffusion-off He/H=0.55 gives %.2f km/s, difference %.2f km/s'
      % (DIFF_TAG, f_match_diff, f_match_diff/2.35482, f_match,
         f_match_diff - f_match))
print('  scan (f [km/s], red, blue, ratio, FWHM, EW):')
for f in (0.0, 10.0, 15.0, 20.0, f_match_diff, 25.0):
    d = broadened_metrics(lamd, excd, f)
    print('    %6.2f  %.3f  %.3f  %.2f  %.3f  %.3f'
          % (f, d['red_depth'], d['blue_depth'], d['red_blue'], d['fwhm_A'],
             d['ew']))

fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2))

# Panel (a): Radial number densities of H and He species under diffusion
p_ion_diff = os.path.join('exhale', DIFF_TAG, 'output', 'Ion_species.txt')
cols_diff = [l for l in open(p_ion_diff) if l.startswith('# columns')][0].split()[2:]
jd = {c: i for i, c in enumerate(cols_diff)}
data_diff = np.loadtxt(p_ion_diff)

rd = data_diff[:, 0]
n_HI_d = data_diff[:, jd['HI']]
n_HII_d = data_diff[:, jd['HII']]
n_H_tot_d = n_HI_d + n_HII_d

n_HeI_d = data_diff[:, jd['HeI']]
n_HeII_d = data_diff[:, jd['HeII']]
n_HeIII_d = data_diff[:, jd['HeIII']]
n_HeTR_d = data_diff[:, jd['HeITR']]
n_He_tot_d = n_HeI_d + n_HeII_d + n_HeIII_d + n_HeTR_d

ax = axs[0]
ax.semilogy(rd, n_H_tot_d, color='royalblue', lw=1.6, label=r'H (total)')
ax.semilogy(rd, n_HI_d, color='royalblue', lw=1.1, ls='--', label=r'H\,\textsc{i}')
ax.semilogy(rd, n_HII_d, color='royalblue', lw=1.0, ls=':', label=r'H\,\textsc{ii}')

ax.semilogy(rd, n_He_tot_d, color='crimson', lw=1.6, label=r'He (total)')
ax.semilogy(rd, n_HeI_d, color='crimson', lw=1.1, ls='--', label=r'He\,\textsc{i}($1^1S$)')
ax.semilogy(rd, n_HeII_d, color='crimson', lw=1.0, ls=':', label=r'He\,\textsc{ii}')
ax.semilogy(rd, np.maximum(n_HeTR_d, 1e-6), color='darkmagenta', lw=1.5, ls='-', label=r'He($2\,^3S$)')

ax.axvline(13.6, color='0.7', ls=':', lw=0.8)
ax.text(13.8, 1e11, r'$R_\star = 13.6\,R_{\rm p}$', fontsize=6.5, color='0.4', rotation=90)

ax.set_xlim(1.0, 20.0)
ax.set_ylim(1e-4, 3e13)
ax.set_xlabel(r'radius $r$ [$R_{\rm p}$]')
ax.set_ylabel(r'number density $n$ [cm$^{-3}$]')
ax.grid(alpha=0.2, which='both')
ax.legend(fontsize=6.2, loc='upper right', framealpha=0.9, ncol=2)
ax.set_title(r'(a) H and He number densities ($K_{zz}=10^9$, $\mathrm{He/H}=2.13$)', fontsize=8)

# Panel (b): Broadened transit profile vs observation
ax2 = axs[1]
ax2.errorbar(o_air, o_flux, yerr=o_fsig, fmt='ko', ms=2.0, lw=0.6, capsize=0,
             zorder=3, label=r'LHS\,1140\,b (2024 GP-corrected)')
ax2.axhline(1.0 - 0.006, color='0.6', ls='--', lw=0.8, label=r'2025 detection limit')
ax2.axhline(1.0, color='0.9', lw=0.6, zorder=0)

ax2.plot(lamd + dlam_air, 1.0 - excd/100, color='C4', lw=1.0, ls=':',
         label=r'diffusion ($K_{zz}=10^{9}$, $\mathrm{He/H}=2.13$), as solved')
ax2.plot(lamd + dlam_air, 1.0 - broaden(lamd, excd, f_match_diff)/100,
         color='C4', lw=1.5,
         label=r'same, $+' + '%.1f' % f_match_diff + r'$\,km\,s$^{-1}$ FWHM Gaussian')
ax2.plot(lam0 + dlam_air, 1.0 - broaden(lam0, exc0, f_match)/100, color='C2',
         lw=1.0, ls='--', dashes=(4, 2),
         label=r'diffusion off ($\mathrm{He/H}=0.55$), $+' + '%.1f' % f_match + r'$\,km\,s$^{-1}$')

ax2.set_xlim(10827, 10831.7); ax2.set_ylim(0.982, 1.006)
ax2.ticklabel_format(axis='x', useOffset=False, style='plain')
ax2.set_xlabel(r'air wavelength [\AA]')
ax2.set_ylabel(r'normalized flux')
ax2.grid(alpha=0.18)
ax2.legend(fontsize=5.6, loc='lower left', framealpha=0.9)
ax2.set_title(r'(b) broadened transit profile vs.\ observation', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_diff_broadened.pdf')
plt.close()
print('wrote lhs1140b_diff_broadened.pdf')

# ====== Figure 7: the same scan applied to the p-winds best fit ============
# The convolution conserves equivalent width, so whether one kernel can match
# depth and width at once is decided by the EW the curve already has at f=0.
LAM0_VAC = 10833.0

PW = (('pwinds_oracle/tspec_authors_turb.txt',
       'p-winds best fit, turbulence on'),
      ('pwinds_oracle/tspec_authors_noturb.txt',
       'p-winds best fit, turbulence off'),
      ('pwinds_oracle/tspec_matched_gj1132.txt',
       'p-winds depth-matched, T=6100 K'))

pw_match = {}
for fn, lab in PW:
    t = np.loadtxt(fn)
    lv, ex = t[:, 0], t[:, 3]
    fm = matched_kernel_curve(lv, ex, frame='vacuum', lam0=LAM0_VAC)
    d0 = broadened_metrics(lv, ex, 0.0, 'vacuum', LAM0_VAC)
    dm = broadened_metrics(lv, ex, fm, 'vacuum', LAM0_VAC)
    pw_match[fn] = fm
    print('%-38s EW(f=0)=%.3f  matched %.2f km/s  red=%.3f blue=%.3f '
          'ratio=%.2f FWHM=%.3f EW=%.3f'
          % (lab, d0['ew'], fm, dm['red_depth'], dm['blue_depth'],
             dm['red_blue'], dm['fwhm_A'], dm['ew']))

fig, ax = plt.subplots(figsize=(4.6, 3.2))
ax.errorbar(o_air, o_flux, yerr=o_fsig, fmt='ko', ms=2.2, lw=0.6, capsize=0,
            zorder=3, label=r'LHS\,1140\,b, 2024 (GP-corrected)')
ax.axhline(1.0 - 0.006, color='0.6', ls='--', lw=0.8,
           label=r'2025 detection limit')
ax.axhline(1.0, color='0.9', lw=0.6, zorder=0)

fn_turb = 'pwinds_oracle/tspec_authors_turb.txt'
t = np.loadtxt(fn_turb)
f_pw = pw_match[fn_turb]
ax.plot(t[:, 0]/AIR + dlam_air, 1.0 - t[:, 3]/100, color='mediumpurple',
        lw=1.5, label=r'p-winds, C26 best fit, as published')
ax.plot(t[:, 0]/AIR + dlam_air,
        1.0 - broaden(t[:, 0], t[:, 3], f_pw, LAM0_VAC)/100,
        color='mediumpurple', lw=1.5, ls='--',
        label=r'same, $+' + '%.1f' % f_pw + r'$\,km\,s$^{-1}$ FWHM Gaussian')
ax.plot(lam0 + dlam_air, 1.0 - broaden(lam0, exc0, f_match)/100, color='C2',
        lw=1.0, label=r'EXHALE, $\mathrm{He/H}=0.55$, $+'
        + '%.1f' % f_match + r'$\,km\,s$^{-1}$')

ax.set_xlim(10827, 10831.7); ax.set_ylim(0.982, 1.006)
ax.ticklabel_format(axis='x', useOffset=False, style='plain')
ax.set_xlabel(r'air wavelength [\AA]')
ax.set_ylabel(r'normalized flux')
ax.grid(alpha=0.18)
ax.legend(fontsize=5.8, loc='lower left', framealpha=0.9)
plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_pwinds_broadened.pdf')
plt.close()
print('wrote lhs1140b_pwinds_broadened.pdf')

# ====== Figure 8: the GJ 699 proxy SED, five times the XUV ================
# Same planet, same code, same compositions; only the stellar spectrum is
# swapped for the second proxy the paper uses.  Panel (a) puts the two GJ 699
# lines nearest the measurement on the data; panel (b) shows how the
# equivalent-width crossing moves between the two SEDs.
SCAN699 = [('heh0p04_gj699', 0.04), ('heh0p06_gj699', 0.06),
           ('solar_gj699', 0.0833333), ('heh0p25_gj699', 0.25),
           ('heh0p55_gj699', 0.55), ('heh1_gj699', 1.0),
           ('heh1000_gj699', 1000.0)]
hh9 = np.array([s[1] for s in SCAN699])
ew9 = np.array([red_ew(s[0]) for s in SCAN699])

sel9 = hh9 <= 1.0
f9 = lambda t: np.interp(np.log10(t), np.log10(hh9[sel9]),
                         np.log10(ew9[sel9])) - np.log10(EW_obs)
cross9 = brentq(f9, hh9[0], 1.0)

fig, axs = plt.subplots(1, 2, figsize=(7.1, 3.0))

a = axs[0]
a.errorbar(o_air, o_flux, yerr=o_fsig, fmt='ko', ms=2.0, lw=0.6, capsize=0,
           zorder=3, label=r'LHS\,1140\,b, 2024 (GP-corrected)')
a.axhline(1.0 - 0.006, color='0.6', ls='--', lw=0.8,
          label=r'2025 detection limit')
a.axhline(1.0, color='0.9', lw=0.6, zorder=0)
c1132 = exhale_curve('solar')
a.plot(c1132[0] + dlam_air, 1.0 - c1132[1]/100, color='C1', lw=0.9, ls=':',
       label=r'GJ\,1132 SED, solar (red $'
       + '%.2f' % fit_metrics(c1132[0], c1132[1], frame='air')['red_depth']
       + r'\%$)')
for tag, col, name in (('heh0p06_gj699', 'C4', r'$\mathrm{He/H}=0.06$ (EW match)'),
                       ('solar_gj699', 'C1', r'solar, $\mathrm{He/H}=0.083$')):
    c = exhale_curve(tag)
    a.plot(c[0] + dlam_air, 1.0 - c[1]/100, color=col, lw=1.4,
           label=r'GJ\,699 SED, ' + name + r' (red $'
           + '%.2f' % fit_metrics(c[0], c[1], frame='air')['red_depth']
           + r'\%$)')
a.set_xlim(10827, 10831.7); a.set_ylim(0.982, 1.006)
a.ticklabel_format(axis='x', useOffset=False, style='plain')
a.set_xlabel(r'air wavelength [\AA]'); a.set_ylabel(r'normalized flux')
a.grid(alpha=0.18); a.legend(fontsize=5.8, loc='lower left', framealpha=0.9)
a.set_title(r'(a) five times the XUV, no turbulence term', fontsize=8)

a = axs[1]
a.loglog(hh, ew_n, 'o-', color='C0', ms=3, lw=1.1, label=r'GJ\,1132 SED')
a.loglog(hh9, ew9, 'o-', color='C3', ms=3, lw=1.1, label=r'GJ\,699 SED')
a.axhspan(EW_obs - EW_err, EW_obs + EW_err, color='0.85', zorder=0,
          label=r'observed')
a.axvline(cross[r'no turbulence'], color='C0', ls=':', lw=0.9)
a.axvline(cross9, color='C3', ls=':', lw=0.9)
a.text(cross9*1.06, 0.93, r'$' + '%.3f' % cross9 + r'$', color='C3',
       fontsize=6.5, transform=a.get_xaxis_transform())
a.text(cross[r'no turbulence']*1.06, 0.93,
       r'$' + '%.3f' % cross[r'no turbulence'] + r'$', color='C0',
       fontsize=6.5, transform=a.get_xaxis_transform())
a.set_xlabel(r'$\mathrm{He/H}$ number ratio')
a.set_ylabel(r'red-pair EW [\%\,\AA]')
a.grid(alpha=0.22, which='both')
a.legend(fontsize=6, loc='lower right')
a.set_title(r'(b) the composition estimate follows the SED', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_gj699.pdf')
plt.close()
print('GJ 699 crossing (no turbulence): He/H = %.4f  H:He = %.2f  %.2fx solar'
      % (cross9, 1/cross9, cross9/0.0833333))
print('GJ 699 matched kernel (He/H=0.06): %.2f km/s'
      % matched_kernel('heh0p06_gj699'))
print('wrote lhs1140b_gj699.pdf')

# ====== Figure 9: the adiabatic-cooling He 2^3S population bump ============
# Schulik & Owen (2025): as adiabatic cooling drops T from ~10^4 to ~10^3 K
# the collisional de-excitation q(T) collapses and the metastable fraction
# rises outward. Panels (a) and (b) locate that rise; panel (c) asks at what
# velocity it sits.
import bump_analysis as BA

BUMP = [('heh0p55', r'GJ\,1132 SED, He/H\,$=0.55$', 'C3'),
        ('solar_gj699', r'GJ\,699 SED, solar He/H', 'C0')]

fig, axb = plt.subplots(1, 3, figsize=(7.1, 2.7))

a = axb[0]
at = a.twinx()
for tag, lab, c in BUMP:
    de = BA.load(tag, False)
    da = BA.load(tag, True)
    a.semilogy(da['r'], np.maximum(da['f3'], 1e-14), color=c, lw=1.3,
               label=lab)
    a.semilogy(de['r'], np.maximum(de['f3'], 1e-14), color=c, lw=1.0, ls=':')
    at.plot(de['r'], de['T'], color=c, lw=0.8, alpha=0.35)
a.set_ylabel(r'$f_3 = n(2\,^3S)/n_{\rm He}$')
a.set_ylim(1e-11, 1e-3)
at.set_ylabel(r'$T$ [K]', fontsize=8); at.set_ylim(0, 6000)
at.tick_params(labelsize=7)
a.legend(fontsize=5.6, loc='lower right')
a.set_title(r'(a) metastable fraction and $T$', fontsize=8)

a = axb[1]
at = a.twinx()
for tag, lab, c in BUMP:
    de = BA.load(tag, False)
    da = BA.load(tag, True)
    a.semilogy(da['r'], np.maximum(da['HeITR'], 1e-4), color=c, lw=1.3)
    a.semilogy(de['r'], np.maximum(de['HeITR'], 1e-4), color=c, lw=1.0,
               ls=':')
    at.plot(de['r'], de['v']/1e5, color=c, lw=0.8, alpha=0.45)
    b = BA.bump_metrics(de)
    a.plot([b['r_peak']], [b['n_peak']], 'o', color=c, ms=3.5)
a.set_ylabel(r'$n(2\,^3S)$ [cm$^{-3}$]'); a.set_ylim(1e-2, 3e2)
at.set_ylabel(r'$v$ [km\,s$^{-1}$]', fontsize=8); at.set_ylim(0, 4.5)
at.tick_params(labelsize=7)
a.set_title(r'(b) population and outflow speed', fontsize=8)

for a in axb[:2]:
    a.set_xlim(1, 20); a.set_xlabel(r'$r$ [$R_{\rm p}$]')
    a.grid(alpha=0.2); a.tick_params(labelsize=7)
    a.yaxis.label.set_size(8); a.xaxis.label.set_size(8)

a = axb[2]
for tag, lab, c in BUMP:
    for adv, ls in ((True, '-'), (False, ':')):
        vd = BA.velocity_distribution(BA.load(tag, adv))
        a.plot(vd['v'], 100*vd['cum'], color=c, lw=1.3 if adv else 1.0,
               ls=ls, label=lab if adv else None)
a.axvline(9.5, color='0.35', ls='--', lw=1.0)
a.text(9.2, 68, r'$\sigma$ demanded by the line', rotation=90,
       fontsize=6, ha='right', va='center', color='0.35')
a.set_xlim(0, 11); a.set_ylim(0, 100)
a.set_xlabel(r'$|v_r|$ [km\,s$^{-1}$]')
a.set_ylabel(r'cumulative $2\,^3S$ column [\%]')
a.grid(alpha=0.2); a.tick_params(labelsize=7)
a.yaxis.label.set_size(8); a.xaxis.label.set_size(8)
a.legend(fontsize=5.6, loc='lower right')
a.set_title(r'(c) the column, in velocity', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_bump.pdf')
plt.close()
for tag, lab, _ in BUMP:
    de, da = BA.load(tag, False), BA.load(tag, True)
    be, ba_ = BA.bump_metrics(de), BA.bump_metrics(da)
    vde, vda = BA.velocity_distribution(de), BA.velocity_distribution(da)
    print('bump %-14s eq: r_peak %.2f T %.0f v %.2f | adv: r(f3 max) %.2f '
          'v %.2f | adv column: median |v| %.2f, >5 km/s %.2f%%'
          % (tag, be['r_peak'], be['T_peak'], be['v_peak'], ba_['r_f3'],
             ba_['v_f3'], vda['q'][0.5], 100*vda['frac'][5.0]))
print('wrote lhs1140b_bump.pdf')

# ====== Figure 9: what K_zz moves, at a fixed reservoir composition ========
# Four converged runs that differ only in the He_Kzz line: same GJ 1132 SED,
# same He/H = 0.55 reservoir, He_diffusion on throughout.  The question the
# figure answers is which parts of the solution the eddy coefficient acts on.
# Rows and provenance: kzz_decision.md sections 3 and 6.1; the operator
# itself: ../docs/binary_diffusion_design.md.
KZZ_RUNS = [('heh0p55_diff_ctrl',    0.0,   r'$K_{zz} = 0$'),
            ('heh0p55_diff_kzz1e8',  1.0e8, r'$K_{zz} = 10^{8}$'),
            ('heh0p55_diff_kzz1e9',  1.0e9, r'$K_{zz} = 10^{9}$ (adopted)'),
            ('heh0p55_diff_kzz1e10', 1.0e10, r'$K_{zz} = 10^{10}$')]
KZZ_HEH = 0.55


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
                head = line.split('columns:', 1)[1].strip()
                # The older schema wrote column names that contain spaces
                # ("F_He=4pi r^2 (F_adv+J)[g/s]") and separated the columns by
                # two or more spaces; the current one writes single-token names
                # separated by one space.  Split on whichever the line uses, so
                # a name is never mistaken for several columns.
                toks = (re.split(r'\s{2,}', head) if re.search(r'\s{2,}', head)
                        else head.split())
                if 'D_eff[cm2/s]' in toks:
                    icol = toks.index('D_eff[cm2/s]')
                break
    return p, icol


def kzz_homopause(d, kzz):
    """Radius where the run's own molecular D_eff equals its K_zz.

    Same construction as exhale/kzz_scan_table.py: the D_eff column of the
    element face-flux table is the stage-resolved binary coefficient the
    operator used, and the eddy term is added to it, so the crossing of
    D_eff with K_zz is the homopause of that run.
    """
    p, icol = _face_flux_file(d)
    if kzz <= 0.0 or p is None:
        return None
    a = np.loadtxt(p, usecols=(1, icol))
    s = a[:, 1] - kzz
    for j in range(len(s) - 1):
        if s[j]*s[j + 1] < 0.0:
            f = s[j]/(s[j] - s[j + 1])
            return a[j, 0] + f*(a[j + 1, 0] - a[j, 0])
    return None


def kzz_profiles(tag):
    """Hydrodynamic, elemental and metastable profiles of one K_zz run."""
    d = os.path.join('exhale', tag)
    hy = np.loadtxt(os.path.join(d, 'output', 'Hydro_ioniz.txt'))
    p = os.path.join(d, 'output', 'Ion_species_adv.txt')
    names = [l for l in open(p) if l.startswith('# columns')][0].split()[2:]
    a = np.loadtxt(p)
    j = {n: i for i, n in enumerate(names)}
    ra = a[:, j['r[Rp]']]
    nHe = sum(a[:, j[s]] for s in ('HeI', 'HeII', 'HeIII', 'HeITR'))
    nH = a[:, j['HI']] + a[:, j['HII']]
    return dict(r=hy[:, 0], n=hy[:, 1], v=hy[:, 2]/1e5, T=hy[:, 4],
                r_adv=ra, ratio=(nHe/np.maximum(nH, 1e-99))/KZZ_HEH,
                tr=a[:, j['HeITR']], d=d)


KZZ_COL = [plt.cm.viridis(x) for x in (0.02, 0.32, 0.55, 0.76)]
fig, axg = plt.subplots(2, 3, figsize=(7.1, 4.3))
ax = axg.ravel()
kzz_prof, kzz_hp = [], []
for (tag, kzz, lab), col in zip(KZZ_RUNS, KZZ_COL):
    P = kzz_profiles(tag)
    kzz_prof.append((tag, kzz, P))
    hp = kzz_homopause(P['d'], kzz)
    kzz_hp.append(hp)
    ax[0].plot(P['r'], P['T'], color=col, lw=1.3, label=lab)
    ax[1].plot(P['r'], P['v'], color=col, lw=1.3)
    ax[2].semilogy(P['r'], P['n'], color=col, lw=1.3)
    ax[3].semilogy(P['r_adv'], np.maximum(P['ratio'], 1e-6), color=col, lw=1.3)
    ax[4].semilogy(P['r_adv'], np.maximum(P['tr'], 1e-6), color=col, lw=1.3)
    ax[5].plot(P['r_adv'], np.maximum(P['ratio'], 1e-6), color=col, lw=1.3)
    if hp is not None:
        ax[5].axvline(hp, color=col, ls=':', lw=1.2)
        ax[5].plot([hp], [np.interp(hp, P['r_adv'], P['ratio'])],
                   'o', color=col, ms=3.5, zorder=5)

ax[0].set_ylabel(r'$T$ [K]')
ax[0].set_title(r'(a) temperature', fontsize=8)
ax[1].set_ylabel(r'$v$ [km\,s$^{-1}$]')
ax[1].set_title(r'(b) velocity', fontsize=8)
ax[2].set_ylabel(r'$\rho/m_{\rm H}$ [cm$^{-3}$]')
ax[2].set_title(r'(c) density', fontsize=8)
ax[3].set_ylabel(r'$(\mathrm{He/H})/(\mathrm{He/H})_0$')
ax[3].set_ylim(1e-4, 3.0)
ax[3].set_title(r'(d) elemental helium, against the reservoir', fontsize=8)
ax[4].set_ylabel(r'$n(2\,^3S)$ [cm$^{-3}$]')
ax[4].set_ylim(1e-4, 1e2)
ax[4].set_title(r'(e) metastable helium', fontsize=8)
for a in ax[:5]:
    a.set_xscale('log')
    a.set_xlim(1, 20)
    a.set_xticks([1, 2, 5, 10, 20])
    a.set_xticklabels([r'1', r'2', r'5', r'10', r'20'])
    a.set_xlabel(r'$r$ [$R_{\rm p}$]')

ax[5].set_xlim(1.0, 1.25)
ax[5].set_ylim(0.0, 1.05)
ax[5].set_xlabel(r'$r$ [$R_{\rm p}$]')
ax[5].set_ylabel(r'$(\mathrm{He/H})/(\mathrm{He/H})_0$')
ax[5].set_title(r'(f) the same, at the base', fontsize=8)
for a in ax:
    a.grid(alpha=0.2)
    a.tick_params(labelsize=7)
    a.yaxis.label.set_size(8)
    a.xaxis.label.set_size(8)
ax[0].legend(fontsize=5.8, loc='upper right', framealpha=0.9)
plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_kzz_profiles.pdf')
plt.close()

for (tag, kzz, P), hp in zip(kzz_prof, kzz_hp):
    mflux = np.interp(20.0, P['r'], P['n']*P['v']*1e5*P['r']**2)
    print('kzz %-22s K_zz=%.0e  homopause %s  T(1.5Rp)=%.0f  v(10Rp)=%.3f  '
          'n(2Rp)=%.3g  rho v r^2(20Rp)=%.3g  (He/H)/HeH: 1.05=%.3f 5=%.3f  '
          'n(2^3S) peak %.3g at r=%.2f'
          % (tag, kzz, ('%.4f' % hp) if hp else 'below the base',
             np.interp(1.5, P['r'], P['T']),
             np.interp(10.0, P['r'], P['v']),
             np.interp(2.0, P['r'], P['n']), mflux,
             np.interp(1.05, P['r_adv'], P['ratio']),
             np.interp(5.0, P['r_adv'], P['ratio']),
             P['tr'].max(), P['r_adv'][int(np.argmax(P['tr']))]))
print('wrote lhs1140b_kzz_profiles.pdf')


# ====== Figure: the composition the line implies, against K_zz ==============
# One point per decade of the eddy coefficient: the He/H reservoir ratio whose
# red-pair equivalent width equals the measured one.  Same solve as
# exhale/heh_diff_scan_table.py -- log-log interpolation of the scanned
# compositions, Brent's method -- but on this file's own red_ew() so the
# figure and the memo's other panels measure the line the same way.
# Provenance of the runs: kzz_decision.md section 6.1.
CROSS_SCANS = [
    (0.0, ['heh4p5_diff_kzz0', 'heh4p6_diff_kzz0', 'heh4p8_diff_kzz0',
           'heh4p9_diff_kzz0', 'heh5_diff_kzz0']),
    (1.0e5, ['heh4p5_diff_kzz1e5', 'heh4p6_diff_kzz1e5',
             'heh4p8_diff_kzz1e5', 'heh4p9_diff_kzz1e5',
             'heh5_diff_kzz1e5']),
    (1.0e6, ['heh2_diff_kzz1e6', 'heh4p3_diff_kzz1e6', 'heh4p4_diff_kzz1e6',
             'heh4p7_diff_kzz1e6', 'heh5_diff_kzz1e6',
             'heh10_diff_kzz1e6']),
    (1.0e7, ['heh3p4_diff_kzz1e7', 'heh3p6_diff_kzz1e7',
             'heh3p8_diff_kzz1e7', 'heh3p9_diff_kzz1e7',
             'heh4_diff_kzz1e7', 'heh4p2_diff_kzz1e7']),
    (1.0e8, ['heh0p55_diff_kzz1e8', 'heh2p7_diff_kzz1e8',
             'heh3_diff_kzz1e8', 'heh3p5_diff_kzz1e8',
             'heh4_diff_kzz1e8']),
    (1.0e9, ['heh0p55_diff_kzz1e9', 'heh1_diff_kzz1e9', 'heh2_diff_kzz1e9',
             'heh2p13_diff_kzz1e9', 'heh4_diff_kzz1e9']),
    (1.0e10, ['heh0p55_diff_kzz1e10', 'heh1_diff_kzz1e10',
              'heh1p4_diff_kzz1e10', 'heh1p5_diff_kzz1e10',
              'heh1p6_diff_kzz1e10']),
    (1.0e11, ['heh1_diff_kzz1e11', 'heh1p2_diff_kzz1e11',
              'heh1p4_diff_kzz1e11', 'heh1p5_diff_kzz1e11']),
]

# The decades below 1e5 are not scanned in composition.  At the fixed
# composition He/H = 5 their equivalent width differs from the K_zz = 0 one
# by far less than the measurement error, which is the measured statement the
# figure draws as a band rather than as points of their own.
FLAT_PROBE = [(1.0e1, 'heh5_diff_kzz1e1'), (1.0e2, 'heh5_diff_kzz1e2'),
              (1.0e3, 'heh5_diff_kzz1e3'), (1.0e4, 'heh5_diff_kzz1e4')]
FLAT_REF = 'heh5_diff_kzz0'
HEH_WELLMIXED = 0.55          # crossing with the operator off, section 6


def heh_at_ew(cases, target):
    """He/H whose red EW equals target, log-log through the scanned cases."""
    h = np.array([key_value(os.path.join('exhale', c, 'input.inp'),
                            'He/H number ratio') for c in cases])
    e = np.array([red_ew(c) for c in cases])
    o = np.argsort(h)
    h, e = h[o], e[o]
    if not (e.min() <= target <= e.max()):
        return np.nan
    return 10.0**brentq(lambda t: np.interp(t, np.log10(h), np.log10(e))
                        - np.log10(target), np.log10(h[0]), np.log10(h[-1]))


def key_value(path, key):
    for line in open(path):
        if line.startswith(key + ':'):
            return float(line.split(':')[1])
    raise ValueError('no ' + key + ' in ' + path)


kx, kc, klo, khi = [], [], [], []
for kzz, cases in CROSS_SCANS:
    c = heh_at_ew(cases, EW_obs)
    lo = heh_at_ew(cases, EW_obs - EW_err)     # less line -> less helium
    hi = heh_at_ew(cases, EW_obs + EW_err)
    kx.append(kzz)
    kc.append(c)
    klo.append(lo)
    khi.append(hi)
    print('cross K_zz=%-8.0e He/H = %.4f  (1 sigma %.4f - %.4f)'
          % (kzz, c, lo, hi))

heh_plateau, plateau_lo, plateau_hi = kc[0], klo[0], khi[0]
ew_ref = red_ew(FLAT_REF)
for kzz, tag in FLAT_PROBE:
    print('flat  K_zz=%-8.0e He/H = 5: EW = %.6f, %+.6f vs K_zz = 0 '
          '(%.3f sigma)' % (kzz, red_ew(tag), red_ew(tag) - ew_ref,
                            (red_ew(tag) - ew_ref)/EW_err))

fig, ax = plt.subplots(figsize=(4.6, 3.2))
X0, X1 = 3.0, 3.0e11
FLAT_X1 = 3.0e5                                # last indistinguishable decade

ax.axhspan(plateau_lo, plateau_hi, xmin=0.0,
           xmax=(np.log10(FLAT_X1) - np.log10(X0))/(np.log10(X1)
                                                    - np.log10(X0)),
           color='0.85', zorder=0)
ax.plot([X0, FLAT_X1], [heh_plateau]*2, color='0.35', ls='--', lw=1.2,
        label=r'$K_{zz}=0$ limit, He/H $= %.2f$' % heh_plateau)
ax.plot([k for k, _ in FLAT_PROBE], [heh_plateau]*len(FLAT_PROBE),
        marker='o', ls='none', mfc='none', mec='0.35', ms=5, mew=1.0,
        label=r'EW within $0.2\sigma$ of $K_{zz}=0$')

ax.axhline(HEH_WELLMIXED, color='C3', ls=':', lw=1.2)
ax.text(6.0, HEH_WELLMIXED*1.08,
        r'well-mixed limit (operator off), He/H $= 0.55$',
        fontsize=6.5, color='C3', va='bottom')

ax.axvline(1.0e9, color='C0', ls='-.', lw=1.0, alpha=0.7)
ax.text(1.0e9/1.5, 0.44, r'adopted', fontsize=6.5, color='C0',
        rotation=90, va='bottom', ha='right')

kxp = np.array(kx[1:])
kcp = np.array(kc[1:])
ax.errorbar(kxp, kcp,
            yerr=[kcp - np.array(klo[1:]), np.array(khi[1:]) - kcp],
            fmt='s', color='C0', ms=4, lw=1.0, capsize=2,
            label=r'EW-matched He/H, $\pm1\sigma$')

pw = np.array([1.0e8, 1.0e11])
ax.plot(pw, kc[5]*(pw/1.0e9)**(-0.136), color='C1', lw=1.1, alpha=0.8,
        label=r'$\propto K_{zz}^{-0.136}$')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(X0, X1)
ax.set_ylim(0.4, 8.0)
ax.set_xticks([1e1, 1e3, 1e5, 1e7, 1e9, 1e11])
ax.set_yticks([0.5, 1, 2, 5])
ax.set_yticklabels([r'0.5', r'1', r'2', r'5'])
ax.set_xlabel(r'$K_{zz}$ [cm$^{2}$\,s$^{-1}$]')
ax.set_ylabel(r'He/H matching the measured EW')
ax.grid(alpha=0.2)
ax.legend(fontsize=6.2, loc='upper right', framealpha=0.9)
plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_heh_vs_kzz.pdf')
plt.close()
print('wrote lhs1140b_heh_vs_kzz.pdf')


# ====== Figure: the flux-closed solution, and the column it stands on =======
# The elemental-flux closure of ../docs/phase_e_flux_closure_design.md
# section 6: the lower atmosphere is a photochemical column (Photochem, with
# the climate step solved) handed over as a profile, and the elemental fluxes
# are iterated to continuity across the matching level, so He/H at the match
# is a solution and not an input.  Runs: exhale/flux_closure/{ref,lo,hi}.
# Record: ../docs/Update_EXHALE.md section 79.
CLOSURE_ARMS = [('flux_closure/ref', 'k00', r'$1.0\times$ start'),
                ('flux_closure/lo', 'k05', r'$0.3\times$ start'),
                ('flux_closure/hi', 'k06', r'$3.0\times$ start')]
CLOSURE_MAIN = CLOSURE_ARMS[0]

for tag, sub, lab in CLOSURE_ARMS:
    c = exhale_curve(tag, sub)
    if c is None:
        print('closure %-18s missing' % tag)
        continue
    m = broadened_metrics(c[0], c[1], 0.0)
    print('closure %-18s %-16s red=%.3f blue=%.3f FWHM=%.4f A  EW=%.4f %%A  '
          'matched kernel %.2f km/s (sigma %.2f)'
          % (tag, lab.replace('$', '').replace('\\', ''), m['red_depth'],
             m['blue_depth'], m['fwhm_A'], m['ew'], matched_kernel(tag, sub),
             matched_kernel(tag, sub)/2.35482))

lamc, excc = exhale_curve(*CLOSURE_MAIN[:2])
f_clos = matched_kernel_curve(lamc, excc)
mc = broadened_metrics(lamc, excc, f_clos)
print('closure at its matched kernel: red=%.3f blue=%.3f ratio=%.2f '
      'FWHM=%.3f EW=%.3f' % (mc['red_depth'], mc['blue_depth'],
                             mc['red_blue'], mc['fwhm_A'], mc['ew']))

# The column handed over: T(p) and the elemental oxygen the cold trap removes.
prof = np.loadtxt('lower_profile/lower_atmosphere_profile.dat')
pnames = [l for l in open('lower_profile/lower_atmosphere_profile.dat')
          if l.startswith('# columns')][0].split()[2:]
pj = {n: i for i, n in enumerate(pnames)}
p_bar, T_col = prof[:, pj['p']], prof[:, pj['T']]
XO, XC, XN = (prof[:, pj[k]] for k in ('X_O', 'X_C', 'X_N'))
qH2O = prof[:, pj['q_H2O']]
jm = int(np.argmin(np.abs(np.log(p_bar) - np.log(1.0e-6))))   # matching level
print('column: %d levels, %.3g to %.3g bar; match at %.3g bar; '
      'O/H %.4g (deep) -> %.4g (match), factor %.0f; C/H %.4g -> %.4g; '
      'N/H %.4g -> %.4g; q_H2O %.3g (deep) -> %.3g (match)'
      % (len(p_bar), p_bar[0], p_bar[-1], p_bar[jm], XO[0], XO[jm],
         XO[0]/XO[jm], XC[0], XC[jm], XN[0], XN[jm], qH2O[0], qH2O[jm]))

fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.9),
                        gridspec_kw=dict(width_ratios=[1.55, 1.0, 1.15]))

a = axs[0]
a.errorbar(o_air, o_flux, yerr=o_fsig, fmt='ko', ms=2.0, lw=0.6, capsize=0,
           zorder=3, label=r'LHS\,1140\,b, 2024 (GP-corrected)')
a.axhline(1.0, color='0.9', lw=0.6, zorder=0)
for (tag, sub, lab), col, ls in zip(CLOSURE_ARMS, ('C0', 'C2', 'C1'),
                                    ('-', '--', ':')):
    c = exhale_curve(tag, sub)
    if c is None:
        continue
    a.plot(c[0] + dlam_air, 1.0 - c[1]/100, color=col, lw=1.1, ls=ls,
           label=r'flux-closed, ' + lab)
a.plot(lamc + dlam_air, 1.0 - broaden(lamc, excc, f_clos)/100, color='C3',
       lw=1.5, label=r'closed solution, $+' + '%.1f' % f_clos
       + r'$\,km\,s$^{-1}$ FWHM')
a.set_xlim(10827, 10831.7); a.set_ylim(0.982, 1.006)
a.ticklabel_format(axis='x', useOffset=False, style='plain')
a.set_xlabel(r'air wavelength [\AA]'); a.set_ylabel(r'normalized flux')
a.grid(alpha=0.18); a.legend(fontsize=5.4, loc='lower right', framealpha=0.9)
a.set_title(r'(a) the line at the flux-closed composition', fontsize=8)

P_LO, P_HI = p_bar.min()/1.4, p_bar.max()*1.4

a = axs[1]
a.semilogy(T_col, p_bar, color='C0', lw=1.4)
a.set_ylim(P_HI, P_LO)
a.set_xlabel(r'temperature [K]'); a.set_ylabel(r'pressure [bar]')
a.set_xlim(150, 600)
a.axhline(1.031, color='0.35', ls='-.', lw=0.9)
a.text(560, 1.031, r'tropopause', fontsize=6.5, color='0.35', ha='right',
       va='bottom')
a.axhline(1.0e-6, color='C3', ls=':', lw=0.9)
a.text(560, 1.0e-6, r'match', fontsize=6.5, color='C3', ha='right',
       va='bottom')
a.grid(alpha=0.18)
a.set_title(r'(b) the solved climate', fontsize=8)

a = axs[2]
a.loglog(np.maximum(XO, 1e-12), p_bar, color='C3', lw=1.4, label=r'O/H')
a.loglog(np.maximum(XC, 1e-12), p_bar, color='C2', lw=1.1, ls='--',
         label=r'C/H')
a.loglog(np.maximum(XN, 1e-12), p_bar, color='C1', lw=1.1, ls=':',
         label=r'N/H')
a.set_ylim(P_HI, P_LO)
a.set_xlim(1e-8, 3e-3)
a.set_xlabel(r'nuclei per H'); a.set_yticklabels([])
a.axhline(1.031, color='0.35', ls='-.', lw=0.9)
a.axhline(1.0e-6, color='C3', ls=':', lw=0.9)
a.grid(alpha=0.18)
a.legend(fontsize=6, loc='upper left', framealpha=0.9)
a.set_title(r'(c) what crosses the cold trap', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_closure.pdf')
plt.close()
print('wrote lhs1140b_closure.pdf')


# ====== Figure: the Knudsen number of the solutions, and where it bites =====
# Collisional validity of a continuum wind solution: ../docs/collisional_
# validity.md, tool ../src/utils/collisional_validity.py (it reads existing
# run directories and changes nothing).  Record: ../docs/Update_EXHALE.md
# section 81.  The three representative LHS 1140 b solutions against the
# HD 209458 b control, on the advection-corrected profiles.
sys.path.insert(0, '../src/utils')
import collisional_validity as CV

KN_CASES = [('exhale/heh0p55', 'C0', '-',
             r'well mixed, He/H\,$=0.55$'),
            ('exhale/heh2p13_diff_kzz1e9', 'C2', '--',
             r'diffusion, He/H\,$=2.13$'),
            ('exhale/flux_closure/hi/k06', 'C1', ':',
             r'flux-closed'),
            ('../backup/phase_d_baseline/new_kzz1e9_d3b', 'C3', '-.',
             r'HD\,209458\,b (control)')]

kn_res = []
for path, col, ls, lab in KN_CASES:
    try:
        kn_res.append((CV.collisional_diagnosis(path, adv=True), col, ls, lab))
    except Exception as exc:                                   # noqa: BLE001
        print('Kn %-38s unavailable (%s)' % (path, exc))

for res, col, ls, lab in kn_res:
    print('Kn %-38s sonic %s  exobase %s  Kn=0.1 at %s  max Kn crit %.3g'
          % (res['case'],
             'none' if res['r_sonic'] is None else '%.3f' % res['r_sonic'],
             'above %.2f' % res['r_top'] if res['r_exobase'] is None
             else '%.3f' % res['r_exobase'],
             '--' if res['r_kn_threshold'] is None
             else '%.3f' % res['r_kn_threshold'],
             res['Kn_crit_max']))

fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.9))

a = axs[0]
# The two radii the memo argues about: where the He 10830 line forms, and
# where the p-winds retrieval puts its isothermal Parker sonic point.
a.axvspan(1.1, 3.0, color='0.85', zorder=0)
a.text(1.8, 2.6, r'He\,10830', fontsize=6.5, color='0.35', ha='center')
a.axvspan(8.0, 9.5, color='C4', alpha=0.16, zorder=0)
a.text(8.7, 2.6, r'p-winds sonic pt.', fontsize=6.5, color='C4',
       ha='center')
for res, col, ls, lab in kn_res:
    a.loglog(res['r'], res['Kn_bulk'], color=col, ls=ls, lw=1.2, label=lab)
    if res['r_sonic'] is not None:
        a.plot([res['r_sonic']],
               [np.interp(res['r_sonic'], res['r'], res['Kn_bulk'])],
               marker='o', ms=4, color=col, mfc='none', zorder=5)
a.axhline(0.1, color='0.4', lw=0.8, ls='-')
a.text(1.15, 0.115, r'$\mathrm{Kn} = 0.1$ (continuum limit)', fontsize=6.5,
       color='0.35')
a.axhline(1.0, color='0.4', lw=0.8, ls='--')
a.text(1.15, 1.15, r'$\mathrm{Kn} = 1$ (exobase)', fontsize=6.5, color='0.35')
a.set_xlim(1.02, 30); a.set_ylim(2.0e-5, 6.0)
a.set_xlabel(r'$r/R_p$'); a.set_ylabel(r'$\mathrm{Kn}_{\rm bulk}$')
a.grid(alpha=0.18, which='both')
a.legend(fontsize=5.6, loc='lower right', framealpha=0.9)
a.set_title(r'(a) bulk Knudsen number', fontsize=8)

a = axs[1]
res0 = kn_res[0][0] if kn_res else None
if res0 is not None:
    SPLIT = [('HI', 'C0', '-', r'H\,\textsc{i}'),
             ('HeI', 'C2', '--', r'He\,\textsc{i}'),
             ('HII', 'C1', ':', r'H\,\textsc{ii}'),
             ('HeII', 'C5', '-.', r'He\,\textsc{ii}'),
             ('e', 'C3', (0, (3, 1, 1, 1)), r'$e^-$')]
    for s, col, ls, lab in SPLIT:
        if s in res0['Kn']:
            a.loglog(res0['r'], res0['Kn'][s], color=col, ls=ls, lw=1.1,
                     label=lab)
    a.loglog(res0['r'], res0['Kn_bulk'], color='k', lw=1.6, alpha=0.7,
             label=r'bulk')
a.axhline(0.1, color='0.4', lw=0.8)
a.axhline(1.0, color='0.4', lw=0.8, ls='--')
a.set_xlim(1.02, 30); a.set_ylim(1.0e-8, 2.0e1)
a.set_xlabel(r'$r/R_p$'); a.set_ylabel(r'$\mathrm{Kn}_s$')
a.grid(alpha=0.18, which='both')
a.legend(fontsize=6, loc='lower right', ncol=2, framealpha=0.9)
a.set_title(r'(b) by species, well-mixed solution', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_knudsen.pdf')
plt.close()
print('wrote lhs1140b_knudsen.pdf')


# ====== Figure: the reservoir ladder under closure, and the XUV grid ========
# (a) The equivalent width of the flux-closed solution against the reservoir
# He/H it was closed at, and where it crosses the measured line.  Runs
# exhale/flux_closure/{ref,heh3,heh5,heh8,heh9p7,heh10p3,heh12}, one arm per
# reservoir, each converged to its own fixed point; the last iterate k of
# each arm is the converged one.  (b) The red-pair depth against the XUV
# scaling, from both models that reproduce the 2024 equivalent width: the
# scalar base at He/H = 2.13 and the flux-closed solution at He/H = 9.71.
# The lower atmosphere is held fixed across the grid, so it isolates the
# wind's response.  Record: ../kzz_decision.md section 8.
CLOSURE_LADDER = ['flux_closure/ref', 'flux_closure/heh3', 'flux_closure/heh5',
                  'flux_closure/heh8', 'flux_closure/heh9p7',
                  'flux_closure/heh10p3', 'flux_closure/heh12']
XUV_LIMIT = 0.6                       # 2025 non-detection, per cent in depth
XUV_GRID = [0.01, 0.10, 0.15, 0.20, 0.25, 0.30, 0.33]
XUV_FAMILIES = [('heh2p13', 'heh2p13_diff_kzz1e9', 2.13,
                 r'scalar base, He/H $= 2.13$', 'C0', 'o', '-'),
                ('closure9p7', 'flux_closure/heh9p7/k03', 9.71,
                 r'flux-closed, He/H $= 9.71$', 'C3', 's', '--')]


def closure_last_k(tag):
    """Converged iterate of a closure arm, and the He/H it returns."""
    rows = [l.split() for l in open(os.path.join('exhale', tag,
                                                 'closure_history.txt'))
            if l.strip() and not l.startswith('#')]
    return 'k%02d' % int(rows[-1][0]), float(rows[-1][11])


def red_depth(tag, sub=''):
    """Red-pair depth [%] of the three-Gaussian fit written by the transit."""
    p = os.path.join('exhale', tag, sub, 'tpm_He10830_metrics.txt')
    if not os.path.isfile(p):
        return np.nan
    for line in open(p):
        if line.startswith('red_depth'):
            return float(line.split()[1])
    return np.nan


lad_h, lad_e = [], []
for tag in CLOSURE_LADDER:
    sub, heh = closure_last_k(tag)
    lad_h.append(heh)
    lad_e.append(red_ew(tag, sub))
    print('ladder %-22s %-4s He/H = %8.4f  EW = %.4f %%A'
          % (tag, sub, heh, lad_e[-1]))
lad_h, lad_e = np.array(lad_h), np.array(lad_e)
o = np.argsort(lad_h)
lad_h, lad_e = lad_h[o], lad_e[o]


def heh_at_closed_ew(target):
    """Reservoir whose flux-closed EW equals target, log-log on the ladder."""
    if not (lad_e.min() <= target <= lad_e.max()):
        return np.nan
    return 10.0**brentq(lambda t: np.interp(t, np.log10(lad_h),
                                            np.log10(lad_e))
                        - np.log10(target), np.log10(lad_h[0]),
                        np.log10(lad_h[-1]))


heh_closed = heh_at_closed_ew(EW_obs)
heh_clo_lo = heh_at_closed_ew(EW_obs - EW_err)
heh_clo_hi = heh_at_closed_ew(EW_obs + EW_err)
print('closure ladder crossing: He/H = %.3f (1 sigma %.3f - %.3f)'
      % (heh_closed, heh_clo_lo, heh_clo_hi))
for i in range(len(lad_h) - 1):
    print('  local EW slope %5.2f - %5.2f: %.2f'
          % (lad_h[i], lad_h[i+1],
             np.log(lad_e[i+1]/lad_e[i])/np.log(lad_h[i+1]/lad_h[i])))

fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.9))

a = axs[0]
a.axhspan(EW_obs - EW_err, EW_obs + EW_err, color='0.85', zorder=0)
a.axhline(EW_obs, color='0.35', lw=1.0, ls='-', zorder=1)
a.text(2.2, EW_obs*1.03, r'measured, $1.108 \pm 0.030$', fontsize=6.5,
       color='0.35', va='bottom')
a.plot(lad_h, lad_e, marker='o', ms=4, lw=1.2, color='C0',
       label=r'flux-closed ladder')
a.plot([heh_closed], [EW_obs], marker='*', ms=11, ls='none', color='C3',
       zorder=4, label=r'crossing, He/H $= %.1f$' % heh_closed)
a.errorbar([heh_closed], [EW_obs],
           xerr=[[heh_closed - heh_clo_lo], [heh_clo_hi - heh_closed]],
           fmt='none', ecolor='C3', lw=1.0, capsize=2, zorder=4)
a.axvline(2.0924, color='C2', ls=':', lw=1.0)
a.text(2.0924*0.95, 0.45, r'closure at the well-mixed reservoir',
       fontsize=6.2, color='C2', rotation=90, va='bottom', ha='right')
a.set_xscale('log'); a.set_yscale('log')
a.set_xlim(1.8, 14.0); a.set_ylim(0.35, 1.6)
a.set_xticks([2, 3, 5, 8, 12])
a.set_xticklabels([r'2', r'3', r'5', r'8', r'12'])
a.set_xticks([], minor=True)
a.set_yticks([0.4, 0.6, 0.8, 1.0, 1.4])
a.set_yticklabels([r'0.4', r'0.6', r'0.8', r'1.0', r'1.4'])
a.set_yticks([], minor=True)
a.set_xlabel(r'reservoir He/H below the match')
a.set_ylabel(r'red-pair $EW$ [\%\,\AA]')
a.grid(alpha=0.2, which='both')
a.legend(fontsize=6.2, loc='lower right', framealpha=0.9)
a.set_title(r'(a) the composition that closes and matches', fontsize=8)

a = axs[1]
a.axhline(XUV_LIMIT, color='0.35', lw=1.0)
a.text(0.0088, XUV_LIMIT*1.12, r'2025 limit, $0.6$\,\%', fontsize=6.5,
       color='0.35', va='bottom', ha='left')
for key, fid, heh, lab, col, mk, ls in XUV_FAMILIES:
    xs, ds = [1.0], [red_depth(*os.path.split(fid)) if '/' in fid
                     else red_depth(fid)]
    for f in XUV_GRID:
        tag = 'xuv%s_%s' % (('%.2f' % f).replace('.', 'p'), key)
        d = red_depth(tag)
        if np.isfinite(d):
            xs.append(f); ds.append(d)
    xs, ds = np.array(xs), np.array(ds)
    o = np.argsort(xs)
    a.plot(xs[o], ds[o], marker=mk, ms=4, lw=1.2, ls=ls, color=col,
           label=lab)
    print('xuv %-12s ' % key + '  '.join('%.2f:%.3f' % (x, d)
                                         for x, d in zip(xs[o], ds[o])))
a.axvline(0.30, color='C2', ls=':', lw=1.0)
a.text(0.30*1.06, 4.0e-4, r'$0.3\times$ fiducial', fontsize=6.2, color='C2',
       rotation=90, va='bottom')
a.set_xscale('log'); a.set_yscale('log')
a.set_xlim(0.008, 1.4); a.set_ylim(2.0e-4, 8.0)
a.set_xticks([0.01, 0.1, 0.3, 1.0])
a.set_xticklabels([r'0.01', r'0.1', r'0.3', r'1'])
a.set_xticks([], minor=True)
a.set_xlabel(r'$F_{\rm XUV}$ / fiducial')
a.set_ylabel(r'red-pair depth [\%]')
a.grid(alpha=0.2, which='both')
a.legend(fontsize=6.2, loc='upper left', framealpha=0.9)
a.set_title(r'(b) the XUV the 2025 non-detection allows', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_closure_ladder.pdf')
plt.close()
print('wrote lhs1140b_closure_ladder.pdf')

# ==== BEGIN thermostat block (docs sec:basemetals, Fig. lhs1140b_thermostat) ==
# Temperature and metastable density for the three solutions that bracket the
# lower-boundary axis: the metal-free scalar base at the reservoir that
# reproduces the line, the flux-closed solution on the photochemical column at
# essentially the same reservoir, and the closed solution at 10.31 -- the rung
# just below the crossing. All three carry H/He element diffusion at
# K_zz = 1e9 cm^2/s. Profiles are the advection-corrected ones the transit
# tool consumes.

THERMO_CASES = [
    ('heh2p13_diff_kzz1e9',
     r'scalar base, no metals, He/H $= 2.13$', 'C0', '-'),
    ('flux_closure/hi/k06',
     r'photochemical base, He/H $= 2.09$', 'C1', '--'),
    ('flux_closure/heh10p3/k01',
     r'photochemical base, He/H $= 10.31$', 'C3', '-.'),
]


def adv_profile(sub):
    """r [R_p], T [K] and n(2^3S) [cm^-3] from a run's *_adv.txt outputs.

    Columns are located by the '# columns' schema header rather than by
    position, so the loader follows the schema the way examples/exhale_io.py
    does.
    """
    def read(path):
        with open(path) as fh:
            fh.readline()
            cols = fh.readline().split()[2:]
        return cols, np.loadtxt(path)

    hc, hd = read(os.path.join('exhale', sub, 'output', 'Hydro_ioniz_adv.txt'))
    ic, idd = read(os.path.join('exhale', sub, 'output', 'Ion_species_adv.txt'))
    return hd[:, 0], hd[:, hc.index('T[K]')], idd[:, ic.index('HeITR')]


fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.9))

for ax in axs:
    ax.axvspan(1.0, 3.0, color='0.90', zorder=0)

for sub, lab, col, ls in THERMO_CASES:
    r, T, ntr = adv_profile(sub)
    axs[0].plot(r, T, ls=ls, lw=1.3, color=col, label=lab)
    axs[1].plot(r, ntr, ls=ls, lw=1.3, color=col, label=lab)
    m = (r >= 1.0) & (r <= 10.0)
    m2 = (r >= 1.0) & (r <= 2.0)
    print('thermostat %-26s Tmax %6.1f K  T(2Rp) %6.1f K  '
          'max n(2^3S) %7.2f  int 1-10 %7.2f  int 1-2 %7.2f'
          % (sub, T.max(), np.interp(2.0, r, T), ntr.max(),
             np.trapz(ntr[m], r[m]), np.trapz(ntr[m2], r[m2])))

a = axs[0]
a.text(1.75, 0.055, r'He\,10830 line-forming region', fontsize=6.2,
       color='0.45', rotation=90, va='bottom', ha='center',
       transform=a.get_xaxis_transform())
a.set_xscale('log')
a.set_xlim(1.0, 10.0)
a.set_ylim(0.0, 6400.0)
a.set_xticks([1, 1.5, 2, 3, 5, 10])
a.set_xticklabels([r'1', r'1.5', r'2', r'3', r'5', r'10'])
a.set_xticks([], minor=True)
a.set_xlabel(r'$r$ [$R_p$]')
a.set_ylabel(r'$T$ [K]')
a.grid(alpha=0.2, which='both')
a.legend(fontsize=6.2, loc='upper right', framealpha=0.9)
a.set_title(r'(a) the wind temperature', fontsize=8)

a = axs[1]
a.text(1.75, 0.055, r'He\,10830 line-forming region', fontsize=6.2,
       color='0.45', rotation=90, va='bottom', ha='center',
       transform=a.get_xaxis_transform())
a.set_xscale('log')
a.set_yscale('log')
a.set_xlim(1.0, 10.0)
a.set_ylim(0.1, 300.0)
a.set_xticks([1, 1.5, 2, 3, 5, 10])
a.set_xticklabels([r'1', r'1.5', r'2', r'3', r'5', r'10'])
a.set_xticks([], minor=True)
a.set_xlabel(r'$r$ [$R_p$]')
a.set_ylabel(r'$n(2\,^3S)$ [cm$^{-3}$]')
a.grid(alpha=0.2, which='both')
a.legend(fontsize=6.2, loc='upper right', framealpha=0.9)
a.set_title(r'(b) the metastable population', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_thermostat.pdf')
plt.close()
print('wrote lhs1140b_thermostat.pdf')
# ==== END thermostat block ===================================================

# ==== BEGIN photochem column block (docs sec:photochem) ======================
# The Photochem solution EXHALE is handed, column by column. Every column of
# the schema is drawn except the two elemental-flux columns, which are
# constants carrying the trial fluxes stated in the header.

PC_PROFILE = 'lower_profile/lower_atmosphere_profile.dat'
PC_TROP_BAR = 1.031          # tropopause, from the header's climate note
PC_MATCH_BAR = 1.0e-6        # matching level EXHALE reads its base from
R_JUP_CM = 6.9911e9
R_EARTH_CM = 6.3725e8


def read_lower_profile(path):
    """Named columns of a lower-atmosphere profile, plus its header notes."""
    head = [l for l in open(path) if l.startswith('#')]
    names = [l for l in head if l.startswith('# columns:')][0]
    names = names.split(':', 1)[1].split()
    tab = np.loadtxt(path)
    return {n: tab[:, i] for i, n in enumerate(names)}, head


def mark_levels(ax):
    """Tropopause and matching level on a pressure axis."""
    ax.axhline(PC_TROP_BAR, color='0.55', lw=0.8, ls=':')
    ax.axhline(PC_MATCH_BAR, color='C2', lw=0.8, ls='--')


def press_axis(ax, lo, hi):
    ax.set_yscale('log')
    ax.set_ylim(hi, lo)                      # low pressure at the top
    ax.set_ylabel(r'$p$ [bar]')
    ax.grid(alpha=0.2, which='both')


PC, PC_HEAD = read_lower_profile(PC_PROFILE)
pc_p = PC['p']
PC_LO, PC_HI = pc_p.min(), pc_p.max()
print('photochem column: %d levels, p %.4g -> %.4g bar, T %.1f -> %.1f K, '
      'Kzz unique %s'
      % (len(pc_p), pc_p[0], pc_p[-1], PC['T'][0], PC['T'][-1],
         np.unique(PC['Kzz'])))


def pc_at(name, p_bar, tab=None):
    """Value of a column at a pressure, linear in log p (deep-to-top table)."""
    t = PC if tab is None else tab
    return np.interp(np.log10(p_bar), np.log10(t['p'])[::-1], t[name][::-1])


for lab, pv in (('deep', pc_p[0]), ('tropopause', PC_TROP_BAR),
                ('1 mbar', 1.0e-3), ('match', PC_MATCH_BAR),
                ('top', pc_p[-1])):
    print('  %-11s p=%9.3e T=%6.1f q_H2=%.4g q_H=%.3e q_H2O=%.3e '
          'q_CH4=%.3e q_NH3=%.3e q_CO=%.3e q_N2=%.3e q_HCN=%.3e q_OH=%.3e '
          'q_C2H2=%.3e q_CO2=%.3e X_He=%.4f X_C=%.4e X_N=%.4e X_O=%.4e'
          % ((lab, pv) + tuple(pc_at(k, pv) for k in
             ('T', 'q_H2', 'q_H', 'q_H2O', 'q_CH4', 'q_NH3', 'q_CO', 'q_N2',
              'q_HCN', 'q_OH', 'q_C2H2', 'q_CO2', 'X_He', 'X_C', 'X_N',
              'X_O'))))

# Elemental carriers at the match, per hydrogen nucleus. The hydrogen-nucleus
# fraction of the gas is 2 q_H2 + q_H, so a carrier with n nuclei of element
# El contributes n q / (2 q_H2 + q_H) to El/H.
f_H_nuclei = 2.0*pc_at('q_H2', PC_MATCH_BAR) + pc_at('q_H', PC_MATCH_BAR)
print('  H-nucleus fraction at the match: %.5f' % f_H_nuclei)
for el, carriers in (('C', (('q_CH4', 1), ('q_CO', 1), ('q_CO2', 1),
                            ('q_HCN', 1), ('q_C2H2', 2))),
                     ('N', (('q_NH3', 1), ('q_N2', 2), ('q_HCN', 1))),
                     ('O', (('q_H2O', 1), ('q_CO', 1), ('q_CO2', 2),
                            ('q_OH', 1)))):
    per_H = {k: n*pc_at(k, PC_MATCH_BAR)/f_H_nuclei for k, n in carriers}
    tot = sum(per_H.values())
    print('  %s: X_%s = %.5e, carriers sum %.5e (%.1f%% of it)'
          % (el, el, pc_at('X_'+el, PC_MATCH_BAR), tot,
             tot/pc_at('X_'+el, PC_MATCH_BAR)*100))
    for k, _ in carriers:
        print('       %-8s %.4e per H  (%5.2f%% of the carriers)'
              % (k, per_H[k], per_H[k]/tot*100))

# ---- page 1: structure, and the carriers that hold the bulk ----
fig, axs = plt.subplots(2, 3, figsize=(7.1, 5.4))

a = axs[0, 0]
a.plot(PC['T'], pc_p, lw=1.3, color='C3')
press_axis(a, PC_LO, PC_HI); mark_levels(a)
a.set_xlabel(r'$T$ [K]')
a.set_xlim(160.0, 560.0)
a.text(540.0, PC_TROP_BAR*0.55, r'tropopause', fontsize=6.0, color='0.45',
       ha='right')
a.text(540.0, PC_MATCH_BAR*0.5, r'match', fontsize=6.0, color='C2',
       ha='right')
a.set_title(r'(a) temperature', fontsize=8)

a = axs[0, 1]
a.plot(PC['Kzz'], pc_p, lw=1.3, color='C0')
press_axis(a, PC_LO, PC_HI); mark_levels(a)
a.set_xscale('log'); a.set_xlim(1.0e7, 1.0e11)
a.set_xlabel(r'$K_{zz}$ [cm$^2$\,s$^{-1}$]')
a.set_title(r'(b) eddy coefficient, constant', fontsize=8)

a = axs[0, 2]
a.plot(PC['n_tot'], pc_p, lw=1.3, color='C0', label=r'$n_{\rm tot}$')
press_axis(a, PC_LO, PC_HI); mark_levels(a)
a.set_xscale('log')
a.set_xlabel(r'$n_{\rm tot}$ [cm$^{-3}$] \ / \ $\rho$ [g\,cm$^{-3}$]')
a.plot(PC['rho'], pc_p, lw=1.3, color='C1', ls='--', label=r'$\rho$')
a.legend(fontsize=6.2, loc='lower left', framealpha=0.9)
a.set_title(r'(c) number and mass density', fontsize=8)

a = axs[1, 0]
a.plot((PC['r']*R_JUP_CM - PC['r'][0]*R_JUP_CM)/1.0e5, pc_p, lw=1.3,
       color='C4')
press_axis(a, PC_LO, PC_HI); mark_levels(a)
a.set_xlabel(r'height above the deep boundary [km]')
a.set_title(r'(d) the radius the column maps to', fontsize=8)
a.text(0.04, 0.06, r'match at $%.4f\,R_{\rm J} = %.3f\,R_\oplus$'
       % (pc_at('r', PC_MATCH_BAR),
          pc_at('r', PC_MATCH_BAR)*R_JUP_CM/R_EARTH_CM),
       fontsize=6.0, transform=a.transAxes)

a = axs[1, 1]
q_He = PC['X_He']*(2.0*PC['q_H2'] + PC['q_H'])
for y, lab, col, ls in ((q_He, r'He', 'C7', '-'),
                        (PC['q_H2'], r'H$_2$', 'C0', '-'),
                        (PC['q_CH4'], r'CH$_4$', 'C1', '--'),
                        (PC['q_NH3'], r'NH$_3$', 'C2', '-.'),
                        (PC['q_H'], r'H', 'C3', ':')):
    a.plot(y, pc_p, lw=1.3, color=col, ls=ls, label=lab)
press_axis(a, PC_LO, PC_HI); mark_levels(a)
a.set_xscale('log'); a.set_xlim(1.0e-22, 5.0)
a.set_xlabel(r'mixing ratio $q$')
a.legend(fontsize=6.0, loc='lower left', ncol=2, framealpha=0.9)
a.set_title(r'(e) the bulk carriers', fontsize=8)

a = axs[1, 2]
for k, lab, col, ls in (('X_He', r'He/H', 'C7', '-'),
                        ('X_C', r'C/H', 'C1', '--'),
                        ('X_N', r'N/H', 'C2', '-.'),
                        ('X_O', r'O/H', 'C0', ':')):
    a.plot(PC[k], pc_p, lw=1.3, color=col, ls=ls, label=lab)
press_axis(a, PC_LO, PC_HI); mark_levels(a)
a.set_xscale('log'); a.set_xlim(1.0e-7, 20.0)
a.set_xlabel(r'nuclei per H nucleus')
a.legend(fontsize=6.0, loc='lower left', ncol=2, framealpha=0.9)
a.set_title(r'(f) elemental ratios', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_photochem_column.pdf')
plt.close()
print('wrote lhs1140b_photochem_column.pdf')

# ---- page 2: the O, C and N systems, and what the reservoir does to them ----
LADDER_ARMS = [('flux_closure/hi/k06', r'He/H $= 2.09$', 'C0', '-'),
               ('flux_closure/heh5/k04', r'He/H $= 5.01$', 'C1', '--'),
               ('flux_closure/heh10p3/k01', r'He/H $= 10.31$', 'C3', '-.')]

fig, axs = plt.subplots(2, 3, figsize=(7.1, 5.4))

a = axs[0, 0]
for k, lab, col, ls in (('q_H2O', r'H$_2$O', 'C0', '-'),
                        ('q_CO', r'CO', 'C1', '--'),
                        ('q_CO2', r'CO$_2$', 'C2', '-.'),
                        ('q_OH', r'OH', 'C3', ':')):
    a.plot(PC[k], pc_p, lw=1.3, color=col, ls=ls, label=lab)
press_axis(a, PC_LO, PC_HI); mark_levels(a)
a.set_xscale('log'); a.set_xlim(1.0e-30, 1.0e-2)
a.set_xlabel(r'mixing ratio $q$')
a.legend(fontsize=6.0, loc='lower left', ncol=2, framealpha=0.9)
a.set_title(r'(a) the oxygen system', fontsize=8)

a = axs[0, 1]
for k, lab, col, ls in (('q_CH4', r'CH$_4$', 'C1', '-'),
                        ('q_CO', r'CO', 'C0', '--'),
                        ('q_CO2', r'CO$_2$', 'C2', '-.'),
                        ('q_HCN', r'HCN', 'C3', ':'),
                        ('q_C2H2', r'C$_2$H$_2$', 'C4', (0, (3, 1, 1, 1)))):
    a.plot(PC[k], pc_p, lw=1.3, color=col, ls=ls, label=lab)
press_axis(a, PC_LO, PC_HI); mark_levels(a)
a.set_xscale('log'); a.set_xlim(1.0e-36, 1.0e-2)
a.set_xlabel(r'mixing ratio $q$')
a.legend(fontsize=6.0, loc='lower left', ncol=2, framealpha=0.9)
a.set_title(r'(b) the carbon system', fontsize=8)

a = axs[0, 2]
for k, lab, col, ls in (('q_NH3', r'NH$_3$', 'C2', '-'),
                        ('q_N2', r'N$_2$', 'C0', '--'),
                        ('q_HCN', r'HCN', 'C3', ':')):
    a.plot(PC[k], pc_p, lw=1.3, color=col, ls=ls, label=lab)
press_axis(a, PC_LO, PC_HI); mark_levels(a)
a.set_xscale('log'); a.set_xlim(1.0e-24, 1.0e-2)
a.set_xlabel(r'mixing ratio $q$')
a.legend(fontsize=6.0, loc='lower left', framealpha=0.9)
a.set_title(r'(c) the nitrogen system', fontsize=8)

arms = []
for sub, lab, col, ls in LADDER_ARMS:
    f = os.path.join('exhale', sub, 'lower_atmosphere_profile.dat')
    if not os.path.isfile(f):
        print('  ladder arm missing: %s' % f)
        continue
    tab, _ = read_lower_profile(f)
    arms.append((tab, lab, col, ls))
    print('  ladder column %-26s T(deep) %6.1f K  q_H2(match) %.5f  '
          'O/H(match) %.4e' % (sub, tab['T'][0],
                               pc_at('q_H2', PC_MATCH_BAR, tab),
                               pc_at('X_O', PC_MATCH_BAR, tab)))

for a, key, xl, xlim, ttl in (
        (axs[1, 0], 'T', r'$T$ [K]', (160.0, 560.0),
         r'(d) the column, by reservoir'),
        (axs[1, 1], 'q_H2', r'$q_{\rm H_2}$', (1.0e-2, 1.0),
         r'(e) H$_2$, by reservoir'),
        (axs[1, 2], 'X_O', r'O/H [nuclei per H]', (1.0e-7, 1.0e-2),
         r'(f) the cold trap, by reservoir')):
    for tab, lab, col, ls in arms:
        a.plot(tab[key], tab['p'], lw=1.3, color=col, ls=ls, label=lab)
    press_axis(a, PC_LO, PC_HI); mark_levels(a)
    if key != 'T':
        a.set_xscale('log')
    a.set_xlim(*xlim)
    a.set_xlabel(xl)
    a.legend(fontsize=6.0, loc='lower left', framealpha=0.9)
    a.set_title(ttl, fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_photochem_column_2.pdf')
plt.close()
print('wrote lhs1140b_photochem_column_2.pdf')
# ==== END photochem column block =============================================

# ==== BEGIN composition profiles block (docs sec:basemetals, Fig. lhs1140b_composition_profiles)
# Elemental number-density fractions and ionization fractions against radius,
# for the same three solutions as the thermostat figure and from the same
# advection-corrected outputs. n_tot follows the code's calc_ntot convention:
# every gas particle counts once (H I, H II, He I, He II, He III and, where
# present, all metal ion stages), electrons excluded, and He 2^3S is NOT added
# because it is an excited level inside He I (bsp_is_excited_level).

COMP_CASES = [
    ('heh2p13_diff_kzz1e9', 2.13,
     r'scalar base, no metals: $2.13$', 'C0', '-'),
    ('flux_closure/hi/k06', 2.0924,
     r'photochemical base: $2.09$', 'C1', '--'),
    ('flux_closure/heh10p3/k01', 10.3116,
     r'photochemical base: $10.31$', 'C3', '-.'),
]
COMP_METALS = ['C', 'O', 'N', 'Mg', 'Si', 'Ca', 'Na', 'K', 'S', 'Fe']


def composition_profile(sub):
    """Nuclei fractions and ionization fractions of one run, against r."""
    def read(path):
        with open(path) as fh:
            fh.readline()
            cols = fh.readline().split()[2:]
        return cols, np.loadtxt(path)

    ic, tab = read(os.path.join('exhale', sub, 'output', 'Ion_species_adv.txt'))
    r = tab[:, 0]
    c = {nm: tab[:, k+1] for k, nm in enumerate(ic[1:])}
    n_H = c['HI'] + c['HII']
    n_He = c['HeI'] + c['HeII'] + c['HeIII']
    n_met = np.zeros_like(n_H)
    for el in COMP_METALS:
        for stage in (el+'I', el+'II', el+'III'):
            if stage in c:
                n_met = n_met + c[stage]
    n_tot = n_H + n_He + n_met
    return dict(r=r, fH=n_H/n_tot, fHe=n_He/n_tot, fmet=n_met/n_tot,
                heh=n_He/n_H, xHII=c['HII']/n_H,
                xHeII=c['HeII']/n_He, xHeIII=c['HeIII']/n_He)


fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.9))
for ax in axs:
    ax.axvspan(1.0, 3.0, color='0.90', zorder=0)

for sub, heh0, lab, col, ls in COMP_CASES:
    P = composition_profile(sub)
    r = P['r']
    axs[0].plot(r, P['fH'], ls=ls, lw=1.3, color=col, label=lab)
    axs[0].plot(r, P['fHe'], ls=ls, lw=1.0, color=col, alpha=0.55)
    axs[1].plot(r, P['heh'], ls=ls, lw=1.3, color=col, label=lab)
    axs[1].axhline(heh0, color=col, ls=':', lw=0.8)
    axs[2].plot(r, P['xHII'], ls=ls, lw=1.3, color=col, label=lab)
    axs[2].plot(r, P['xHeII'], ls=ls, lw=1.0, color=col, alpha=0.55)
    axs[2].plot(r, P['xHeIII'], ls=ls, lw=0.8, color=col, alpha=0.35)
    at = lambda y, q: np.interp(q, r, y)
    print('composition %-26s @2Rp  nH/ntot %.4f  nHe/ntot %.4f  metals %.3e  '
          'He/H %.4f (reservoir %.4f)  x(HII) %.4f  x(HeII) %.4f  x(HeIII) %.3e'
          % (sub, at(P['fH'], 2.0), at(P['fHe'], 2.0), at(P['fmet'], 2.0),
             at(P['heh'], 2.0), heh0, at(P['xHII'], 2.0), at(P['xHeII'], 2.0),
             at(P['xHeIII'], 2.0)))
    print('     metal share of n_tot: max %.3e over 1-30 Rp'
          % P['fmet'][r >= 1.0].max())

for a in axs:
    a.set_xscale('log')
    a.set_xlim(1.0, 10.0)
    a.set_xticks([1, 1.5, 2, 3, 5, 10])
    a.set_xticklabels([r'1', r'1.5', r'2', r'3', r'5', r'10'])
    a.set_xticks([], minor=True)
    a.set_xlabel(r'$r$ [$R_p$]')
    a.grid(alpha=0.2, which='both')
    a.text(1.75, 0.055, r'He\,10830 line-forming region', fontsize=6.0,
           color='0.45', rotation=90, va='bottom', ha='center',
           transform=a.get_xaxis_transform())

a = axs[0]
a.set_ylim(0.0, 1.0)
a.set_ylabel(r'nuclei fraction of $n_{\rm tot}$')
a.text(3.15, 0.115, r'heavy: $n_{\rm H}/n_{\rm tot}$', fontsize=6.2)
a.text(3.15, 0.045, r'light: $n_{\rm He}/n_{\rm tot}$', fontsize=6.2)
a.legend(fontsize=6.0, loc='upper left', framealpha=0.9,
         title=r'reservoir He/H', title_fontsize=6.0)
a.set_title(r'(a) elemental fractions', fontsize=8)

a = axs[1]
a.set_yscale('log')
a.set_ylim(0.1, 60.0)
a.set_ylabel(r'He/H by nuclei')
a.axhline(1.0, color='0.35', lw=0.8)
a.text(1.06, 1.10, r'He $=$ H', fontsize=6.2, color='0.35')
a.text(1.06, 0.125, r'dotted: the reservoir each run was given', fontsize=6.0)
a.set_title(r'(b) helium against hydrogen', fontsize=8)

a = axs[2]
a.set_yscale('log')
a.set_ylim(1.0e-5, 3.0)
a.set_ylabel(r'ionization fraction')
a.text(3.05, 6.6e-5, r'heavy: $x({\rm H\,II})$', fontsize=6.0)
a.text(3.05, 3.4e-5, r'light: $x({\rm He\,II})$', fontsize=6.0)
a.text(3.05, 1.75e-5, r'faint: $x({\rm He\,III})$', fontsize=6.0)
a.set_title(r'(c) ionization state', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_composition_profiles.pdf')
plt.close()
print('wrote lhs1140b_composition_profiles.pdf')
# ==== END composition profiles block =========================================
