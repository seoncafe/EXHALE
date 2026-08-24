#!/usr/bin/env python3
"""Figures for docs/lhs1140b_exhale_vs_pwinds.tex.

Run from LHS1140b/. Writes PDF (vector) into ../docs/figures/.
"""
import csv, os, sys
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
