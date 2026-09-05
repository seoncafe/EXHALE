#!/usr/bin/env python3
"""Figures for docs/lhs1140b_exhale_vs_pwinds.tex.

Run from LHS1140b/. Writes PDF (vector) into ../docs/figures/.

STALE (P48; Update_EXHALE section 137): the lhs1140b_*.pdf figures currently in
../docs/figures/ were made while this script read the profile files' GHOST rows
as solution cells. It now reads through exhale_io.loadtxt_cells, so re-running
it moves the He 10830 depths and equivalent widths (0.2-3.8 per cent on the
depths measured elsewhere) and trims one point from each end of every radial
curve. Nothing was regenerated; regeneration awaits instruction. Affected
files: ../docs/figures/README_STALE_P48.md.
"""
import csv, os, re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq

sys.path.insert(0, '..')
# The EXHALE profile files carry two GHOST rows at each end -- a fixed base
# state below, a zero-gradient / WENO3 extrapolation above -- and they are not
# solution cells. exhale_io.loadtxt_cells drops them; every read of a
# Hydro_ioniz / Ion_species file below goes through it. The p-winds oracle
# profiles, the tpm_*.txt transit curves, the SED, the lower-atmosphere
# profiles and output/element_flux_profile.txt (written over the faces
# j = 1..N-1) have no ghost rows and keep plain np.loadtxt.
sys.path.insert(0, os.path.join('..', 'examples'))
from exhale_io import loadtxt_cells

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
EW_digitized = np.trapz(o_dep[mo], o_vac[mo])
EW_digitized_err = np.sqrt(np.sum((o_sig[mo]*np.median(np.diff(o_vac)))**2))

# Every crossing quoted in the memo, and every one in the results.txt of the
# run directories the crossings are read from, is solved against the digitized
# value as it is printed there: 1.108 +/- 0.030 %A.  Use the same two numbers
# here so that the figures and the tables cross at the same composition.
EW_obs, EW_err = 1.108, 0.030

def curve_file(path):
    """Air wavelength and excess absorption [%] of one transit file."""
    if not os.path.isfile(path):
        return None
    s = np.loadtxt(path)
    return s[:, 0], (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100   # air, %


def exhale_curve(tag, sub=''):
    return curve_file(os.path.join('exhale', tag, sub, 'tpm_He10830.txt'))

def red_ew(tag, sub=''):
    c = exhale_curve(tag, sub)
    if c is None:
        return np.nan
    lam = c[0]*AIR
    m = (lam >= EW_LO) & (lam <= EW_HI)
    return np.trapz(c[1][m], lam[m])

# The arm the memo calls the EW-matched solution: the run solved AT the
# well-mixed crossing He/H = 0.4132, on the current binary
# (crossings_j96/results.txt Sect. 6, row W1132p).
EWMATCH_TAG = 'crossings_j96/cf_W1132_P'
EWMATCH_HEH = 0.4132

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

# Three compositions of the well-mixed ladder of record on the current binary
# (crossings_j96 lineage W): the two ends of the ladder and the solved arm
# nearest the crossing 0.4132.  Each arm's run directory was copied into
# exhale/refresh_j96/fig4style and synthesized twice, nominal and with
# EXHALE_TRANSIT_TURB=1 into tpm_turb/, so both curves are on one physics.
# Drawn deepest first, so that the shallowest arm is not hidden under the
# other two: over this range the profile hardly moves with composition.
for fn, fn_turb, col, wid, lab in (
        ('exhale/refresh_j96/fig4style/wm_heh0p46/tpm_He10830.txt',
         'exhale/refresh_j96/fig4style/wm_heh0p46/tpm_turb/tpm_He10830.txt',
         'C1', 2.6, r'$\mathrm{He/H}=0.46$'),
        ('exhale/refresh_j96/fig4style/wm_heh0p415/tpm_He10830.txt',
         'exhale/refresh_j96/fig4style/wm_heh0p415/tpm_turb/tpm_He10830.txt',
         'C2', 1.5, r'$\mathrm{He/H}=0.415$ (EW match)'),
        ('exhale/refresh_j96/fig4style/wm_heh0p40/tpm_He10830.txt',
         'exhale/refresh_j96/fig4style/wm_heh0p40/tpm_turb/tpm_He10830.txt',
         'C0', 0.8, r'$\mathrm{He/H}=0.40$')):
    c = curve_file(fn)
    if c is not None:
        axs[1].plot(c[0] + dlam_air, 1.0 - c[1]/100, color=col,
                    lw=0.7*wid, ls=':', alpha=0.8)
    ct = curve_file(fn_turb)
    if ct is not None:
        axs[1].plot(ct[0] + dlam_air, 1.0 - ct[1]/100, color=col, lw=wid,
                    label=lab)
axs[1].plot([], [], color='0.4', lw=0.9, ls=':', label=r'(dotted: no turbulence)')
axs[1].legend(fontsize=6, loc='lower left', framealpha=0.9)
axs[1].set_title(r'(b) EXHALE, wind solved; solid: turbulence on', fontsize=8)
plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_fig4style.pdf')
plt.close()

# ============ Figure 2: equivalent width vs He/H ===========================
# Two ladders, read the same way, both on the current binary.  The wide scan
# of the memo's Sect. sec:scan is the one that carries the turbulence pair;
# every one of its seven arms was re-solved and synthesized twice
# (exhale/refresh_j96/scan).  The crossing of record is the ten-arm ladder
# over He/H = 0.40-0.46, exhale/crossings_j96 (lineage W, results.txt Sect. 1),
# whose spacing brackets the measurement four times more closely.


def red_ew_file(path):
    """Red-pair EW [%A] of one synthesized transit file."""
    lam, exc = curve_file(path)
    lam = lam*AIR
    m = (lam >= EW_LO) & (lam <= EW_HI)
    return np.trapz(exc[m], lam[m])


def ew_crossing(h, e, target):
    """He/H at which the ladder (h, e) reaches target, log-log."""
    o = np.argsort(h)
    h, e = np.asarray(h)[o], np.asarray(e)[o]
    return 10.0**brentq(lambda t: np.interp(t, np.log10(h), np.log10(e))
                        - np.log10(target), np.log10(h[0]), np.log10(h[-1]))


# the wide scan over the range that holds the crossing, re-solved on the
# current binary and synthesized twice per arm (exhale/refresh_j96/scan)
_SC = 'exhale/refresh_j96/scan/%s/%stpm_He10830.txt'
SCAN_EARLIER = [(h, _SC % (t, ''), _SC % (t, 'tpm_turb/')) for h, t in (
    (0.0833, 'solar'), (0.25, 'heh0p25'), (0.50, 'heh0p5'),
    (0.55, 'heh0p55'), (0.60, 'heh0p6'), (0.70, 'heh0p7'),
    (1.00, 'heh1'))]

# the ladder of record, current binary (crossings_j96 lineage W, solved on
# the composition-renormalized SvS85 split of Update_EXHALE Sect. 96).  The
# two arms below 0.40 that directory also carries are left out here: they are
# below the +1 sigma end and their equivalent width is not monotone.
LADDER_WM = [('crossings_j96/wm_heh0p40', 0.400),
             ('crossings_j96/wm_heh0p405', 0.405),
             ('crossings_j96/wm_heh0p41', 0.410),
             ('crossings_j96/wm_heh0p415', 0.415),
             ('crossings_j96/wm_heh0p42', 0.420),
             ('crossings_j96/wm_heh0p425', 0.425),
             ('crossings_j96/wm_heh0p43', 0.430),
             ('crossings_j96/wm_heh0p44', 0.440),
             ('crossings_j96/wm_heh0p45', 0.450),
             ('crossings_j96/wm_heh0p46', 0.460)]

hh = np.array([c[0] for c in SCAN_EARLIER])
ew_n = np.array([red_ew_file(c[1]) for c in SCAN_EARLIER])
ew_t = np.array([red_ew_file(c[2]) for c in SCAN_EARLIER])
hw = np.array([c[1] for c in LADDER_WM])
ew_w = np.array([red_ew(c[0]) for c in LADDER_WM])

cross = {}
for y, lab in ((ew_n, r'no turbulence'), (ew_t, r'turbulence on')):
    cross[lab] = ew_crossing(hh, y, EW_obs)
wm_cross = ew_crossing(hw, ew_w, EW_obs)
wm_lo = ew_crossing(hw, ew_w, EW_obs - EW_err)
wm_hi = ew_crossing(hw, ew_w, EW_obs + EW_err)

fig, ax = plt.subplots(figsize=(3.5, 2.8))
ax.axhspan(EW_obs - EW_err, EW_obs + EW_err, color='0.85', zorder=0,
           label=r'observed')
ax.axvspan(wm_lo, wm_hi, color='C2', alpha=0.13, zorder=0)
for y, lab, col, ls in ((ew_n, r'no turbulence', 'C3', '-'),
                        (ew_t, r'turbulence on', 'C0', '--')):
    ax.loglog(hh, y, marker='o', ls=ls, color=col, ms=3, lw=1.0, alpha=0.85,
              label=r'wide scan, ' + lab)
    ax.axvline(cross[lab], color=col, ls=':', lw=0.8)
ax.loglog(hw, ew_w, 's-', color='C2', ms=3.5, lw=1.5,
          label=r'ladder of record')
ax.axvline(wm_cross, color='C2', ls=':', lw=1.2)
ax.plot([wm_cross], [EW_obs], marker='*', ms=10, ls='none', color='C2',
        zorder=5)
ax.text(wm_cross*0.95, 0.26, r'$' + '%.4f' % wm_cross + r'$', color='C2',
        fontsize=6.5, rotation=90, va='bottom', ha='right')
ax.text(cross[r'turbulence on']*1.06, 0.26,
        r'$' + '%.3f' % cross[r'no turbulence'] + r'$ / $'
        + '%.3f' % cross[r'turbulence on'] + r'$', color='0.35',
        fontsize=6.5, rotation=90, va='bottom', ha='left')
ax.set_xlim(0.075, 1.2)
ax.set_ylim(0.24, 1.9)
ax.set_xticks([0.1, 0.2, 0.4, 0.7, 1.0])
ax.set_xticklabels([r'0.1', r'0.2', r'0.4', r'0.7', r'1'])
ax.set_xticks([], minor=True)
ax.set_yticks([0.3, 0.5, 0.8, 1.1, 1.6])
ax.set_yticklabels([r'0.3', r'0.5', r'0.8', r'1.1', r'1.6'])
ax.set_yticks([], minor=True)
ax.set_xlabel(r'$\mathrm{He/H}$ number ratio')
ax.set_ylabel(r'red-pair EW [\%\,\AA]')
ax.grid(alpha=0.22, which='both')
ax.legend(fontsize=5.8, loc='upper left', framealpha=0.9)
plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_ew_vs_heh.pdf')
plt.close()
print('well-mixed crossing (current binary): He/H = %.4f '
      '(1 sigma %.4f - %.4f), H:He = %.2f, %.1fx solar'
      % (wm_cross, wm_lo, wm_hi, 1/wm_cross, wm_cross/0.0833))
for lab, x in cross.items():
    print('  wide scan (%s): He/H = %.4f' % (lab, x))
print('  wide scan, turbulence moves the crossing by %.1f per cent'
      % (100*abs(cross[r'turbulence on']/cross[r'no turbulence'] - 1)))

# ============ Figure 3: the structural difference ==========================
lab_i = [l for l in open('exhale/heh1000/output/Ion_species.txt')
         if l.startswith('# columns')][0].split()[2:]
k = {n: j for j, n in enumerate(lab_i)}
ion = loadtxt_cells('exhale/heh1000/output/Ion_species.txt')
hyd = loadtxt_cells('exhale/heh1000/output/Hydro_ioniz.txt')
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

print(f'EW target = {EW_obs:.4f} +/- {EW_err:.4f}  '
      f'(digitized {EW_digitized:.4f} +/- {EW_digitized_err:.4f})')
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

# the run solved at the well-mixed crossing 0.4132 on the current binary;
# no turbulence run exists for it.
c = exhale_curve(EWMATCH_TAG)
ax.plot(c[0] + dlam_air, 1.0 - c[1]/100, color='C2', lw=1.5,
        label=r'EXHALE, $\mathrm{He/H}=0.4132$ (EW-matched)')

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
_EWD = os.path.join('exhale', EWMATCH_TAG, 'output')
hy = loadtxt_cells(os.path.join(_EWD, 'Hydro_ioniz.txt'))
lab_i = [l for l in open(os.path.join(_EWD, 'Ion_species.txt'))
         if l.startswith('# columns')][0].split()[2:]
ki = {n: j for j, n in enumerate(lab_i)}
io = loadtxt_cells(os.path.join(_EWD, 'Ion_species.txt'))
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
print('EW-matched solution (He/H=0.4132): T_max %.0f K at r=%.2f; '
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

lam0, exc0 = exhale_curve(EWMATCH_TAG)
f_match = matched_kernel(EWMATCH_TAG)
print('matched kernel (He/H=0.4132, no turbulence): %.2f km/s '
      '(sigma %.2f km/s)' % (f_match, f_match/2.35482))
for tag, sub, lab in (('crossings_j96/wm_heh0p40', '', 'He/H=0.40'),
                      ('crossings_j96/wm_heh0p46', '', 'He/H=0.46')):
    print('  same measurement, %-26s %.2f km/s' % (lab, matched_kernel(tag, sub)))
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
       label=r'EXHALE, $\mathrm{He/H}=0.4132$, as solved')
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
# The EW crossing moves from He/H = 0.4132 to 1.6108 once binary H/He element
# diffusion is on at He_Kzz = 1e9 (exhale/crossings_j96/results.txt Sect. 3).
# Repeat the kernel measurement on the EW-matched diffusive solution and
# compare the matched kernel with the diffusion-off one.  Both arms are the
# runs solved AT their own crossing (crossings_j96 Sect. 6).
DIFF_TAG = 'crossings_j96/cf_K1e9_P'
DIFF_HEH = 1.6108
lamd, excd = exhale_curve(DIFF_TAG)
f_match_diff = matched_kernel_curve(lamd, excd)
print('matched kernel (%s, He/H=%.2f, He_Kzz=1e9): %.2f km/s (sigma %.2f), '
      'diffusion-off He/H=0.4132 gives %.2f km/s, difference %.2f km/s'
      % (DIFF_TAG, DIFF_HEH, f_match_diff, f_match_diff/2.35482, f_match,
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
data_diff = loadtxt_cells(p_ion_diff)

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
ax.text(13.8, 1e-2, r'$R_\star = 13.6\,R_{\rm p}$', fontsize=6.5, color='0.4', rotation=90)

ax.set_xlim(1.0, 20.0)
ax.set_ylim(1e-4, 3e13)
ax.set_xlabel(r'radius $r$ [$R_{\rm p}$]')
ax.set_ylabel(r'number density $n$ [cm$^{-3}$]')
ax.grid(alpha=0.2, which='both')
ax.legend(fontsize=6.2, loc='upper right', framealpha=0.9, ncol=2)
ax.set_title(r'(a) H and He number densities ($K_{zz}=10^9$, $\mathrm{He/H}=1.6108$)', fontsize=8)

# Panel (b): Broadened transit profile vs observation
ax2 = axs[1]
ax2.errorbar(o_air, o_flux, yerr=o_fsig, fmt='ko', ms=2.0, lw=0.6, capsize=0,
             zorder=3, label=r'LHS\,1140\,b (2024 GP-corrected)')
ax2.axhline(1.0 - 0.006, color='0.6', ls='--', lw=0.8, label=r'2025 detection limit')
ax2.axhline(1.0, color='0.9', lw=0.6, zorder=0)

ax2.plot(lamd + dlam_air, 1.0 - excd/100, color='C4', lw=1.0, ls=':',
         label=r'diffusion ($K_{zz}=10^{9}$, $\mathrm{He/H}=1.6108$), as solved')
ax2.plot(lamd + dlam_air, 1.0 - broaden(lamd, excd, f_match_diff)/100,
         color='C4', lw=1.5,
         label=r'same, $+' + '%.1f' % f_match_diff + r'$\,km\,s$^{-1}$ FWHM Gaussian')
ax2.plot(lam0 + dlam_air, 1.0 - broaden(lam0, exc0, f_match)/100, color='C2',
         lw=1.0, ls='--', dashes=(4, 2),
         label=r'diffusion off ($\mathrm{He/H}=0.4132$), $+' + '%.1f' % f_match + r'$\,km\,s$^{-1}$')

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
        lw=1.0, label=r'EXHALE, $\mathrm{He/H}=0.425$, $+'
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
# Panel (b) reads both crossings from ladders solved on the current binary
# (Update_EXHALE Sect. 96): the GJ 1132 well-mixed one of Fig. ew
# (crossings_j96 lineage W) and the GJ 699 one of crossings_j96 Sect. 2
# (lineage G, of record).  Panel (a) draws the arm of that same ladder
# nearest the crossing 0.0483, g_heh0p048, and the two solar curves it
# carries for contrast are crossings_j96/{sol1132,sol699}, the runs of
# exhale/misc_gm25/{solar,gj699_solar} re-solved on the current binary, so
# that every curve in the figure stands on the same physics.
LADDER_699 = [('crossings_j96/g_heh0p042', 0.042),
              ('crossings_j96/g_heh0p044', 0.044),
              ('crossings_j96/g_heh0p045', 0.045),
              ('crossings_j96/g_heh0p046', 0.046),
              ('crossings_j96/g_heh0p047', 0.047),
              ('crossings_j96/g_heh0p048', 0.048),
              ('crossings_j96/g_heh0p049', 0.049),
              ('crossings_j96/g_heh0p050', 0.050),
              ('crossings_j96/g_heh0p052', 0.052)]
hh9 = np.array([c[1] for c in LADDER_699])
ew9 = np.array([red_ew(c[0]) for c in LADDER_699])
cross9 = ew_crossing(hh9, ew9, EW_obs)
cross9_lo = ew_crossing(hh9, ew9, EW_obs - EW_err)
cross9_hi = ew_crossing(hh9, ew9, EW_obs + EW_err)

fig, axs = plt.subplots(1, 2, figsize=(7.1, 3.0))

a = axs[0]
a.errorbar(o_air, o_flux, yerr=o_fsig, fmt='ko', ms=2.0, lw=0.6, capsize=0,
           zorder=3, label=r'LHS\,1140\,b, 2024 (GP-corrected)')
a.axhline(1.0 - 0.006, color='0.6', ls='--', lw=0.8,
          label=r'2025 detection limit')
a.axhline(1.0, color='0.9', lw=0.6, zorder=0)
c1132 = exhale_curve('crossings_j96/sol1132')
a.plot(c1132[0] + dlam_air, 1.0 - c1132[1]/100, color='C1', lw=0.9, ls=':',
       label=r'GJ\,1132 SED, solar (red $'
       + '%.2f' % fit_metrics(c1132[0], c1132[1], frame='air')['red_depth']
       + r'\%$)')
# zorder puts the shallower arm on top: the two overlap over the core.
for tag, col, zo, name in (('crossings_j96/g_heh0p048', 'C4', 2.4,
                            r'$\mathrm{He/H}=0.048$ (EW match)'),
                           ('crossings_j96/sol699', 'C1', 2.2,
                            r'solar, $\mathrm{He/H}=0.083$')):
    c = exhale_curve(tag)
    a.plot(c[0] + dlam_air, 1.0 - c[1]/100, color=col, lw=1.4, zorder=zo,
           label=r'GJ\,699 SED, ' + name + r' (red $'
           + '%.2f' % fit_metrics(c[0], c[1], frame='air')['red_depth']
           + r'\%$)')
a.set_xlim(10827, 10831.7); a.set_ylim(0.982, 1.006)
a.ticklabel_format(axis='x', useOffset=False, style='plain')
a.set_xlabel(r'air wavelength [\AA]'); a.set_ylabel(r'normalized flux')
a.grid(alpha=0.18); a.legend(fontsize=5.8, loc='lower left', framealpha=0.9)
a.set_title(r'(a) five times the XUV, no turbulence term', fontsize=8)

a = axs[1]
a.axhspan(EW_obs - EW_err, EW_obs + EW_err, color='0.85', zorder=0,
          label=r'observed')
a.plot(hw, ew_w, 's-', color='C0', ms=3.5, lw=1.3, label=r'GJ\,1132 SED')
a.plot(hh9, ew9, 'o-', color='C3', ms=3.5, lw=1.3, label=r'GJ\,699 SED')
a.axvline(wm_cross, color='C0', ls=':', lw=0.9)
a.axvline(cross9, color='C3', ls=':', lw=0.9)
a.plot([wm_cross, cross9], [EW_obs]*2, marker='*', ms=9, ls='none',
       color='0.25', zorder=5)
a.text(cross9*1.07, 0.93, r'$' + '%.4f' % cross9 + r'$', color='C3',
       fontsize=6.5, transform=a.get_xaxis_transform())
a.text(wm_cross*0.94, 0.93, r'$' + '%.4f' % wm_cross + r'$', color='C0',
       fontsize=6.5, ha='right', transform=a.get_xaxis_transform())
a.set_xscale('log')
a.set_xlim(0.035, 0.62)
a.set_ylim(0.95, 1.23)
a.set_xticks([0.04, 0.06, 0.1, 0.2, 0.4])
a.set_xticklabels([r'0.04', r'0.06', r'0.1', r'0.2', r'0.4'])
a.set_xticks([], minor=True)
a.set_xlabel(r'$\mathrm{He/H}$ number ratio')
a.set_ylabel(r'red-pair EW [\%\,\AA]')
a.grid(alpha=0.22, which='both')
a.legend(fontsize=6, loc='lower right')
a.set_title(r'(b) the composition estimate follows the SED', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_gj699.pdf')
plt.close()
print('GJ 699 crossing (no turbulence): He/H = %.4f (1 sigma %.4f - %.4f)  '
      'H:He = %.2f  %.2fx solar; GJ 1132 crossing %.4f, ratio %.1f'
      % (cross9, cross9_lo, cross9_hi, 1/cross9, cross9/0.0833333,
         wm_cross, wm_cross/cross9))
print('GJ 699 matched kernel (He/H=0.06, retired coefficient): %.2f km/s'
      % matched_kernel('heh0p06_gj699'))
print('GJ 699 matched kernel (He/H=0.048, current binary): %.2f km/s'
      % matched_kernel('crossings_j96/g_heh0p048'))
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
a.axvline(f_match/2.35482, color='0.35', ls='--', lw=1.0)
a.text(f_match/2.35482 - 0.3, 68, r'$\sigma$ demanded by the line', rotation=90,
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
# The four arms are exhale/refresh_j96/misc/kzz*, the runs of
# exhale/heh0p55_diff_* re-solved on the current binary.  The K_zz = 0 arm no
# longer carries the defect of ../docs/Update_EXHALE.md Sect. 88.4(a): Sect. 96
# removed the divergent division, and its metastable curve is drawn for its
# value.
KZZ_RUNS = [('refresh_j96/misc/kzz0',    0.0,   r'$K_{zz} = 0$'),
            ('refresh_j96/misc/kzz1e8',  1.0e8, r'$K_{zz} = 10^{8}$'),
            ('refresh_j96/misc/kzz1e9',  1.0e9, r'$K_{zz} = 10^{9}$ (adopted)'),
            ('refresh_j96/misc/kzz1e10', 1.0e10, r'$K_{zz} = 10^{10}$')]
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
    hy = loadtxt_cells(os.path.join(d, 'output', 'Hydro_ioniz.txt'))
    p = os.path.join(d, 'output', 'Ion_species_adv.txt')
    names = [l for l in open(p) if l.startswith('# columns')][0].split()[2:]
    a = loadtxt_cells(p)
    j = {n: i for i, n in enumerate(names)}
    ra = a[:, j['r[Rp]']]
    # He I already contains the 2^3S metastable (it is a level of He I, not
    # a separate species -- species_table.f90), so the elemental sum must
    # NOT add the HeITR column on top of it.
    nHe = sum(a[:, j[s]] for s in ('HeI', 'HeII', 'HeIII'))
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
# Provenance of the runs: exhale/crossings_j96/results.txt section 4, every
# decade re-bracketed on the current binary (the composition-renormalized
# SvS85 split of Update_EXHALE Sect. 96); K_zz = 1e9 is section 3 of the
# same file.
CROSS_SCANS = [
    (0.0, ['crossings_j96/k0_heh3p50', 'crossings_j96/k0_heh3p64',
           'crossings_j96/k0_heh3p70', 'crossings_j96/k0_heh3p75',
           'crossings_j96/k0_heh3p79', 'crossings_j96/k0_heh3p85',
           'crossings_j96/k0_heh3p93']),
    (1.0e5, ['crossings_j96/kz1e5_heh3p35', 'crossings_j96/kz1e5_heh3p50',
             'crossings_j96/kz1e5_heh3p64', 'crossings_j96/kz1e5_heh3p79',
             'crossings_j96/kz1e5_heh3p93']),
    (1.0e6, ['crossings_j96/kz1e6_heh3p19', 'crossings_j96/kz1e6_heh3p32',
             'crossings_j96/kz1e6_heh3p46', 'crossings_j96/kz1e6_heh3p60',
             'crossings_j96/kz1e6_heh3p74']),
    (1.0e7, ['crossings_j96/kz1e7_heh2p70', 'crossings_j96/kz1e7_heh2p82',
             'crossings_j96/kz1e7_heh2p94', 'crossings_j96/kz1e7_heh3p06',
             'crossings_j96/kz1e7_heh3p18']),
    (1.0e8, ['crossings_j96/kz1e8_heh2p05', 'crossings_j96/kz1e8_heh2p15',
             'crossings_j96/kz1e8_heh2p23', 'crossings_j96/kz1e8_heh2p32',
             'crossings_j96/kz1e8_heh2p41']),
    (1.0e9, ['crossings_j96/k9_heh1p50', 'crossings_j96/k9_heh1p55',
             'crossings_j96/k9_heh1p575', 'crossings_j96/k9_heh1p60',
             'crossings_j96/k9_heh1p625', 'crossings_j96/k9_heh1p65',
             'crossings_j96/k9_heh1p70', 'crossings_j96/k9_heh1p75']),
    (1.0e10, ['crossings_j96/kz1e10_heh1p06', 'crossings_j96/kz1e10_heh1p10',
              'crossings_j96/kz1e10_heh1p15', 'crossings_j96/kz1e10_heh1p20',
              'crossings_j96/kz1e10_heh1p25', 'crossings_j96/kz1e10_heh1p29']),
    (1.0e11, ['crossings_j96/kz1e11_heh0p795', 'crossings_j96/kz1e11_heh0p83',
              'crossings_j96/kz1e11_heh0p865', 'crossings_j96/kz1e11_heh0p90',
              'crossings_j96/kz1e11_heh0p93']),
]

# The decades below 1e5 are not scanned in composition.  At the fixed
# composition He/H = 3.8075, the K_zz = 0 crossing, their equivalent width
# differs from the K_zz = 0 one by far less than the measurement error, which
# is the measured statement the figure draws as a band rather than as points
# of their own.
FLAT_PROBE = [(1.0e1, 'crossings_j96/probe1e1_heh3p8075'),
              (1.0e2, 'crossings_j96/probe1e2_heh3p8075'),
              (1.0e3, 'crossings_j96/probe1e3_heh3p8075'),
              (1.0e4, 'crossings_j96/probe1e4_heh3p8075')]
FLAT_REF = 'crossings_j96/probe0_heh3p8075'
HEH_WELLMIXED = 0.4132        # crossing with the operator off, Fig. ew


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
    print('flat  K_zz=%-8.0e He/H = %.4f: EW = %.6f, %+.3e vs K_zz = 0 '
          '(%.3f sigma)'
          % (kzz, key_value(os.path.join('exhale', tag, 'input.inp'),
                            'He/H number ratio'),
             red_ew(tag), red_ew(tag) - ew_ref,
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
ax.text(6.0, HEH_WELLMIXED*1.06,
        r'well-mixed limit (operator off), He/H $= '
        + '%.4f' % HEH_WELLMIXED + r'$',
        fontsize=6.5, color='C3', va='bottom')

ax.axvline(1.0e9, color='C0', ls='-.', lw=1.0, alpha=0.7)
ax.text(1.0e9/1.5, 0.50, r'adopted', fontsize=6.5, color='C0',
        rotation=90, va='bottom', ha='right')

kxp = np.array(kx[1:])
kcp = np.array(kc[1:])
ax.errorbar(kxp, kcp,
            yerr=[kcp - np.array(klo[1:]), np.array(khi[1:]) - kcp],
            fmt='s', color='C0', ms=4, lw=1.0, capsize=2,
            label=r'EW-matched He/H, $\pm1\sigma$')

pw = np.array([1.0e8, 1.0e11])
ax.plot(pw, kc[5]*(pw/1.0e9)**(-0.145), color='C1', lw=1.1, alpha=0.8,
        label=r'$\propto K_{zz}^{-0.145}$')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(X0, X1)
ax.set_ylim(0.38, 8.0)
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
# is a solution and not an input.  Runs: exhale/refresh_j96/closure/{rw_ref,
# rw_lo,rw_hi}, the converged iterates of the closure re-run end to end,
# re-solved
# to their own fixed point on the current Penning coefficient with the
# photochemical columns held fixed.
# Record: ../docs/Update_EXHALE.md section 79; exhale/misc_gm25/results.txt.
CLOSURE_ARMS = [('refresh_j96/closure/rw_ref', '', r'$1.0\times$ start'),
                ('refresh_j96/closure/rw_lo', '', r'$0.3\times$ start'),
                ('refresh_j96/closure/rw_hi', '', r'$3.0\times$ start')]
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

# The same four runs Table tab:validity reports, all on the current binary.
KN_CASES = [('exhale/refresh_j96/scan/heh0p55', 'C0', '-',
             r'well mixed, He/H\,$=0.55$'),
            ('exhale/refresh_j96/thermostat/ms2p13_r1e4', 'C2', '--',
             r'diffusion, He/H\,$=2.13$'),
            ('exhale/refresh_j96/thermostat/fc2p09', 'C1', ':',
             r'flux-closed'),
            ('exhale/refresh_j96/validity/hd209_control', 'C3', '-.',
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
# exhale/crossings_j96: nineteen reservoirs on the current binary
# (Update_EXHALE Sect. 96), every one closed in full and then re-solved on
# its own output with the profile held fixed, so that arms which stopped at
# different k are read at their common wind fixed point (the rw_ directories,
# which is where every number the memo quotes for this ladder is read).  The
# 8.20 and 8.25 arms straddle the measurement, so the crossing is
# interpolated inside a measured bracket, and so are both edges of its
# 1 sigma interval.  Record: exhale/crossings_j96/results.txt.
# exhale/flux_closure/heh12 is from the earlier ladder: the same reservoir
# reached by a 1.50x jump from the 8.0 arm instead of a 1.08x step,
# converging to the same base and the same matching-level composition but to
# a hotter, more extended outer region, and is plotted apart as discarded
# (docs/Update_EXHALE.md Sect. 83); its line is re-solved on the current
# binary in crossings_j96/L12p0jump.  (b) The
# red-pair depth against the XUV scaling, from both models that reproduce the
# 2024 equivalent width: the scalar base at He/H = 2.13 and the flux-closed
# solution at He/H = 9.71.  The lower atmosphere is held fixed across the
# grid, so it isolates the wind's response.  Record: ../kzz_decision.md
# section 8.
CLOSURE_LADDER = [('crossings_j96/rwL_7p0', 7.00),
                  ('crossings_j96/rwL_7p4', 7.40),
                  ('crossings_j96/rwL_7p5', 7.50),
                  ('crossings_j96/rwL_7p7', 7.70),
                  ('crossings_j96/rw_8p0', 8.00),
                  ('crossings_j96/rw_8p1', 8.10),
                  ('crossings_j96/rw_8p2', 8.20),
                  ('crossings_j96/rw_8p25', 8.25),
                  ('crossings_j96/rw_8p3', 8.30),
                  ('crossings_j96/rw_8p45', 8.45),
                  ('crossings_j96/rw_8p5', 8.50),
                  ('crossings_j96/rw_8p6', 8.60),
                  ('crossings_j96/rw_8p7', 8.70),
                  ('crossings_j96/rw_8p75', 8.75),
                  ('crossings_j96/rw_8p9', 8.90),
                  ('crossings_j96/rw_9p0', 9.00),
                  ('crossings_j96/rw_9p1', 9.10),
                  ('crossings_j96/rw_9p3', 9.30),
                  ('crossings_j96/rw_9p5', 9.50)]
CLOSURE_DISCARDED = 'flux_closure/heh12'   # seeded across a 1.50x jump
# He/H at the match is chemistry and stands; the line is read off the same
# arm re-solved on the current binary.
CLOSURE_DISCARDED_LINE = 'crossings_j96/L12p0jump'
XUV_LIMIT = 0.6                       # 2025 non-detection, per cent in depth
# Both XUV families are exhale/xuvkzz_gm25/{S,C}*, every arm re-solved on the
# current Penning coefficient (exhale/xuvkzz_gm25/results.txt).  The scaled
# points are keyed by the fraction; the closed family has no 0.30 arm.
XUV_GRID = [0.01, 0.10, 0.15, 0.20, 0.25, 0.30, 0.33]
XUV_FAMILIES = [('S', 'refresh_j96/xuvkzz/S1p00', 2.13,
                 r'scalar base, He/H $= 2.13$', 'C0', 'o', '-'),
                ('C', 'refresh_j96/xuvkzz/C1p00', 9.71,
                 r'flux-closed, He/H $= 9.71$', 'C3', 's', '--')]
# Where each model's depth crosses the 0.6 per cent limit, read linearly in
# both variables on the scaled grid alone
# (exhale/refresh_j96/xuvkzz/results.txt).
XUV_CROSSINGS = [(0.326, 'C0'), (0.368, 'C3')]


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
for tag, heh in CLOSURE_LADDER:
    lad_h.append(heh)
    lad_e.append(red_ew(tag))
    print('ladder %-28s He/H = %8.4f  EW = %.4f %%A'
          % (tag, heh, lad_e[-1]))
lad_h, lad_e = np.array(lad_h), np.array(lad_e)
o = np.argsort(lad_h)
lad_h, lad_e = lad_h[o], lad_e[o]

dsc_sub, dsc_h = closure_last_k(CLOSURE_DISCARDED)
dsc_e = red_ew(CLOSURE_DISCARDED_LINE)
print('discarded %-18s %-4s He/H = %8.4f  EW = %.4f %%A'
      % (CLOSURE_DISCARDED, dsc_sub, dsc_h, dsc_e))


def heh_at_closed_ew(target):
    """Reservoir whose flux-closed EW equals target, log-log on the ladder.

    Piecewise linear in log-log between the scanned rungs.  The measured EW
    and both edges of its 1 sigma interval fall between solved rungs, so
    nothing here is extrapolated; the continuation is kept for the case of a
    target outside the ladder.
    """
    lx, ly, lt = np.log10(lad_h), np.log10(lad_e), np.log10(target)
    if ly[0] <= lt <= ly[-1]:
        return 10.0**brentq(lambda t: np.interp(t, lx, ly) - lt,
                            lx[0], lx[-1])
    i, j = (-2, -1) if lt > ly[-1] else (0, 1)
    return 10.0**(lx[i] + (lt - ly[i])*(lx[j] - lx[i])/(ly[j] - ly[i]))


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
a.text(2.2, EW_obs*1.035, r'measured, $1.108 \pm 0.030$', fontsize=6.5,
       color='0.35', va='bottom')
a.plot(lad_h, lad_e, marker='o', ms=4, lw=1.2, color='C0',
       label=r'flux-closed ladder')
a.plot([dsc_h], [dsc_e], marker='o', ms=6, ls='none', mfc='none',
       mec='0.55', mew=1.2, zorder=3,
       label=r'discarded, seeded across a jump')
a.plot([heh_closed], [EW_obs], marker='*', ms=11, ls='none', color='C3',
       zorder=4, label=r'crossing, He/H $= %.2f$' % heh_closed)
a.errorbar([heh_closed], [EW_obs],
           xerr=[[heh_closed - heh_clo_lo], [heh_clo_hi - heh_closed]],
           fmt='none', ecolor='C3', lw=1.0, capsize=2, zorder=4)
a.axvline(2.0924, color='C2', ls=':', lw=1.0)
a.text(2.0924*0.95, 1.01, r'closure at the well-mixed reservoir',
       fontsize=6.2, color='C2', rotation=90, va='bottom', ha='right')
a.set_xscale('log'); a.set_yscale('log')
a.set_xlim(1.8, 16.0); a.set_ylim(1.00, 1.45)
a.set_xticks([2, 3, 5, 8, 9, 11, 14])
a.set_xticklabels([r'2', r'3', r'5', r'8', r'9', r'11', r'14'])
a.set_xticks([], minor=True)
a.set_yticks([1.0, 1.1, 1.2, 1.3, 1.4])
a.set_yticklabels([r'1.0', r'1.1', r'1.2', r'1.3', r'1.4'])
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
    xs, ds = [1.0], [red_depth(fid)]
    for f in XUV_GRID:
        tag = 'refresh_j96/xuvkzz/%s%s' % (key, ('%.2f' % f).replace('.', 'p'))
        d = red_depth(tag)
        if np.isfinite(d):
            xs.append(f); ds.append(d)
    xs, ds = np.array(xs), np.array(ds)
    o = np.argsort(xs)
    a.plot(xs[o], ds[o], marker=mk, ms=4, lw=1.2, ls=ls, color=col,
           label=lab)
    print('xuv %-12s ' % key + '  '.join('%.2f:%.3f' % (x, d)
                                         for x, d in zip(xs[o], ds[o])))
for xc, cc in XUV_CROSSINGS:
    a.axvline(xc, color=cc, ls=':', lw=1.0)
a.text(0.326*0.94, 4.0e-4, r'$0.33\times$', fontsize=6.2, color='C0',
       rotation=90, va='bottom', ha='right')
a.text(0.368*1.06, 4.0e-4, r'$0.37\times$', fontsize=6.2, color='C3',
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
# essentially the same reservoir, and the closed solution at 10.31, a rung of
# the closure ladder (the crossing later moved to 11.73, straddled by the
# 11.11 and 12.01 rungs; 10.31 is kept here because the three curves are the
# ones the thermostat argument was measured on). All three carry H/He
# element diffusion at
# K_zz = 1e9 cm^2/s. Profiles are the advection-corrected ones the transit
# tool consumes.

THERMO_CASES = [
    ('refresh_j96/thermostat/ms2p13_r1e4',
     r'scalar base, no metals, He/H $= 2.13$', 'C0', '-'),
    ('refresh_j96/thermostat/fc2p09',
     r'photochemical base, He/H $= 2.09$', 'C1', '--'),
    ('refresh_j96/thermostat/fc10p3',
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
        return cols, loadtxt_cells(path)

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
    ('refresh_j96/thermostat/ms2p13_r1e4', 2.13,
     r'scalar base, no metals: $2.13$', 'C0', '-'),
    ('refresh_j96/thermostat/fc2p09', 2.0924,
     r'photochemical base: $2.09$', 'C1', '--'),
    ('refresh_j96/thermostat/fc10p3', 10.3116,
     r'photochemical base: $10.31$', 'C3', '-.'),
]
COMP_METALS = ['C', 'O', 'N', 'Mg', 'Si', 'Ca', 'Na', 'K', 'S', 'Fe']


def composition_profile(sub):
    """Nuclei fractions and ionization fractions of one run, against r."""
    def read(path):
        with open(path) as fh:
            fh.readline()
            cols = fh.readline().split()[2:]
        return cols, loadtxt_cells(path)

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

# ==== Figure: density and hydrogen ionization, EXHALE against p-winds =======
# The plan's Phase F item 1 asks for five profiles at the prescribed
# H:He = 1e-3 -- density, velocity, temperature, ionization and metastable
# helium.  The three-panel structure figure above carries velocity,
# temperature and the metastable population; this one carries the remaining
# two, from the same pair of solutions (exhale/heh1000, He/H = 1000, against
# the matched_gj1132 oracle).  The dotted curve in panel (b) is
# photoionization equilibrium evaluated on p-winds' own density and
# temperature, so the gap between it and the solid p-winds curve is how far
# out of equilibrium the imposed Parker flow holds hydrogen.
M_H = 1.6726219e-24
H_PLANCK, C_LIGHT, EV_ERG = 6.62607015e-27, 2.99792458e10, 1.602176634e-12
R_JUP = 6.9911e9


def h_photoionization_rate(sed_file):
    """Unattenuated H I photoionization rate [1/s] from an SED at the planet.

    Verner et al. (1996) ground-state cross section, integrated against the
    photon flux below 911.267 A.  Reproduces the rate p-winds computes
    internally for the same file to 0.1 %.
    """
    d = np.loadtxt(sed_file)
    lam, flx = d[:, 0], d[:, 1]                       # [A], [erg/s/cm2/A]
    m = (lam > 0.0) & (lam <= 911.267)
    lam, flx = lam[m], flx[m]
    e_erg = H_PLANCK*C_LIGHT/(lam*1e-8)
    e_ev = e_erg/EV_ERG
    e0, s0, ya, pp = 4.298e-1, 5.475e4, 3.288e1, 2.963
    x = e_ev/e0
    fy = (x - 1.0)**2*x**(0.5*pp - 5.5)*(1.0 + np.sqrt(x/ya))**(-pp)
    sig = np.where(e_ev >= 13.6, s0*fy*1e-18, 0.0)
    return np.trapz(flx/e_erg*sig, lam)


hyd_e = loadtxt_cells('exhale/heh1000/output/Hydro_ioniz.txt')
lab_i = [l for l in open('exhale/heh1000/output/Ion_species.txt')
         if l.startswith('# columns')][0].split()[2:]
ki = {n: j for j, n in enumerate(lab_i)}
ion_e = loadtxt_cells('exhale/heh1000/output/Ion_species.txt')
pw = np.loadtxt('pwinds_oracle/profile_matched_gj1132.txt')
r_p, v_p, rho_p, fhii_p = pw[:, 0], pw[:, 1], pw[:, 2], pw[:, 3]

r_e = hyd_e[:, 0]
rho_e = hyd_e[:, 1]*M_H                                # [g/cm3]
nH_e = ion_e[:, ki['HI']] + ion_e[:, ki['HII']]
xhii_e = ion_e[:, ki['HII']]/np.maximum(nH_e, 1e-300)

# photoionization equilibrium on the p-winds structure: with electrons from
# hydrogen alone (p-winds' own closure) phi (1-x) = alpha n_H x^2
T_PW = 6100.0
H_FRAC = 1.0e-3/(1.0 + 1.0e-3)
phi_H = h_photoionization_rate('sed/lhs1140_sed_gj1132_at_b.txt')
alpha_B = 2.59e-13*(T_PW/1e4)**(-0.7)                  # p-winds' case B
nH_p = H_FRAC*rho_p/(M_H*(H_FRAC + 4.0*(1.0 - H_FRAC)))
q = alpha_B*nH_p
xeq_p = (-phi_H + np.sqrt(phi_H**2 + 4.0*phi_H*q))/(2.0*q)

fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.6))
ax[0].semilogy(r_e, rho_e, 'C3', lw=1.3, label=r'EXHALE')
ax[0].semilogy(r_p, rho_p, 'C0', lw=1.3, label=r'p-winds')
ax[0].set_ylabel(r'$\rho$ [g\,cm$^{-3}$]')
ax[0].set_ylim(1e-20, 1e-11)
ax[0].set_title(r'(a) mass density', fontsize=8)
ax[1].semilogy(r_e, np.maximum(xhii_e, 1e-4), 'C3', lw=1.3, label=r'EXHALE')
ax[1].semilogy(r_p, fhii_p, 'C0', lw=1.3, label=r'p-winds')
ax[1].semilogy(r_p, xeq_p, 'C0', lw=1.0, ls=':',
               label=r'p-winds structure, photoion.\ equilibrium')
ax[1].set_ylabel(r'$x({\rm H\,II})$')
ax[1].set_ylim(1e-3, 2.0)
ax[1].set_title(r'(b) hydrogen ionization', fontsize=8)
for a in ax:
    a.set_xlim(1, 20)
    a.set_xlabel(r'$r$ [$R_{\rm p}$]')
    a.axvline(13.6, color='0.7', ls=':', lw=0.8)
    a.grid(alpha=0.2)
    a.legend(fontsize=6)
plt.tight_layout()
plt.savefig(f'{OUT}/lhs1140b_structure_rho_ion.pdf')
plt.close()

_at = lambda x, y, rr: np.interp(rr, x, y)
print('EXHALE (He/H = 1000) against p-winds matched_gj1132:')
print('  phi_H(unattenuated) = %.3e /s  ->  t_ion = %.3e s = %.0f d'
      % (phi_H, 1.0/phi_H, 1.0/phi_H/86400.0))
for rr in (1.2, 2.0, 5.0, 10.0, 20.0):
    print('  r = %5.1f Rp: rho %9.3e / %9.3e = %6.2f ; '
          'x(HII) %7.4f / %7.4f (equilibrium on the p-winds structure %7.4f)'
          % (rr, _at(r_e, rho_e, rr), _at(r_p, rho_p, rr),
             _at(r_e, rho_e, rr)/_at(r_p, rho_p, rr),
             _at(ion_e[:, 0], xhii_e, rr), _at(r_p, fhii_p, rr),
             _at(r_p, xeq_p, rr)))
print('  base cell: rho %9.3e / %9.3e = %6.1f (the two lower boundary '
      'conditions, not a profile difference)' % (rho_e[0], rho_p[0],
                                                 rho_e[0]/rho_p[0]))
sel = (r_p >= 1.05)
lr = np.log10(np.interp(r_p[sel], r_e, rho_e)/rho_p[sel])
print('  1.05-20 Rp, log10 rho(EXHALE/p-winds): %+.2f at r = %.2f to '
      '%+.2f at r = %.2f'
      % (lr.min(), r_p[sel][np.argmin(lr)], lr.max(), r_p[sel][np.argmax(lr)]))
m3 = (ion_e[:, 0] >= 1.0) & (ion_e[:, 0] <= 20.0)
Rp_cm = 0.157692*R_JUP
col_e = np.trapz(ion_e[m3, ki['HeITR']], ion_e[m3, 0]*Rp_cm)
col_p = np.trapz(pw[:, 5], r_p*Rp_cm)
print('  radial He 2^3S column 1-20 Rp: EXHALE %.3e, p-winds %.3e, ratio %.3f'
      % (col_e, col_p, col_e/col_p))
print('wrote lhs1140b_structure_rho_ion.pdf')
# ==== END density / ionization block ========================================
