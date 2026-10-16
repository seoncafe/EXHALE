#!/usr/bin/env python3
"""Stellar power the wind absorbs, by band, against what the column's two
models take at their top (Q1, part 3).

The spectrum is photochem_work/flux_photochem.txt of k05 (the SED at the
planet, erg cm^-2 s^-1 nm^-1; the wind reads the same SED file with dayside
dilution 1).  The wind's optical depth along the radial (substellar) beam,
from the face R0 to the top, is

  tau(lambda) = sum_X sigma_X(lambda) N_X,   N_X = int n_X dr  (Ion_species.txt)

with sigma for H, He and H2 the total photoabsorption cross sections of
photochem_clima_data (the ones Photochem itself attenuates with), He II and
C, N, O neutral and singly ionized outer-shell fits of Verner et al. (1996)
(references/verner_photo/photo.dat).  The H2 photoabsorption of the table is
smoothed over the Lyman and Werner lines and is therefore NOT used for the
line band 91.2-120.1 nm: there the absorbed power is the code's own Solomon
pumping, k_LW n_H2 / <p_diss> <h nu> (Lyman_Werner.txt; <p_diss> = 0.135 of
Draine & Bertoldi 1996 Table 2 and <h nu> = 11.7354 eV as lyman_werner.f90
states them).  Lyman alpha (121.0-122.2 nm) is listed with the H I column
only; its core is scattered, not absorbed, by the wind.

Usage: python3 stellar_bands_absorbed_in_wind.py <k05 dir> <xsections dir> <verner photo.dat>
"""
import os
import sys
import numpy as np
import h5py

run, xdir, vfile = sys.argv[1:4]
R_J = 7.1492e9
HC_eVnm = 1239.84198
EV = 1.602176634e-12

wl, F = np.loadtxt(os.path.join(run, 'photochem_work', 'flux_photochem.txt'),
                   skiprows=1, unpack=True)       # nm, erg cm^-2 s^-1 nm^-1
res = {}
for line in open(os.path.join(run, 'EXHALE_resolved.out')):
    t = line.split()
    if len(t) >= 2 and not line.startswith('#'):
        res[t[0]] = t[1]
Rp = float(res['planet_radius_RJ']) * R_J
R0 = float(res['base_face_radius_cm'])
hdr = [l for l in open(os.path.join(run, 'output', 'Ion_species.txt'))
       if l.startswith('# columns')][0].split()[2:]
ion = np.loadtxt(os.path.join(run, 'output', 'Ion_species.txt'))[2:-2]
r = ion[:, 0] * Rp
f = np.empty(r.size + 1)
f[0] = R0
f[1:-1] = 0.5 * (r[1:] + r[:-1])
f[-1] = 2 * r[-1] - f[-2]
dr = np.diff(f)
N = {k: np.sum(ion[:, hdr.index(k)] * dr) for k in hdr[1:]}
print('radial columns above R0 [cm^-2]: ' + ', '.join(
    '%s %.3e' % (k, N[k]) for k in ('HI', 'HII', 'HeI', 'HeII', 'HeITR', 'H2', 'H3p',
                                   'CI', 'NI', 'OI', 'CII', 'NII', 'OII')))


def tab_xs(sp):
    with h5py.File(os.path.join(xdir, sp + '.h5'), 'r') as h:
        w = np.array(h['wavelengths'])
        a = np.array(h['photoabsorption'])
    return np.interp(wl, w, a, left=0.0, right=0.0)


vt = np.loadtxt(vfile)


def verner(Z, Ne):
    row = vt[(vt[:, 0] == Z) & (vt[:, 1] == Ne)][0]
    Eth, Emax, E0, s0, ya, P, yw, y0, y1 = row[2:]
    E = HC_eVnm / wl
    xx = E / E0 - y0
    y = np.sqrt(xx**2 + y1**2)
    Fy = ((xx - 1)**2 + yw**2) * y**(0.5 * P - 5.5) * (1 + np.sqrt(y / ya))**(-P)
    return np.where(E >= Eth, s0 * Fy * 1e-18, 0.0)


sig = {'HI': tab_xs('H'), 'HeI': tab_xs('He'), 'H2': tab_xs('H2'),
       'HeII': verner(2, 1), 'CI': verner(6, 6), 'NI': verner(7, 7),
       'OI': verner(8, 8), 'CII': verner(6, 5), 'NII': verner(7, 6),
       'OII': verner(8, 7)}
lw_lines = (wl >= 91.2) & (wl < 120.1)
tau = np.zeros_like(wl)
# the H and H2 tables carry the Lyman series and the Lyman/Werner lines
# smoothed; longward of 91.2 nm neither is a continuum, so both are dropped
# there (the H2 lines are the code's pumping below, Ly alpha is separate)
for k, s in sig.items():
    if k in ('H2', 'HI'):
        s = np.where(wl >= 91.2, 0.0, s)
    tau += s * N[k]
absorbed = F * (1 - np.exp(-np.minimum(tau, 700)))

# Lyman-Werner line absorption from the code's own pumping
lw = np.loadtxt(os.path.join(run, 'output', 'Lyman_Werner.txt'))[2:-2]
vol = 4 * np.pi / 3 * (f[1:]**3 - f[:-1]**3)
diss = np.sum(lw[:, 6] * lw[:, 3] * vol) / (4 * np.pi * R0**2)   # cm^-2 s^-1
P_LW_abs = diss / 0.135 * 11.7354 * EV
F_LW = np.trapz(F[lw_lines], wl[lw_lines])

bands = [('X-ray 0.9-10 nm', 0.9, 10), ('EUV 10-50.4 nm (He I edge)', 10, 50.4),
         ('EUV 50.4-91.2 nm', 50.4, 91.2), ('LW 91.2-100 nm', 91.2, 100),
         ('LW 100-120.1 nm', 100, 120.1), ('Ly alpha 121.0-122.2 nm', 121.0, 122.2),
         ('FUV 120.1-200 nm (excl. Ly a)', 120.1, 200), ('NUV 200-400 nm', 200, 400),
         ('400-3000 nm', 400, 3000)]
print('\n%-32s %12s %12s %10s' % ('band', 'F_top', 'abs. wind', 'share'))
rows = []
for name, a, b in bands:
    m = (wl >= a) & (wl < b)
    if 'FUV' in name:
        m &= ~((wl >= 121.0) & (wl < 122.2))
    Ft = np.trapz(F[m], wl[m]) if m.sum() > 1 else 0.0
    Fa = np.trapz(absorbed[m], wl[m]) if m.sum() > 1 else 0.0
    if name.startswith('LW'):
        Fa = np.nan
    rows.append((name, Ft, Fa))
    print('%-32s %12.4e %12.4e %10.4f' % (name, Ft, Fa, Fa / Ft if Ft > 0 else 0))
print('total SED 0.9-3000 nm: %.4e erg cm^-2 s^-1' % np.trapz(F, wl))
print('LW band 91.2-120.1 nm: F_top %.4e (code F_LW_star 1.0897), absorbed by Solomon pumping '
      '(k_LW n_H2 / 0.135 x 11.7354 eV) %.4e erg cm^-2 s^-1 (share %.3f); dissociations %.3e cm^-2 s^-1'
      % (F_LW, P_LW_abs, P_LW_abs / F_LW, diss))
m = wl < 91.2
print('below 91.2 nm: F_top %.4e, absorbed in the wind %.4e (share %.4f), '
      'transmitted %.4e' % (np.trapz(F[m], wl[m]), np.trapz(absorbed[m], wl[m]),
                            np.trapz(absorbed[m], wl[m]) / np.trapz(F[m], wl[m]),
                            np.trapz((F - absorbed)[m], wl[m])))
m = wl < 100
print('below 100 nm: F_top %.4e, absorbed (continua; LW lines 91.2-100 excluded) %.4e'
      % (np.trapz(F[m], wl[m]), np.trapz(absorbed[m], wl[m])))
for lam in (1, 2, 5, 10, 20, 30, 50, 60, 80, 90):
    i = np.argmin(abs(wl - lam))
    print('  tau_wind(%5.1f nm) = %.3e  (He I %.2e, H2 %.2e, H I %.2e)'
          % (wl[i], tau[i], sig['HeI'][i] * N['HeI'], sig['H2'][i] * N['H2'], sig['HI'][i] * N['HI']))
# Lyman alpha core: line-centre tau of the H I column (Doppler at 1000 K)
s0 = 0.4162 * 0.02654 / (np.sqrt(np.pi) * (np.sqrt(2 * 1.380649e-16 * 1000 / 1.6735e-24) / 1215.67e-8))
print('Ly alpha line-centre tau of N(H I) at 1000 K: %.3e' % (s0 * N['HI']))
# heating the wind deposits, for comparison
hb = np.loadtxt(os.path.join(run, 'output', 'Heating_breakdown.txt'))[2:-2]
# Ly alpha: Voigt line of the wind's H I column (T = 1000 K), the share of the
# stellar Ly alpha band (121.0-122.2 nm, SED as given) with tau > 1, which the
# wind SCATTERS (no absorption; the 1-D wind does not transport it)
from scipy.special import voigt_profile
lam0 = 121.567
m = (wl >= 121.0) & (wl < 122.2)
nu0 = 2.99792458e10 / (lam0 * 1e-7)
dnu = 2.99792458e10 / (wl[m] * 1e-7) - nu0
bD = np.sqrt(2 * 1.380649e-16 * 1000 / 1.6735e-24)
sigG = nu0 * bD / 2.99792458e10 / np.sqrt(2)
gamL = 6.265e8 / (4 * np.pi)
sig_lya = 0.02654 * 0.4162 * voigt_profile(dnu, sigG, gamL)
tl = sig_lya * N['HI']
print('Ly alpha band %.4e erg cm^-2 s^-1; flux at tau_HI > 1: %.4e (share %.3f), '
      'F(1-exp(-tau)): %.4e' % (np.trapz(F[m], wl[m]), np.trapz(F[m] * (tl > 1), wl[m]),
                                np.trapz(F[m] * (tl > 1), wl[m]) / np.trapz(F[m], wl[m]),
                                np.trapz(F[m] * (1 - np.exp(-tl)), wl[m])))

# the column under the wind: where the transmitted X-rays (and, for Photochem,
# the unattenuated XUV) are absorbed, with the photoabsorption tables of He,
# H2, CH4, NH3, H2O on the k05 profile; columns N_X(p) = x_X p/(m g)
cols = ('p r T n_tot rho Kzz q_H2 q_H X_He X_C X_N X_O F_H F_He q_H2O q_CO '
        'q_CO2 q_CH4 q_NH3 q_HCN q_N2 q_OH q_C2H2').split()
pr = np.loadtxt(os.path.join(run, 'lower_atmosphere_profile.dat'))
pr = pr[np.argsort(pr[:, 0])]
pb = pr[:, 0]
mb = pr[:, cols.index('rho')] / pr[:, cols.index('n_tot')]
g = 6.674e-8 * float(res['planet_mass_MJ']) * 1.89813e30 / (pr[:, cols.index('r')] * R_J)**2
qq = {k: pr[:, cols.index('q_' + k)] for k in ('H2', 'CH4', 'NH3', 'H2O')}
qq['He'] = 1 - sum(pr[:, cols.index(c)] for c in cols if c.startswith('q_'))
Ncol = {k: np.concatenate([[0.0], np.cumsum(0.5 * (v[1:] / (mb[1:] * g[1:]) + v[:-1] / (mb[:-1] * g[:-1]))
                                              * np.diff(pb) * 1e6)]) for k, v in qq.items()}
sigc = {'He': tab_xs('He'), 'H2': tab_xs('H2'), 'CH4': tab_xs('CH4'),
        'NH3': tab_xs('NH3'), 'H2O': tab_xs('H2O')}
for lab, FF, lo, hi in (('transmitted by the wind', F - absorbed, 0.9, 91.2),
                        ('unattenuated (what Photochem takes)', F, 0.9, 91.2),
                        ('unattenuated (what Photochem takes)', F, 91.2, 120.1)):
    m = (wl >= lo) & (wl < hi)
    print('column absorption of the %s %.1f-%.1f nm (%.4e erg cm^-2 s^-1 at normal incidence):'
          % (lab, lo, hi, np.trapz(FF[m], wl[m])))
    prev = np.trapz(FF[m], wl[m])
    for pcut in (1.0e-6 * np.exp(0.5), 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0):
        i = np.searchsorted(pb, pcut)
        t = sum(sigc[k][m] * Ncol[k][i] for k in sigc)
        rem = np.trapz(FF[m] * np.exp(-np.minimum(t, 700)), wl[m])
        print('   p <= %.3e bar: absorbed %.4e, remaining %.4e' % (pcut, prev - rem if pcut == 1.0e-6 * np.exp(0.5) else np.trapz(FF[m], wl[m]) - rem, rem))
print('wind heating total (Heating_breakdown col4) per face area: %.4e erg cm^-2 s^-1'
      % (np.sum(hb[:, 3] * vol) / (4 * np.pi * R0**2)))
