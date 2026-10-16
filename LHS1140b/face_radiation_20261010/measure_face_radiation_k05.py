#!/usr/bin/env python3
"""Measured effect of the wind's face optical depth on the LHS 1140 b k05
column (md/matching_face_radiation_python_20261010.md section 4).

Reads (never writes) the stored k05 run: the SED at the planet
(photochem_work/flux_photochem.txt), the column profile
(lower_atmosphere_profile.dat) and the face file built from the stored wind
outputs by build_wind_face_optical_depth.py (the stand-in).  Uses the
adapter's own functions (src/utils/photochem_to_lower_profile.py), so what
is measured is what the adapter applies.

  1. band powers before and after the wind, for Photochem and the climate
     model, under the two geometries (per_model, wind);
  2. the photoabsorptions the attenuation removes from the column, by
     species, with Photochem's own cross sections on the k05 profile;
  3. the X-ray heating of the column: deposit profile, heat share, total;
  4. the climate model's absorbed stellar flux under the two geometries
     (clima's own radiative transfer, on copies of the k05 inputs).
Figures to docs/figures/matching_face_*.png.

Usage: measure_face_radiation_k05.py <k05 dir> <face file> <scratch dir>
"""
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(TREE, 'src', 'utils'))
import photochem_to_lower_profile as ad          # noqa: E402

import matplotlib                                 # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt                   # noqa: E402

run, face_path, scratch = sys.argv[1:4]
FIG = os.path.join(TREE, 'docs', 'figures')
R_J = 7.1492e9
G = 6.674e-8
M_J = 1.89813e30
EV = 1.602176634e-12
HC_EV_NM = 1239.841984

wl, F0 = np.loadtxt(os.path.join(run, 'photochem_work', 'flux_photochem.txt'),
                    skiprows=1, unpack=True)
face = ad.read_wind_face_optical_depth(face_path)
lines = ad.h2_lyman_werner_line_depth(face['layers'])
geos = {'per_model': ad.face_geometry('per_model', None, 48.0),
        'wind': ad.face_geometry('wind', face['dilution'])}
out = open(os.path.join(HERE, 'k05_face_radiation_measured.txt'), 'w')


def say(s=''):
    print(s)
    out.write(s + '\n')


say('# MEASURED on k05 with the stand-in face file %s (sha256 %s)'
    % (os.path.basename(face_path), face['sha'][:16]))
say('# H2 layers: ' + ', '.join('%.3e cm^-2 at %.0f K' % l
                                 for l in face['layers']))

# ---- 1. band powers --------------------------------------------------- #
bands = (('X-ray 0.9-10 nm', 0.0, 10.0), ('EUV 10-80.4 nm', 10.0, 80.38),
         ('80.4-91.2 nm', 80.38, 91.2), ('LW 91.2-100 nm', 91.2, 100.0),
         ('LW 100-120.1 nm', 100.0, 120.1))
trans = {}
for gname, g in geos.items():
    for model in ('photochem', 'climate'):
        trans[(gname, model)] = ad.wind_face_transmission(
            wl, face, g['%s_mu' % model], g['%s_beam' % model], lines)
t_line_only = ad.wind_face_transmission(
    wl, dict(face, tau=np.zeros_like(face['tau'])), [1.0], [1.0], lines)
t_cont_only = ad.wind_face_transmission(wl, face, [1.0], [1.0], None)
say('\n## 1. Band powers at normal incidence [erg cm^-2 s^-1]: the star, and'
    ' what crosses the wind for each model (energy-weighted over its slants)')
say('%-18s %10s | %-28s | %-28s' % ('band', 'star', 'per_model: Photochem  climate',
                                     'wind: Photochem  climate'))
for name, lo, hi in bands:
    e0 = ad.band_power(wl, F0, lo, hi)
    v = [ad.band_power(wl, F0*trans[k], lo, hi) for k in
         (('per_model', 'photochem'), ('per_model', 'climate'),
          ('wind', 'photochem'), ('wind', 'climate'))]
    say('%-18s %10.4e | %13.4e %13.4e | %13.4e %13.4e' % ((name, e0) + tuple(v)))
for name, lo, hi in (('LW 91.2-120.1 nm', 91.2, 120.1),):
    say('%s at mu = 1: lines only %.4f, continuum only %.4f, both %.4f of %.4e'
        % (name, ad.band_power(wl, F0*t_line_only, lo, hi)/ad.band_power(wl, F0, lo, hi),
           ad.band_power(wl, F0*t_cont_only, lo, hi)/ad.band_power(wl, F0, lo, hi),
           ad.band_power(wl, F0*trans[('wind', 'photochem')], lo, hi)
           / ad.band_power(wl, F0, lo, hi), ad.band_power(wl, F0, lo, hi)))
for gname, g in geos.items():
    say('geometry %-9s: horizontal energy factor Photochem %.4f (zenith %.2f deg,'
        ' diurnal %s), climate %.4f'
        % (gname, g['photochem_energy_factor'], g['photochem_zenith_deg'],
           g['photochem_diurnal_fac'], g['climate_energy_factor']))

# ---- the k05 column ---------------------------------------------------- #
res = {}
for line in open(os.path.join(run, 'EXHALE_resolved.out')):
    t = line.split()
    if len(t) >= 2 and not line.startswith('#'):
        res[t[0]] = t[1]
cols = None
for line in open(os.path.join(run, 'lower_atmosphere_profile.dat')):
    if line.startswith('# columns:'):
        cols = line.split()[2:]
pr = np.loadtxt(os.path.join(run, 'lower_atmosphere_profile.dat'))
pr = pr[np.argsort(-pr[:, 0])]                 # deep -> shallow (top last)
p_bar = pr[:, cols.index('p')]
r_cm = pr[:, cols.index('r')]*R_J
n_tot = pr[:, cols.index('n_tot')]
rho = pr[:, cols.index('rho')]
T = pr[:, cols.index('T')]
q = {c[2:]: pr[:, cols.index(c)] for c in cols if c.startswith('q_')}
q['He'] = 1.0 - sum(q.values())
z = r_cm - r_cm[-1]                            # ascending, top = face
mix = {sp: q[sp] for sp in ad.XRAY_COLUMN_ABSORBERS if sp in q}
xs = {sp: ad.photoabsorption_table(sp) for sp in
      ('He', 'H2', 'H', 'CH4', 'NH3', 'H2O', 'CO', 'CO2', 'N2', 'HCN', 'C2H2')}
xs = {k: v for k, v in xs.items() if v is not None}
say('\n# column: %d levels %.3e-%.3e bar; at the match q_H2 %.4f q_He %.4f'
    ' q_CH4 %.3e q_NH3 %.3e T %.1f K n %.3e'
    % (p_bar.size, p_bar[-1], p_bar[0], q['H2'][-1], q['He'][-1],
       q['CH4'][-1], q['NH3'][-1], T[-1], n_tot[-1]))

# ---- 2. photoabsorptions removed --------------------------------------- #
say('\n## 2. Photons absorbed in the column below 120.1 nm [cm^-2 s^-1],'
    ' Photochem geometry, Photochem cross sections (quantum yield 1)')
lam_m = wl < ad.WIND_FACE_TAU_LAMBDA_MAX_NM
lam = wl[lam_m]
dlam = np.diff(ad.sample_cells(wl))[lam_m]
phot = F0[lam_m]*dlam/(HC_EV_NM/lam*EV)          # photons cm^-2 s^-1 per cell
dz = np.diff(z)
Ncol = {}
for sp in xs:
    if sp in q:
        n = q[sp]*n_tot
        N = np.empty_like(n)
        N[-1] = 0.0
        for j in range(n.size - 2, -1, -1):
            N[j] = N[j + 1] + 0.5*(n[j] + n[j + 1])*dz[j]
        Ncol[sp] = N
sig = {sp: np.interp(lam, xs[sp][0], xs[sp][1], left=0.0, right=0.0)
       for sp in Ncol}
tau_col = sum(sig[sp][None, :]*Ncol[sp][:, None] for sp in Ncol)  # (nz, nlam)
i_half = int(np.argmin(np.abs(p_bar - 1e-6*np.exp(0.5))))


def absorbed_by_species(tphot, mu, beam, k_layer=0):
    """photons absorbed above level k_layer by species (column total for
    k_layer = 0), summed over the beams."""
    outd = {}
    for m, b in zip(mu, beam):
        trans_k = np.exp(-np.minimum(tau_col[k_layer]/m, 700.0))
        absorbed = b*tphot*(1.0 - trans_k)          # photons per cell
        tot = np.where(tau_col[k_layer] > 0, tau_col[k_layer], 1.0)
        for sp in Ncol:
            outd[sp] = outd.get(sp, 0.0) + float(np.sum(
                absorbed*sig[sp]*Ncol[sp][k_layer]/tot))
    return outd


for gname in ('per_model', 'wind'):
    g = geos[gname]
    t = trans[(gname, 'photochem')][lam_m]
    for lab, k in (('whole column', 0), ('above p_match e^0.5', i_half)):
        a0 = absorbed_by_species(phot, g['photochem_mu'], g['photochem_beam'], k)
        a1 = absorbed_by_species(phot*t, g['photochem_mu'], g['photochem_beam'], k)
        say('%-9s %-20s ' % (gname, lab) + '; '.join(
            '%s %.3e -> %.3e (%.1f%% removed)' % (sp, a0[sp], a1[sp],
                                                100*(1 - a1[sp]/a0[sp]) if a0[sp] > 0 else 0)
            for sp in ('He', 'H2', 'CH4', 'NH3', 'H2O', 'C2H2', 'HCN') if sp in a0))
    for lo, hi in ((0, 80.38), (80.38, 91.2), (91.2, 120.1)):
        mm = (lam >= lo) & (lam < hi)
        e_star = float(np.sum((F0[lam_m]*dlam)[mm]))*g['photochem_energy_factor']
        e_after = float(np.sum((F0[lam_m]*dlam*t)[mm]))*g['photochem_energy_factor']
        say('%-9s Photochem %.2f-%.1f nm: %.4e -> %.4e erg cm^-2 s^-1 on the'
            ' horizontal (%.4e photons -> %.4e)'
            % (gname, lo, hi, e_star, e_after,
               float(np.sum(phot[mm]))*g['photochem_energy_factor'],
               float(np.sum((phot*t)[mm]))*g['photochem_energy_factor']))

# ---- 3. X-ray heating -------------------------------------------------- #
say('\n## 3. Photoionizing radiation the wind transmits (< 80.38 nm),'
    ' deposited in the column')
x_top = {'He': q['He'][-1], 'H2': q['H2'][-1], 'H': q['H'][-1]}
E_max = HC_EV_NM/wl.min() + 1.0
eta = ad.photoelectron_heat_share(x_top['He'], x_top['H2'], x_top['H'],
                                  T[-1], E_max)
sol = eta.record
import electron_degradation_spencer_fano as sf   # noqa: E402
for E0 in (30, 100, 300, 1000, 1300):
    r = sol.record(E0)
    say('  Spencer-Fano at E0 = %4d eV (x_He %.3f, x_H2 %.3f, x_e %.0e, T_e'
        ' %.0f K): heat share %.4f = f_heat %.4f + rot %.4f + vib %.4f;'
        ' ionization %.4f, He excitation %.4f, H2 electronic %.4f,'
        ' H2 dissociation (chemical) %.4f'
        % (E0, x_top['He'], x_top['H2'], ad.XRAY_COLUMN_ELECTRON_FRACTION,
           T[-1], eta(E0),
           (r['heat_coulomb'] + r['heat_recoil'] + r['heat_thermalized']
            + r['heat_fragments'])/E0, r['E_rot_H2']/E0, r['E_vib_H2']/E0,
           (r['E_ion_HI'] + r['E_ion_HeI'] + r['E_ion_He23S'] + r['E_ion_H2']
            + r['E_ion_H2_diss'])/E0,
           (r['E_exc_He_singlet'] + r['E_exc_He_triplet'])/E0,
           (r['E_exc_H2_singlet'] + r['E_exc_H2_triplet_bound'])/E0,
           r['E_diss_H2_chem']/E0))
xray = {}
for gname in ('per_model', 'wind'):
    g = geos[gname]
    qa_t = np.zeros_like(z)
    qh_t = np.zeros_like(z)
    rec_t = dict(incident=0.0, absorbed=0.0, heat=0.0, by_absorber={})
    for m, b in zip(g['climate_mu'], g['climate_beam']):
        F_m = F0*ad.wind_face_transmission(wl, face, [m], [1.0])
        qa, qh, rm = ad.xray_deposit_in_column(z, n_tot, mix, wl, F_m, [m], [b],
                                               eta, xs=xs)
        qa_t += qa
        qh_t += qh
        for k in ('incident', 'absorbed', 'heat'):
            rec_t[k] += rm[k]
        for k, v in rm['by_absorber'].items():
            rec_t['by_absorber'][k] = rec_t['by_absorber'].get(k, 0.0) + v
    w = np.gradient(z)
    cum_from_top = np.cumsum((qh_t*w)[::-1])        # top -> deep
    p_desc = p_bar[::-1]
    def above(pc):
        return float(np.sum((qh_t*w)[p_bar <= pc]))
    xray[gname] = (qa_t, qh_t, rec_t)
    say('%-9s (climate beams, energy factor %.3f): incident %.4e, absorbed %.4e'
        ' (%s), heat %.4e erg cm^-2 s^-1 (share %.4f)'
        % (gname, g['climate_energy_factor'], rec_t['incident'], rec_t['absorbed'],
           ', '.join('%s %.2f' % (k, v/rec_t['absorbed']) for k, v in
                     sorted(rec_t['by_absorber'].items(), key=lambda kv: -kv[1])[:4]),
           rec_t['heat'], rec_t['heat']/rec_t['absorbed']))
    say('          heat above p_match e^0.5 %.4e, above 1e-5 bar %.4e, above'
        ' 1e-4 bar %.4e, above 1e-3 bar %.4e; half of it above %.3e bar'
        % (above(1e-6*np.exp(0.5)), above(1e-5), above(1e-4), above(1e-3),
           float(np.interp(0.5*cum_from_top[-1], cum_from_top, p_desc))))

# ---- 4. the climate model's absorbed stellar flux ---------------------- #
say('\n## 4. clima, absorbed stellar flux at the top for the two geometries'
    ' (k05 climate inputs copied to scratch, unattenuated spectrum)')
cd = os.path.join(scratch, 'clima_geometry')
os.makedirs(cd, exist_ok=True)
for f in ('climate_settings.yaml', 'climate_species.yaml', 'flux_photochem.txt'):
    shutil.copy(os.path.join(run, 'photochem_work', f), os.path.join(cd, f))
from photochem.clima import AdiabatClimate      # noqa: E402
import radiative_convective_column as rcc       # noqa: E402
mixing = rcc.deep_molecular_partition(dict(H=1.0, He=2.09, C=2.77588e-4,
                                           N=8.18465e-5, O=6.06178e-4))
isr = {}
for gname, g in geos.items():
    c = AdiabatClimate(os.path.join(cd, 'climate_species.yaml'),
                       os.path.join(cd, 'climate_settings.yaml'),
                       os.path.join(cd, 'flux_photochem.txt'))
    if g['climate_zenith_u'] is not None:
        c.rad.zenith_u = np.asarray(g['climate_zenith_u'], dtype=float)
    if g['climate_photon_scale'] is not None:
        c.rad.photon_scale_factor = g['climate_photon_scale']
    P_i = np.array([mixing.get(sp, 0.0)*20.0e6 for sp in c.species_names])
    isr[gname] = c.TOA_fluxes(400.0, P_i)
    say('%-9s TOA (absorbed stellar, OLR) at T_deep 400 K: %.4e, %.4e erg cm^-2'
        ' s^-1; bolometric %.4e' % (gname, isr[gname][0], isr[gname][1],
                                    c.rad.bolometric_flux()*1e3))
say('ratio of absorbed stellar flux wind/per_model: %.4f'
    % (isr['wind'][0]/isr['per_model'][0]))
out.close()

# ---- figures ----------------------------------------------------------- #
plt.rcParams.update({'font.size': 9})
fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.6))
m = wl < 121.0
ax[0].semilogy(wl[m], np.maximum(face['tau'][np.clip(np.searchsorted(
    face['lo'], wl[m], side='right') - 1, 0, None)], 1e-6), 'k-', lw=1,
    label=r'$\tau_{\rm cont}$ (radial)')
nu, tl = lines
lam_l = 2.99792458e17/nu
ax[0].semilogy(lam_l[::20], np.maximum(tl[::20], 1e-6), '-', color='C0',
               lw=0.3, alpha=0.6, label=r'H$_2$ Lyman-Werner lines')
ax[0].set_xlabel(r'$\lambda$ [nm]')
ax[0].set_ylabel(r'$\tau_{\rm wind}$ at the face')
ax[0].set_ylim(1e-4, 1e4)
ax[0].legend(loc='upper right', fontsize=7)
ax[0].set_title(r'k05 stand-in: radial optical depth of the wind')
mm = wl < 125.0
ax[1].semilogy(wl[mm], F0[mm], 'k-', lw=0.8, label='star (at the planet)')
ax[1].semilogy(wl[mm], np.maximum(F0[mm]*trans[('wind', 'photochem')][mm], 1e-8),
               '-', color='C3', lw=0.8, label=r'below the wind, $\mu=1$')
ax[1].set_xlabel(r'$\lambda$ [nm]')
ax[1].set_ylabel(r'$F_\lambda$ [erg cm$^{-2}$ s$^{-1}$ nm$^{-1}$]')
ax[1].set_ylim(1e-6, 1e2)
ax[1].legend(loc='lower right', fontsize=7)
ax[1].set_title('spectrum the column takes')
fig.tight_layout()
fig.savefig(os.path.join(FIG, 'matching_face_tau_and_spectrum_k05.png'), dpi=150)

fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.6))
for gname, ls in (('wind', '-'), ('per_model', '--')):
    qa, qh, rec = xray[gname]
    w = np.gradient(z)
    lnp = np.log(p_bar)
    dlnp = np.abs(np.gradient(lnp))
    ax[0].loglog(qh*w/dlnp, p_bar, ls, color='C3',
                 label=r'heat, ' + gname.replace('_', ' '))
    ax[0].loglog(qa*w/dlnp, p_bar, ls, color='C0',
                 label=r'absorbed, ' + gname.replace('_', ' '))
ax[0].axhline(1e-6*np.exp(0.5), color='0.5', lw=0.6)
ax[0].set_ylim(1e-1, 1e-6)
ax[0].set_xlim(1e-4, 10)
ax[0].set_xlabel(r'$dF/d\ln p$ [erg cm$^{-2}$ s$^{-1}$]')
ax[0].set_ylabel(r'$p$ [bar]')
ax[0].legend(fontsize=7)
ax[0].set_title('transmitted X-rays in the k05 column')
Eg = np.linspace(20, 1300, 300)
ax[1].plot(Eg, eta(Eg), 'k-', label='heat share (used)')
rr = [sol.record(E) for E in Eg]
ax[1].plot(Eg, [(r['heat_coulomb'] + r['heat_recoil'] + r['heat_thermalized']
                 + r['heat_fragments'])/E for r, E in zip(rr, Eg)], '--',
           label=r'$f_{\rm heat}$ (elastic, Coulomb, fragments)')
ax[1].plot(Eg, [(r['E_ion_HI'] + r['E_ion_HeI'] + r['E_ion_He23S'] + r['E_ion_H2']
                 + r['E_ion_H2_diss'])/E for r, E in zip(rr, Eg)], ':',
           label='ionization')
ax[1].plot(Eg, [(r['E_exc_He_singlet'] + r['E_exc_He_triplet'])/E
                for r, E in zip(rr, Eg)], '-.', label='He excitation')
ax[1].set_xlabel(r'photoelectron energy $E_0$ [eV]')
ax[1].set_ylabel('share of $E_0$')
ax[1].set_ylim(0, 0.6)
ax[1].legend(fontsize=7)
ax[1].set_title(r'Spencer-Fano, He/H$_2$ column at the match')
fig.tight_layout()
fig.savefig(os.path.join(FIG, 'matching_face_xray_heating_k05.png'), dpi=150)
print('figures written')
