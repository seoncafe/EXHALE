#!/usr/bin/env python3
"""He 10830 excess absorption against the H:He ratio (LHS 1140b).

This reproduces the central argument of Cherubim et al. (2026): at the low
XUV flux of LHS 1140b, hydrogen at more than about 1% absorbs the photons
below 911 A that ionize ground-state helium, so metastable helium is not
produced and the observed line cannot be matched.  Their conclusion,
H:He <~ 1e-3, follows from that dependence.

Everything except H:He is held at the oracle's reference point
(run_pwinds.py, case matched_gj1132): the v23 GJ 1132 spectrum,
Mdot = 2.03e8 g/s, T = 6100 K -- the temperature inside the paper's
explored range at which our SED reproduces the observed depth.

Outputs
  tspec_hhe_<tag>.txt   vacuum wavelength [A], transmission, excess [%],
                        instrument-convolved excess [%]
  scan_hhe.txt          summary table: H:He, n_He3 peak, metrics
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '..', '..', '..', 'p-winds')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '..', '..')))
import astropy.units as u
import astropy.constants as c
from p_winds import parker, hydrogen, helium, transit, lines, tools
from he_line_metrics import fit_metrics

# reference point (../system_parameters.md, pwinds_oracle/README.md)
R_pl_earth, M_pl_earth = 1.730, 5.60
M_star, R_star_sun, a_au, inc_deg = 0.1844, 0.2159, 0.0946, 89.96
T_wind = 5160.0
SED = '../sed/lhs1140_sed_gj1132_at_b.txt'

GRID = [(m, h) for m in (1.0e8, 1.45e8, 2.03e8)
              for h in (1.0e-3, 10**-3.5, 1.0e-4)]

R_pl = (R_pl_earth*c.R_earth/c.R_jup).decompose().value
M_pl = (M_pl_earth*c.M_earth/c.M_jup).decompose().value
R_pl_m = (R_pl_earth*c.R_earth).to(u.m).value
R_star_m = (R_star_sun*c.R_sun).to(u.m).value
b_impact = (a_au*u.au*np.cos(np.radians(inc_deg))
            / (R_star_sun*c.R_sun)).decompose().value

RESOLVING_POWER = 68000.0
AIR_TO_VAC = 10832.057/10829.09114


def convolve_lsf(wl_m, spec):
    wl0 = np.mean(wl_m)
    sig = wl0/RESOLVING_POWER/(2.0*np.sqrt(2.0*np.log(2.0)))
    dl = np.median(np.diff(wl_m))
    x = np.arange(-5*sig, 5*sig + dl, dl)
    k = np.exp(-0.5*(x/sig)**2); k /= k.sum()
    return np.convolve(spec, k, mode='same')


spectrum = tools.make_spectrum_from_file(
    SED, units={'wavelength': u.angstrom,
                'flux': u.erg/u.s/u.cm**2/u.angstrom}, skiprows=0)

r = np.logspace(0, np.log10(20.0), 200)
l0, l1, l2, f0, f1, f2, a_ij = lines.he_3_properties()
wl_grid = np.linspace(1.0826e-6, 1.0838e-6, 400)
m_He = 4.002602*1.66053907e-27
m_h = c.m_p.to(u.g).value

grid, t_depth, r_from_planet = transit.draw_transit(
    planet_to_star_ratio=R_pl_m/R_star_m, planet_physical_radius=R_pl_m,
    impact_parameter=b_impact, grid_size=121, supersampling=10)

rows = []
print("# grid at T = 5160 K, paper-matched SED")
print(f"{'H:He':>10} {'mu':>6} {'f_ion(top)':>11} {'n_He3max':>11} "
      f"{'red[%]':>8} {'blue[%]':>8} {'red/blue':>9} {'FWHM[A]':>8}")

for mdot, hhe in GRID:
    h_frac = hhe/(1.0 + hhe)
    f_ion, mu_bar = hydrogen.ion_fraction(
        r, R_pl, T_wind, h_frac, mdot, M_pl,
        mean_molecular_weight_0=4.0*(1 - h_frac) + h_frac,
        star_mass=M_star, semimajor_axis=a_au,
        spectrum_at_planet=spectrum, relax_solution=True, return_mu=True)
    vs = parker.sound_speed(T_wind, mu_bar)
    rs = parker.radius_sonic_point(M_pl, vs)
    rhos = parker.density_sonic_point(mdot, rs, vs)
    v_arr, rho_arr = parker.structure(r*R_pl/rs)
    _, f_he_3 = helium.population_fraction(
        r, v_arr, rho_arr, f_ion, R_pl, T_wind, h_frac, vs, rs, rhos,
        spectrum_at_planet=spectrum,
        initial_state=np.array([1.0, 0.0]), relax_solution=True)

    he_frac = 1.0 - h_frac
    n_he3 = (rho_arr*rhos/(m_h*(h_frac + 4.0*he_frac)))*he_frac*f_he_3

    t_spec = transit.radiative_transfer_2d(
        grid, r_from_planet, r*R_pl_m, n_he3*1e6, v_arr*vs*1e3,
        np.array([l0, l1, l2]), np.array([f0, f1, f2]),
        np.array([a_ij]*3), wl_grid, T_wind, m_He,
        wind_broadening_method='average')

    excess = (1.0 - t_spec/np.max(t_spec))*100.0
    exc_conv = convolve_lsf(wl_grid, excess)
    wl_vac = wl_grid*1e10*AIR_TO_VAC

    tag = f"m{mdot:.2e}_h{hhe:g}".replace('.', 'p').replace('-', 'm')
    np.savetxt(f'tspec_grid5160_{tag}.txt',
               np.column_stack([wl_vac, t_spec, excess, exc_conv]),
               header=(f'LHS 1140b p-winds, H:He = {hhe:g}, '
                       f'Mdot = {mdot:.3e} g/s, T = {T_wind:.0f} K\n'
                       'vacuum wavelength [A]  transmission  excess[%]  '
                       'excess convolved R=68000 [%]'))

    m = fit_metrics(wl_vac, exc_conv)
    print(f"{hhe:10.4g} {mu_bar:6.3f} {f_ion[-1]:11.3f} {n_he3.max():11.3e} "
          f"{m['red_depth']:8.3f} {m['blue_depth']:8.3f} "
          f"{m['red_blue']:9.2f} {m['fwhm_A']:8.3f}")
    rows.append([hhe, mu_bar, f_ion[-1], n_he3.max(), m['red_depth'],
                 m['blue_depth'], m['red_blue'], m['fwhm_A']])

print('# grid done')
