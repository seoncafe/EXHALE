#!/usr/bin/env python3
"""Mdot scan around the published retrieval point (Phase A3 supplement).

With the Phase-A2 spectrum (fxuv33; the paper's exact v23 SED is
unrecoverable), the forward model at the published point overshoots the
observed depths.  This scan maps how the blended-red depth, blue depth and
their ratio move with Mdot at T = 5160 K, H:He = 1e-3, to locate where OUR
spectrum reproduces the 2024 observation (red 1.24%, blue 0.25%,
red/blue amplitude ratio 6.7).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '..', '..', '..', 'p-winds')))
import astropy.units as u
import astropy.constants as c
from p_winds import parker, hydrogen, helium, transit, lines, tools

R_pl_earth, M_pl_earth = 1.730, 5.60
M_star, R_star_sun, a_au, inc_deg = 0.1844, 0.2159, 0.0946, 89.96
T_wind, h_he = 5160.0, 1.0e-3
h_frac = h_he/(1.0 + h_he)

R_pl = (R_pl_earth*c.R_earth/c.R_jup).decompose().value
M_pl = (M_pl_earth*c.M_earth/c.M_jup).decompose().value
R_pl_m = (R_pl_earth*c.R_earth).to(u.m).value
R_star_m = (R_star_sun*c.R_sun).to(u.m).value
b_impact = (a_au*u.au*np.cos(np.radians(inc_deg))/(R_star_sun*c.R_sun)).decompose().value

RESOLVING_POWER = 68000.0
AIR_TO_VAC = 10832.057/10829.09114

def convolve_lsf(wl_m, spec):
    wl0 = np.mean(wl_m); fwhm = wl0/RESOLVING_POWER
    sig = fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
    dl = np.median(np.diff(wl_m))
    x = np.arange(-5*sig, 5*sig + dl, dl)
    k = np.exp(-0.5*(x/sig)**2); k /= k.sum()
    return np.convolve(spec, k, mode='same')

spectrum = tools.make_spectrum_from_file(
    '../sed/lhs1140_sed_gj1132_fxuv33_at_b.txt',
    units={'wavelength': u.angstrom, 'flux': u.erg/u.s/u.cm**2/u.angstrom},
    skiprows=5)

r = np.logspace(0, np.log10(20.0), 200)
l0, l1, l2, f0, f1, f2, a_ij = lines.he_3_properties()
wl_grid = np.linspace(1.0826e-6, 1.0838e-6, 400)
m_He = 4.002602*1.66053907e-27
m_h = c.m_p.to(u.g).value

grid, t_depth, r_from_planet = transit.draw_transit(
    planet_to_star_ratio=R_pl_m/R_star_m, planet_physical_radius=R_pl_m,
    impact_parameter=b_impact, grid_size=121, supersampling=10)

print("# Mdot[g/s]  red_conv[%]  blue_conv[%]  red/blue  n_He3max[cm-3]")
rows = []
for mdot in [5.0e7, 7.0e7, 1.0e8, 1.4e8, 2.03e8]:
    f_ion, mu_bar = hydrogen.ion_fraction(
        r, R_pl, T_wind, h_frac, mdot, M_pl,
        mean_molecular_weight_0=4.0*(1-h_frac)+h_frac,
        star_mass=M_star, semimajor_axis=a_au,
        spectrum_at_planet=spectrum, relax_solution=True, return_mu=True)
    vs = parker.sound_speed(T_wind, mu_bar)
    rs = parker.radius_sonic_point(M_pl, vs)
    rhos = parker.density_sonic_point(mdot, rs, vs)
    v_array, rho_array = parker.structure(r*R_pl/rs)
    f_he_1, f_he_3 = helium.population_fraction(
        r, v_array, rho_array, f_ion, R_pl, T_wind, h_frac, vs, rs, rhos,
        spectrum_at_planet=spectrum,
        initial_state=np.array([1.0, 0.0]), relax_solution=True)
    he_frac = 1.0 - h_frac
    n_he_3 = (rho_array*rhos/(m_h*(h_frac+4.0*he_frac)))*he_frac*f_he_3
    t_spec = transit.radiative_transfer_2d(
        grid, r_from_planet, r*R_pl_m, n_he_3*1e6, v_array*vs*1e3,
        np.array([l0, l1, l2]), np.array([f0, f1, f2]),
        np.array([a_ij]*3), wl_grid, T_wind, m_He,
        wind_broadening_method='average')
    excess = (1.0 - t_spec/np.max(t_spec))*100.0
    ec = convolve_lsf(wl_grid, excess)
    wlA = wl_grid*1e10*AIR_TO_VAC
    red = ec[(wlA > 10832.8) & (wlA < 10833.7)].max()
    blu = ec[(wlA > 10831.6) & (wlA < 10832.5)].max()
    print(f"{mdot:10.3e}  {red:8.3f}  {blu:8.3f}  {red/blu:6.2f}  {n_he_3.max():.3e}")
    rows.append([mdot, red, blu, red/blu, n_he_3.max()])
np.savetxt('scan_mdot_fxuv33.txt', np.array(rows),
           header='Mdot[g/s] red_conv[%] blue_conv[%] red/blue n_He3max[cm-3]'
                  '\nT=5160K H:He=1e-3 SED=fxuv33; targets red 1.24 blue 0.25 ratio 6.7')
