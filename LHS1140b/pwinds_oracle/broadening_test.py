#!/usr/bin/env python3
"""p-winds oracle for LHS 1140b (Phase A3).

Forward-models the He 10830 transmission spectrum at the published retrieval
point of Cherubim et al. (2026): Mdot = 2.03e8 g/s, T_wind = 5160 K,
H:He = 1e-3, using the Phase-A2 spectra and the frozen p-winds clone at
../../../p-winds (commit bd01d2f, 2024-08-16).

Every parameter is from ../system_parameters.md; the SEDs are the v23-built
(paper-era Mega-MUSCLES) products of ../sed.  Outputs, one set per case:

  tspec_<case>.txt   vacuum wavelength [A], transmission, excess absorption
                     [%], instrument-convolved (R = 68,000) excess [%]
  profile_<case>.txt r [R_p], v [km/s], rho [g/cm3], f_HII, f_He3 [triplet
                     fraction], n_He3 [cm^-3]

Comparison targets (2024 detection): blended-red excess 1.24 +0.22/-0.23 %,
blue 0.25 +0.14/-0.12 %, FWHM 0.86 A.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '..', '..', '..', 'p-winds')))
import astropy.units as u
import astropy.constants as c
from p_winds import parker, hydrogen, helium, transit, lines, tools

# ---------------- system parameters (system_parameters.md) ----------------
R_pl_earth = 1.730          # [R_Earth]
M_pl_earth = 5.60           # [M_Earth]
M_star     = 0.1844         # [M_sun]
R_star_sun = 0.2159         # [R_sun]
a_au       = 0.0946         # [AU]
inc_deg    = 89.96          # [deg]

# retrieval point (the oracle target)
T_wind   = 5160.0           # [K]
mdot     = 2.03e8           # [g/s]
h_he     = 1.0e-3           # H:He number ratio
h_frac   = h_he/(1.0 + h_he)   # H fraction of nuclei

R_pl = (R_pl_earth*c.R_earth/c.R_jup).decompose().value   # [R_jup]
M_pl = (M_pl_earth*c.M_earth/c.M_jup).decompose().value   # [M_jup]
R_pl_m = (R_pl_earth*c.R_earth).to(u.m).value
R_star_m = (R_star_sun*c.R_sun).to(u.m).value
b_impact = (a_au*u.au*np.cos(np.radians(inc_deg))/(R_star_sun*c.R_sun)) \
    .decompose().value

# case -> (SED file, T_wind [K]).  'published' = the paper's best-fit
# parameter point (Mdot 2.03e8, T 5160 K, H:He 1e-3); both v23-built SEDs.
# 'matched_gj1132' is the T inside the paper's explored 5160-6400 K range at
# which our reconstruction reproduces the observed red depth (README.md).
CASES = {'published_gj1132': ('../sed/lhs1140_sed_gj1132_at_b.txt', 5160.0),
         'matched_gj1132':   ('../sed/lhs1140_sed_gj1132_at_b.txt', 6100.0),
         'published_gj699':  ('../sed/lhs1140_sed_gj699_at_b.txt', 5160.0)}

# instrument (WINERED): Gaussian LSF, R = 68,000
RESOLVING_POWER = 68000.0
AIR_TO_VAC = 10832.057/10829.09114   # from the triplet line 0 (vacuum/air)

def convolve_lsf(wl_m, spec):
    """Gaussian instrument convolution at R = RESOLVING_POWER."""
    wl0 = np.mean(wl_m)
    fwhm = wl0/RESOLVING_POWER
    sig = fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
    dl = np.median(np.diff(wl_m))
    x = np.arange(-5*sig, 5*sig + dl, dl)
    k = np.exp(-0.5*(x/sig)**2); k /= k.sum()
    return np.convolve(spec, k, mode='same')


# ---- broadening test: matched_gj1132 physics, three RT treatments ----
import sys
sys.path.insert(0, '../..')
from he_line_metrics import fit_metrics

case, (sedfile, T_wind) = 'matched_gj1132', CASES['matched_gj1132']
sedfile, T_wind = CASES['matched_gj1132']
spectrum = tools.make_spectrum_from_file(
    sedfile, units={'wavelength': u.angstrom,
                    'flux': u.erg/u.s/u.cm**2/u.angstrom}, skiprows=0)
r = np.logspace(0, np.log10(20.0), 200)
f_ion, mu_bar = hydrogen.ion_fraction(
    r, R_pl, T_wind, h_frac, mdot, M_pl,
    mean_molecular_weight_0=4.0*(1-h_frac) + 1.0*h_frac,
    star_mass=M_star, semimajor_axis=a_au,
    spectrum_at_planet=spectrum, relax_solution=True, return_mu=True)
vs = parker.sound_speed(T_wind, mu_bar)
rs = parker.radius_sonic_point(M_pl, vs)
rhos = parker.density_sonic_point(mdot, rs, vs)
r_array = r*R_pl/rs
v_array, rho_array = parker.structure(r_array)
f_he_1, f_he_3 = helium.population_fraction(
    r, v_array, rho_array, f_ion,
    R_pl, T_wind, h_frac, vs, rs, rhos, spectrum_at_planet=spectrum,
    initial_state=np.array([1.0, 0.0]), relax_solution=True)
m_h = c.m_p.to(u.g).value
he_frac = 1.0 - h_frac
n_nuclei = rho_array*rhos/(m_h*(h_frac + 4.0*he_frac))
n_he_3 = he_frac*n_nuclei*f_he_3

grid, t_depth, r_from_planet = transit.draw_transit(
    planet_to_star_ratio=R_pl_m/R_star_m,
    planet_physical_radius=R_pl_m,
    impact_parameter=b_impact,
    grid_size=121, supersampling=10)
l0, l1, l2, f0, f1, f2, a_ij = lines.he_3_properties()
wl_grid = np.linspace(1.0826e-6, 1.0838e-6, 400)
m_He = 4.002602*1.66053907e-27
r_m = r*R_pl_m
n_he3_m3 = n_he_3*1e6
v_ms = v_array*vs*1e3

print(f"matched_gj1132 physics (T={T_wind:.0f} K, Mdot={mdot:.2e}, H:He={h_frac/(1-h_frac):.4g})")
for method, turb, tag in (('average', False, 'average'),
                          ('formal',  False, 'formal'),
                          ('formal',  True,  'formal_turb')):
    t_spec = transit.radiative_transfer_2d(
        grid, r_from_planet, r_m, n_he3_m3, v_ms,
        np.array([l0, l1, l2]), np.array([f0, f1, f2]),
        np.array([a_ij, a_ij, a_ij]),
        wl_grid, T_wind, m_He, wind_broadening_method=method,
        turbulence_broadening=turb)
    excess = (1.0 - t_spec/np.max(t_spec))*100.0
    excess_conv = convolve_lsf(wl_grid, excess)
    wl_vac_A = wl_grid*1e10*AIR_TO_VAC
    m = fit_metrics(wl_vac_A, excess_conv, frame='vacuum')
    print(f"  {tag:12s} red {m['red_depth']:6.3f} %  blue {m['blue_depth']:6.3f} %"
          f"  r/b {m['red_blue']:5.2f}  FWHM {m['fwhm_A']:.3f} A")
    np.savetxt(f'tspec_broadening_{tag}.txt',
               np.column_stack([wl_vac_A, t_spec, excess, excess_conv]),
               header='vacuum wavelength [A]  transmission  excess[%]  '
                      'excess convolved R=68000 [%]  ('+tag+')')
print('  observed     red  1.254 %  blue  0.198 %  r/b  6.35  FWHM 0.841 A')
