#!/usr/bin/env python3
"""p-winds forward model and likelihood for LHS 1140 b He 10830, with the
non-thermal broadening promoted to a free parameter.

The forward model is the one released by the authors of Cherubim et al. (2026)
(Zenodo 15723779, ``lhs1140b_2024B_pwinds.py``), reproduced line for line:
same SED file and normalization (const-res v23 GJ 1132 x Lx-ratio 0.59), same
system parameters, same altitude grid, same 5-phase transit averaging, same
fixed bulk line-of-sight velocity v_wind = 2.26 km/s, same instrumental
profile (4.4 km/s FWHM Gaussian) and the same released spectrum with its GP
component added back.  ``log_likelihood`` is theirs, including the
``log(2 pi sigma^2)`` term, so the numbers here are directly comparable with
the grid stored in ``lhs1140b_grid_likelihood_results_1132_v26e3.pickle``.

The one addition is ``v_nt``: an isotropic non-thermal velocity dispersion
(1-D Gaussian sigma, m/s) added in quadrature to the Doppler width of the
line, i.e. to the opacity, not to the emergent spectrum.

How it is applied.  With ``wind_broadening_method='average'`` p-winds builds
the Gaussian width of the Voigt profile as

    sigma_v^2 = k T / m + v_wind_broad^2 + v_turb^2,
    v_turb^2  = (5/6) k T / m     (Lampon et al. 2020, the authors' switch),

and in that method the temperature argument of ``radiative_transfer_2d``
enters nowhere else.  Passing ``turbulence_broadening=False`` together with an
effective temperature

    T_rt = (11/6) T + m_He v_nt^2 / k_B

therefore reproduces the authors' width exactly at v_nt = 0 and adds v_nt in
quadrature otherwise.  ``verify_reproduction.py`` checks both statements
numerically.  The atmospheric model always sees the true T.
"""

import os
import numpy as np
import astropy.constants as c
import astropy.units as u
from astropy.io import fits
from astropy.convolution import convolve
from p_winds import parker, hydrogen, helium, transit, lines

ZENODO = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/LHS1140b_zenodo/'

# --- constants, exactly as in the authors' script -------------------------
pc2m = 3.0857e16
au2m = 1.496e11
Rjup = 7.1492e7
d_1140 = 14.96 * pc2m
d_1132 = 12.613 * pc2m
smax_1140 = 0.0946           # au

R_pl = 0.154                 # Jupiter radii
M_pl = 0.0176                # Jupiter masses
planet_to_star_ratio = 0.072092
impact_parameter = 0.23

m_h = c.m_p.to(u.g).value            # g
m_He = 4 * 1.67262192369e-27         # kg
k_B = 1.380649e-23                   # J/K

r = np.logspace(0, np.log10(10), 100)          # planetary radii
initial_f_ion = 0.0
initial_f_he = np.array([1.0, 0.0])
relax_solution = True
exact_phi = True
sample_phases = np.linspace(-0.50, 0.50, 5)
n_samples = len(sample_phases)
transit_grid_size = 100
supersampling = 5
V_WIND = 2.26e3                      # m/s, fixed by the authors

w0, w1, w2, f0, f1, f2, a_ij = lines.he_3_properties()
w_array = np.array([w0, w1, w2])
f_array = np.array([f0, f1, f2])
a_array = np.array([a_ij, a_ij, a_ij])

# --- host spectrum: const-res v23 GJ 1132, Lx ratio 0.59 ------------------
Lx_ratio_1132 = 0.59
norm_constant_1132 = Lx_ratio_1132 * (d_1140 / (smax_1140 * au2m)) ** 2
_spec = fits.getdata(
    os.path.join(ZENODO,
                 'hlsp_muscles_multi_multi_gj1132_broadband_v23_const-res-sed.fits'), 1)
host_spectrum = {'wavelength': _spec['WAVELENGTH'],
                 'flux_lambda': _spec['FLUX'] * norm_constant_1132,
                 'wavelength_unit': u.angstrom,
                 'flux_unit': u.erg / u.s / u.cm ** 2 / u.angstrom}


# --- observed spectrum ----------------------------------------------------
def load_data():
    """Released WINERED 2024B in-transit spectrum, prepared as the authors do.

    Returns (x [m, air], y, yerr) ready for ``log_likelihood``.
    """
    import pickle
    with open(os.path.join(ZENODO, 'lhs1140b_trans_spec_GP_IT_vac.pickle'),
              'rb') as fh:
        he_spec = pickle.load(fh)
    wl_obs = he_spec['wavelength']
    s = 1e4 / wl_obs
    n = (1 + 0.0000834254 + 0.02406147 / (130 - s ** 2)
         + 0.00015998 / (38.9 - s ** 2))
    wl_obs = wl_obs / n                      # vacuum -> air
    f_obs = he_spec['flux'] * 0.01 + 1
    u_obs = he_spec['error'] * 0.01
    f_obs = f_obs + he_spec['GP'] * 0.01
    return wl_obs / 1e10, f_obs, u_obs


def instrument_kernel(wl_obs_m):
    """Authors' instrumental profile: 4.4 km/s FWHM Gaussian on the data grid."""
    wl_obs = wl_obs_m * 1e10
    width_v = 4.4                                     # km/s FWHM
    sigma_wl = (width_v / (2 * (2 * np.log(2)) ** 0.5)
                / c.c.to(u.km / u.s).value * np.mean(wl_obs))
    x = wl_obs
    mu = np.mean(wl_obs)
    return (1 / sigma_wl / (2 * np.pi) ** 0.5
            * np.exp(-0.5 * (x - mu) ** 2 / sigma_wl ** 2))


# --- transit ray tracing: parameter independent, computed once ------------
_R_pl_physical = R_pl * Rjup
_F_MAPS = None
_T_DEPTHS = None
_R_MAPS = None


def _ray_tracing():
    global _F_MAPS, _T_DEPTHS, _R_MAPS
    if _F_MAPS is None:
        fm, td, rm = [], [], []
        for i in range(n_samples):
            a, b, cc = transit.draw_transit(
                planet_to_star_ratio, impact_parameter=impact_parameter,
                supersampling=supersampling, phase=sample_phases[i],
                planet_physical_radius=_R_pl_physical,
                grid_size=transit_grid_size)
            fm.append(a)
            td.append(b)
            rm.append(cc)
        _F_MAPS, _T_DEPTHS, _R_MAPS = fm, td, rm
    return _F_MAPS, _T_DEPTHS, _R_MAPS


# --- forward model --------------------------------------------------------
def atmospheric_model(log_m_dot, log_T, h_fraction):
    """Authors' ``atmospheric_model``: returns n(He 2^3S) [cm^-3] and v [km/s]."""
    m_dot = 10 ** log_m_dot
    T = 10 ** log_T
    he_fraction = 1 - h_fraction
    he_h_fraction = he_fraction / h_fraction
    mean_f_ion = 0.90
    mu_0 = (1 + 4 * he_h_fraction) / (1 + he_h_fraction + mean_f_ion)

    f_r, mu_bar = hydrogen.ion_fraction(
        r, R_pl, T, h_fraction, m_dot, M_pl, mu_0,
        spectrum_at_planet=host_spectrum, initial_f_ion=initial_f_ion,
        relax_solution=relax_solution, exact_phi=exact_phi, return_mu=True)

    vs = parker.sound_speed(T, mu_bar)
    rs = parker.radius_sonic_point(M_pl, vs)
    rhos = parker.density_sonic_point(m_dot, rs, vs)
    v_array, rho_array = parker.structure(r * R_pl / rs)

    f_he_1, f_he_3 = helium.population_fraction(
        r, v_array, rho_array, f_r, R_pl, T, h_fraction, vs, rs, rhos,
        spectrum_at_planet=host_spectrum, initial_state=initial_f_he,
        relax_solution=relax_solution)

    n_he = rho_array * rhos * he_fraction / (h_fraction + 4 * he_fraction) / m_h
    return f_he_3 * n_he, v_array * vs


def rt_temperature(T, v_nt):
    """Effective temperature that carries thermal + Lampon turbulence + v_nt."""
    return (11.0 / 6.0) * T + m_He * v_nt ** 2 / k_B


def transmission_model(wavelength_array, n_he_3, T, v_array, v_nt=0.0,
                       v_wind=V_WIND, broadening='average',
                       authors_turbulence=None):
    """5-phase averaged transmission spectrum.

    ``v_nt`` is folded into the Doppler width through ``rt_temperature``.  Pass
    ``authors_turbulence=True/False`` to bypass that and call p-winds with the
    true temperature and its own ``turbulence_broadening`` switch (used by the
    verification script).
    """
    r_SI = r * _R_pl_physical
    v_SI = v_array * 1000
    n_he_3_SI = n_he_3 * 1e6
    f_maps, t_depths, r_maps = _ray_tracing()

    if authors_turbulence is None:
        temp_rt = rt_temperature(T, v_nt)
        turb = False
    else:
        temp_rt = T
        turb = authors_turbulence

    spectra = []
    for i in range(n_samples):
        spec = transit.radiative_transfer_2d(
            f_maps[i], r_maps[i], r_SI, n_he_3_SI, v_SI, w_array, f_array,
            a_array, wavelength_array, temp_rt, m_He,
            bulk_los_velocity=v_wind, wind_broadening_method=broadening,
            turbulence_broadening=turb)
        spectra.append(spec + t_depths[i])
    return np.mean(np.array(spectra), axis=0)


def cascading_model(theta, wavelength_array, kernel):
    """(log_mdot, log_T, h_fraction[, v_nt[, v_wind]]) -> convolved spectrum.

    ``v_nt`` defaults to 0 (the authors' broadening) and ``v_wind`` to the
    2.26 km/s the authors hold fixed in their released grid search.
    """
    log_m_dot, log_T, h_fraction = theta[0], theta[1], theta[2]
    v_nt = theta[3] if len(theta) > 3 else 0.0
    v_wind = theta[4] if len(theta) > 4 else V_WIND
    try:
        n_he_3, v = atmospheric_model(log_m_dot, log_T, h_fraction)
        t_spec = transmission_model(wavelength_array, n_he_3, 10 ** log_T, v,
                                    v_nt=v_nt, v_wind=v_wind)
        return convolve(t_spec, kernel, boundary='extend')
    except Exception:
        return np.full_like(wavelength_array, 1.0)


def log_likelihood(theta, x, y, yerr, kernel):
    """Authors' Gaussian log-likelihood, constant term included."""
    try:
        model = cascading_model(theta, x, kernel)
        sigma2 = yerr ** 2
        return -0.5 * np.sum((y - model) ** 2 / sigma2
                             + np.log(2 * np.pi * sigma2))
    except RuntimeError:
        return -np.inf


def chi2_from_loglike(loglike, yerr):
    """chi^2 corresponding to the authors' log-likelihood definition."""
    return -2.0 * loglike - np.sum(np.log(2 * np.pi * yerr ** 2))
