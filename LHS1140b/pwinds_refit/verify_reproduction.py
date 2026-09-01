#!/usr/bin/env python3
"""Validation of the refit forward model against the authors' released grid.

Three checks:

1. The likelihood at the authors' grid argmax reproduces the value stored in
   their ``lhs1140b_grid_likelihood_results_1132_v26e3.pickle``.
2. The same for a handful of other grid nodes spread over the cube.
3. The ``v_nt = 0`` effective-temperature substitution is numerically
   identical to calling p-winds with ``turbulence_broadening=True``, and the
   substitution reproduces a hand-computed Doppler width for v_nt > 0.

Run: /opt/miniconda3/bin/python3 verify_reproduction.py
"""
import pickle
import numpy as np
import refit_lib as R

GRID = R.ZENODO + 'lhs1140b_grid_likelihood_results_1132_v26e3.pickle'

x, y, yerr = R.load_data()
kernel = R.instrument_kernel(x)

with open(GRID, 'rb') as fh:
    g = pickle.load(fh)
lv = g['likelihood_values']
log_mdot_values = g['log_mdot_values']
temp_values = g['temp_values']
h_values = g['h_fraction_values']

print('=== check 1/2: authors\' grid nodes ===')
print(f'{"i,j,k":>10s} {"log Mdot":>9s} {"T [K]":>8s} {"h_frac":>10s} '
      f'{"authors lnL":>13s} {"ours":>13s} {"diff":>10s}')
nodes = [np.unravel_index(np.nanargmax(lv), lv.shape),
         (8, 6, 9), (7, 6, 9), (9, 6, 9), (8, 5, 9), (8, 7, 9),
         (8, 6, 8), (8, 6, 10), (10, 10, 5), (5, 3, 12)]
seen = set()
for idx in nodes:
    if idx in seen:
        continue
    seen.add(idx)
    i, j, k = idx
    theta = (log_mdot_values[i], np.log10(temp_values[j]), h_values[k])
    ours = R.log_likelihood(theta, x, y, yerr, kernel)
    ref = lv[i, j, k]
    print(f'{str(idx):>10s} {theta[0]:9.4f} {temp_values[j]:8.1f} '
          f'{h_values[k]:10.3e} {ref:13.4f} {ours:13.4f} {ours - ref:10.2e}')

print()
print('=== check 3: the v_nt substitution ===')
best = (np.log10(g['best_fit']['mdot']), np.log10(g['best_fit']['temp']),
        g['best_fit']['h_fraction'])
T = 10 ** best[1]
n3, v = R.atmospheric_model(*best)
wl = x

s_auth = R.transmission_model(wl, n3, T, v, authors_turbulence=True)
s_sub = R.transmission_model(wl, n3, T, v, v_nt=0.0)
print(f'authors turbulence=True vs v_nt=0 substitution: '
      f'max |diff| = {np.max(np.abs(s_auth - s_sub)):.3e}')

s_off = R.transmission_model(wl, n3, T, v, authors_turbulence=False)
print(f'turbulence=False (no substitution)             : '
      f'max |diff| vs authors = {np.max(np.abs(s_auth - s_off)):.3e}'
      '   (must be non-zero: the switch does something)')

# Doppler width bookkeeping
kT_m = R.k_B * T / R.m_He
for v_nt in (0.0, 5e3, 10e3):
    T_rt = R.rt_temperature(T, v_nt)
    sig_sub = np.sqrt(R.k_B * T_rt / R.m_He)
    sig_want = np.sqrt(11.0 / 6.0 * kT_m + v_nt ** 2)
    print(f'v_nt = {v_nt/1e3:5.1f} km/s : sigma(sub) = {sig_sub/1e3:7.4f} '
          f'km/s, sigma(expected) = {sig_want/1e3:7.4f} km/s, '
          f'rel diff {abs(sig_sub/sig_want - 1):.2e}')
print('(the wind term v_wind_broad enters both sides identically and is not '
      'shown here)')
