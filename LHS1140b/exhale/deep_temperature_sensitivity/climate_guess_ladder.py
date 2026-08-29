#!/usr/bin/env python3
"""The climate solve alone, from a ladder of initial guesses.

`radiative_convective_column.solve_radiative_convective_column` is called at
the LHS 1140 b configuration of `nh_refusal_diagnosis/scan_reservoir_pc090.sh`
and nothing else, so what varies between rows is only the initial guess handed
to `AdiabatClimate.surface_temperature_bg_gas`.  The stellar flux file is
written by the adapter's own `write_flux_file`, so the spectrum the solve sees
is the one the stored columns were made with.

The columns are the deep boundary temperature, the tropopause the
pseudoadiabat meets the stratosphere at, and the water the solution carries
above it, each to full precision: the spread down a column is the
reproducibility of the climate step by itself, before any chemistry has run.

usage: climate_guess_ladder.py [He/H]
"""
import argparse
import os
import sys

import numpy as np

EX = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00'
sys.path.insert(0, os.path.join(EX, 'src', 'utils'))
import lower_profile_schema as sch                           # noqa: E402
import photochem_to_lower_profile as ad                      # noqa: E402
import radiative_convective_column as rcc                    # noqa: E402

SOLAR = dict(ad.SOLAR_ABUNDANCES)
GUESSES = [400.0, 400.000001, 400.0001, 400.01, 400.1, 401.0, 410.0, 450.0,
           350.0, 300.0]

HEH = float(sys.argv[1]) if len(sys.argv) > 1 else 9.0
here = os.path.dirname(os.path.abspath(__file__))
work = os.path.join(here, 'climate_work_heh%s' % str(HEH).replace('.', 'p'))
os.makedirs(work, exist_ok=True)

ab = dict(SOLAR)
ab.pop('S', None)
ab['He'] = HEH

flux_args = argparse.Namespace(
    stellar_flux=os.path.join(EX, 'LHS1140b', 'sed',
                              'lhs1140_sed_gj1132_at_b.txt'),
    wavelength_unit='A', flux_at_planet=True, r_star=None, a_orb=None)
flux_file = os.path.join(work, 'stellar_flux.txt')
ad.write_flux_file(flux_args, flux_file)

print('# LHS 1140 b, He/H = %g, deep boundary 20 bar, 60 climate layers' % HEH)
print('# %-12s %-24s %-24s %-24s %-24s'
      % ('T_guess[K]', 'T_deep[K]', 'P_trop[bar]', 'T_trop[K]', 'f_H2O_trop'))
rows = []
for g in GUESSES:
    try:
        sol = rcc.solve_radiative_convective_column(
            work, flux_file, 0.0176220*sch.MJ, 0.157692*sch.RJ, ab, 20.0,
            p_top_dyn=1.0e-2, nlayers=60, t_deep_guess=g, verbose=False)
    except BaseException as exc:
        print('  %-12g %s' % (g, str(exc).splitlines()[0][:70]))
        sys.stdout.flush()
        continue
    rows.append((g, sol))
    print('  %-12.6f %-24.17g %-24.17g %-24.17g %-24.17g'
          % (g, sol['T_deep'], sol['P_trop_bar'], sol['T_trop'],
             sol['f_H2O_trop']))
    sys.stdout.flush()

if rows:
    td = np.array([s['T_deep'] for _, s in rows])
    pt = np.array([s['P_trop_bar'] for _, s in rows])
    ft = np.array([s['f_H2O_trop'] for _, s in rows])
    print('#')
    print('# T_deep     spread %.3e K   (%.3e relative)'
          % (td.max()-td.min(), (td.max()-td.min())/td.mean()))
    print('# P_trop     spread %.3e bar (%.3e relative)'
          % (pt.max()-pt.min(), (pt.max()-pt.min())/pt.mean()))
    print('# f_H2O_trop spread %.3e     (%.3e relative)'
          % (ft.max()-ft.min(), (ft.max()-ft.min())/ft.mean()))
