"""Shared setup for the background-pressure bracket diagnosis.

Builds the same `AdiabatClimate` object the LHS 1140 b closure ladder builds
(`src/utils/radiative_convective_column.py` through
`photochem_to_lower_profile.py --climate --climate-p-deep 20.0`), so the
residual measured here is the one `make_profile_bg_gas` sees.
"""
import os
import sys

import numpy as np

EX = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00'
sys.path.insert(0, os.path.join(EX, 'src', 'utils'))

import radiative_convective_column as rcc          # noqa: E402
from photochem_to_lower_profile import SOLAR_ABUNDANCES  # noqa: E402

MJ = 1.898e30      # placeholder, replaced below from the schema
FLUX = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'flux_photochem.txt')
P_DEEP_BAR = 20.0
P_TOP_DYN = 1.0e-2
NLAYERS = 60
T_DEEP_GUESS = 400.0
T_TROP_GUESS = 120.0

import lower_profile_schema as sch                 # noqa: E402
MP_G = 0.0176220*sch.MJ
RP_CM = 0.157692*sch.RJ


def build(heh, workdir, verbose=False):
    """The climate object, the surface partial pressures and the background
    gas index, for reservoir He/H = `heh`."""
    from photochem.clima import AdiabatClimate
    os.makedirs(workdir, exist_ok=True)
    ab = dict(SOLAR_ABUNDANCES)
    ab['He'] = float(heh)
    mixing = rcc.deep_molecular_partition(ab)
    bg = rcc.background_gas(mixing)
    sp_file, st_file = rcc.write_climate_files(workdir, MP_G, RP_CM, 0.0,
                                               NLAYERS)
    c = AdiabatClimate(sp_file, st_file, FLUX)
    c.solve_for_T_trop = True
    c.T_trop = T_TROP_GUESS
    c.P_top = P_TOP_DYN
    c.verbose = verbose
    names = list(c.species_names)
    P_deep = P_DEEP_BAR*rcc.BAR
    P_i = np.array([mixing.get(sp, 0.0)*P_deep for sp in names])
    return c, P_i, bg, names.index(bg), P_deep, mixing


# --- the patched Fortran algorithm, transcribed ---------------------------
NSCAN = 49
LOG_SPAN = 12.0
PTOL = 1.0e-10
MAXIT = 100


def residual(c, P_i, ind, P_surf, T_surf, x):
    """`pressure_residual` of the patch: place 10**x in the background slot,
    build the profile, return (f, valid, message)."""
    P = P_i.copy()
    P[ind] = 10.0**x
    try:
        c.make_profile(T_surf, P)
    except Exception as e:
        return float('nan'), False, str(e)
    return (c.P_surf - P_surf)/P_surf, True, ''


def scan(c, P_i, ind, P_surf, T_surf, nscan=NSCAN, span=LOG_SPAN):
    """Return the full scan table [(x, f, valid, message), ...]."""
    out = []
    x0 = np.log10(P_surf)
    for i in range(nscan):
        x = x0 - span + span*i/(nscan - 1)
        f, ok, msg = residual(c, P_i, ind, P_surf, T_surf, x)
        out.append((x, f, ok, msg))
    return out


def bracket_from_scan(rows, ptol=PTOL):
    """Replay the patch's bracketing decision on a scan table.

    Returns (status, detail).  status is 'root', 'bracket' or 'none'.
    """
    prev_valid = False
    xp = fp = None
    for x, f, ok in [(r[0], r[1], r[2]) for r in rows]:
        if not ok:
            prev_valid = False
            continue
        if abs(f) <= ptol:
            return 'root', (x, f)
        if prev_valid and f*fp < 0.0:
            return 'bracket', (xp, x, fp, f)
        xp, fp, prev_valid = x, f, True
    return 'none', None
