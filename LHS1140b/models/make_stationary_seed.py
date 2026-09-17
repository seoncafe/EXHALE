#!/usr/bin/env python3
"""Rebuild an EXHALE restart seed with the velocity a stationary wind carries.

Why this exists
---------------
An archived EXHALE state can be stationary in the wind and not stationary in
the base layer.  On LHS 1140 b the archived solution of
`LHS1140b/examples/scalar_base_kzz1e9` carries v = -538 cm/s in its first
cell and a mass flux rho v r^2 that stands 1.8e4 above the wind's value
there (`docs/lhs1140b_stationary_L1_20260913.md` sections 1 and 4.1).  The
stationary velocity of that layer is fixed by continuity alone,

    v(r) = F_wind / (rho(r) r^2),

and on this planet it is 0.0297 cm/s at the base, Mach 2.7e-7
(`docs/lhs1140b_stationary_L3_20260913.md` section 2, table T1).  A seed
carrying the archived thermodynamics and THIS velocity is the state item L5a
of `docs/PLAN_20260913_lhs_stationary.md` hands to the partitioned
stationary solve.

What sets the loaded state, and what therefore has to be edited
---------------------------------------------------------------
`load_IC` builds the primitive state as `W = (rho, v, p)` with

* `rho` from `calc_rho` of the SPECIES DENSITIES of `Ion_species_IC.txt`
  (`src/modules/files_IO/load_IC.f90` near line 799), not from the density
  column of `Hydro_ioniz_IC.txt`, which is read into a scratch variable and
  discarded (lines 348 and 351);
* `v` and `p` from `Hydro_ioniz_IC.txt`;
* the temperature NOT from the T column at all: T follows as p over the
  particle count the composition carries.

So a seed that is to hold a prescribed density must scale the species
columns, and one that is to hold a prescribed temperature must set the
pressure to the particle count times that temperature.  Editing the density
and temperature columns of the hydrodynamic file alone changes nothing but
the file: MEASURED on this state, a seed whose ln p and ln T were filtered
in cells 1 to 30 loaded with a density identical to the unfiltered one in
every cell and an odd-even component of ln T of 4.82e-2 against the
unfiltered 4.73e-2.  This script therefore writes both files and keeps them
consistent: the species of a cell are scaled by one factor, which carries
the density and the particle count together, and the pressure is then set
from that particle count and the temperature the seed is to hold.

What the script writes
----------------------
`<dst>/Hydro_ioniz_IC.txt` and `<dst>/Ion_species_IC.txt`.

`--base-anchor A` scales the whole column, species and pressure together, by
A before anything else.  That sets the pressure the lower boundary sees: the
base face velocity the characteristic condition writes is
`v_b = v_i + (p_b - p_i)/(rho_i c_i)` (`base_boundary.f90` line 352), so at
Mach 1e-7 a column whose base face pressure misses the reservoir's by a few
percent enters the domain at hundreds of cm/s whatever its interior velocity
is.  On the LHS 1140 b archived state that offset is 2.92 percent and the
inflow is -1858 cm/s; A = 0.971649 takes it to +0.12 cm/s.  A is calibrated
against the code's own boundary condition, not derived here: two step-free
post-processing passes (`Do only PP: True` with `CFL: 1.0e-12`) at two values
of A give the base ghost velocity, which is linear in A, and one secant step
lands it on the stationary velocity.  Apply it BEFORE `--odd-even` at your
peril: the filter moves rho at cell 1 by percent and undoes the anchor.

`--odd-even` removes the odd-even component of ln rho and ln T over the
first cells.  That component reaches 4.8e-2 in ln rho at cell 3 of the
LHS 1140 b seed and, in a layer at Mach 2.5e-6, enters the interface mass
flux through its acoustic part as d/Mach times rho v
(`docs/lhs1140b_stationary_L1_20260913.md` section 4).  The filter is the
one L1 used: four passes of a half-odd-even subtraction, at full strength
over cells 1 to `--oe-cells` and tapered linearly to zero `--oe-taper`
cells further out.

Usage
-----
    make_stationary_seed.py <src_output_dir> <dst_output_dir> [options]

`<src_output_dir>` holds `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` (or
the `Hydro_ioniz.txt` / `Ion_species.txt` of a finished run, which are
accepted under those names too).
"""

import argparse
import os
import sys

import numpy as np

# Grid ghost rows per end, EXHALE's Ng (src/modules/init/define_grid.f90).
NGHOST = 2

# Columns of Hydro_ioniz*.txt, schema 2.
C_R, C_RHO, C_V, C_P, C_T = 0, 1, 2, 3, 4

K_B = 1.380649e-16      # erg K^-1, src/modules/init/parameters.f90


def read_table(path):
    """Return (header lines, data array) of an EXHALE profile file."""
    header, rows = [], []
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith('#'):
                header.append(line)
            elif line.strip():
                rows.append(line.rstrip('\n'))
    data = np.array([[float(x) for x in row.split()] for row in rows])
    return header, data


def write_table(path, header, data):
    with open(path, 'w') as fh:
        for line in header:
            fh.write(line)
        for row in data:
            fh.write(' ' + '  '.join('%.17E' % x for x in row) + '\n')


def pick(src_dir, stem):
    """The restart name if it is there, else the finished-run name."""
    for name in (stem + '_IC.txt', stem + '.txt'):
        path = os.path.join(src_dir, name)
        if os.path.exists(path):
            return path
    raise SystemExit('no %s_IC.txt or %s.txt in %s' % (stem, stem, src_dir))


def far_wind_flux(r, flux, r_far):
    """The wind's mass flux: the median of rho v r^2 outside r_far."""
    far = r >= r_far
    if far.sum() < 3:
        raise SystemExit('fewer than 3 cells outside r = %g R_p' % r_far)
    return float(np.median(flux[far])), int(far.sum())


def matching_radius(r, flux, f_wind, band):
    """The innermost cell from which the archived flux stays within `band`.

    Not the first crossing: the base layer of an unrelaxed state crosses the
    wind's value accidentally while oscillating around it by decades.  The
    radius wanted is the one beyond which the archived state IS the wind.
    """
    off = np.abs(flux / f_wind - 1.0) > band
    if not off.any():
        return 0, r[0]
    last = int(np.max(np.where(off)[0]))
    if last + 1 >= len(r):
        raise SystemExit('the archived flux is outside the band in the '
                         'outermost cell; no matching radius exists')
    return last + 1, float(r[last + 1])


def kill_odd_even(x, weight, npass):
    """Subtract the odd-even component of ln x, weighted cell by cell."""
    lx = np.log(x).copy()
    for _ in range(npass):
        odd = np.zeros_like(lx)
        odd[1:-1] = lx[1:-1] - 0.5 * (lx[:-2] + lx[2:])
        lx = lx - 0.5 * odd * weight
    return np.exp(lx)


def odd_even_amplitude(x, lo, hi):
    lx = np.log(x)
    odd = lx[1:-1] - 0.5 * (lx[:-2] + lx[2:])
    return float(np.max(np.abs(odd[lo - 1:hi - 1])))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('src_dir', help='directory holding the archived state')
    ap.add_argument('dst_dir', help='directory the seed is written into')
    ap.add_argument('--base-anchor', type=float, default=1.0,
                    help='scale the whole column, species and pressure '
                         'together, by this factor before anything else, so '
                         'that its base face meets the reservoir (default 1, '
                         'no change)')
    ap.add_argument('--r-far', type=float, default=2.0,
                    help='radius [R_p] beyond which the state is the wind, '
                         'over which F_wind is the median (default 2.0)')
    ap.add_argument('--band', type=float, default=0.10,
                    help='fractional band around F_wind that defines the '
                         'matching radius (default 0.10)')
    ap.add_argument('--blend', type=int, default=10,
                    help='cells over which the stationary velocity is '
                         'blended linearly into the archived one, ending at '
                         'the matching radius (default 10)')
    ap.add_argument('--odd-even', action='store_true',
                    help='also remove the odd-even component of ln rho and '
                         'ln T in the first cells')
    ap.add_argument('--oe-cells', type=int, default=30,
                    help='cells filtered at full strength (default 30)')
    ap.add_argument('--oe-taper', type=int, default=15,
                    help='cells over which the filter tapers to zero '
                         '(default 15)')
    ap.add_argument('--oe-passes', type=int, default=4,
                    help='passes of the half-odd-even subtraction '
                         '(default 4)')
    args = ap.parse_args(argv)

    hydro_path = pick(args.src_dir, 'Hydro_ioniz')
    ion_path = pick(args.src_dir, 'Ion_species')
    header, d = read_table(hydro_path)
    ion_header, ion = read_table(ion_path)
    nrow = d.shape[0]
    if ion.shape[0] != nrow:
        raise SystemExit('the two files have different row counts')
    ncell = nrow - 2 * NGHOST
    inner = slice(NGHOST, NGHOST + ncell)      # physical cells 1..ncell

    r = d[:, C_R]
    rho = d[:, C_RHO].copy()                   # m_H cm^-3, sum n_s m_s
    v = d[:, C_V].copy()
    p = d[:, C_P].copy()
    T = d[:, C_T].copy()
    n_part = p / (K_B * T)                     # cm^-3, gas plus electrons

    note = []
    if args.base_anchor != 1.0:
        rho *= args.base_anchor
        p *= args.base_anchor
        n_part *= args.base_anchor
        ion[:, 1:] *= args.base_anchor
        note.append('the whole column scaled by %.6f so that its base face '
                    'meets the reservoir' % args.base_anchor)
        print('base anchor: the column scaled by %.6f' % args.base_anchor)
    rho_new, T_new = rho.copy(), T.copy()

    if args.odd_even:
        w = np.zeros(nrow)
        for i in range(NGHOST, NGHOST + ncell):
            j = i - NGHOST + 1
            if j <= args.oe_cells:
                w[i] = 1.0
            else:
                w[i] = max(0.0, (args.oe_cells + args.oe_taper - j)
                           / float(args.oe_taper))
        rho_new = kill_odd_even(rho, w, args.oe_passes)
        T_new = kill_odd_even(T, w, args.oe_passes)
        lo, hi = NGHOST, NGHOST + args.oe_cells
        print('odd-even filter over cells 1 to %d: max |drho/rho| %.3e, '
              'max |dT/T| %.3e'
              % (args.oe_cells,
                 np.max(np.abs(rho_new[lo:hi] / rho[lo:hi] - 1.0)),
                 np.max(np.abs(T_new[lo:hi] / T[lo:hi] - 1.0))))
        print('odd-even of ln rho over those cells %.3e -> %.3e; of ln T '
              '%.3e -> %.3e'
              % (odd_even_amplitude(rho, NGHOST, NGHOST + args.oe_cells),
                 odd_even_amplitude(rho_new, NGHOST, NGHOST + args.oe_cells),
                 odd_even_amplitude(T, NGHOST, NGHOST + args.oe_cells),
                 odd_even_amplitude(T_new, NGHOST, NGHOST + args.oe_cells)))
        note.append('the odd-even component of ln rho and ln T removed over '
                    'cells 1 to %d, tapered to zero by cell %d, the density '
                    'through the species columns'
                    % (args.oe_cells, args.oe_cells + args.oe_taper))

    # The species of a cell all carry the same factor, so the density, the
    # particle count and every abundance ratio move together.
    scale = rho_new / rho
    ion[:, 1:] *= scale[:, None]
    n_part_new = n_part * scale
    p_new = n_part_new * K_B * T_new

    flux = rho_new * v * r ** 2
    f_wind, nfar = far_wind_flux(r[inner], flux[inner], args.r_far)
    k_c, r_c = matching_radius(r[inner], flux[inner], f_wind, args.band)
    cell_c = k_c + 1
    print('F_wind = %.6e [mH cm/s R_p^2], median of %d cells at r >= %g R_p'
          % (f_wind, nfar, args.r_far))
    print('matching radius: cell %d, r = %.5f R_p (the innermost cell from '
          'which |rho v r^2 / F_wind - 1| <= %.2f holds outward)'
          % (cell_c, r_c, args.band))

    v_stat = f_wind / (rho_new * r ** 2)
    v_new = v.copy()
    # Below the blend window: the stationary velocity.  Inside it: a linear
    # ramp to the archived one, which the window's outer end already equals
    # to within `band`.  Above: the archived wind, untouched.
    i_c = NGHOST + k_c
    i_0 = max(NGHOST, i_c - args.blend)
    v_new[:i_0] = v_stat[:i_0]
    for i in range(i_0, i_c + 1):
        s = (i - i_0) / float(max(1, i_c - i_0))
        v_new[i] = (1.0 - s) * v_stat[i] + s * v[i]
    # The inner ghosts carry the same continuity law as the cells they
    # border, so the ghost rows and cell 1 hold one mass flux.
    v_new[:NGHOST] = v_stat[:NGHOST]

    print('v [cm/s] cells 1 to 8, archived : '
          + ' '.join('%9.4g' % x for x in v[NGHOST:NGHOST + 8]))
    print('v [cm/s] cells 1 to 8, seed     : '
          + ' '.join('%9.4g' % x for x in v_new[NGHOST:NGHOST + 8]))
    fn = rho_new * v_new * r ** 2 / f_wind
    print('rho v r^2 / F_wind, seed, cells 1, 2, %d, %d, %d: %.4f %.4f '
          '%.4f %.4f %.4f'
          % (cell_c - args.blend, cell_c, ncell,
             fn[NGHOST], fn[NGHOST + 1], fn[i_0], fn[i_c],
             fn[NGHOST + ncell - 1]))

    d[:, C_RHO] = rho_new
    d[:, C_V] = v_new
    d[:, C_P] = p_new
    d[:, C_T] = T_new

    note.insert(0, 'v replaced by F_wind/(rho r^2) below cell %d '
                   '(r = %.5f R_p) with a %d-cell linear blend into the '
                   'archived velocity, F_wind = %.6e mH cm/s R_p^2 measured '
                   'as the median of rho v r^2 at r >= %g R_p'
                   % (cell_c, r_c, args.blend, f_wind, args.r_far))
    tail = ('# stationary seed: ' + '; '.join(note)
            + '; source ' + os.path.abspath(hydro_path)
            + '; made by LHS1140b/models/make_stationary_seed.py '
              '(PLAN_20260913_lhs_stationary item L5a)\n')

    os.makedirs(args.dst_dir, exist_ok=True)
    write_table(os.path.join(args.dst_dir, 'Hydro_ioniz_IC.txt'),
                list(header) + [tail], d)
    write_table(os.path.join(args.dst_dir, 'Ion_species_IC.txt'),
                list(ion_header) + [tail], ion)
    print('wrote %s/{Hydro_ioniz_IC.txt,Ion_species_IC.txt}' % args.dst_dir)
    return 0


if __name__ == '__main__':
    sys.exit(main())
