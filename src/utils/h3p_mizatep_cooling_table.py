#!/usr/bin/env python3
"""Build src/modules/lower_atmosphere/h3p_mizatep_cooling_table.f90 from the
ExoMol cooling-function file of the MiZATeP H3+ line list.

Usage
-----
    python3 h3p_mizatep_cooling_table.py <1H3_p__MiZATeP.cf> [<out.f90>]

The default output is the module path above, relative to the EXHALE root
(this script's ../..).

WHAT THE INPUT IS.  `1H3_p__MiZATeP.cf` from the ExoMol database
(https://www.exomol.com/db/H3_p/1H3_p/MiZATeP/1H3_p__MiZATeP.cf, retrieved
2026-10-10): the LTE cooling function W(T) of H3+ from the MiZATeP line list
of Mizus, Alijah, Zobov, Lodi, Kyuberis, Yurchenko, Tennyson & Polyansky
(2017, MNRAS 468, 1717), "total energy emitted by a single molecule in one
second per unit solid angle" (their section 5), in erg s^-1 sr^-1
molecule^-1 on a 1 K grid from 1 to 5000 K.  Two columns, T [K] and W; the
entries below 8 K are NAN.

WHICH ROWS ARE KEPT, AND WHY 20 K.  The file divides the emission sum by the
partition function of the companion `.pf` file, which is printed to four
decimals; below about 16 K that rounding is the dominant error (the .pf gives
Q(8 K) = 0.0001 against the summed 6.08e-5).  MEASURED against the full
MiZATeP emission sum divided by the summed Q (md/h3p_cooling_cause_20261010.md,
its miz_power_calc.npz): .cf/W = 0.608 at 8 K, 1.115 at 9 K, 1.058 at 10 K,
0.993 at 13 K, 0.997 at 15 K, 1.0005 at 17 K, 1.0001 at 20 K and within
1e-4 of 1 from 25 K up.  The table therefore starts at 20 K, which is also
the lowest temperature of Mizus et al. (2017) Table 7, and runs to 5000 K,
the upper limit of the paper ("at least 99 per cent complete for
temperatures up to 5000 K").  The values are copied verbatim (nine
significant figures, erg s^-1 sr^-1); the Fortran converts to W.
"""
import math
import os
import sys

T_FIRST = 20
T_LAST = 5000


def read_cf(path):
    rows = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) != 2:
                continue
            rows[float(parts[0])] = parts[1]
    return rows


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    cf_path = sys.argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        here, '..', 'modules', 'lower_atmosphere',
        'h3p_mizatep_cooling_table.f90')
    rows = read_cf(cf_path)
    values = []
    w_prev = 0.0
    for t in range(T_FIRST, T_LAST + 1):
        text = rows[float(t)]
        w = float(text)
        if not (math.isfinite(w) and w > 0.0):
            sys.exit(f'non-positive or nonfinite W at {t} K: {text}')
        if not w > w_prev:
            sys.exit(f'W not strictly increasing at {t} K')
        w_prev = w
        values.append(text.replace('E', 'd'))
    n = len(values)

    L = []
    L.append('      module h3p_mizatep_cooling_table')
    L.append('      ! The LTE cooling function W(T) of H3+ from the MiZATeP line list')
    L.append('      ! (Mizus et al. 2017, MNRAS 468, 1717): the power one molecule with')
    L.append('      ! an LTE internal population at T emits into one steradian,')
    L.append('      ! in erg s^-1 sr^-1 molecule^-1, at T = %d, %d, ..., %d K.'
             % (T_FIRST, T_FIRST + 1, T_LAST))
    L.append('      !')
    L.append('      ! GENERATED FILE.  Rebuild with')
    L.append('      !     python3 src/utils/h3p_mizatep_cooling_table.py <1H3_p__MiZATeP.cf>')
    L.append('      ! Edit the script, not this file.  Source: the ExoMol file')
    L.append('      ! 1H3_p__MiZATeP.cf (retrieved 2026-10-10), copied verbatim; the')
    L.append('      ! script says why the rows below 20 K are not used.  The one')
    L.append('      ! reader of this table is h3p_lte_emission_mizatep (h3p_cooling.f90).')
    L.append('      implicit none')
    L.append('      private')
    L.append('      public :: h3p_cf_t_first, h3p_cf_t_last, h3p_cf_nodes, h3p_cf_w_erg')
    L.append('      integer, parameter :: h3p_cf_t_first = %d   ! K, first node' % T_FIRST)
    L.append('      integer, parameter :: h3p_cf_t_last  = %d ! K, last node' % T_LAST)
    L.append('      integer, parameter :: h3p_cf_nodes   = %d' % n)
    L.append('      ! W(T) at T = h3p_cf_t_first + i - 1 K  [erg s^-1 sr^-1 molecule^-1]')
    L.append('      real*8, parameter :: h3p_cf_w_erg(h3p_cf_nodes) = [ &')
    per = 4
    for i in range(0, n, per):
        chunk = ', '.join(values[i:i + per])
        tail = ' ]' if i + per >= n else ', &'
        L.append('           ' + chunk + tail)
    L.append('      end module h3p_mizatep_cooling_table')
    with open(out_path, 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    print(f'wrote {out_path}: {n} nodes, {T_FIRST}-{T_LAST} K')


if __name__ == '__main__':
    main()
