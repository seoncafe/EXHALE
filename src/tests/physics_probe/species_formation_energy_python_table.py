#!/usr/bin/env python3
"""The Python formation-energy table against the Fortran one.

`src/utils/formation_energy_flux_diagnostic.py` carries its own `EPS_EV`
dictionary, because it reads output files and never links the code.  Two
tables of the same physical quantity are two chances to disagree, so this
check compares them entry by entry against
`species_formation_energy_table.dat`, which the Fortran driver
`species_formation_energy_table.f90` writes from the production table.

Run the Fortran driver first: this script reads the file it leaves in
$EXHALE_TEST_WORK, or in the current directory when that is unset.

One line per assertion:  PASS|FAIL <name> measured= reference= tol=
Exit status is nonzero if any assertion failed.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src', 'utils'))

# The names of the .dat file, and the key each one carries in the
# diagnostic's EPS_EV.  The metal stages are written as metal<column> and
# are compared through the diagnostic's own ion labels.
BASE_NAMES = {
    'HI': 'HI', 'HII': 'HII', 'HeI': 'HeI', 'HeII': 'HeII',
    'HeIII': 'HeIII', 'HeITR': 'HeITR', 'H2': 'H2', 'H2p': 'H2p',
    'H3p': 'H3p', 'HeHp': 'HeHp',
}
# f_sp column of each metal stage, in the canonical order of
# species_table.f90 (columns 7 to 33).
METAL_COLUMN = {
    7: 'CI', 8: 'CII', 9: 'CIII', 10: 'OI', 11: 'OII', 12: 'OIII',
    13: 'NI', 14: 'NII', 15: 'NIII', 16: 'MgI', 17: 'MgII', 18: 'MgIII',
    19: 'SiI', 20: 'SiII', 21: 'SiIII', 22: 'CaI', 23: 'CaII',
    24: 'CaIII', 25: 'NaI', 26: 'NaII', 27: 'KI', 28: 'KII', 29: 'SI',
    30: 'SII', 31: 'FeI', 32: 'FeII', 33: 'FeIII',
}
# Entries the Fortran table carries that the diagnostic does not.  Reported
# rather than asserted: adding them is a change to a file this test does not
# own, and the values and their sources are printed so it can be made.
NOT_IN_DIAGNOSTIC = ('OH', 'H2O', 'CO', 'O1D')

# Absolute tolerance [eV].  The two tables are transcriptions of the same
# published constants, so what is allowed is the decimal truncation of the
# Python literals, not a physical difference.
TOL_EV = 1.0e-6

n_fail = 0


def verdict(name, ok, measured, reference, tol):
    global n_fail
    tag = 'PASS'
    if not ok:
        tag = 'FAIL'
        n_fail += 1
    print('%s %s measured=%.15e reference=%.15e tol=%s'
          % (tag, name, measured, reference, tol))


def read_table():
    work = os.environ.get('EXHALE_TEST_WORK', '')
    candidates = []
    if work:
        candidates.append(os.path.join(work, 'species_formation_energy_table.dat'))
    candidates.append('species_formation_energy_table.dat')
    candidates.append(os.path.join(ROOT, 'build', 'tests', 'physics_probe',
                                   'species_formation_energy_table.dat'))
    for path in candidates:
        if os.path.exists(path):
            table = {}
            with open(path) as fh:
                for line in fh:
                    if line.startswith('#') or not line.strip():
                        continue
                    key, value = line.split()
                    table[key] = float(value)
            return path, table
    return None, {}


def main():
    global n_fail
    path, fortran = read_table()
    if not fortran:
        print('FAIL species_formation_energy_dat measured=0 reference=1 '
              'tol=0  (run species_formation_energy_table.x first)')
        return 1
    print('# Fortran table read from %s' % path)

    from formation_energy_flux_diagnostic import EPS_EV as python_table

    for dat_name, py_name in sorted(BASE_NAMES.items()):
        if dat_name not in fortran:
            verdict('python_table_' + dat_name, False, 0.0, 0.0, 'present')
            continue
        if py_name not in python_table:
            verdict('python_table_' + dat_name, False, 0.0, 0.0, 'present')
            continue
        f = fortran[dat_name]
        p = python_table[py_name]
        verdict('python_table_' + dat_name, abs(f - p) <= TOL_EV, p, f,
                '%.1e absolute' % TOL_EV)

    for column, py_name in sorted(METAL_COLUMN.items()):
        dat_name = 'metal%02d' % column
        if dat_name not in fortran or py_name not in python_table:
            verdict('python_table_' + py_name, False, 0.0, 0.0, 'present')
            continue
        f = fortran[dat_name]
        p = python_table[py_name]
        verdict('python_table_' + py_name, abs(f - p) <= TOL_EV, p, f,
                '%.1e absolute' % TOL_EV)

    missing = [n for n in NOT_IN_DIAGNOSTIC if n in fortran]
    if missing:
        print('# Fortran entries the diagnostic table does not carry, with')
        print('# the values it would need (0 K formation energies from the')
        print('# NIST-JANAF Shomate table of oxygen_rates.f90, and the O(1D)')
        print('# excitation energy of the same module):')
        for name in missing:
            print('#   %-5s %14.6f eV' % (name, fortran[name]))

    if n_fail:
        print('species_formation_energy_python_table: %d assertion(s) failed'
              % n_fail)
        return 1
    print('species_formation_energy_python_table: every assertion passed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
