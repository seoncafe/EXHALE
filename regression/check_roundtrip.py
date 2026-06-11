"""Compare an Ion_species IC file with the ATES_DUMP_IC re-dump.

Usage: check_roundtrip.py <reference.txt> <dump.txt> {schema2|legacy}

schema2: every species column (H/He + all metal ions) of the dump must match
         the reference to round-off (the loader divides by rho*n0 and the
         writer multiplies back, so agreement is ~1e-14 relative).
legacy:  H/He(+HeITR) columns must match the reference; metal ion columns
         must instead equal the neutral-from-abundance fallback (neutral
         stage > 0 everywhere, higher stages exactly 0).
"""
import sys
import numpy as np

ref_f, dump_f, mode = sys.argv[1], sys.argv[2], sys.argv[3]
ref  = np.loadtxt(ref_f)
dump = np.loadtxt(dump_f)

NHHE = 7           # r + HI HII HeI HeII HeIII HeITR
RTOL = 1e-12

def rel(a, b):
    scale = np.maximum(np.abs(b), np.abs(b).max() * 1e-30 + 1e-300)
    return np.abs(a - b) / scale

fail = 0

# H/He(+TR) columns must always round-trip
for c in range(1, NHHE):
    bad = rel(dump[:, c], ref[:, c]) > RTOL
    if bad.any():
        print(f'  FAIL H/He col {c}: max rel diff {rel(dump[:,c],ref[:,c]).max():.2e}')
        fail = 1

if mode == 'schema2':
    if dump.shape[1] != ref.shape[1]:
        print(f'  FAIL: column count {dump.shape[1]} != {ref.shape[1]}'); fail = 1
    for c in range(NHHE, min(dump.shape[1], ref.shape[1])):
        bad = rel(dump[:, c], ref[:, c]) > RTOL
        if bad.any():
            print(f'  FAIL metal col {c}: max rel diff {rel(dump[:,c],ref[:,c]).max():.2e}')
            fail = 1
    if not fail:
        print(f'  PASS: all {ref.shape[1]-1} species columns restored (rtol {RTOL})')
elif mode == 'legacy':
    # Metal columns: neutral stages of elements with non-zero abundance must
    # be positive; non-neutral stages must be exactly zero. Canonical layout:
    # cols 7.. : CI CII CIII OI OII OIII NI NII NIII MgI MgII MgIII SiI SiII
    #            SiIII CaI CaII CaIII NaI NaII KI KII SI SII FeI FeII FeIII
    neutral_cols = [7, 10, 13, 16, 19, 22, 25, 27, 29, 31]
    ion_cols     = [c for c in range(7, 34) if c not in neutral_cols]
    for c in ion_cols:
        if c < dump.shape[1] and np.any(dump[:, c] != 0.0):
            print(f'  FAIL legacy fallback: ionized metal col {c} non-zero'); fail = 1
    # At least one neutral metal column should be positive (abundances set)
    pos = [c for c in neutral_cols if c < dump.shape[1] and np.all(dump[:, c] > 0)]
    if not pos:
        print('  FAIL legacy fallback: no neutral metal columns populated'); fail = 1
    if not fail:
        print(f'  PASS: H/He restored; metals at neutral fallback ({len(pos)} elements)')
else:
    print('unknown mode'); fail = 1

sys.exit(fail)
