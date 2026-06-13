#!/usr/bin/env python3
"""P0 oracle runner / checker for the Wind-AE Fortran port.

For each fixture case (inputs/ + defs.h + windsoln_ref.csv):
  1. install the case's defs.h and inputs into the Wind-AE C tree,
  2. rebuild and run ./bin/relaxed_ae,
  3. compare the produced windsoln.csv against windsoln_ref.csv.

Usage:
  run_oracle.py [--windae <wind_ae pkg dir>] [case ...]
  run_oracle.py --diff A.csv B.csv [rtol]   # plain CSV compare only

Comparison: header lines must match textually (modulo float formatting);
data columns are compared with per-column rtol (default 0 = byte-equal
after parsing, i.e. exact float equality).
"""
import sys, os, shutil, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WINDAE = '/nfs/mocafe/kiseon/RT_Codes/Exoplanetary_Atmospheres/wind-ae-main/wind_ae'
CASES = ['mc09_fiducial', 'hd209458b', 'hd189733b']

def read_windsoln(fn):
    header, rows = [], []
    for line in open(fn):
        if line.startswith('#'):
            header.append(line.rstrip('\n'))
        elif line.strip():
            rows.append([float(x) for x in line.split(',')])
    return header, np.array(rows)

def compare(fa, fb, rtol=0.0, atol=0.0):
    ha, da = read_windsoln(fa)
    hb, db = read_windsoln(fb)
    ok = True
    if len(ha) != len(hb):
        print(f'  header line count differs: {len(ha)} vs {len(hb)}'); ok = False
    if da.shape != db.shape:
        print(f'  data shape differs: {da.shape} vs {db.shape}'); return False
    if rtol == 0.0 and atol == 0.0:
        same = np.array_equal(da, db)
        if not same:
            bad = ~np.isclose(da, db, rtol=1e-14, atol=0)
            print(f'  exact mismatch in {bad.sum()} entries; '
                  f'max rel dev = {np.nanmax(np.abs(da-db)/np.maximum(np.abs(db),1e-300)):.3e}')
            ok = False
    else:
        bad = ~np.isclose(da, db, rtol=rtol, atol=atol)
        if bad.any():
            rel = np.abs(da-db)/np.maximum(np.abs(db), 1e-300)
            print(f'  {bad.sum()} entries exceed rtol={rtol}; max rel dev = {np.nanmax(rel):.3e}')
            ok = False
    return ok

def run_case(case, rtol):
    d = os.path.join(HERE, case)
    print(f'=== {case}')
    shutil.copy(os.path.join(d, 'defs.h'), os.path.join(WINDAE, 'src/defs.h'))
    for f in os.listdir(os.path.join(d, 'inputs')):
        shutil.copy(os.path.join(d, 'inputs', f), os.path.join(WINDAE, 'inputs', f))
    subprocess.run(['make', 'compile'], cwd=os.path.join(WINDAE, 'src'),
                   capture_output=True, check=True)
    r = subprocess.run(['./bin/relaxed_ae'], cwd=WINDAE, capture_output=True, text=True)
    if r.returncode != 0:
        print('  relaxed_ae FAILED:', r.stdout[-500:], r.stderr[-200:]); return False
    ok = compare(os.path.join(WINDAE, 'saves/windsoln.csv'),
                 os.path.join(d, 'windsoln_ref.csv'), rtol=rtol)
    print('  PASS' if ok else '  FAIL')
    return ok

if __name__ == '__main__':
    args = sys.argv[1:]
    if args and args[0] == '--diff':
        rtol = float(args[3]) if len(args) > 3 else 0.0
        sys.exit(0 if compare(args[1], args[2], rtol=rtol) else 1)
    rtol = 0.0
    if '--rtol' in args:
        i = args.index('--rtol'); rtol = float(args[i+1]); del args[i:i+2]
    cases = args or CASES
    allok = all(run_case(c, rtol) for c in cases)
    print('ALL PASS' if allok else 'FAILURES PRESENT')
    sys.exit(0 if allok else 1)
