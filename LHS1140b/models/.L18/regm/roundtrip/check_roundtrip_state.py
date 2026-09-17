#!/usr/bin/env python3
"""Restart round-trip identity: write -> read -> write must be the identity.

Compares a reference Ion_species/Hydro_ioniz pair with the pair the reloaded
run dumped (EXHALE_DUMP_IC=1, which writes the loaded state before the first
equilibrium sweep). Three statements:

  * every species column of Ion_species matches to round-off, matched BY
    LABEL from the '# columns' header rather than by position, so a schema
    change is reported as a missing column and not as a silent mismatch;
  * rho, v and p of Hydro_ioniz match to round-off (T is derived from the
    composition and is allowed the particle count's own round-off);
  * the '# coupling:' line -- what the run was produced under -- survives.

This replaces check_roundtrip.py, which compared by column POSITION (so a
schema change reads as a numeric mismatch instead of as a missing column) and
whose case had no readable IC to run against. Driven by roundtrip_check.sh in
this directory.
Usage: check_roundtrip_state.py ref_ion dump_ion ref_hyd dump_hyd
"""
import sys

import numpy as np

RTOL = 1.0e-12


def read_columns(path):
    labels, coupling = None, None
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                break
            body = line.lstrip('#').strip()
            if body.lower().startswith('columns'):
                labels = body[len('columns'):].split()
            elif body.lower().startswith('coupling:'):
                coupling = body
    if labels is None:
        raise SystemExit('%s: no "# columns" header' % path)
    data = np.loadtxt(path, comments='#')
    return {lab: data[:, i] for i, lab in enumerate(labels)}, labels, coupling


def rel(a, b):
    scale = np.maximum(np.abs(b), np.abs(b).max() * 1.0e-30 + 1.0e-300)
    return np.abs(a - b) / scale


def compare(ref_f, dump_f, only=None):
    ref, ref_lab, ref_cpl = read_columns(ref_f)
    dmp, dmp_lab, dmp_cpl = read_columns(dump_f)
    bad = 0
    missing = [c for c in ref_lab if c not in dmp]
    if missing:
        print('  FAIL %s: columns absent from the dump: %s'
              % (dump_f, ' '.join(missing)))
        bad += 1
    for c in ref_lab:
        if c in missing or c == 'r[Rp]':
            continue
        if only is not None and c not in only:
            continue
        d = rel(dmp[c], ref[c]).max()
        if d > RTOL:
            print('  FAIL %-14s max rel diff %.3e' % (c, d))
            bad += 1
    return bad, ref_cpl, dmp_cpl


def main():
    ref_ion, dump_ion, ref_hyd, dump_hyd = sys.argv[1:5]
    bad = 0
    n, _, _ = compare(ref_ion, dump_ion)
    bad += n
    # The '# coupling:' line is written on Hydro_ioniz.txt only.
    n, cpl_ref, cpl_dump = compare(ref_hyd, dump_hyd,
                                   only={'rho[mH/cm3]', 'v[cm/s]', 'p[cgs]'})
    bad += n
    if cpl_ref is None:
        print('  FAIL reference carries no "# coupling:" header')
        bad += 1
    elif cpl_dump != cpl_ref:
        print('  FAIL coupling header not preserved:\n    ref  %s\n    dump %s'
              % (cpl_ref, cpl_dump))
        bad += 1
    else:
        print('  ok   coupling header preserved: %s' % cpl_ref)
    if bad:
        print('\n%d round-trip check(s) FAILED' % bad)
        return 1
    print('\nrestart round trip is the identity to %.0e' % RTOL)
    return 0


if __name__ == '__main__':
    sys.exit(main())
