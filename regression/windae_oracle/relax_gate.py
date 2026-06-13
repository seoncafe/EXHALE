#!/usr/bin/env python3
"""P1 relaxation gate: diff the Fortran relax_ae_f output against the C
windsoln_ref first M=1501 rows. Run relax_ae_f first (see plan)."""
import sys, numpy as np
def load(fn, nmax=None):
    rows=[[float(x) for x in l.split(',')] for l in open(fn)
          if not l.startswith('#') and l.strip()]
    a=np.array(rows); return a[:nmax] if nmax else a
def gate(fort_csv, ref_csv, m=1501, rtol=1e-6):
    F=load(fort_csv); W=load(ref_csv, m)
    names=['r','rho','v','T','Ys_HI','Ys_HeI','Ncol_HI','Ncol_HeI','q','z']
    worst=0.0; wv=''
    for i in range(10):
        rel=np.abs(F[:,i]-W[:,i])/np.maximum(np.abs(W[:,i]),1e-300)
        if rel.max()>worst: worst=rel.max(); wv=names[i]
    ok = worst < rtol
    print(f'{"PASS" if ok else "FAIL"}: max rel={worst:.3e} ({wv}); rtol={rtol}')
    return ok
if __name__=='__main__':
    sys.exit(0 if gate(sys.argv[1], sys.argv[2],
                       int(sys.argv[3]) if len(sys.argv)>3 else 1501) else 1)
