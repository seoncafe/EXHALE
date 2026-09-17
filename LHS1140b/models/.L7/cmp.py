import sys, numpy as np
def load(fn):
    cols=None; rows=[]
    for line in open(fn):
        if line.startswith('#'):
            t=line[1:].split()
            if t and t[0]=='columns': cols=t[1:]
            continue
        rows.append([float(x) for x in line.split()])
    return cols, np.array(rows)
ca,A=load(sys.argv[1]); cb,B=load(sys.argv[2])
print(f"{sys.argv[1].split('/')[-1]}  rows {A.shape} -> {B.shape}")
for i,c in enumerate(ca):
    j=cb.index(c) if c in cb else None
    if j is None: print(f"  {c}: absent in second"); continue
    a=A[:,i]; b=B[:,j]
    den=np.maximum(np.abs(a),np.abs(b))
    d=np.where(den>0, np.abs(a-b)/np.where(den>0,den,1), 0.0)
    k=int(np.argmax(d))
    print(f"  {c:8s} max|rel| {d.max():.3e} at row {k}   (a={a[k]:.6e} b={b[k]:.6e})")
extra=[c for c in cb if c not in ca]
for c in extra:
    j=cb.index(c); b=B[:,j]
    print(f"  {c:8s} NEW  min {b.min():.6e}  max {b.max():.6e}")
