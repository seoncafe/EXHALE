import numpy as np, sys, hashlib
fn = sys.argv[1]
raw = open(fn,'rb').read()
off = 0
hdr = np.frombuffer(raw, dtype=np.int32, count=6, offset=off); off += 24
nY, nfs1, nfs2, N, Ng, nvar = [int(x) for x in hdr]
print("header: nY=%d f_sp=(%d,%d) N=%d Ng=%d nvar_jac=%d" % (nY,nfs1,nfs2,N,Ng,nvar))
def rd(n):
    global off
    a = np.frombuffer(raw, dtype=np.float64, count=n, offset=off); off += 8*n
    return a
dtau = rd(1)[0]; print("dtau =", repr(dtau))
Y = rd(nY); fsp = rd(nfs1*nfs2); F = rd(nY); Fjac = rd(nY); D = rd(nY); Drow = rd(nY)
print("bytes consumed %d of %d" % (off, len(raw)))
for nm,a in (("Y",Y),("F",F),("Fjac",Fjac),("D",D),("Drow",Drow)):
    print("%-5s md5=%s min=%.6E max=%.6E norm2=%.6E" % (nm, hashlib.md5(a.tobytes()).hexdigest(), a.min(), a.max(), np.linalg.norm(a)))
print("F == Fjac bitwise:", np.array_equal(F.view(np.int64), Fjac.view(np.int64)))
print("D == Drow bitwise:", np.array_equal(D.view(np.int64), Drow.view(np.int64)))
for k,nm in enumerate(("mass","momentum","energy")):
    r = Drow[k::nvar]; f = F[k::nvar]; s = np.abs(f/r)
    j = int(np.argmax(s))
    print("row %-8s Drow: min %.6E at cell %d, max %.6E at cell %d ; |F/Drow| max %.6E at cell %d" %
          (nm, r.min(), int(np.argmin(r))+1, r.max(), int(np.argmax(r))+1, s[j], j+1))
    print("      cells 1,2,3 : Drow %.6E %.6E %.6E ; F %.6E %.6E %.6E" % (r[0],r[1],r[2],f[0],f[1],f[2]))
