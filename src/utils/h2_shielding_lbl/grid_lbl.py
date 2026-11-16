import sys, os, time
sys.path.insert(0,'.')
import numpy as np, lbl
T = float(sys.argv[1]); wlmax = float(sys.argv[2])
N = np.concatenate(([0.0], 10.0**np.arange(17.0, 22.01, 0.25)))
t0=time.time()
kov, kno, ll = lbl.shielding_curve(T, N, wl_max=wlmax, keep_frac=1.0-1e-7)
np.savez('lbl_T%04d_w%04d.npz' % (int(T), int(wlmax)),
         N=N[1:], f_ov=(kov/kov[0])[1:], f_no=(kno/kno[0])[1:],
         k_unatt=kov[0], nlines=len(ll['nu0']))
print('T=%g wlmax=%g  %.0f s' % (T, wlmax, time.time()-t0))
for i,n in enumerate(N[1:]):
    print('%10.3e %11.4e %11.4e' % (n, (kov/kov[0])[i+1], (kno/kno[0])[i+1]))
