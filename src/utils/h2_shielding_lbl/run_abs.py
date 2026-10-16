import sys, os, time; sys.path.insert(0,'.')
import numpy as np, lbl_table as LT
T = float(sys.argv[1])
NCOL = 39
logN = np.linspace(12.0, np.log10(5.0e21), NCOL)
N = np.concatenate(([0.0], 10.0**logN))
t0=time.time()
sp, ps = LT.absolute_curves(T, N)
np.savez('abs_T%04d.npz'%int(T), logN=logN, sigma_pump=sp[1:], p_single=ps[1:],
         sigma_pump0=sp[0], p_single0=ps[0])
print('T=%g  %.0fs  sigma_pump(0)=%.4e  p_single(0)=%.4f  sigma_pump(5e21)=%.4e'
      % (T, time.time()-t0, sp[0], ps[0], sp[-1]))
