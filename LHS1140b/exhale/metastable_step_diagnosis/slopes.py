#!/usr/bin/env python3
"""Local logarithmic slope d ln EW / d ln(He/H) of the stored scan runs,
   for converting an equivalent-width change into a crossing shift."""
import numpy as np, os, glob
B = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/exhale'
A = 10832.057/10829.09114; LO, HI = 10832.60, 10834.20
def ew(d):
    p = os.path.join(d, 'tpm_He10830.txt')
    if not os.path.exists(p): return None
    a = np.loadtxt(p); lam = a[:,0]*A
    dep = (a[:,2].max()-a[:,2])/a[:,2].max()*100.0
    m = (lam >= LO) & (lam <= HI)
    return np.trapz(dep[m], lam[m])
grp = {'well mixed, diffusion off':
         [('heh0p25',0.25),('heh0p5',0.5),('heh0p55',0.55),('heh0p6',0.6),
          ('heh0p7',0.7),('heh1',1.0)],
       'K_zz = 1e9':
         [('heh2_diff_kzz1e9',2.0),('heh2p13_diff_kzz1e9',2.13),
          ('heh3_diff_kzz1e9',3.0),('heh4_diff_kzz1e9',4.0),
          ('heh5_diff_kzz1e9',5.0),('heh8_diff_kzz1e9',8.0),
          ('heh9p7_diff_kzz1e9',9.7)]}
for k, v in grp.items():
    print('##', k); pts = []
    for tag, x in v:
        e = ew(os.path.join(B, tag))
        if e is not None: pts.append((x, e)); print('  %-24s He/H=%-6g EW=%.5f'%(tag,x,e))
    for i in range(1, len(pts)):
        s = np.log(pts[i][1]/pts[i-1][1])/np.log(pts[i][0]/pts[i-1][0])
        print('    slope %g-%g : %.4f'%(pts[i-1][0], pts[i][0], s))
print('## flux closure arms')
for d in sorted(glob.glob(B+'/flux_closure/*/k0*')):
    e = ew(d)
    if e is not None: print('  %-42s EW=%.5f'%(d.replace(B+'/',''), e))
