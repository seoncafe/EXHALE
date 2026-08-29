import numpy as np, os, sys
SC=sys.argv[1]
EX='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/exhale'
def prof(D,f):
    h=open(os.path.join(D,'Ion_species.txt')).readlines()[1].split()[2:]
    a=np.loadtxt(os.path.join(D,f)); return a[:,0],a[:,h.index('HeITR')]
for d in ['heh0p55_diff_ctrl','heh0p55_diff_kzz1e10']:
    for f in ['Ion_species.txt','Ion_species_adv.txt']:
        r0,n0=prof(os.path.join(EX,d,'output'),f)
        r1,n1=prof(os.path.join(SC,'rep/LHS1140b/exhale',d+'_t1','output'),f)
        rel=np.abs(n1/np.maximum(n0,1e-300)-1)
        print('%-22s %-20s  max|dn3/n3| = %.3e at r=%.4f ; median %.2e'%(d,f,rel.max(),r0[rel.argmax()],np.median(rel)))
        # steps in rerun
        st=[(r1[i],n1[i]/n1[i-1]) for i in range(1,len(n1)) if n1[i-1]>0 and n1[i]>0 and r1[i]>1.05 and not (0.9<n1[i]/n1[i-1]<1.12)]
        st0=[(r0[i],n0[i]/n0[i-1]) for i in range(1,len(n0)) if n0[i-1]>0 and n0[i]>0 and r0[i]>1.05 and not (0.9<n0[i]/n0[i-1]<1.12)]
        print('    stored steps r>1.05 (first 6): ',' '.join('%.4f(%.3f)'%x for x in st0[:6]),' n=%d'%len(st0))
        print('    rerun  steps r>1.05 (first 6): ',' '.join('%.4f(%.3f)'%x for x in st[:6]),' n=%d'%len(st))
