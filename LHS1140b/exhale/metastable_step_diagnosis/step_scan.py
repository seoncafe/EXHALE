import numpy as np, os, sys
base='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/exhale/'
runs=['heh0p55_diff_ctrl','heh0p55_diff_kzz1e8','heh0p55_diff_kzz1e9','heh0p55_diff_kzz1e10',
      'heh2p13_diff_kzz1e9','flux_closure/heh11p1/k01','flux_closure/hi/k06']
def cols(p):
    return open(p).readlines()[1].split()[2:]
for run in runs:
    D=os.path.join(base,run,'output')
    if not os.path.isdir(D): print(run,'MISSING'); continue
    h=cols(os.path.join(D,'Ion_species.txt')); itr=h.index('HeITR')
    eq=np.loadtxt(os.path.join(D,'Ion_species.txt'))
    ad=np.loadtxt(os.path.join(D,'Ion_species_adv.txt'))
    r=eq[:,0]
    def steps(a,label):
        v=a[:,itr]; out=[]
        for i in range(1,len(v)):
            if v[i-1]>0 and v[i]>0:
                rat=v[i]/v[i-1]
                if rat<0.9 or rat>1.12: out.append((r[i],rat,v[i-1],v[i]))
        print('  %-4s steps=%d'%(label,len(out)), ' '.join('r=%.4f(%.3f)'%(x[0],x[1]) for x in out[:20]))
    print(run)
    steps(eq,'eq'); steps(ad,'adv')
