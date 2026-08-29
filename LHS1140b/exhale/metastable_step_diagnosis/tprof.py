import numpy as np, os, sys
base='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/exhale/'
run=sys.argv[1]
D=os.path.join(base,run,'output')
hy=np.loadtxt(os.path.join(D,'Hydro_ioniz.txt'))
hh=open(os.path.join(D,'Hydro_ioniz.txt')).readlines()[1].split()[2:]
iT=hh.index('T[K]')
r=hy[:,0]; T=hy[:,iT]
h=open(os.path.join(D,'Ion_species.txt')).readlines()[1].split()[2:]; itr=h.index('HeITR')
eq=np.loadtxt(os.path.join(D,'Ion_species.txt')); ad=np.loadtxt(os.path.join(D,'Ion_species_adv.txt'))
# find eq steps
v=eq[:,itr]
print('run',run,'  Tmax=%.3f at r=%.4f'%(T.max(),r[T.argmax()]))
for i in range(1,len(v)):
    if v[i-1]>0 and v[i]>0 and (v[i]/v[i-1]<0.9 or v[i]/v[i-1]>1.12) and r[i]>1.05:
        print('  EQ  step i=%d r=%.4f rat=%.4f  T=%.2f -> %.2f'%(i,r[i],v[i]/v[i-1],T[i-1],T[i]))
w=ad[:,itr]
for i in range(1,len(w)):
    if w[i-1]>0 and w[i]>0 and (w[i]/w[i-1]<0.9 or w[i]/w[i-1]>1.12) and r[i]>1.05:
        print('  ADV step i=%d r=%.4f rat=%.4f  T=%.2f -> %.2f'%(i,r[i],w[i]/w[i-1],T[i-1],T[i]))
