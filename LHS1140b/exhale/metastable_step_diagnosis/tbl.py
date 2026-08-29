import numpy as np, os, sys
base='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/exhale/'
run=sys.argv[1]; r0=float(sys.argv[2]); r1=float(sys.argv[3])
D=os.path.join(base,run,'output')
h=open(os.path.join(D,'Ion_species.txt')).readlines()[1].split()[2:]
itr=h.index('HeITR'); i1=h.index('HeI'); i2=h.index('HeII'); ih=h.index('HI'); ihii=h.index('HII')
eq=np.loadtxt(os.path.join(D,'Ion_species.txt')); ad=np.loadtxt(os.path.join(D,'Ion_species_adv.txt'))
hy=np.loadtxt(os.path.join(D,'Hydro_ioniz.txt'))
hh=open(os.path.join(D,'Hydro_ioniz.txt')).readlines()[1].split()[2:]
r=eq[:,0]; sel=(r>=r0)&(r<=r1)
print('cols hydro:',hh)
iT=hh.index('T[K]') if 'T[K]' in hh else 1
iv=[k for k,c in enumerate(hh) if c.startswith('v')]
print('%-5s %-8s %11s %8s %11s %8s %11s %11s %11s'%('i','r','TR_eq','rat','TR_adv','rat','HeI_eq','HeII_eq','v'))
idx=np.where(sel)[0]
for i in idx:
    re=eq[i,itr]/eq[i-1,itr] if i>0 else 0
    ra=ad[i,itr]/ad[i-1,itr] if i>0 else 0
    print('%-5d %-8.4f %11.5e %8.4f %11.5e %8.4f %11.4e %11.4e %11.4e'%(i,r[i],eq[i,itr],re,ad[i,itr],ra,eq[i,i1],eq[i,i2],hy[i,iv[0]] if iv else 0))
