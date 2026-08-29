import numpy as np, os, sys
base='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/exhale/'
run=sys.argv[1]; lo=int(sys.argv[2]); hi=int(sys.argv[3])
D=os.path.join(base,run,'output')
h=open(os.path.join(D,'Ion_species.txt')).readlines()[1].split()[2:]
idx={c:k for k,c in enumerate(h)}
eq=np.loadtxt(os.path.join(D,'Ion_species.txt')); ad=np.loadtxt(os.path.join(D,'Ion_species_adv.txt'))
hy=np.loadtxt(os.path.join(D,'Hydro_ioniz.txt')); hya=np.loadtxt(os.path.join(D,'Hydro_ioniz_adv.txt'))
hh=open(os.path.join(D,'Hydro_ioniz.txt')).readlines()[1].split()[2:]
iT=hh.index('T[K]'); iv=hh.index('v[cm/s]'); irho=hh.index('rho[mH/cm3]')
r=eq[:,0]
print('%-4s %-8s %9s %9s %10s %10s %10s %10s %10s %10s %10s'%('i','r','T_eq','T_adv','n3_eq','n3_adv','HeI_adv','HeII_adv','HI_adv','HII_adv','v'))
for i in range(lo,hi+1):
    print('%-4d %-8.4f %9.2f %9.2f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.3e'%(
      i,r[i],hy[i,iT],hya[i,iT],eq[i,idx['HeITR']],ad[i,idx['HeITR']],
      ad[i,idx['HeI']],ad[i,idx['HeII']],ad[i,idx['HI']],ad[i,idx['HII']],hy[i,iv]))
