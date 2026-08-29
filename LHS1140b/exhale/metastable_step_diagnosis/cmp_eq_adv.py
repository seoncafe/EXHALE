import numpy as np, sys, os
base='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/exhale/'
run=sys.argv[1] if len(sys.argv)>1 else 'heh0p55_diff_kzz1e10'
D=os.path.join(base,run,'output')
def load(f):
    return np.loadtxt(os.path.join(D,f))
eq=load('Ion_species.txt'); ad=load('Ion_species_adv.txt')
hy=load('Hydro_ioniz.txt'); hya=load('Hydro_ioniz_adv.txt')
hdr=open(os.path.join(D,'Ion_species.txt')).readlines()[1].split()[2:]
itr=hdr.index('HeITR'); ihe1=hdr.index('HeI'); ihe2=hdr.index('HeII'); ihi=hdr.index('HI'); ihii=hdr.index('HII')
r=eq[:,0]
print('%-8s %-9s %12s %12s %8s | %12s %12s %8s'%('i','r','TR_eq','TR_adv','adv/eq','HeI_eq','HeI_adv','adv/eq'))
sel=(r>1.25)&(r<1.62)
for i in np.where(sel)[0]:
    print('%-8d %-9.4f %12.5e %12.5e %8.4f | %12.5e %12.5e %8.4f'%(i,r[i],eq[i,itr],ad[i,itr],ad[i,itr]/eq[i,itr],eq[i,ihe1],ad[i,ihe1],ad[i,ihe1]/eq[i,ihe1]))
