import numpy as np, os, sys
base='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/exhale/'
kb_eV=8.617333262e-5
def Q31(T):
    return np.where(T<=4.0e3, 1.9e-9*(3.0e2/T)**0.07, 9.1e-9*(3.0e2/T)**0.50)
def q31a(T):
    ups=2.847*np.exp(-1.252e-5*T)-1.953*np.exp(-3.558e-4*T)
    return 2.10e-8*np.sqrt(13.60/(kb_eV*T))*np.exp(-0.80/(kb_eV*T))*ups/3.0
def q31b(T):
    ups=1.185*np.exp(-4.749e-6*T)-0.9131*np.exp(-1.669e-4*T)
    return 2.10e-8*np.sqrt(13.60/(kb_eV*T))*np.exp(-1.40/(kb_eV*T))*ups/3.0
A31=1.272e-4
run=sys.argv[1]; lo=int(sys.argv[2]); hi=int(sys.argv[3])
adv = len(sys.argv)>4 and sys.argv[4]=='adv'
D=os.path.join(base,run,'output')
suf='_adv' if adv else ''
h=open(os.path.join(D,'Ion_species.txt')).readlines()[1].split()[2:]
itr=h.index('HeITR'); ihi=h.index('HI'); ihii=h.index('HII'); ihe2=h.index('HeII'); ihe3=h.index('HeIII')
sp=np.loadtxt(os.path.join(D,'Ion_species%s.txt'%suf))
hy=np.loadtxt(os.path.join(D,'Hydro_ioniz%s.txt'%suf))
hh=open(os.path.join(D,'Hydro_ioniz.txt')).readlines()[1].split()[2:]; iT=hh.index('T[K]')
hb=np.loadtxt(os.path.join(D,'Heating_breakdown.txt'))
r=sp[:,0]; T=hy[:,iT]; nHI=sp[:,ihi]; n3=sp[:,itr]
ne=hb[:,2]; hePen=hb[:,13]; heTR=hb[:,7]
print('%-4s %-8s %8s %10s %10s %10s %10s %8s %8s %8s'%('i','r','T','n3','Q31*nHI','A31','q3x*ne','Pen_frac','n3rat','pred'))
for i in range(lo,hi+1):
    Lp=Q31(T[i])*nHI[i]; Lq=(q31a(T[i])+q31b(T[i]))*ne[i]
    Ltot_known=Lp+A31+Lq
    f=Lp/Ltot_known
    print('%-4d %-8.4f %8.2f %10.4e %10.4e %10.4e %10.4e %8.4f %8.4f %8.4f'%(
        i,r[i],T[i],n3[i],Lp,A31,Lq,f, n3[i]/n3[i-1] if i>0 else 0, 0))
# jump prediction at the crossing
