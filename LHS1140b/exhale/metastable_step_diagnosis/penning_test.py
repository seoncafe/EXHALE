"""Counterfactual test: is the He 2^3S step exactly the Taylor+2025 Penning
   branch switch at T = 4000 K?  Rebuild n3 with Q31 forced onto the branch the
   neighbouring cell used, and compare with the smooth extrapolation."""
import numpy as np, os, sys
base='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/exhale/'
kb_eV=8.617333262e-5; A31=1.272e-4
def Q31_lo(T): return 1.9e-9*(3.0e2/T)**0.07
def Q31_hi(T): return 9.1e-9*(3.0e2/T)**0.50
def Q31(T):    return np.where(T<=4.0e3, Q31_lo(T), Q31_hi(T))
def q31a(T):
    u=2.847*np.exp(-1.252e-5*T)-1.953*np.exp(-3.558e-4*T)
    return 2.10e-8*np.sqrt(13.60/(kb_eV*T))*np.exp(-0.80/(kb_eV*T))*u/3.0
def q31b(T):
    u=1.185*np.exp(-4.749e-6*T)-0.9131*np.exp(-1.669e-4*T)
    return 2.10e-8*np.sqrt(13.60/(kb_eV*T))*np.exp(-1.40/(kb_eV*T))*u/3.0
def run(name, which='eq'):
    D=os.path.join(base,name,'output'); suf='' if which=='eq' else '_adv'
    h=open(os.path.join(D,'Ion_species.txt')).readlines()[1].split()[2:]
    itr=h.index('HeITR'); ihi=h.index('HI')
    sp=np.loadtxt(os.path.join(D,'Ion_species%s.txt'%suf))
    hy=np.loadtxt(os.path.join(D,'Hydro_ioniz%s.txt'%suf))
    hh=open(os.path.join(D,'Hydro_ioniz.txt')).readlines()[1].split()[2:]
    ne=np.loadtxt(os.path.join(D,'Heating_breakdown.txt'))[:,2]
    r=sp[:,0]; T=hy[:,hh.index('T[K]')]; nHI=sp[:,ihi]; n3=sp[:,itr]
    print('== %s (%s)  Tmax=%.1f K'%(name,which,T.max()))
    out=[]
    for i in range(1,len(r)):
        if r[i]<1.05: continue
        if (T[i-1]-4000.0)*(T[i]-4000.0)>0: continue     # no 4000 K crossing
        rat=n3[i]/n3[i-1]
        up = T[i]>4000.0
        # counterfactual: put cell i on the SAME Penning branch as cell i-1
        Lk  = Q31(T[i])*nHI[i]  + A31 + (q31a(T[i])+q31b(T[i]))*ne[i]
        Qcf = Q31_lo(T[i]) if up else Q31_hi(T[i])
        Lcf = Qcf*nHI[i]        + A31 + (q31a(T[i])+q31b(T[i]))*ne[i]
        n3cf= n3[i]*Lk/Lcf
        # smooth extrapolation of the pre-step trend (log-linear on last 2 cells)
        ext = n3[i-1]*(n3[i-1]/n3[i-2])
        print('  cross at i=%d r=%.4f  T %.2f -> %.2f  n3 %.4f -> %.4f (rat %.4f)'%(i,r[i],T[i-1],T[i],n3[i-1],n3[i],rat))
        print('     Q31 branch jump factor           = %.4f'%(Q31_hi(T[i])/Q31_lo(T[i])))
        print('     Penning share of known loss      = %.4f'%(Q31(T[i])*nHI[i]/Lk))
        print('     n3 rebuilt on neighbour branch   = %.4f'%n3cf)
        print('     n3 smooth extrapolation          = %.4f   -> residual %.2f%%'%(ext,100*(n3cf/ext-1)))
if __name__=='__main__':
    for nm in ['flux_closure/heh10p3/k01','heh0p55','heh0p55_diff_kzz1e10','heh2p13_diff_kzz1e9','flux_closure/heh11p1/k01',
               'heh0p55_diff_kzz1e9','heh0p55_diff_kzz1e8','heh0p55_diff_ctrl','flux_closure/hi/k06']:
        for w in ('eq','adv'):
            try: run(nm,w)
            except Exception as e: print(nm,w,'ERR',e)
