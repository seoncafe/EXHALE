import numpy as np
pi=np.pi
# --- code grid (set_energy_vectors, PL branch) ---
def grid(e_sub_low, NlTR):
    e=[]
    if NlTR>0:
        e += [e_sub_low*(13.6/e_sub_low)**((j-1.0)/NlTR) for j in range(1,NlTR+1)]
    for (emin,emax,n) in [(13.6,24.6,50),(24.6,54.4,50),(54.4,123.98,50),(123.98,1240.0,50)]:
        e += [emin*(emax/emin)**((j-1.0)/n) for j in range(1,n+1)]
    e=np.array(e); de=np.empty_like(e)
    de[0]=0.5*(e[1]-e[0]); de[1:-1]=0.5*(e[2:]-e[:-2]); de[-1]=0.5*(e[-1]-e[-2])
    return e,de
def sigma_H(E,Z=1.0):
    E0=13.6*Z*Z; E=np.atleast_1d(E).astype(float); s=np.full_like(E,6.3/(Z*Z))
    m=E>E0; eps=np.sqrt(E[m]/E0-1.0)
    s[m]=6.3/(Z*Z)*(E0/E[m])**4*np.exp(4-4*np.arctan(eps)/eps)/(1-np.exp(-2*pi/eps))
    s[E<0.99999*E0]=0.0; return s
def Jinc(E, PLind=-1.0, e_low=13.6, e_mid=123.98, e_top=1240.0, J_EUV=1.0, J_X=0.0):
    P1=PLind+1
    if abs(P1)<1e-12: Jn=J_EUV/np.log(e_mid/e_low); JXn=J_X/np.log(e_top/e_mid)
    else: Jn=J_EUV*P1/(e_mid**P1-e_low**P1); JXn=J_X*P1/(e_top**P1-e_mid**P1)
    E=np.atleast_1d(E); return np.where(E<e_mid, Jn*E**PLind, JXn*E**PLind)
# exact reference
Ef=np.logspace(np.log10(13.6),np.log10(1240.0),400001)
for (lab,esub,ntr) in [("floor 4.80 (He23S)",4.80,20),("floor 5.139 (Na I)",5.139,20),("floor 13.6",13.6,0)]:
    e,de=grid(esub,ntr)
    for Z,name in [(1.0,'HI'),(2.0,'HeII')]:
        code=np.sum(Jinc(e)*sigma_H(e,Z)/e*de)
        exact=np.trapz(Jinc(Ef)*sigma_H(Ef,Z)/Ef,Ef)
        print(f"{lab:22s} P_{name:5s} code/exact = {code/exact:.4f}")
    # heating H I
    code=np.sum(Jinc(e)*sigma_H(e)*(1-13.6/np.maximum(e,13.6))*de); exact=np.trapz(Jinc(Ef)*sigma_H(Ef)*(1-13.6/Ef),Ef)
    print(f"{lab:22s} Heat_HI code/exact = {code/exact:.4f}")
    # bin at 13.6
    k=np.argmin(abs(e-13.6)); print(f"   e_v at threshold={e[k]:.3f}, de_v={de[k]:.4f}, spacing below={e[k]-e[k-1] if k>0 else float('nan'):.4f}, above={e[k+1]-e[k]:.4f}")
# --- B2: sub-Lyman field vs photosphere (wasp_full) ---
h=4.135667696e-15; c=2.99792458e10; kB=8.617333262e-5; hp_erg=6.62607015e-27
Teff=6459.0; Rs=1.458*6.957e10; a=0.02544*1.495978707e13
J_EUV=10**30.42/(4*pi*a*a)
def planck_per_eV(E):  # pi B_nu (R/a)^2 / h   [erg cm^-2 s^-1 eV^-1]
    nu=E/h; Bnu=2*hp_erg*nu**3/c**2/(np.exp(E/(kB*Teff))-1.0)
    return pi*Bnu*(Rs/a)**2/h
print("\nwasp_full: J_EUV(planet) = %.3e erg/cm2/s"%J_EUV)
for E in [4.8,6.0,8.0,10.0,13.0]:
    pl=planck_per_eV(E); code=Jinc(E,J_EUV=J_EUV)[0]
    print(f"E={E:5.1f} eV  Planck {pl:.3e}  code PL {code:.3e}  ratio {pl/code:.3g}")
# P(He23S) over 4.78-13.6 with code's wing-A cross section
def vfky(E,Eth,E0,s0,ya,P,yw,y0,y1):
    E=np.atleast_1d(E); x=E/E0-y0; z=np.sqrt(x*x+y1*y1); Q=5.5-0.5*P
    s=s0*((x-1)**2+yw**2)*z**(-Q)*(1+np.sqrt(z/ya))**(-P); s[E<Eth]=0; return s
Eb=np.linspace(4.78,13.6,200001); sig=vfky(Eb,4.78,2.645,20.8,1e12,3.42,2.681,1.956,2.603)*1e-18
xi=0.5
P_code=xi*np.trapz(Jinc(Eb,J_EUV=J_EUV)*sig/(Eb*1.602176634e-12),Eb)
P_pl  =xi*np.trapz(planck_per_eV(Eb)*sig/(Eb*1.602176634e-12),Eb)
print(f"P(He23S) 4.78-13.6 eV: code PL {P_code:.3e} s^-1, Planck {P_pl:.3e} s^-1, ratio {P_pl/P_code:.0f}; sigma(4.8eV)={vfky(4.8,4.78,2.645,20.8,1e12,3.42,2.681,1.956,2.603)[0]:.2f} Mb")
# total flux on grid vs nominal
e,de=grid(4.80,20); print("integrated grid flux / nominal J_EUV (LX off):", np.sum(Jinc(e,J_EUV=J_EUV)*de)/J_EUV)
