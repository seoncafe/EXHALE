#!/usr/bin/env python3
"""Reproduce the paper's own p-winds model at its published best fit, using
the authors' released script settings (LHS1140b_zenodo/lhs1140b_2024B_pwinds.py):
  - const-res v23 SED x Lx-ratio 0.59 (their exact normalization)
  - R_pl = 0.154 R_J, a = 0.09382 AU, b = 0.23, planet/star = 0.072092
  - turbulence_broadening = True  (Lampon+2020 term -- the extra width we had
    earlier attributed to an undocumented ~13 km/s kernel)
  - 5-phase transit averaging over [-0.5, 0.5]
  - v_wind = 2.26 km/s, NOT applied inside the model here (see below)

Parameter vector: the MCMC medians (Mdot 2.03e8 g/s, T 5160 K, H:He 1.01e-3),
which is the vector the paper's Fig. 4 legend states for the plotted curve.

Result at the MCMC medians:
  turbulence ON  -> red 1.376 %, FWHM 0.575 A   (paper Fig 4: 1.323 %, ~0.71)
  turbulence OFF -> red 1.447 %, FWHM 0.458 A
So the authors' own two switches, not a mystery kernel, produce the plotted
line shape; the reproduction matches the paper to ~3 % in depth.

The residual width (ours ~0.58 A vs the paper's ~0.71 A) is NOT explained by
the choice of parameter vector: the grid argmax (1.83e8 g/s, 5131.6 K,
8.338e-4) and these medians give FWHM values differing by < 0.01 A. The cause
of the width residual remains unidentified.

Bulk line-of-sight velocity convention: v_wind is applied at plot time, not
inside the model, so that the files written here follow the same convention as
every other model file in this study (tspec_published_*, the EXHALE tpm_*
outputs). This is exact for wind_broadening_method='average', where
bulk_los_velocity acts as a rigid shift: running with v_wind and with 0 gives
red/blue depths and FWHM identical to the third decimal, and only the fitted
line shift changes (+0.080 A vs -0.002 A).
"""
import numpy as np, astropy.constants as c, astropy.units as u, sys
from astropy.io import fits
from astropy.convolution import convolve
from p_winds import parker, hydrogen, helium, transit, lines
sys.path.insert(0,'/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00')
from he_line_metrics import fit_metrics

au2m=1.496e11; pc2m=3.0857e16; Rjup=7.1492e7
d_1140=14.96*pc2m; d_1132=12.613*pc2m; smax_1140=0.0946
M_sun=1.98847e30; M_star=0.18*M_sun
R_pl=0.154; M_pl=0.0176; a_pl=0.09382
planet_to_star_ratio=0.072092; impact_parameter=0.23
m_h=c.m_p.to(u.g).value; m_He=4*1.67262192369e-27
r=np.logspace(0,np.log10(10),100)
w0,w1,w2,f0,f1,f2,a_ij=lines.he_3_properties()
w_array=np.array([w0,w1,w2]); f_array=np.array([f0,f1,f2]); a_array=np.array([a_ij]*3)
sample_phases=np.linspace(-0.5,0.5,5); n_samples=5

Lx_ratio_1132=0.59
norm=Lx_ratio_1132*(d_1140/(smax_1140*au2m))**2
ZEN='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/LHS1140b_zenodo/'
spec=fits.getdata(ZEN+'hlsp_muscles_multi_multi_gj1132_broadband_v23_const-res-sed.fits',1)
host={'wavelength':spec['WAVELENGTH'],'flux_lambda':spec['FLUX']*norm,
      'wavelength_unit':u.angstrom,'flux_unit':u.erg/u.s/u.cm**2/u.angstrom}

def atmo(log_mdot,log_T,hfrac):
    m_dot=10**log_mdot; T=10**log_T; he=1-hfrac; hh=he/hfrac
    mu0=(1+4*hh)/(1+hh+0.90)
    f_r,mu_bar=hydrogen.ion_fraction(r,R_pl,T,hfrac,m_dot,M_pl,mu0,
        spectrum_at_planet=host,initial_f_ion=0.0,relax_solution=True,
        exact_phi=True,return_mu=True)
    vs=parker.sound_speed(T,mu_bar); rs=parker.radius_sonic_point(M_pl,vs)
    rhos=parker.density_sonic_point(m_dot,rs,vs)
    v_arr,rho_arr=parker.structure(r*R_pl/rs)
    f1_,f3_=helium.population_fraction(r,v_arr,rho_arr,f_r,R_pl,T,hfrac,vs,rs,rhos,
        spectrum_at_planet=host,initial_state=np.array([1.,0.]),relax_solution=True)
    n_he=rho_arr*rhos*he/(hfrac+4*he)/m_h
    return f3_*n_he, v_arr*vs

def tmodel(wl, v_wind, n3, log_T, v_arr, turb, broad='average'):
    R=R_pl*Rjup; rSI=r*R; vSI=v_arr*1000; n3SI=n3*1e6
    specs=[]
    for i in range(n_samples):
        fm,td,rm=transit.draw_transit(planet_to_star_ratio,impact_parameter=impact_parameter,
            supersampling=5,phase=sample_phases[i],planet_physical_radius=R,grid_size=100)
        s=transit.radiative_transfer_2d(fm,rm,rSI,n3SI,vSI,w_array,f_array,a_array,
            wl,10**log_T,m_He,bulk_los_velocity=v_wind,wind_broadening_method=broad,
            turbulence_broadening=turb)
        specs.append(s+td)
    return np.mean(specs,axis=0)

# instrumental profile on observed grid; use a fine grid instead
wl=np.linspace(1.0827e-6,1.0839e-6,600)
# v_wind = 2.26 km/s is applied at plot time (see the module docstring), so the
# model itself is evaluated with zero bulk line-of-sight velocity.
v_wind=0.0
# MCMC medians -- the vector the paper's Fig. 4 legend states.
# (grid argmax, for reference: Mdot 1.83e8 g/s, T 5131.6 K, H:He 8.338e-4)
best=(np.log10(2.03e8), np.log10(5160.0), 1.01e-3)
n3,v=atmo(*best)
AIR=10832.057/10829.09114
for turb,tag in ((True,'turbulence ON (authors)'),(False,'turbulence OFF (our oracle)')):
    ts=tmodel(wl,v_wind,n3,best[1],v,turb)
    exc=(1-ts/ts.max())*100
    # instrument convolution R=68000
    lam0=np.mean(wl)*1e10; fwhm=lam0/68000; sig=fwhm/2.3548
    dl=np.median(np.diff(wl))*1e10; x=np.arange(-5*sig,5*sig+dl,dl)
    k=np.exp(-0.5*(x/sig)**2); k/=k.sum()
    excc=np.convolve(exc,k,mode='same')
    lam_vac=wl*1e10*AIR
    m=fit_metrics(lam_vac,excc,frame='vacuum')
    print(f'{tag:28s} red {m["red_depth"]:.3f}%  blue {m["blue_depth"]:.3f}%  r/b {m["red_blue"]:.2f}  FWHM {m["fwhm_A"]:.3f} A')
    out='tspec_authors_turb.txt' if turb else 'tspec_authors_noturb.txt'
    np.savetxt(out, np.column_stack([lam_vac, ts, exc, excc]),
               header=('LHS 1140b at the MCMC medians quoted by the paper Fig. 4 '
                       'legend\n(Mdot 2.03e8 g/s, T 5160 K, H:He 1.01e-3), the '
                       'authors own settings:\nconst-res v23 x A=0.59, '
                       'turbulence_broadening=%s, 5-phase average\n'
                       'no bulk line-of-sight velocity applied here: '
                       'v_wind = 2.26 km/s is applied at plot time\n'
                       'vacuum wavelength [A]  transmission  excess[%%]  '
                       'excess convolved R=68000 [%%]') % turb)
    print('    -> '+out)
print('paper Fig4 purple (digitized):  red 1.323%  ~  FWHM 0.708 A')
print('observed:                       red 1.254%  blue 0.198%  r/b 6.35  FWHM 0.841 A')
