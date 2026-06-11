#!/usr/bin/env python3
"""Compare WASP-121b cooling: coronal Fe II (baseline) vs density-dependent
Fe II multilevel-SE (new). Reads the two Cooling_breakdown.txt dumps and the
Hydro_ioniz.txt temperature profiles, and reports the Fe II suppression vs
radius plus the change in the dominant coolant across 1.15-1.4 Rp (Huang Fig 10).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 1-indexed columns in Cooling_breakdown.txt
R, T, NE, TOT = 0, 1, 2, 3
MgII, CaII, NaI, FeI, FeII = 18, 24, 26, 32, 33   # 0-indexed metal columns

base = np.loadtxt('output_coronalFeII/Cooling_breakdown.txt')
new  = np.loadtxt('output/Cooling_breakdown.txt')

# align on radius (both runs share the same grid by construction)
rb, rn = base[:, R], new[:, R]
assert base.shape == new.shape and np.allclose(rb, rn), 'grid mismatch'
r = rn

feII_b, feII_n = base[:, FeII], new[:, FeII]
tot_b,  tot_n  = base[:, TOT],  new[:, TOT]
mgII_n = new[:, MgII]
Tb, Tn = base[:, T], new[:, T]

# suppression factor (avoid /0)
supp = np.where(feII_n > 0, feII_b / np.maximum(feII_n, 1e-300), np.nan)

def at(rad):
    i = np.argmin(np.abs(r - rad))
    return i

print('# r/Rp   T_base   T_new    FeII_base   FeII_new   suppress   '
      'FeII/tot_b  FeII/tot_n  MgII/tot_n')
for rad in [1.00, 1.02, 1.05, 1.10, 1.15, 1.20, 1.30, 1.40]:
    i = at(rad)
    fb = feII_b[i] / tot_b[i] if tot_b[i] > 0 else 0
    fn = feII_n[i] / tot_n[i] if tot_n[i] > 0 else 0
    mn = mgII_n[i] / tot_n[i] if tot_n[i] > 0 else 0
    print(f'{r[i]:6.3f} {Tb[i]:8.1f} {Tn[i]:8.1f}  {feII_b[i]:.3e}  '
          f'{feII_n[i]:.3e}  {supp[i]:8.1f}  {fb:9.4f}  {fn:9.4f}  {mn:9.4f}')

# dominant coolant in the 1.15-1.4 Rp band (new run)
band = (r >= 1.15) & (r <= 1.40)
labels = ['CI','CII','CIII','OI','OII','OIII','NI','NII','NIII','MgI','MgII',
          'MgIII','SiI','SiII','SiIII','CaI','CaII','CaIII','NaI','NaII','KI',
          'KII','SI','SII','FeI','FeII','FeIII']
print('\n# Band-integrated (1.15-1.4 Rp) metal-line cooling share, NEW run:')
metal_new = new[:, 8:8+27]
band_sum = metal_new[band].sum(axis=0)
order = np.argsort(band_sum)[::-1]
for k in order[:6]:
    print(f'  {labels[k]:5s} {band_sum[k]/band_sum.sum()*100:6.2f}%')

# base totals
print(f'\n# Total cooling at base r~1.0:  baseline {tot_b[at(1.0)]:.3e}  '
      f'new {tot_n[at(1.0)]:.3e}')

# ---- figure ----
fig, ax = plt.subplots(2, 2, figsize=(11, 8))

ax[0,0].semilogy(r, feII_b, 'r-',  label=r'Fe\,II coronal (baseline)')
ax[0,0].semilogy(r, feII_n, 'b-',  label=r'Fe\,II density-dep.\ (SE)')
ax[0,0].semilogy(r, mgII_n, 'g--', label=r'Mg\,II (new run)')
ax[0,0].set_xlabel(r'$r/R_p$'); ax[0,0].set_ylabel(r'cooling [erg cm$^{-3}$ s$^{-1}$]')
ax[0,0].set_xlim(1.0, 1.5); ax[0,0].legend(fontsize=8); ax[0,0].set_title('Fe II line cooling vs radius')

ax[0,1].semilogy(r, supp, 'k-')
ax[0,1].set_xlabel(r'$r/R_p$'); ax[0,1].set_ylabel(r'$\Lambda_{\rm coronal}/\Lambda_{\rm SE}$')
ax[0,1].set_xlim(1.0, 1.5); ax[0,1].set_title('Fe II suppression factor (coronal / SE)')
ax[0,1].axhline(1.0, color='grey', lw=0.6, ls=':')

ax[1,0].plot(r, Tb, 'r-', label='baseline (coronal)')
ax[1,0].plot(r, Tn, 'b-', label='new (density-dep.)')
ax[1,0].set_xlabel(r'$r/R_p$'); ax[1,0].set_ylabel(r'$T$ [K]')
ax[1,0].set_xlim(1.0, 1.5); ax[1,0].legend(fontsize=8); ax[1,0].set_title('Temperature profile')

ax[1,1].plot(r, feII_b/np.maximum(tot_b,1e-300), 'r-', label=r'Fe\,II/total (baseline)')
ax[1,1].plot(r, feII_n/np.maximum(tot_n,1e-300), 'b-', label=r'Fe\,II/total (new)')
ax[1,1].plot(r, mgII_n/np.maximum(tot_n,1e-300), 'g--', label=r'Mg\,II/total (new)')
ax[1,1].set_xlabel(r'$r/R_p$'); ax[1,1].set_ylabel('fraction of total cooling')
ax[1,1].set_xlim(1.0, 1.5); ax[1,1].set_ylim(0, 1.0); ax[1,1].legend(fontsize=8)
ax[1,1].set_title('Fractional cooling budget')

fig.suptitle('WASP-121b: density-dependent Fe II line cooling (multilevel SE)')
fig.tight_layout()
fig.savefig('FeII_density_correction_compare.png', dpi=130)
print('\nwrote FeII_density_correction_compare.png')
