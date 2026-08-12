#!/usr/bin/env python3
"""Build examples/EXHALE_analysis.ipynb from named cells via nbformat.
Run:  python3 build_notebook.py   (then optionally
      jupyter nbconvert --to notebook --execute --inplace EXHALE_analysis.ipynb)
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))
code = lambda s: cells.append(nbf.v4.new_code_cell(s))

md("""# EXHALE: analysis examples

This notebook reads converged EXHALE runs with the `exhale_io` helper and
reproduces the profile and parameter-comparison figures of the user manual
(`docs/EXHALE_user_manual.tex`). Run it from the `examples/` directory.
""")

code("""import os, numpy as np
import matplotlib.pyplot as plt
import exhale_io as aio
ROOT = '..'
W = os.path.join(ROOT, 'WASP-121b')

def load_tut(sub):
    d = os.path.join(ROOT, 'examples', sub)
    adv = os.path.exists(os.path.join(d, 'output', 'Hydro_ioniz_adv.txt'))
    return aio.load_run(os.path.join(d,'output'), os.path.join(d,'input.inp'), adv=adv)
""")

md("## 1. A converged run: the generic hot-Jupiter tutorial\n"
   "`exhale_io.load_run` returns named, cgs-unit arrays (radius in $R_p$).")
code("""t = load_tut('tutorial')
print('Mdot = %.3f Mp/Gyr' % aio.mdot_Mp_per_Gyr(t))
fig, ax = plt.subplots(2, 2, figsize=(9,7))
ax[0,0].semilogy(t.r, t.n);  ax[0,0].set_ylabel(r'$n$ [cm$^{-3}$]')
ax[0,1].plot(t.r, t.v_kms);  ax[0,1].set_ylabel(r'$v$ [km s$^{-1}$]')
ax[1,0].plot(t.r, t.T);      ax[1,0].set_ylabel(r'$T$ [K]')
xh = t.x_ion(['HI','HII'])
ax[1,1].plot(t.r, xh['HI'], label='H I'); ax[1,1].plot(t.r, xh['HII'], label='H II')
ax[1,1].set_ylabel('H ionization fraction'); ax[1,1].legend()
for a in ax.flat: a.set_xlabel(r'$r$ [$R_p$]'); a.grid(alpha=0.3)
plt.tight_layout(); plt.show()
""")

md("## 2. Effect of metal-line cooling (metals on vs off)\n"
   "Same planet with and without `metals.inp`.")
code("""on, off = load_tut('tutorial'), load_tut('tutorial_nometals')
fig, ax = plt.subplots(1, 2, figsize=(10,4))
ax[0].plot(off.r, off.T, '--', label='metals off'); ax[0].plot(on.r, on.T, label='metals on')
ax[0].set_ylabel('T [K]'); ax[0].legend()
xo, xn = off.x_ion(['HI','HII']), on.x_ion(['HI','HII'])
ax[1].plot(off.r, xo['HII'], '--', label='H II (off)'); ax[1].plot(on.r, xn['HII'], label='H II (on)')
ax[1].set_ylabel('H II fraction'); ax[1].legend()
for a in ax: a.set_xlabel(r'$r$ [$R_p$]'); a.grid(alpha=0.3)
plt.tight_layout(); plt.show()
""")

md("## 3. Cooling breakdown by channel (WASP-121b, metals on)\n"
   "From `Cooling_breakdown.txt`: which coolant dominates at each altitude.")
code("""c = aio.load_cooling(os.path.join(W,'output','Cooling_breakdown.txt'))
fig, ax = plt.subplots(figsize=(8,5))
ax.semilogy(c['r'], c['cool_total'], 'k', lw=2, label='total')
for nm in ['coex_HI','rec','FeII','MgII','CII','OII','CaII']:
    if nm in c['chan']:
        ax.semilogy(c['r'], np.clip(c['chan'][nm],1e-30,None), label=nm)
ax.set_ylim(c['cool_total'].max()*1e-4, c['cool_total'].max()*2)
ax.set_xlabel(r'$r$ [$R_p$]'); ax.set_ylabel(r'cooling [erg cm$^{-3}$ s$^{-1}$]')
ax.legend(ncol=2, fontsize=8); ax.grid(alpha=0.3); plt.show()
""")

md("## 4. Spherical vs Roche (RLOF) and the mass-loss rate\n"
   "Removing `Domain mode: Spherical` switches on the Roche tidal potential.")
code("""A = aio.load_run(os.path.join(W,'output_caseA2058'), os.path.join(W,'input.inp.caseA2058'))
B = aio.load_run(os.path.join(W,'output'), os.path.join(W,'input.inp'))
print('Mdot: spherical A = %.3f, Roche B = %.3f Mp/Gyr'
      % (aio.mdot_Mp_per_Gyr(A), aio.mdot_Mp_per_Gyr(B)))
fig, ax = plt.subplots(1, 3, figsize=(13,4))
for r_,lab in [(A,'spherical (A)'),(B,'Roche (B)')]:
    ax[0].plot(r_.r, r_.T, label=lab); ax[1].plot(r_.r, r_.v_kms); ax[2].semilogy(r_.r, r_.n)
ax[0].set_ylabel('T [K]'); ax[0].legend()
ax[1].set_ylabel(r'$v$ [km s$^{-1}$]'); ax[2].set_ylabel(r'$n$ [cm$^{-3}$]')
for a in ax: a.set_xlabel(r'$r$ [$R_p$]'); a.set_xlim(1,2); a.grid(alpha=0.3)
plt.tight_layout(); plt.show()
""")

md("## 5. Transmission: metal transit radii vs Huang et al. (2023)\n"
   "Convert TPM line depths $h$ to $R_{\\rm eff}/R_\\star=\\sqrt{(R_p/R_\\star)^2+h}$. "
   "Run `EXHALE_transit.py` (spherical/triaxial) for the depths; values below are Case A.")
code("""Rp, Rstar = 2.058*aio.RJ, 1.458*6.957e10
td = (Rp/Rstar)**2; reff = lambda h: np.sqrt(td + h/100.0)
lines = ['Mg II (4A)','Ca II K','Na D2']
model  = [reff(2.820), reff(7.343), reff(1.046)]
huangA = [0.182, 0.199, 0.152]; huangD = [0.302, 0.278, 0.147]
x = np.arange(len(lines)); w = 0.27
fig, ax = plt.subplots(figsize=(7,4.2))
ax.bar(x-w, model, w, label='this work (A)'); ax.bar(x, huangA, w, label='Huang A')
ax.bar(x+w, huangD, w, label='Huang D')
ax.set_xticks(x); ax.set_xticklabels(lines); ax.set_ylabel(r'$R_{\\rm eff}/R_\\star$')
ax.legend(); ax.grid(alpha=0.3, axis='y'); plt.show()
""")

nb['cells'] = cells
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3',
                             'language': 'python'}
with open('EXHALE_analysis.ipynb', 'w') as f:
    nbf.write(nb, f)
print('wrote EXHALE_analysis.ipynb (%d cells)' % len(cells))
