#!/usr/bin/env python3
"""Compare the WASP-121b _adv metal ionization between the frozen-eq
post-process (pp_metals=1) and the re-solve (pp_metals=2).

Mode 1 freezes the metal ion densities at their converged equilibrium values
when building the advected (_adv) profiles. Mode 2 re-solves the metal
ionization balance at the advection-corrected H/He and the post-process
temperature, so the metal stage split tracks the (more ionized, slightly
different-T) advected wind rather than the equilibrium structure.

This reads:
  output_ppm1_dd/Ion_species_adv.txt   (mode 1, density-dependent FeII)
  output/Ion_species_adv.txt           (mode 2, re-solve)
  output/Ion_species.txt               (the frozen equilibrium split)
  {dir}/Hydro_ioniz_adv.txt            (advected T)
  output/Hydro_ioniz.txt               (equilibrium T)
and reports the singly-ionized fraction of each element vs radius, plus a
figure. n_e here is the H/He electron density written to the file; the metal
contribution to n_e is small but is what the re-solve self-consistently uses.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Ion_species[_adv].txt: col0=r, 1=HI,2=HII,3=HeI,4=HeII,5=HeIII,6=HeITR,
# then metals mion 1..27 at col 6+im (0-indexed).  mion: 1=CI..3=CIII,
# 4..6=O, 7..9=N, 10..12=Mg, 13..15=Si, 16..18=Ca, 19..20=Na, 21..22=K,
# 23..24=S, 25..27=Fe.
def col(im):           # 0-indexed adv/eq column for canonical mion index im
    return 6 + im

# element -> (neutral mion, list of stage mion indices, singly-ionized mion)
ELEM = {
    'C':  (1,  [1, 2, 3],   2),
    'O':  (4,  [4, 5, 6],   5),
    'N':  (7,  [7, 8, 9],   8),
    'Mg': (10, [10, 11, 12], 11),
    'Ca': (16, [16, 17, 18], 17),
    'Na': (19, [19, 20],     20),
    'Fe': (25, [25, 26, 27], 26),
}

# Mode 1 (frozen) and mode 2 (re-solve) _adv profiles are preserved in their
# own dirs; the canonical output/ is kept at the mode-1 default. The equilibrium
# split and equilibrium T are mode-independent (byte-identical), so they are read
# from output/.
m1 = np.loadtxt('output_ppm1_dd/Ion_species_adv.txt')       # frozen eq metals
m2 = np.loadtxt('output_ppm2_resolve/Ion_species_adv.txt')  # re-solved metals
eq = np.loadtxt('output/Ion_species.txt')                   # equilibrium split
hd_a = np.loadtxt('output_ppm2_resolve/Hydro_ioniz_adv.txt')  # advected T (mode 2)
hd_e = np.loadtxt('output/Hydro_ioniz.txt')                 # equilibrium T

r  = m2[:, 0]
assert np.allclose(r, m1[:, 0]) and np.allclose(r, eq[:, 0]), 'grid mismatch'
Ta = hd_a[:, 4]
Te = hd_e[:, 4]
neH = m2[:, 2] + m2[:, 4] + 2.0*m2[:, 5]   # H/He electrons (HII+HeII+2HeIII)

def ion_frac(arr, elem):
    n0, stages, ising = ELEM[elem]
    tot = sum(arr[:, col(s)] for s in stages)
    tot = np.maximum(tot, 1e-300)
    return arr[:, col(ising)] / tot

def at(rad):
    return int(np.argmin(np.abs(r - rad)))

rads = [1.00, 1.02, 1.05, 1.10, 1.15, 1.20, 1.30, 1.40, 1.48]

print('# Singly-ionized fraction X+/X_tot in the _adv wind: '
      'eq(frozen) vs re-solve, per element')
print('# also T_eq, T_adv [K] and the H/He electron density n_eH [cm^-3]')
hdr = '  r/Rp   T_eq    T_adv   ' + ' '.join(f'{e:>14s}' for e in ELEM)
print(hdr)
print('        [K]     [K]     ' +
      ' '.join(f'{"eq->resolve":>14s}' for _ in ELEM))
for rad in rads:
    i = at(rad)
    cells = []
    for e in ELEM:
        f1 = ion_frac(m1, e)[i]
        f2 = ion_frac(m2, e)[i]
        cells.append(f'{f1:5.3f}->{f2:5.3f}')
    print(f'{r[i]:6.3f} {Te[i]:7.0f} {Ta[i]:7.0f}  ' +
          ' '.join(f'{c:>14s}' for c in cells))

# largest relative change in X+ density across the wind, per element
print('\n# max |Delta n(X+)| / n(X+) (re-solve vs frozen), wind 1.0-1.48 Rp')
band = (r >= 1.0) & (r <= 1.48)
for e in ELEM:
    _, _, ising = ELEM[e]
    n1 = m1[:, col(ising)]
    n2 = m2[:, col(ising)]
    rel = np.where(n1 > 0, np.abs(n2 - n1)/np.maximum(n1, 1e-300), 0.0)
    j = np.argmax(np.where(band, rel, 0.0))
    print(f'  {e:3s}  max rel change {rel[j]*100:6.2f}%  at r/Rp={r[j]:.3f} '
          f'(n+ {n1[j]:.3e} -> {n2[j]:.3e})')

# ---- figure: singly-ionized fraction vs radius for the line-cooling metals ----
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
panel = [('Fe', ax[0, 0]), ('Mg', ax[0, 1]), ('Ca', ax[1, 0]), ('Na', ax[1, 1])]
for e, a in panel:
    a.plot(r, ion_frac(m1, e), 'r-', label='frozen eq (mode 1)')
    a.plot(r, ion_frac(m2, e), 'b-', label='re-solve (mode 2)')
    a.set_xlabel(r'$r/R_p$')
    a.set_ylabel(e + r' II / ' + e + r' total')
    a.set_xlim(1.0, 1.48)
    a.set_ylim(0, 1.02)
    a.legend(fontsize=8)
    a.set_title(e + r' singly-ionized fraction (\_adv wind)')
fig.suptitle(r'WASP-121b: metal ionization re-solve in the advected wind '
             r'(pp\_metals 1 vs 2)')
fig.tight_layout()
fig.savefig('metal_resolve_compare.png', dpi=130)
print('\nwrote metal_resolve_compare.png')
