#!/usr/bin/env python3
"""Downward radiative power of the wind through the matching face (Q1, step 3).

Reads the converged k05 state (read only) and writes, for every radiative
channel of output/Cooling_breakdown.txt:
  P_tot  = sum_i Lambda_i 4 pi r_i^2 dr_i / (4 pi R0^2)            [erg cm^-2 s^-1]
  P_down = sum_i Lambda_i 4 pi r_i^2 dr_i W(r_i) / (4 pi R0^2)
where W(r) = (1 - sqrt(1 - (R0/r)^2))/2 is the fraction of the isotropic
emission of a thin shell at r whose rays meet the sphere of radius R0
(solid angle of the sphere seen from r, over 4 pi).  The emission is taken
optically thin in the wind (the H3+ self-absorption is measured separately).

Usage: python3 wind_downward_emission.py <k05 run directory> <out prefix>
"""
import sys
import os
import numpy as np

run = sys.argv[1]
out = sys.argv[2]
R_J = 7.1492e9
res = {}
for line in open(os.path.join(run, 'EXHALE_resolved.out')):
    if line.startswith('#') or not line.strip():
        continue
    t = line.split()
    if len(t) >= 2:
        res[t[0]] = t[1]
Rp = float(res['planet_radius_RJ']) * R_J
R0 = float(res['base_face_radius_cm'])

cb = np.loadtxt(os.path.join(run, 'output', 'Cooling_breakdown.txt'))[2:-2]
hdr = [l for l in open(os.path.join(run, 'output', 'Cooling_breakdown.txt'))
       if l.startswith('# col1')][0]
names = ['r', 'T', 'ne', 'cool_total', 'reco', 'coio', 'coex_HI', 'coex_HeI',
         'coex_HeII', 'brem', 'H3p_IR', 'H2_IR', 'H2O_IR', 'CO_IR',
         'CI', 'CII', 'CIII', 'OI', 'OII', 'OIII', 'NI', 'NII', 'NIII', 'MgI',
         'MgII', 'MgIII', 'SiI', 'SiII', 'SiIII', 'CaI', 'CaII', 'CaIII', 'NaI',
         'NaII', 'KI', 'KII', 'SI', 'SII', 'FeI', 'FeII', 'FeIII']
assert cb.shape[1] == len(names), (cb.shape, len(names))
r = cb[:, 0] * Rp
# faces: the base face R0, then midpoints of the centres, last face mirrored
f = np.empty(r.size + 1)
f[0] = R0
f[1:-1] = 0.5 * (r[1:] + r[:-1])
f[-1] = 2 * r[-1] - f[-2]
vol = 4 * np.pi / 3 * (f[1:]**3 - f[:-1]**3)
area = 4 * np.pi * R0**2
W = 0.5 * (1 - np.sqrt(np.clip(1 - (R0 / r)**2, 0, None)))

hy = np.loadtxt(os.path.join(run, 'output', 'Hydro_ioniz.txt'))[2:-2]
ion_hdr = [l for l in open(os.path.join(run, 'output', 'Ion_species.txt'))
           if l.startswith('# columns')][0].split()[2:]
ion = np.loadtxt(os.path.join(run, 'output', 'Ion_species.txt'))[2:-2]
nH3p = ion[:, ion_hdr.index('H3p')]
nH2 = ion[:, ion_hdr.index('H2')]
print('R0 = %.6e cm, Rp = %.6e cm, R0/Rp = %.6f, first centre %.6f Rp'
      % (R0, Rp, R0 / Rp, r[0] / Rp))
print('check: cool column of Hydro_ioniz vs channel sum max rel diff %.2e'
      % np.max(np.abs(hy[:, 6] - cb[:, 3]) / np.maximum(np.abs(hy[:, 6]), 1e-300)))

rows = []
for k, nm in enumerate(names[3:], start=3):
    L = cb[:, k]
    Ptot = np.sum(L * vol) / area
    Pdn = np.sum(L * vol * W) / area
    m12 = r < 1.2 * Rp
    P12 = np.sum((L * vol)[m12]) / area
    if Ptot != 0.0 or Pdn != 0.0:
        rows.append((nm, Ptot, Pdn, P12))
with open(out + '_channels.txt', 'w') as fo:
    fo.write('# channel  P_tot  P_down(thin)  P_tot(r<1.2Rp)  [erg cm^-2 s^-1 of the face R0]\n')
    for row in rows:
        fo.write('%-10s % .5e % .5e % .5e\n' % row)
        print('%-10s P_tot % .4e  P_down % .4e  P(<1.2Rp) % .4e' % row)

# cumulative H3+ emission and the radii of its quartiles
L = cb[:, names.index('H3p_IR')] * vol
c = np.cumsum(L) / np.sum(L)
for q in (0.1, 0.25, 0.5, 0.75, 0.9):
    print('H3+ emission quantile %.2f at r = %.4f Rp' % (q, np.interp(q, c, r / Rp)))
cd = np.cumsum(L * W) / np.sum(L * W)
for q in (0.25, 0.5, 0.75):
    print('H3+ DOWNWARD emission quantile %.2f at r = %.4f Rp' % (q, np.interp(q, cd, r / Rp)))
# emission-weighted temperature of the downward H3+ power
print('downward-power-weighted T(H3+) = %.1f K; T at 1.0,1.04,1.08,1.12,1.2 Rp:'
      % (np.sum(L * W * cb[:, 1]) / np.sum(L * W)),
      np.interp([1.0, 1.04, 1.08, 1.12, 1.2], r / Rp, cb[:, 1]))
# H3+ and H2 radial columns from the face outward, and from each shell down
dr = f[1:] - f[:-1]
N3 = np.sum(nH3p * dr)
N2 = np.sum(nH2 * dr)
print('radial column above R0: N(H3+) = %.3e cm^-2, N(H2) = %.3e cm^-2' % (N3, N2))
N3_below = np.cumsum(nH3p * dr) - 0.5 * nH3p * dr
print('emission-weighted H3+ column between the emitter and the face (radial): %.3e cm^-2'
      % (np.sum(L * W * N3_below) / np.sum(L * W)))
np.savetxt(out + '_profile.txt',
           np.column_stack([r / Rp, cb[:, 1], hy[:, 3], nH2, nH3p,
                            cb[:, names.index('H3p_IR')], vol, W]),
           header='r[Rp] T[K] p[cgs] n_H2 n_H3p[cm^-3] Lambda_H3p[erg cm^-3 s^-1] '
                  'shell_volume[cm^3] W_down; face R0 = %.9e cm' % R0)
