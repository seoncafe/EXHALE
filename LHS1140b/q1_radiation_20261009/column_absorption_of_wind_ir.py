#!/usr/bin/env python3
"""Where the wind's downward H3+ infrared emission is absorbed in the column.

Line by line, on a grid uniform in ln(nu) (step dv = 0.1 km/s, 4 points per
Doppler parameter of CH4 at 185 K), 1000-10000 cm^-1:

EMISSION.  The downward H3+ power of every wind cell (k05_wind_profile.txt,
written by wind_downward_emission.py: Lambda_H3+ of the code times the shell
volume times W(r) = (1 - sqrt(1-(R0/r)^2))/2, per unit face area) is spread
over the lines of the ExoMol MiZATeP list with the LTE line shape at the cell
temperature (option --shape lte), or put entirely into the nu2 fundamental
(upper (v1,v2) = (0,1), lower (0,0); option --shape fundamental: the limit of
a vibrational population far below LTE, Miller et al. 2013 Table 6 gives s =
0.005-0.3 at n_H2 = 1e8-1e10 cm^-3 where most of the downward power forms).
Each line is a Gaussian of the Doppler width at its power-weighted
temperature (the wind speed there is below 1 m/s, no shift).  The total is
the code's H3+ power, so only the spectral SHAPE comes from the line list.

ANGLES.  The direction cosine at the face of every ray from a cell that meets
the sphere R0 (exact, thin wind: cos alpha uniform in [sqrt(1-(R0/r)^2), 1],
mu_face = sqrt(1 - (r sin alpha / R0)^2)), weighted by the cell power, in 40
bins; the column is plane parallel (its thickness is below 1e-2 R_p) and the
slant optical depth is tau/mu.

COLUMN.  The k05 lower-atmosphere profile (p, T, composition) on its own
levels, hydrostatic columns dN_X = x_X dp/(m g), m = rho/n_tot.  HITRAN
(hitran.org line-by-line, main isotopologues, downloaded 2026-10-09) lines
of H2O, CH4, NH3, C2H2, HCN, CO, H2 (quadrupole) with S(T) from the HITRAN
partition sums Q(T) (hitran.org/data/Q), Voigt profile, Lorentz width
gamma_air (296/T)^n_air p (the air widths stand for the H2/He background;
H2/He widths of CH4 and NH3 are 0.8-1.3 times the air ones, which moves only
the levels below ~0.05 bar where the Lorentz width exceeds the Doppler
width), far-wing cut at 25 cm^-1.  H2-H2 and H2-He CIA from the climate
model's own tables (photochem_clima_data, Molliere et al. 2019).

OUTPUT: the fraction of the downward power absorbed between consecutive
levels, and the transmitted remainder below the deepest level.

Usage: python3 column_absorption_of_wind_ir.py <k05 dir> <exomol dir>
                <hitran dir> <cia dir> <out prefix> [--shape lte|fundamental]
"""
import os
import sys
import numpy as np
from scipy.special import voigt_profile
from scipy.signal import oaconvolve
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3p_linelist import H3pLines, C2  # noqa: E402

KB = 1.380649e-16
AMU = 1.66053907e-24
CL = 2.99792458e10
G = 6.674e-8
R_J = 7.1492e9

run, exo, hit, cia_dir, out = sys.argv[1:6]
shape = 'lte'
if '--shape' in sys.argv:
    shape = sys.argv[sys.argv.index('--shape') + 1]

# ------------------------------------------------------------------ wind
prof = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'k05_wind_profile.txt'))
hdr = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'k05_wind_profile.txt')).readline()
R0 = float(hdr.split('R0 =')[1].split()[0])
rRp, Tw, _, nH2w, nH3w, Lam, vol, W = prof.T
Rp = R0      # the face is at r = 1.000000 R_p (measured)
r = rRp * Rp
Pcell = Lam * vol * W / (4 * np.pi * R0**2)          # erg cm^-2 s^-1
use = Pcell > 1e-12 * Pcell.sum()
r, Tw, Pcell, nH3w = r[use], Tw[use], Pcell[use], nH3w[use]
P_down = Pcell.sum()
print('downward H3+ power (code Lambda): %.4f erg cm^-2 s^-1' % P_down)

# arrival direction cosines at the face
mu_edges = np.linspace(0.0, 1.0, 41)
w_mu = np.zeros(40)
for ri, Pi in zip(r, Pcell):
    x = np.sqrt(max(1 - (R0 / ri)**2, 0.0))
    ca = x + (1 - x) * (np.arange(400) + 0.5) / 400
    mu = np.sqrt(np.clip(1 - (ri / R0)**2 * (1 - ca**2), 0, 1))
    w_mu += Pi * np.histogram(mu, bins=mu_edges)[0] / 400
w_mu /= w_mu.sum()
mu_c = 0.5 * (mu_edges[1:] + mu_edges[:-1])
print('mean mu at the face %.3f; mean 1/mu %.3f'
      % (np.sum(w_mu * mu_c), np.sum(w_mu / mu_c)))

# line powers
L = H3pLines(exo)
if shape == 'fundamental':
    sv = np.loadtxt(os.path.join(exo, 'states_v.txt'), dtype=int)
    v1 = np.full(sv[:, 0].max() + 1, -9)
    v2 = np.full_like(v1, -9)
    v1[sv[:, 0]] = sv[:, 1]
    v2[sv[:, 0]] = sv[:, 2]
    keep_f = ((v1[L.iu] == 0) & (v2[L.iu] == 1)
              & (v1[L.il] == 0) & (v2[L.il] == 0))
    print('fundamental lines: %d' % keep_f.sum())
Tg_edges = np.array([0, 600, 900, 1200, 1400, 1600, 1800, 2000, 2200, 2400,
                     2600, 2800, 3000, 3300, 3600, 4000, 1e9])
Pj = np.zeros(L.nu.size)
PjT = np.zeros(L.nu.size)
ratio_num = 0.0
for a, b in zip(Tg_edges[:-1], Tg_edges[1:]):
    m = (Tw >= a) & (Tw < b)
    if not m.any():
        continue
    Pg = Pcell[m].sum()
    Tg = np.sum(Pcell[m] * Tw[m]) / Pg
    lp = L.line_power(Tg)
    if shape == 'fundamental':
        lp = np.where(keep_f, lp, 0.0)
    Pj += Pg * lp / lp.sum()
    PjT += Pg * lp / lp.sum() * Tg
Teff = PjT / np.maximum(Pj, 1e-300)
o = np.argsort(Pj)[::-1]
c = np.cumsum(Pj[o]) / Pj.sum()
nkeep = np.searchsorted(c, 0.9995) + 1
sel = o[:nkeep]
print('H3+ lines carrying 99.95%% of the power: %d' % nkeep)
nuL, PL, TL = L.nu[sel], Pj[sel] / Pj[sel].sum() * P_down, Teff[sel]

# wind self-absorption: line-centre tau of each line through the H3+ column
# between the emitting shell and the face, populations LTE at the absorber T
dr = np.gradient(r)
N3 = np.cumsum(nH3w * dr) - 0.5 * nH3w * dr
N_eff = np.sum(Pcell * N3) / P_down
T_abs = np.sum(Pcell * np.cumsum(nH3w * dr * Tw) / np.maximum(np.cumsum(nH3w * dr), 1e-300)) / P_down
sig0 = L.tau_center_per_column(T_abs)[sel]
inv_mu = np.sum(w_mu / mu_c)
tau0 = sig0 * N_eff * inv_mu
# Doppler-profile transmission of a line emitted behind tau0 (upper bound:
# all of the column between emitter and face, slant 1/mu mean)
xg = np.linspace(-4, 4, 81)
phi = np.exp(-xg**2)
phi /= phi.sum()
trans_w = np.array([np.sum(phi * np.exp(-t * np.exp(-xg**2))) for t in tau0[:20000]])
lost = np.sum(PL[:20000] * (1 - trans_w)) / P_down
print('wind self-absorption: N(H3+) emitter->face (power weighted, radial) %.3e cm^-2, '
      'absorber T %.0f K, max line-centre tau %.3e, absorbed share (upper bound) %.3e'
      % (N_eff, T_abs, tau0.max(), lost))

# ------------------------------------------------------------------ column
cols = ('p r T n_tot rho Kzz q_H2 q_H X_He X_C X_N X_O F_H F_He q_H2O q_CO '
        'q_CO2 q_CH4 q_NH3 q_HCN q_N2 q_OH q_C2H2').split()
pr = np.loadtxt(os.path.join(run, 'lower_atmosphere_profile.dat'))
pr = pr[np.argsort(pr[:, 0])]                 # top (p_match) first
pb = pr[:, 0]
Tc = pr[:, cols.index('T')]
mbar = pr[:, cols.index('rho')] / pr[:, cols.index('n_tot')]
rc = pr[:, cols.index('r')] * R_J
Mp = float([l for l in open(os.path.join(run, 'EXHALE_resolved.out'))
            if l.startswith('planet_mass_MJ')][0].split()[1]) * 1.89813e30
gc = G * Mp / rc**2
xsp = {'H2O': 'q_H2O', 'CH4': 'q_CH4', 'NH3': 'q_NH3', 'C2H2': 'q_C2H2',
       'HCN': 'q_HCN', 'CO': 'q_CO', 'H2': 'q_H2'}
mass = {'H2O': 18.01, 'CH4': 16.03, 'NH3': 17.03, 'C2H2': 26.02,
        'HCN': 27.01, 'CO': 28.0, 'H2': 2.016}
qid = {'H2O': 1, 'CH4': 32, 'NH3': 45, 'C2H2': 76, 'HCN': 70, 'CO': 26,
       'H2': 103}
x = {s: pr[:, cols.index(v)] for s, v in xsp.items()}
q_known = sum(x[s] for s in xsp) + pr[:, cols.index('q_H')] \
    + pr[:, cols.index('q_CO2')] + pr[:, cols.index('q_N2')] \
    + pr[:, cols.index('q_OH')]
x_He = np.clip(1 - q_known, 0, 1)
# layers between consecutive profile levels, merged below 3e-6 bar to
# Delta ln p >= 0.35 (the warmed layer keeps the profile's own levels)
keep_lv = [0]
for i in range(1, pb.size):
    if pb[i] <= 3e-6 or np.log(pb[i] / pb[keep_lv[-1]]) >= 0.35 or i == pb.size - 1:
        keep_lv.append(i)
keep_lv = np.array(keep_lv)
f_lev = lambda y: y / (mbar * gc)
def column_between(y):
    """int y dp/(m g) between the kept levels (trapezoid on the fine levels)."""
    fl = f_lev(y)
    seg = 0.5 * (fl[1:] + fl[:-1]) * (pb[1:] - pb[:-1]) * 1e6
    cs = np.concatenate([[0.0], np.cumsum(seg)])
    return np.diff(cs[keep_lv])
dN = {s: column_between(x[s]) for s in xsp}
Ntot = column_between(np.ones_like(pb))
TN = column_between(Tc)
lay_T = TN / Ntot
lay_p = np.sqrt(pb[keep_lv][1:] * pb[keep_lv][:-1])
nden = pb * 1e6 / (KB * Tc)
dCIA = {'H2-H2': column_between(x['H2']**2 * nden),
        'H2-He': column_between(x['H2'] * x_He * nden)}
pb = pb[keep_lv]
nlay = lay_p.size
print('column: %d layers, p %.3e - %.3e bar; N_CH4 above 1 mbar %.3e, total %.3e'
      % (nlay, pb[0], pb[-1], dN['CH4'][lay_p < 1e-3].sum(), dN['CH4'].sum()))


def read_par(path):
    a = np.genfromtxt(path, delimiter=[3, 12, 10, 10, 5, 5, 10, 4],
                      usecols=(1, 2, 4, 6, 7), dtype=float)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4]


def qtips(i):
    t = np.loadtxt(os.path.join(hit, 'q%d.txt' % i))
    return lambda T: np.interp(T, t[:, 0], t[:, 1])


ABS = {}
for s in xsp:
    nu, S, ga, Epp, na = read_par(os.path.join(hit, s + '.par'))
    Q = qtips(qid[s])
    ok = np.isfinite(Epp) & (S > 0)
    ABS[s] = dict(nu=nu[ok], S=S[ok], ga=ga[ok], Epp=Epp[ok],
                  na=np.where(np.isfinite(na[ok]), na[ok], 0.7), Q=Q)
    print('%s: %d lines' % (s, ok.sum()))


def S_T(d, T):
    T0 = 296.0
    return (d['S'] * d['Q'](T0) / d['Q'](T) * np.exp(-C2 * d['Epp'] * (1 / T - 1 / T0))
            * (1 - np.exp(-C2 * d['nu'] / T)) / (1 - np.exp(-C2 * d['nu'] / T0)))


cia = {}
for k in ('H2-H2', 'H2-He'):
    with h5py.File(os.path.join(cia_dir, k + '.h5'), 'r') as h:
        cia[k] = (1e4 / np.array(h['wavelengths']), np.array(h['T']),
                  np.array(h['log10xs']))


def cia_xs(k, nu, T):
    nug, Tg, lx = cia[k]
    o = np.argsort(nug)
    iT = np.clip(np.searchsorted(Tg, T) - 1, 0, Tg.size - 2)
    f = np.clip((T - Tg[iT]) / (Tg[iT + 1] - Tg[iT]), 0, 1)
    l = (1 - f) * lx[:, iT] + f * lx[:, iT + 1]
    return 10**np.interp(nu, nug[o], l[o], left=-99, right=-99)


# ------------------------------------------------------------------ spectral sweep
dv = 0.1e5 / CL                       # d ln nu
nu_lo, nu_hi = 1000.0, 10000.0
if '--nu' in sys.argv:
    k = sys.argv.index('--nu')
    nu_lo, nu_hi = float(sys.argv[k + 1]), float(sys.argv[k + 2])
nchunk_pts = 2**15
lnmin = np.log(nu_lo)
ntot = int(np.log(nu_hi / nu_lo) / dv)
WING = 25.0
absorbed = np.zeros(nlay)
absorbed_by = {s: np.zeros(nlay) for s in list(xsp) + ['CIA']}
emitted = 0.0
band_edges = np.array([1000, 1500, 2000, 2400, 2700, 3000, 3300, 3600, 4500,
                       5500, 6500, 10000.])
band_emit = np.zeros(band_edges.size - 1)
band_trans_warm = np.zeros(band_edges.size - 1)
i_warm = np.searchsorted(lay_p, 1e-6 * np.exp(0.5))   # layers above p_match e^0.5
band_abs_warm = np.zeros(band_edges.size - 1)
band_trans_bot = np.zeros(band_edges.size - 1)
for c0 in range(0, ntot, nchunk_pts):
    n = min(nchunk_pts, ntot - c0)
    nu = np.exp(lnmin + (c0 + np.arange(n)) * dv)
    a, b = nu[0], nu[-1]
    # emission on the chunk
    m = (nuL >= a - 1) & (nuL <= b + 1)
    E = np.zeros(n)
    if m.any():
        for Tlo, Thi in zip(Tg_edges[:-1], Tg_edges[1:]):
            mm = m & (TL >= Tlo) & (TL < Thi)
            if not mm.any():
                continue
            Tm = np.sum(PL[mm] * TL[mm]) / PL[mm].sum()
            sig_ln = np.sqrt(KB * Tm / (3.02293 * AMU)) / CL   # Gaussian std in ln nu
            idx = (np.log(nuL[mm]) - lnmin) / dv - c0
            h = np.zeros(n)
            i0 = np.floor(idx).astype(int)
            fr = idx - i0
            for ii, ww in ((i0, 1 - fr), (i0 + 1, fr)):
                ok = (ii >= 0) & (ii < n)
                np.add.at(h, ii[ok], PL[mm][ok] * ww[ok])
            k = np.arange(-int(6 * sig_ln / dv) - 1, int(6 * sig_ln / dv) + 2)
            ker = np.exp(-0.5 * (k * dv / sig_ln)**2)
            E += np.convolve(h, ker / ker.sum(), mode='same')
    if E.sum() <= 0:
        continue
    emitted += E.sum()
    # optical depth of every layer on the chunk
    nu_c = np.sqrt(a * b)
    dnu_c = nu_c * dv
    kpad = int(WING / dnu_c) + 1
    tau = np.zeros((nlay, n))
    tau_s = {}
    for s, d in ABS.items():
        ms = (d['nu'] >= a - WING) & (d['nu'] <= b + WING)
        if not ms.any():
            continue
        dd = {k: v[ms] for k, v in d.items() if k != 'Q'}
        dd['Q'] = d['Q']
        idx = (np.log(dd['nu']) - lnmin) / dv - c0 + kpad
        i0 = np.floor(idx).astype(int)
        fr = idx - i0
        ts = np.zeros((nlay, n))
        okd = (i0 >= 0) & (i0 + 1 < n + 2 * kpad)
        i0, fr = i0[okd], fr[okd]
        dd = {k: (v[okd] if k != 'Q' else v) for k, v in dd.items()}
        if i0.size == 0:
            continue
        for j in range(nlay):
            if dN[s][j] <= 0:
                continue
            Sj = S_T(dd, lay_T[j]) * dN[s][j]
            # Doppler line-centre tau S N /(sqrt(pi) sigma_D sqrt 2) < 1e-7
            if Sj.max() / (nu_c * 2.5e-6) < 1e-7:
                continue
            h = (np.bincount(i0, weights=Sj * (1 - fr), minlength=n + 2 * kpad)
                 + np.bincount(i0 + 1, weights=Sj * fr, minlength=n + 2 * kpad))[:n + 2 * kpad]
            sigD = nu_c * np.sqrt(KB * lay_T[j] / (mass[s] * AMU)) / CL
            gam = np.sum(Sj * dd['ga'] * (296.0 / lay_T[j])**dd['na']) / Sj.sum() * lay_p[j] / 1.01325
            half = min(kpad, int((8 * sigD + 400 * gam) / dnu_c) + 2)
            kk = np.arange(-half, half + 1) * dnu_c
            ker = voigt_profile(kk, sigD, gam) * dnu_c
            conv = oaconvolve(h, ker, mode='same')[kpad:kpad + n] / dnu_c
            ts[j] = conv
        tau += ts
        tau_s[s] = ts
    tc = np.zeros((nlay, n))
    for j in range(nlay):
        tc[j] = (cia_xs('H2-H2', nu, lay_T[j]) * dCIA['H2-H2'][j]
                 + cia_xs('H2-He', nu, lay_T[j]) * dCIA['H2-He'][j])
    tau += tc
    tau_s['CIA'] = tc
    cum = np.vstack([np.zeros(n), np.cumsum(tau, axis=0)])     # at levels
    Tr = np.zeros_like(cum)
    for wm, mc in zip(w_mu, mu_c):
        Tr += wm * np.exp(-cum / mc)
    dabs = (Tr[:-1] - Tr[1:]) * E[None, :]                     # (nlay, n)
    absorbed += dabs.sum(axis=1)
    tot = np.maximum(tau, 1e-300)
    for s, ts in tau_s.items():
        absorbed_by[s] += (dabs * ts / tot).sum(axis=1)
    bi = np.clip(np.searchsorted(band_edges, nu) - 1, 0, band_edges.size - 2)
    band_emit += np.bincount(bi, weights=E, minlength=band_edges.size - 1)
    band_abs_warm += np.bincount(bi, weights=(Tr[0] - Tr[i_warm]) * E, minlength=band_edges.size - 1)
    band_trans_bot += np.bincount(bi, weights=Tr[-1] * E, minlength=band_edges.size - 1)
    print('chunk %.1f-%.1f cm^-1: emitted %.3e, absorbed above 1e-3 bar %.3e, below %.3e, transmitted %.3e'
          % (a, b, E.sum(), dabs[lay_p < 1e-3].sum(), dabs[lay_p >= 1e-3].sum(), (Tr[-1] * E).sum()),
          flush=True)

print('emitted on the grid %.4f of %.4f' % (emitted, P_down))
f_abs = absorbed / emitted
with open(out + '_layers.txt', 'w') as fo:
    fo.write('# shape %s; downward H3+ power %.5e erg cm^-2 s^-1 (code Lambda), on grid %.5e\n'
             % (shape, P_down, emitted))
    fo.write('# columns: p_top[bar] p_bot[bar] T_layer[K] absorbed_fraction cumulative '
             + ' '.join('by_' + s for s in absorbed_by) + '\n')
    cumf = np.cumsum(f_abs)
    for j in range(nlay):
        fo.write('%.5e %.5e %7.2f %.5e %.5e ' % (pb[j], pb[j + 1], lay_T[j], f_abs[j], cumf[j])
                 + ' '.join('%.3e' % (absorbed_by[s][j] / emitted) for s in absorbed_by) + '\n')
    fo.write('# transmitted below %.3e bar: %.5e\n' % (pb[-1], 1 - cumf[-1]))
    fo.write('# bands [cm^-1]: emitted share, absorbed above p_match*e^0.5 (fraction of the band), transmitted below the deepest level\n')
    for k in range(band_edges.size - 1):
        fo.write('#   %6.0f-%6.0f  %.4f  %.3e  %.3e\n'
                 % (band_edges[k], band_edges[k + 1], band_emit[k] / emitted,
                    band_abs_warm[k] / max(band_emit[k], 1e-300),
                    band_trans_bot[k] / max(band_emit[k], 1e-300)))


def share(p1, p2):
    # a layer belongs to the interval that holds its lower boundary
    m = (pb[1:] > p1 * 1.0000001) & (pb[1:] <= p2 * 1.0000001)
    return f_abs[m].sum()


for lab, p1, p2 in (('warmed layer p_match..p_match e^0.5', 1e-6, 1e-6 * np.exp(0.5)),
                    ('p_match e^0.5 .. 1e-3 bar', 1e-6 * np.exp(0.5), 1e-3),
                    ('1e-3 .. 0.1 bar', 1e-3, 0.1),
                    ('0.1 .. 1 bar', 0.1, 1.0),
                    ('1 bar .. deepest', 1.0, pb[-1] * 1.01)):
    print('%-38s %.4e  (%.4e erg cm^-2 s^-1)' % (lab, share(p1, p2), share(p1, p2) * P_down))
print('transmitted below %.2f bar: %.4e' % (pb[-1], 1 - f_abs.sum()))
