"""Line-by-line H2 Lyman-Werner shielding, with and without line overlap.

Independent of both CLOUDY and the Meudon PDR code: it uses only the H2 UV
line list and level data that ship with the Meudon PDR code (Abgrall, Roueff
& Drira 2000), LTE level populations at a fixed T, and a normally incident
flat-F_lambda beam attenuated through a slab.

  tau(nu, N) = N * sum_lines x_lower * sigma_line(nu)      (overlap included)
  tau_i(nu, N) = N * x_lower(i) * sigma_i(nu)              (one line at a time)

  k_diss(N) = sum_i x_i (pi e^2/m c) f_i p_diss,i * int phi_i(nu) 4pi J_nu/(h nu) dnu
"""
import os, sys
import numpy as np
from scipy.special import wofz

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'PDR71', 'data')
CLIGHT = 2.99792458e10
KB = 1.380649e-16
M_H2 = 2.0 * 1.67353284e-24
PI_E2_MC = 0.0265400            # pi e^2 / (m_e c)  [cm^2 Hz]
HPL = 6.62607015e-27

STATES = ['B', 'Cp', 'Cm', 'Bp', 'Dp', 'Dm']


def read_X_levels():
    g, E_K, v, J = [], [], [], []
    for line in open(os.path.join(DATA, 'Levels', 'level_h2.dat')):
        if line.startswith(('#', '@')):
            continue
        p = line.split()
        if len(p) < 6:
            continue
        g.append(float(p[1])); E_K.append(float(p[2]))
        v.append(int(p[4])); J.append(int(p[5]))
    return np.array(g), np.array(E_K), np.array(v), np.array(J)


def read_upper(state):
    f = open(os.path.join(DATA, 'Levels', 'lev_H2_%s.dat' % state))
    rows = []
    for line in f:
        if line.startswith('#'):
            continue
        p = line.split()
        if len(p) == 3:            # "Nb v_max J_max"
            continue
        if len(p) < 6:
            continue
        rows.append((int(p[0]), int(p[1]), int(p[2]), float(p[3]),
                     float(p[4]), float(p[5])))
    a = np.array(rows)
    return dict(idx=a[:, 0].astype(int), vu=a[:, 1].astype(int),
                Ju=a[:, 2].astype(int), Ecm=a[:, 3], Gam=a[:, 4], pdis=a[:, 5])


def read_lines(state):
    f = open(os.path.join(DATA, 'UVdata', 'UV_H2_%s.dat' % state))
    rows = []
    for line in f:
        if line.startswith('#'):
            continue
        p = line.split()
        if len(p) == 1:
            continue
        if len(p) < 4:
            continue
        rows.append((int(p[0]), int(p[1]), float(p[2]), float(p[3])))
    a = np.array(rows)
    return a[:, 0].astype(int), a[:, 1].astype(int), a[:, 2], a[:, 3]


def build_linelist(T, wl_min=911.75, wl_max=1200.0):
    gX, EX, vX, JX = read_X_levels()
    w = gX * np.exp(-EX / T)
    x = w / w.sum()
    out = []
    for st in STATES:
        up = read_upper(st)
        lu, ll, Aij, EEcm = read_lines(st)
        umap = {i: k for k, i in enumerate(up['idx'])}
        for k in range(len(lu)):
            ku = umap.get(lu[k])
            if ku is None:
                continue
            il = ll[k] - 1
            if il >= len(gX):
                continue
            wl_A = 1.0e8 / EEcm[k]
            if not (wl_min <= wl_A <= wl_max):
                continue
            Ju = up['Ju'][ku]; Jl = JX[il]
            gu_over_gl = (2.0 * Ju + 1.0) / (2.0 * Jl + 1.0)
            f_lu = Aij[k] * gu_over_gl * wl_A ** 2 / 6.6702e15
            out.append((CLIGHT / (wl_A * 1e-8), f_lu, x[il],
                        up['Gam'][ku], up['pdis'][ku]))
    a = np.array(out)
    return dict(nu0=a[:, 0], f=a[:, 1], xl=a[:, 2], Gam=a[:, 3], pdis=a[:, 4])


def _profiles(nu, nu0, dnuD, a_voigt, Gam, core_hw):
    """Voigt profile of one line on the whole grid (Lorentz wing + Voigt core)."""
    phi = (Gam / (4.0 * np.pi ** 2)) / ((nu - nu0) ** 2 + (Gam / (4.0 * np.pi)) ** 2)
    i0 = np.searchsorted(nu, nu0 - core_hw * dnuD)
    i1 = np.searchsorted(nu, nu0 + core_hw * dnuD)
    z = (nu[i0:i1] - nu0) / dnuD + 1j * a_voigt
    phi[i0:i1] = wofz(z).real / (np.sqrt(np.pi) * dnuD)
    return phi


def shielding_curve(T, Ncols, nu_pts=400_000, wl_min=911.75, wl_max=1200.0,
                    core_hw=300.0, keep_frac=1.0 - 1e-6, verbose=False):
    """Returns (k_overlap, k_nooverlap) [arbitrary but common units] for each N."""
    ll = build_linelist(T, wl_min, wl_max)
    nu0, f, xl, Gam, pdis = ll['nu0'], ll['f'], ll['xl'], ll['Gam'], ll['pdis']
    w0 = xl * f * pdis
    order = np.argsort(w0)[::-1]
    c = np.cumsum(w0[order]) / w0.sum()
    nkeep = int(np.searchsorted(c, keep_frac) + 1)
    keep = order[:nkeep]
    nu0, f, xl, Gam, pdis = (a[keep] for a in (nu0, f, xl, Gam, pdis))
    nl = len(nu0)
    if verbose:
        print('  lines kept: %d of %d' % (nl, len(w0)), flush=True)
    b = np.sqrt(2.0 * KB * T / M_H2)
    dnuD = nu0 * b / CLIGHT
    a_voigt = Gam / (4.0 * np.pi * dnuD)

    nu = np.linspace(CLIGHT / (wl_max * 1e-8), CLIGHT / (wl_min * 1e-8), nu_pts)
    dnu = nu[1] - nu[0]
    phot = 1.0 / nu ** 3          # flat F_lambda -> photon number per unit nu

    sig = np.zeros(nu_pts)
    for i in range(nl):
        sig += PI_E2_MC * f[i] * xl[i] * _profiles(nu, nu0[i], dnuD[i],
                                                   a_voigt[i], Gam[i], core_hw)
    Ncols = np.atleast_1d(np.asarray(Ncols, dtype=float))
    att = np.exp(-np.outer(Ncols, sig))
    k_ov = np.zeros(len(Ncols))
    k_no = np.zeros(len(Ncols))
    for i in range(nl):
        phi = _profiles(nu, nu0[i], dnuD[i], a_voigt[i], Gam[i], core_hw)
        w = xl[i] * f[i] * pdis[i] * phi * phot * dnu
        k_ov += att @ w
        s_i = PI_E2_MC * f[i] * xl[i] * phi
        for j, N in enumerate(Ncols):
            k_no[j] += np.sum(w * np.exp(-N * s_i))
        if verbose and i % 200 == 0:
            print('  line %d/%d' % (i, nl), flush=True)
    return k_ov, k_no, ll
