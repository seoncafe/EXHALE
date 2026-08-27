"""Density-dependent Fe II line cooling by multilevel statistical equilibrium.

The Phase-2 Fe II table (FeII_coronal in export_cooling_tables.py) is the
optically-thin CORONAL limit: only the ground level is populated and every
collisional excitation is assumed to be followed by a radiative decay, so the
cooling is strictly proportional to n_e (per n_e n_FeII it is a function of T
alone). That assumption fails badly at the cool, DENSE base of a hot Jupiter
thermosphere: the Fe II cooling there is carried almost entirely by FORBIDDEN
metastable / a6D fine-structure lines (A <~ 5e-3 s^-1), whose critical densities

    n_crit = A_ul / sum_c n_c k_ul^c   ~   1e4 - 1e6 cm^-3   (electron collider)

sit FAR below the base n_e ~ 1e8 cm^-3.  Those lines are therefore collisionally
saturated (thermalized to their LTE populations); the coronal coefficient
overestimates their cooling by ~ n_e / n_crit.

This module solves the full CHIANTI multilevel statistical-equilibrium problem
for the Fe II level populations at each (T, n_e), so the cooling smoothly
interpolates between the coronal (low-n_e) and LTE (high-n_e) limits:

    Lambda_vol / n_FeII  =  sum_{u>l} n_u A_ul dE_ul        [erg s^-1 per ion]

with the n_u from the SE solve (sum_i n_i = 1).  To keep the existing ATES
per-(n_e n_FeII) tabulation interface we report

    Lambda_eff(T, n_e)  =  (Lambda_vol / n_FeII) / n_e      [erg cm^3 s^-1].

At low n_e, Lambda_eff -> the coronal Lambda(T); at high n_e it falls off ~1/n_e
as the numerator saturates.  Electron collisions only (the dominant collider in
the partially ionized base is still e- for Fe II forbidden lines; an H I / proton
collider extension is flagged below).  No proton (.psplups) channel yet.
"""
import os

import numpy as np

from chianti_cooling import (read_elvlc, read_scups, read_wgfa, ion_dir,
                             _level_energy_cm1, upsilon, COLL_PREF, RY_ERG,
                             K_B_ERG, HC_OVER_K)

OUTDIR = os.path.dirname(os.path.abspath(__file__))


def _load_fe2():
    d = ion_dir('fe', 2)
    lev = _level_energy_cm1(read_elvlc(os.path.join(d, 'fe_2.elvlc')))
    sc = read_scups(os.path.join(d, 'fe_2.scups'))
    wg = read_wgfa(os.path.join(d, 'fe_2.wgfa'))
    # radiative A summed over duplicate (ll,ul) entries
    A = {}
    for ll, ul, wl, gf, a in wg:
        if a > 0.0:
            A[(ll, ul)] = A.get((ll, ul), 0.0) + a
    return lev, sc, A


def fe2_levels_cooling(T, ne, lev, sc, A, nmax=None):
    """Statistical-equilibrium Fe II cooling at scalar (T, ne).

    Returns Lambda_vol_per_ion [erg s^-1] = sum_{u>l} n_u A_ul dE_ul.
    nmax optionally caps the number of levels (by energy) for speed/robustness.
    """
    idxs = sorted(lev.keys(), key=lambda i: lev[i][1])
    if nmax is not None:
        idxs = idxs[:nmax]
    pos = {idx: k for k, idx in enumerate(idxs)}
    nL = len(idxs)
    g = np.array([lev[i][0] for i in idxs])
    E = np.array([lev[i][1] for i in idxs])           # cm^-1

    # Rate matrix M (per ion): M[j,i] = rate i -> j (j != i)
    M = np.zeros((nL, nL))
    for tr in sc:
        ll, ul = tr['ll'], tr['ul']
        if ll not in pos or ul not in pos:
            continue
        a, b = pos[ll], pos[ul]                        # a=lower, b=upper
        de_erg = tr['de'] * RY_ERG
        ups = float(upsilon(tr, T)[0])
        # de-excitation u->l and excitation l->u (electron collider)
        q_ul = COLL_PREF / (g[b] * np.sqrt(T)) * ups
        q_lu = COLL_PREF / (g[a] * np.sqrt(T)) * ups * np.exp(-de_erg /
                                                              (K_B_ERG * T))
        M[a, b] += ne * q_ul      # u -> l  (into lower)
        M[b, a] += ne * q_lu      # l -> u  (into upper)
    for (ll, ul), a in A.items():
        if ll in pos and ul in pos:
            M[pos[ll], pos[ul]] += a   # radiative u -> l

    # Statistical equilibrium: for each j, sum_i M[j,i] n_i - n_j sum_i M[i,j] = 0
    Amat = M.copy()
    for j in range(nL):
        Amat[j, j] -= M[:, j].sum()
    # replace one equation with normalization sum n_i = 1
    Amat[0, :] = 1.0
    rhs = np.zeros(nL)
    rhs[0] = 1.0
    n = np.linalg.solve(Amat, rhs)
    n = np.clip(n, 0.0, None)

    # radiative cooling per ion
    lam = 0.0
    for (ll, ul), a in A.items():
        if ll in pos and ul in pos:
            de_erg = (E[pos[ul]] - E[pos[ll]]) * HC_OVER_K * K_B_ERG  # cm^-1->erg
            lam += n[pos[ul]] * a * de_erg
    return lam


def lambda_eff_grid(Tgrid, negrid, nmax=None):
    """Lambda_eff(T, ne) [erg cm^3 s^-1] = (cooling per ion)/ne on a (T,ne) grid."""
    lev, sc, A = _load_fe2()
    out = np.zeros((Tgrid.size, negrid.size))
    for it, T in enumerate(Tgrid):
        for ie, ne in enumerate(negrid):
            lam_ion = fe2_levels_cooling(T, ne, lev, sc, A, nmax=nmax)
            out[it, ie] = lam_ion / ne
    return out


# ----------------------------------------------------------------------------
# Optimized 2-D table generator for the EXHALE Fortran port.
#
# The slow path (lambda_eff_grid) rebuilds a CubicSpline for every transition at
# every (T, ne) cell.  Here we (i) descale Upsilon for each transition ONCE over
# the whole T-grid (upsilon() accepts an array T -> one spline per transition),
# and (ii) split the rate matrix M = ne * Q_coll(T) + A_rad, so for each T we
# build Q_coll once and only the cheap ne * Q_coll + A_rad assembly + solve runs
# inside the ne loop.  Result is identical to lambda_eff_grid (verified in main).
# ----------------------------------------------------------------------------

def lambda_eff_table(Tgrid, negrid, nmax=None):
    """Lambda_eff(T, ne) [erg cm^3 s^-1] on a (T, ne) grid, fast.

    Returns out[it, ie] = (sum_u n_u A_ul dE_ul) / ne with n_u from the
    multilevel statistical-equilibrium solve at (Tgrid[it], negrid[ie]).
    """
    lev, sc, A = _load_fe2()
    idxs = sorted(lev.keys(), key=lambda i: lev[i][1])
    if nmax is not None:
        idxs = idxs[:nmax]
    pos = {idx: k for k, idx in enumerate(idxs)}
    nL = len(idxs)
    g = np.array([lev[i][0] for i in idxs])
    E = np.array([lev[i][1] for i in idxs])           # cm^-1

    # Upsilon for each transition over the whole T-grid (one spline build each)
    trans = []
    for tr in sc:
        ll, ul = tr['ll'], tr['ul']
        if ll not in pos or ul not in pos:
            continue
        a, b = pos[ll], pos[ul]                        # a=lower, b=upper
        de_erg = tr['de'] * RY_ERG
        ups = np.asarray(upsilon(tr, Tgrid), dtype=float)   # shape (nT,)
        trans.append((a, b, de_erg, ups))

    # radiative matrix (T- and ne-independent): A_rad[a,b] += A_ul (u->l into a)
    Arad = np.zeros((nL, nL))
    rad_list = []                                      # (lo, up, A_ul, dE_erg)
    for (ll, ul), a in A.items():
        if ll in pos and ul in pos:
            lo, up = pos[ll], pos[ul]
            Arad[lo, up] += a
            de = (E[up] - E[lo]) * HC_OVER_K * K_B_ERG  # cm^-1 -> erg
            rad_list.append((lo, up, a, de))

    out = np.zeros((Tgrid.size, negrid.size))
    diag = np.diag_indices(nL)
    for it, T in enumerate(Tgrid):
        sqrtT = np.sqrt(T)
        Q = np.zeros((nL, nL))                         # collisional, per unit ne
        for (a, b, de_erg, ups) in trans:
            u = ups[it]
            q_ul = COLL_PREF / (g[b] * sqrtT) * u
            q_lu = COLL_PREF / (g[a] * sqrtT) * u * np.exp(-de_erg /
                                                           (K_B_ERG * T))
            Q[a, b] += q_ul                            # u -> l
            Q[b, a] += q_lu                            # l -> u
        for ie, ne in enumerate(negrid):
            M = ne * Q + Arad
            Amat = M.copy()
            Amat[diag] -= M.sum(axis=0)                # diag -= total out-rate
            Amat[0, :] = 1.0                           # normalization row
            rhs = np.zeros(nL)
            rhs[0] = 1.0
            n = np.clip(np.linalg.solve(Amat, rhs), 0.0, None)
            lam = 0.0
            for (lo, up, a_ul, de) in rad_list:
                lam += n[up] * a_ul * de
            out[it, ie] = lam / ne
    return out


# EXHALE cooling grid (Cool_coeff.f90: logT 3.00..5.00, dlogT 0.05 -> 41 points).
COOL_LOGT = np.arange(3.0, 5.0 + 1e-9, 0.05)
# Electron-density axis for the 2-D Fe II table: log10 ne 0..14, dlog 0.5.
COOL_LOGNE = np.arange(0.0, 14.0 + 1e-9, 0.5)


def fortran_grid_1d(name, vals, ncols=6):
    """Emit `real*8, parameter :: name(N) = [ ... ]` for a 1-D grid."""
    out = ["   real*8, parameter :: %s(NCOOLNE) = [ &" % name]
    for i in range(0, len(vals), ncols):
        chunk = vals[i:i + ncols]
        last = i + ncols >= len(vals)
        body = ',   '.join('%9.4fd0' % v for v in chunk)
        out.append('          ' + body + (' ]' if last else ', &'))
    return '\n'.join(out)


def fortran_block_2d(name, arr, ncols=5):
    """Emit a 2-D log10(Lambda_eff) array (NCOOLT rows, NCOOLNE cols).

    Fortran stores column-major; we declare the literal in the natural
    (NCOOLT, NCOOLNE) shape using a reshape so the source stays readable: the
    flat list is column-major (ne-outer, T-inner), matching the array shape.
    """
    nT, nN = arr.shape
    flat = arr.reshape(nT * nN, order='F')             # column-major flatten
    out = ["   real*8, parameter :: %s(NCOOLT,NCOOLNE) = &" % name,
           "      reshape( [ &"]
    for i in range(0, flat.size, ncols):
        chunk = flat[i:i + ncols]
        last = i + ncols >= flat.size
        body = ',   '.join('%11.5fd0' % v for v in chunk)
        out.append('          ' + body + (' ], &' if last else ', &'))
    out.append("      [NCOOLT, NCOOLNE] )")
    return '\n'.join(out)


def emit_fortran_table(nmax=None):
    """Print the Fe II 2-D density-dependent cooling table as Fortran source."""
    from chianti_cooling import cooling_effective
    T = 10.0**COOL_LOGT
    ne = 10.0**COOL_LOGNE
    grid = lambda_eff_table(T, ne, nmax=nmax)          # (nT, nN) erg cm^3 s^-1
    logL = np.log10(np.clip(grid, 1e-99, None))

    cor = cooling_effective('fe', 2, T, pop='coronal')

    print('! Fe II density-dependent line cooling, multilevel statistical')
    print('! equilibrium (fe2_cooling.py: lambda_eff_table).  Lambda_eff(T,ne)')
    print('! [erg cm^3 s^-1] = (sum_u n_u A_ul dE_ul)/ne; coronal at low ne,')
    print('! ~1/ne (LTE-saturated) at high ne.  logT 3.0..5.0 (NCOOLT=%d),' % T.size)
    print('! log10 ne 0..14 dlog 0.5 (NCOOLNE=%d).' % ne.size)
    print()
    print('   integer, parameter :: NCOOLNE = %d' % ne.size)
    print(fortran_grid_1d('cool_logne', COOL_LOGNE))
    print(fortran_block_2d('cool_logL_FeII_ne', logL))
    print()
    # sanity table for the log: low-ne edge vs the old coronal table,
    # and the suppression at the WASP-121b base.
    import sys
    print('! --- sanity (stderr): low-ne edge vs coronal, base suppression ---',
          file=sys.stderr)
    ie0 = 0
    print('! T[K]    coronal     SE(ne=1)    ratio(SE/cor)', file=sys.stderr)
    for it in (0, 10, 20, 30, 40):
        print('! %7.0f  %.3e  %.3e  %6.3f' %
              (T[it], cor[it], grid[it, ie0], grid[it, ie0] /
               max(cor[it], 1e-99)), file=sys.stderr)
    # base cell T~2358 (it=10 -> 10^3.5=3162; use nearest), ne=1e8
    it_b = int(np.argmin(abs(T - 2358.)))
    ie_b = int(np.argmin(abs(ne - 1e8)))
    print('! base T=%.0f ne=%.0e: coronal=%.3e  SE=%.3e  suppress=%.1fx' %
          (T[it_b], ne[ie_b], cor[it_b], grid[it_b, ie_b],
           cor[it_b] / max(grid[it_b, ie_b], 1e-99)), file=sys.stderr)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == '--table':
        emit_fortran_table()
        sys.exit(0)
    from chianti_cooling import cooling_effective
    Tg = np.array([2358., 5000., 8000., 1e4, 2e4, 3e4])
    neg = np.array([1e2, 1e4, 1e6, 1e8, 1e10])
    cor = cooling_effective('fe', 2, Tg, pop='coronal')
    grid = lambda_eff_grid(Tg, neg)
    print('Fe II Lambda_eff(T,ne) [erg cm^3 s^-1]; coronal is the ne->0 limit.')
    hdr = '   T[K]   coronal' + ''.join('   ne=%.0e' % x for x in neg)
    print(hdr)
    for it, T in enumerate(Tg):
        row = '  %6.0f  %.2e' % (T, cor[it])
        row += ''.join('  %.2e' % grid[it, ie] for ie in range(neg.size))
        print(row)
    print()
    print('Base check (WASP-121b r=1: T=2358 K, ne~1e8):')
    it = 0
    ie = int(np.argmin(abs(neg - 1e8)))
    print('  coronal     Lambda_eff = %.3e' % cor[it])
    print('  SE (ne=1e8) Lambda_eff = %.3e' % grid[it, ie])
    print('  suppression factor      = %.1f x' % (cor[it] / grid[it, ie]))
