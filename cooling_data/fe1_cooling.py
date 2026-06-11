"""Fe I collisional line-cooling builder (Huang et al. 2023, Section 2.5).

CHIANTI v11 has no Fe I model atom, so - following Huang - Fe I line cooling is
built from NIST oscillator strengths + the Van Regemorter approximation, with
the lower levels populated in Boltzmann equilibrium over the metastable
manifold. Inputs are the auditable raw NIST dumps from fetch_fe1_nist.py:
    fe1_nist_lines.tsv    permitted (E1) lines: Aki, g_i, g_k, E_i, E_k
    fe1_nist_levels.tsv   every level (g, E) for the partition function

Cooling coefficient per (n_e * n_FeI_total), optically thin [erg cm^3 s^-1]:

    Lambda(T) = sum_lines  f_l(T) * q_lu(T) * dE_lu ,

    f_l(T)  = g_l exp(-E_l/kT) / U(T)        Boltzmann lower-level fraction
    U(T)    = sum_{even levels, E<e_cut} g exp(-E/kT)
    q_lu    = 2.16 a^-1.68 exp(-a) T^-3/2 f_lu   (Van Regemorter compact,
              Huang Eq. 10; a = dE/kT), f_lu from Aki + g's + lambda.

Metastable manifold = EVEN-parity levels below the lowest odd level (z7D deg at
19351 cm^-1): below that energy there is no lower level for an E1 decay, so the
even levels are genuinely metastable and hold the Fe I population. (A pure
energy cut, as used for Fe II where the cut conveniently sits at the first odd
level, would here wrongly include the radiatively-depopulated z7D odd levels.)
Upper levels of the cooling lines are odd (E1) and assumed to decay radiatively.
"""
import csv
import io
import os

import numpy as np

from chianti_cooling import (van_regemorter_clu_compact, EV_ERG, K_B_ERG,
                             HC_OVER_K)

OUTDIR = os.path.dirname(os.path.abspath(__file__))
CM1_TO_EV = 1.0 / 8065.543937   # eV per cm^-1
F_CONST = 1.4992e-16            # f_lu = F_CONST * (g_k/g_i) * Aki * lambda_A^2

# Metastable-manifold cutoff: the lowest odd-parity level (z7D deg, 19351 cm^-1
# = 2.40 eV). Below it every Fe I level is even -> truly metastable.
ECUT_CM1_DEFAULT = 19351.0


def _num(s):
    s = s.strip().strip('"').replace('=', '').replace('[', '').replace(']', '')
    if s in ('', '*'):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _read_tsv(path):
    rows = [ln for ln in open(path) if not ln.startswith('#')]
    return list(csv.DictReader(io.StringIO(''.join(rows)), delimiter='\t'))


def load_levels(path=None):
    """Even-parity Fe I levels as (E_cm1, g) arrays (the metastable manifold)."""
    path = path or os.path.join(OUTDIR, 'fe1_nist_levels.tsv')
    E, g = [], []
    for r in _read_tsv(path):
        term = (r.get('Term') or '').strip().strip('"')
        Ei = _num(r.get('Level (cm-1)', ''))
        gi = _num(r.get('g', ''))
        if Ei is None or gi is None:
            continue
        if '*' in term:           # odd parity: radiatively connected, skip
            continue
        E.append(Ei)
        g.append(gi)
    return np.array(E), np.array(g)


def load_lines(path=None):
    """Permitted (E1) Fe I lines with an even, low lower level.

    Returns arrays Ei, dE_eV, gi, gk, Aki for lines that have a transition
    probability, blank Type (E1), and an even-parity lower term.
    """
    path = path or os.path.join(OUTDIR, 'fe1_nist_lines.tsv')
    Ei, Ek, gi, gk, Aki = [], [], [], [], []
    for r in _read_tsv(path):
        if (r.get('Type') or '').strip() != '':      # keep only E1 (blank)
            continue
        a = _num(r.get('Aki(s^-1)', ''))
        ei = _num(r.get('Ei(cm-1)', ''))
        ek = _num(r.get('Ek(cm-1)', ''))
        g_i = _num(r.get('g_i', ''))
        g_k = _num(r.get('g_k', ''))
        ti = (r.get('term_i') or '').strip().strip('"')
        if None in (a, ei, ek, g_i, g_k) or ek <= ei:
            continue
        if '*' in ti:             # lower level odd -> not metastable, skip
            continue
        Ei.append(ei); Ek.append(ek); gi.append(g_i); gk.append(g_k); Aki.append(a)
    Ei = np.array(Ei); Ek = np.array(Ek)
    gi = np.array(gi); gk = np.array(gk); Aki = np.array(Aki)
    lam_A = 1.0e8 / (Ek - Ei)                          # vacuum, from dE
    dE_eV = (Ek - Ei) * CM1_TO_EV
    f_lu = F_CONST * (gk / gi) * Aki * lam_A**2
    return Ei, dE_eV, gi, f_lu


def cooling_FeI(T, e_cut_cm1=ECUT_CM1_DEFAULT):
    """Fe I line cooling per (n_e * n_FeI) [erg cm^3 s^-1]; T array [K]."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    Elev, glev = load_levels()
    sel = Elev < e_cut_cm1
    Elev, glev = Elev[sel], glev[sel]
    # partition function over the even metastable manifold, U(T)
    U = np.zeros_like(T)
    for e, g in zip(Elev, glev):
        U += g * np.exp(-e * HC_OVER_K / T)

    Ei, dE_eV, gi, f_lu = load_lines()
    keep = Ei < e_cut_cm1
    Ei, dE_eV, gi, f_lu = Ei[keep], dE_eV[keep], gi[keep], f_lu[keep]

    total = np.zeros_like(T)
    for k in range(Ei.size):
        f_l = gi[k] * np.exp(-Ei[k] * HC_OVER_K / T) / U
        q_lu = van_regemorter_clu_compact(T, dE_eV[k], f_lu[k])
        total += f_l * q_lu * (dE_eV[k] * EV_ERG)
    return total, int(Ei.size), int(Elev.size)


# 41-point ATES cooling grid (Cool_coeff.f90: logT 3.00..5.00, dlogT 0.05).
COOL_LOGT = np.arange(3.0, 5.0 + 1e-9, 0.05)


def fortran_block(name, logL):
    out = [f"   real*8, parameter :: {name}(NCOOLT) = [ &"]
    for i in range(0, len(logL), 5):
        chunk = logL[i:i + 5]
        last = i + 5 >= len(logL)
        body = ',   '.join(f'{v:11.5f}d0' for v in chunk)
        out.append('          ' + body + (' ]' if last else ', &'))
    return '\n'.join(out)


if __name__ == '__main__':
    T = 10.0**COOL_LOGT
    for ecut, lab in [(19351.0, '2.40 eV (1st odd, default)'),
                      (25900.0, '3.21 eV (1st quintet-odd)'),
                      (38459.0, '4.77 eV (Fe II analog)')]:
        L, nline, nlev = cooling_FeI(T, ecut)
        i1, i2 = np.argmin(abs(T - 5000.)), np.argmin(abs(T - 1e4))
        i3 = np.argmin(abs(T - 2e4))
        print(f"e_cut={lab:28s} nlines={nline:4d} nlev={nlev:3d}  "
              f"L(5e3)={L[i1]:.3e} L(1e4)={L[i2]:.3e} L(2e4)={L[i3]:.3e}")

    L, nline, nlev = cooling_FeI(T, ECUT_CM1_DEFAULT)
    logL = np.log10(np.clip(L, 1e-99, None))
    print(f"\n# Fe I line cooling, NIST f-values + Van Regemorter (Huang Fig 5)")
    print(f"# manifold: {nlev} even levels < {ECUT_CM1_DEFAULT:.0f} cm^-1; "
          f"{nline} permitted lines")
    print(fortran_block('cool_logL_FeI', logL))
