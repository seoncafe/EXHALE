#!/usr/bin/env python3
"""Spencer-Fano degradation of photoelectrons in an H / He / H2 / e gas.

Deposition record of a primary electron of energy E0 (eV) in a gas of atomic
hydrogen H(1s), helium He(1^1S) and He(2^3S), molecular hydrogen H2(X, v=0,
thermal J), and thermal electrons, for every primary energy on a uniform grid
at once.

Method (Spencer & Fano 1954; adjoint or "yield" form as in Dalgarno, Yan &
Liu 1999, ApJS 125, 237, section 4, and Xu & McCray 1991, ApJ 375, 190):
for every recorded quantity q (a count or an energy), the expectation
Y_q(E) over the whole slowing-down history of an electron of energy E obeys

    R(E) Y_q(E) = sum_c nu_c(E) [ d_{c,q}(E) + Y_q(E - W_c(E)) + Y_q(s_c) ]
                + sum_t nu_t^ion(E) [ d_{t,q} + int_0^{(E-I_t)/2} P_t(eps|E)
                                       ( Y_q(eps) + Y_q(E - I_t - eps) ) deps ]
                + (L(E)/dE) [ dE l_q(E) + Y_q(E - dE) ]

with nu = n sigma (cm^-1), R the sum of all rates, W_c the energy lost in a
discrete collision c (negative for a superelastic one), s_c the energy of an
electron the event releases besides the scattered one (only the He photon
recycling of the Furlanetto & Stoever limit uses it), P_t the normalized
secondary-electron spectrum of an ionization of target t, and L(E) (eV/cm)
the continuous losses: Coulomb drag on thermal electrons, elastic recoil on
neutrals and the net H2 rotational excitation, with l_q their split. The
grid is uniform in E (bin dE); Y at an energy between two nodes is the linear
interpolation, which conserves energy exactly, so that the closure
sum_q(energy buckets) = E is a test of the bookkeeping, and the resolution
error shows in the values (tested by dE convergence). Electrons below the
floor energy E_floor (= 1.5 k T_e, at least one bin) are thermal; their
remaining energy is heat ("heat_thermalized").

Every discrete loss is kept channel by channel with its own threshold; only
the Coulomb, elastic-recoil and rotational losses (energy steps 0.04-0.2 eV
or continuous) are treated as a drift.

Cross sections (all READ from the downloaded data bases or published
tables; see electron_degradation_data/README.md):
  H   CCC data base (Bray et al. 2002, J. Phys. B 35, R117; files e-H_XSEC_LS):
      elastic 1S.1S (used as the momentum-transfer cross section, an upper
      bound), 1s -> nl for n = 2..10 l-resolved, total ionization TICS.1S.
      n > 10: n^-3 continuation of the n = 10 total (Bethe scaling of the
      oscillator strength, as Dalgarno et al. 1999 section 2.2 do for He).
  He  CCC data base (Fursa & Bray 1995, PRA 52, 1279 theory; files e-He_sig):
      elastic, 1^1S -> n^{1,3}L for n <= 5, TICS.s1S; from He(2^3S):
      every final state n <= 5 including the superelastic 2^3S -> 1^1S, and
      TICS.t2S. n > 5: n^-3 continuation of the n = 5 total of each spin.
  H2  MCCC data base (Scarlett et al. 2021, ADNDT 137, 101361; 2023 for the
      rotational files): X(v=0) -> every bound vibrational level of B, B',
      C, D, EF, GK, H, I, J, a, c, d, e, g, h, i, j; dissociative excitation
      of every state including b 3Sigma_u+; rotational N -> N +- 2 (v = 0,
      N = 0..9); total ionization TICS.
      Yoon et al. 2008 (JPCRD 37, 913): momentum transfer (Table 3),
      v = 0 -> 1 (Table 5), partial ionization H2+ / H+ (Table 13),
      dissociative attachment (Table 14).
      Janev, Reiter & Samm 2003 (Juelich 4105) eq. (74): v = 0 -> 2.
  Secondary-electron spectra below 100 eV primaries: the binary-encounter-
      dipole (BED) singly differential cross section, Kim & Rudd 1994 (PRA
      50, 3954) eq. (44), with the dipole oscillator strengths df/dw of
      their Table I for H and H2 and of Kim, Johnson & Rudd 2000 (PRA 61,
      034702) eq. (6) for He; only its shape in the secondary energy is
      used, the total is the CCC/MCCC one. He(2^3S): binary-encounter-
      Bethe form, Kim & Rudd 1994 eqs. (52)-(54) with Q = 1, U = B.
      The Opal, Peterson & Beaty (1971) form p ~ 1/(1+(eps/ebar)^2.1) is
      kept as an option (spectrum='opal') for the reference limits whose
      authors used it (Furlanetto & Stoever 2010; Dalgarno et al. 1999).
  Coulomb loss: Schunk & Hays (1971) as written by Furlanetto & Stoever
      (2010) eq. (5) with (12): dE/dx = -2 pi e^4 n_e ln Lambda / E,
      ln Lambda = ln(4E/zeta_e), zeta_e = 7.40e-11 (n_e/cm^-3) eV, times
      the thermal factor ((E - E_e)/(E - 0.53 E_e))^2.36, E_e = k T_e, of
      Swartz, Nisbet & Green (1971) as printed by Dalgarno et al. eq. (2).
      Option coulomb='swartz': Dalgarno et al. eq. (2) itself.

Usage as a module:
    gas = Gas(n_HI=..., n_HeI=..., n_He23S=..., n_H2=..., n_e=..., T_e=...)
    sol = solve(gas)                 # Solution on the grid
    rec = sol.record(25.0)           # dict of the deposition record
CLI: python3 electron_degradation_spencer_fano.py --nH 1 --nHe 0.1 --xe 1e-2
"""
import argparse
import math
import os
import re
import time
from dataclasses import dataclass, field

import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'electron_degradation_data')
A0SQ = 2.8002852e-17          # a0^2 in cm^2 (CODATA a0 = 5.29177e-9 cm)
E4 = (1.43996448e-7)**2       # e^4 in eV^2 cm^2 (e^2 = 1.43996e-7 eV cm)
KB_EV = 8.617333262e-5        # eV/K
ME_AMU = 5.48579909e-4        # electron mass in u
RYD = 13.605693               # eV
I_H = 13.598434               # H ionization energy, eV (NIST ASD)
I_HE = 24.587389              # He ionization energy, eV (NIST ASD)
E_HE23S = 19.819615           # He 2^3S, eV (NIST ASD, local csv)
I_H2 = 15.4254                # H2 adiabatic ionization energy, eV (NIST
                              # Chemistry WebBook; MCCC TICS header 15.43)
D0_H2 = 4.47713               # H2 dissociation energy from v = 0, J = 0, eV:
                              # the MCCC b 3Sigma_u+ DE threshold (header),
                              # 4.4781 eV in the spectroscopic literature
I_H2_DI = 18.076              # H2 -> H+ + H(1s) + e threshold = D0 + I_H
EA_H = 0.754                  # H- electron affinity, eV
H2_B_RADDISS = 0.27           # B 1Sigma_u+ radiative dissociation fraction
H2_C_RADDISS = 0.036          # C 1Pi_u   (Dalgarno et al. 1999 section 5.3,
                              # citing Stephens & Dalgarno 1972, Abgrall
                              # et al. 1992, 1997)
H2_B_FRAGMENT_KE = 5.4        # mean fragment kinetic energy of b 3Sigma_u+
                              # dissociation at the vertical transition
                              # (Dalgarno et al. 1999 section 5.3)

# He I levels (eV), statistically averaged over J, READ from NIST ASD
# (references/nist_asd/He_I_lines_500_21000A.csv); 5G not in the file,
# taken equal to 5F (hydrogenic for L >= 3).
HE_LEVELS = {
    't2S': 19.819615, 's2S': 20.615775, 't2P': 20.964094, 's2P': 21.218023,
    't3S': 22.718466, 's3S': 22.920317, 't3P': 23.007075, 't3D': 23.073652,
    's3D': 23.074075, 's3P': 23.087019, 't4S': 23.593959, 's4S': 23.673571,
    't4P': 23.707891, 't4D': 23.736091, 's4D': 23.736335, 't4F': 23.737007,
    's4F': 23.737010, 's4P': 23.742070, 't5S': 23.971971, 's5S': 24.011215,
    't5P': 24.028225, 't5D': 24.042662, 's5D': 24.042803, 't5F': 24.043154,
    's5F': 24.043156, 's5P': 24.045800, 't5G': 24.043154, 's5G': 24.043156,
}

# BED oscillator strengths df/dw = sum_k c_k y^k, y = B/(W + B):
# H and H2 Kim & Rudd (1994) Table I; He Kim, Johnson & Rudd (2000) eq. (6).
BED = {
    'H':  dict(B=13.6057, U=13.6057, N=1,
               c={2: -2.2473e-2, 3: 1.1775, 4: -4.6264e-1, 5: 8.9064e-2}),
    'He': dict(B=24.587, U=39.51, N=2,
               c={3: 8.24012, 4: -10.4769, 5: 3.96496, 6: -0.0445976}),
    'H2': dict(B=15.43, U=25.68, N=2,
               c={3: 1.1262, 4: 6.3982, 5: -7.8055, 6: 2.1440}),
}
OPAL_EBAR = {'H': 8.0, 'He': 15.8, 'H2': 8.3, 'He23S': 2.8}


# --------------------------------------------------------------------------
# data readers
# --------------------------------------------------------------------------
def _read_table(path):
    """(E, sigma[cm^2]) from a CCC or MCCC file; units from the header."""
    E, S = [], []
    unit = A0SQ
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                if 'cm^2' in line:
                    unit = 1.0
                continue
            p = line.split()
            if len(p) < 2:
                continue
            try:
                e, s = float(p[0]), float(p[1])
            except ValueError:
                continue
            E.append(e)
            S.append(max(s, 0.0))
    E = np.array(E)
    S = np.array(S) * unit
    o = np.argsort(E, kind='stable')
    return E[o], S[o]


def _header_value(path, key):
    with open(path) as f:
        for line in f:
            if line.startswith('#') and key in line:
                m = re.search(r'([-+0-9.Ee]+)\s*(eV|Hartrees)?\s*$',
                              line.strip())
                if m:
                    return float(m.group(1))
    return None


def _on_grid(E, S, grid, threshold=None, extrapolate='hold'):
    """sigma on the grid, linear in E, zero below the first point or the
    threshold; above the last point held (unused: data reach >= 100 eV)."""
    out = np.interp(grid, E, S, left=0.0, right=S[-1] if
                    extrapolate == 'hold' else 0.0)
    if threshold is not None:
        out[grid < threshold] = 0.0
    out[grid < E[0]] = 0.0
    return out


YOON_MT = np.array([  # Table 3, eV, 1e-16 cm^2
    (0.001, 7.25), (0.01, 7.7), (0.02, 8.14), (0.05, 9.27), (0.1, 10.47),
    (0.2, 12.02), (0.3, 13.0), (0.5, 14.52), (0.7, 15.66), (1, 17.01),
    (1.2, 17.7), (1.5, 18.15), (1.8, 18.22), (2, 18.11), (2.5, 17.4),
    (3, 16.28), (4, 13.7), (5, 11.59), (6, 10.0), (7, 8.59), (8, 7.48),
    (9, 6.58), (10, 5.78), (12, 4.396), (15, 3.275), (18, 2.529),
    (20, 2.154), (25, 1.476), (30, 1.077), (40, 0.636), (50, 0.417),
    (60, 0.311), (70, 0.243), (80, 0.2), (90, 0.169), (100, 0.149)])
YOON_V01 = np.array([  # Table 5, eV, 1e-16 cm^2 (no entry between 7 and 15)
    (0.516, 0.0), (0.55, 0.01), (0.6, 0.017), (0.66, 0.027), (0.74, 0.038),
    (0.8, 0.047), (0.85, 0.059), (1, 0.094), (1.25, 0.172), (1.48, 0.24),
    (2, 0.367), (2.46, 0.442), (3, 0.486), (3.8, 0.461), (4.6, 0.38),
    (5, 0.343), (5.3, 0.324), (5.6, 0.302), (6.1, 0.27), (7, 0.21),
    (15, 0.0565), (20, 0.0332), (25, 0.0242), (30, 0.017), (40, 0.0112),
    (50, 0.00968), (60, 0.00817), (80, 0.00609), (100, 0.00451)])
YOON_PARTIAL = np.array([  # Table 13: E, H2+ (1e-17), H+ (1e-18) cm^2
    (30, 6.42, 0.86), (35, 7.42, 1.76), (40, 8.12, 2.87), (45, 8.39, 4.08),
    (50, 8.59, 4.82), (55, 8.74, 5.72), (60, 8.82, 6.25), (65, 8.8, 6.82),
    (70, 8.79, 7.05), (80, 8.63, 7.39), (100, 8.24, 7.59)])
YOON_DA = np.array([  # Table 14, eV, 1e-16 cm^2
    (3.56, 9.70e-7), (3.63, 1.90e-6), (3.69, 8.09e-6), (3.75, 1.80e-5),
    (3.76, 2.07e-5), (3.81, 1.94e-5), (3.88, 1.29e-5), (3.94, 8.38e-6),
    (4.0, 6.21e-6), (4.06, 3.88e-6), (4.13, 2.80e-6), (4.25, 2.00e-6),
    (4.375, 9.70e-7), (7.5, 0.0), (7.67, 1.94e-5), (7.84, 2.64e-5),
    (8.08, 3.76e-5), (8.3, 4.90e-5), (8.59, 6.51e-5), (8.87, 8.15e-5),
    (9.17, 9.74e-5), (9.48, 1.10e-4), (9.78, 1.19e-4), (9.99, 1.22e-4),
    (10.2, 1.23e-4), (10.4, 1.23e-4), (10.6, 1.20e-4), (10.8, 1.16e-4),
    (11.0, 1.12e-4), (11.2, 1.08e-4), (11.4, 1.02e-4), (11.6, 9.29e-5),
    (11.8, 8.27e-5), (12.1, 6.80e-5), (12.4, 5.57e-5), (12.6, 5.17e-5),
    (12.9, 5.01e-5), (13.2, 6.92e-5), (13.3, 9.03e-5), (13.5, 1.48e-4),
    (13.7, 1.80e-4), (13.8, 1.95e-4), (13.9, 1.95e-4), (14.0, 1.89e-4),
    (14.1, 1.72e-4), (14.2, 1.55e-4), (14.3, 1.39e-4), (14.4, 1.26e-4),
    (14.5, 1.12e-4), (14.6, 1.02e-4), (14.7, 9.33e-5), (14.8, 8.88e-5),
    (14.9, 8.31e-5), (15.2, 8.03e-5), (15.4, 7.89e-5), (15.5, 7.85e-5),
    (15.6, 7.86e-5), (15.9, 7.97e-5)])


def janev_vib02(E):
    """sigma(v=0->2), Janev, Reiter & Samm (2003) eq. (74), xi = 0.628,
    Delta E = 1.003 eV (cm^2)."""
    dE, xi = 1.003, 0.628
    x = np.asarray(E, float) / dE
    out = np.zeros_like(x)
    m = x > 1.0
    out[m] = 5.78 * xi / dE**4 / x[m]**2 * (1 - 1 / x[m])**(6.11 / dE**0.5)
    return out * 1e-16


# --------------------------------------------------------------------------
# hydrogen case-B cascade: probability that an nl excitation ends in 2p
# --------------------------------------------------------------------------
_LYA_CACHE = {}


def hydrogen_radial_dipole(n1, l1, n2, l2):
    """<n1 l1 | r | n2 l2> in a0 for hydrogen (numerical, sympy R_nl)."""
    from sympy import symbols, lambdify
    from sympy.physics.hydrogen import R_nl
    from scipy.integrate import quad
    r = symbols('r', positive=True)
    f = lambdify(r, R_nl(n1, l1, r, 1) * R_nl(n2, l2, r, 1) * r**3, 'math')
    rmax = 8.0 * max(n1, n2)**2 + 40.0
    val = quad(f, 0, rmax, limit=400, epsabs=1e-12, epsrel=1e-10)[0]
    return val


def hydrogen_A(n, l, n2, l2):
    """Einstein A (s^-1), hydrogen nl -> n2 l2, infinite nuclear mass:
    A = (4/3) alpha^3 omega^3 max(l,l2)/(2l+1) |<r>|^2 / tau0 (a.u.)."""
    alpha, tau0 = 7.2973525693e-3, 2.4188843265857e-17
    w = 0.5 * (1.0 / n2**2 - 1.0 / n**2)
    R = hydrogen_radial_dipole(n, l, n2, l2)
    return 4.0 / 3.0 * alpha**3 * w**3 * max(l, l2) / (2 * l + 1) * R**2 / tau0


def lya_probability_case_b(nmax=10):
    """P(nl -> ... -> 2p) with every np -> 1s Lyman photon reabsorbed in
    place (case B), from hydrogenic Einstein coefficients; P(2s) = 0,
    P(2p) = 1. Returns dict {(n, l): P}."""
    if nmax in _LYA_CACHE:
        return _LYA_CACHE[nmax]
    P = {(2, 0): 0.0, (2, 1): 1.0}
    for n in range(3, nmax + 1):
        for l in range(n):
            num = den = 0.0
            for n2 in range(2, n):
                for l2 in (l - 1, l + 1):
                    if l2 < 0 or l2 >= n2:
                        continue
                    a = hydrogen_A(n, l, n2, l2)
                    num += a * P[(n2, l2)]
                    den += a
            P[(n, l)] = num / den if den > 0 else 0.0
    _LYA_CACHE[nmax] = P
    return P


# --------------------------------------------------------------------------
# gas and options
# --------------------------------------------------------------------------
@dataclass
class Gas:
    n_HI: float = 0.0
    n_HeI: float = 0.0          # He ground singlet 1^1S
    n_He23S: float = 0.0        # He metastable 2^3S
    n_H2: float = 0.0
    n_e: float = 0.0            # thermal electrons
    T_e: float = 1.0e4          # K
    T_rot: float = None         # H2 rotational temperature (default T_e)


@dataclass
class Options:
    dE: float = 0.02
    E_max: float = 60.0
    coulomb: str = 'schunk_hays'     # or 'swartz' (Dalgarno et al. eq. 2)
    E_e_swartz: float = None         # E_e of the Swartz factor (default kT_e)
    spectrum: str = 'bed'            # or 'opal'
    neutral_recoil: bool = True
    he_photon_recycling: bool = False   # FS10 section 4 assumption
    level_set: str = 'full'          # 'full' or 'fs10' (n <= 4 + 1.43 n=4)
    h2_rot_pop: dict = None          # {J: fraction}; default thermal T_rot
    he23s_targets: bool = True
    superelastic_iterations: int = 30
    floor: float = None              # eV; default max(1.5 k T_e, dE)
    coulomb_discrete_fraction: float = 0.0
    # > 0: the Coulomb drag is applied as discrete losses of this fraction
    # of the electron energy, at the rate that gives the same mean loss
    # (Habing & Goldsmith 1971; Shull 1979; Furlanetto & Stoever 2010
    # eq. 6, f = 0.05), wherever f E >= 2 dE; only to reproduce that
    # Monte Carlo convention, the physical treatment is the drift (0)


# record keys --------------------------------------------------------------
ENERGY_KEYS = [
    'heat_coulomb', 'heat_recoil', 'heat_thermalized', 'heat_fragments',
    'E_ion_HI', 'E_ion_HeI', 'E_ion_He23S', 'E_ion_H2', 'E_ion_H2_diss',
    'E_ion_HI_by_He_photons',
    'E_exc_H', 'E_exc_He_singlet', 'E_exc_He_triplet', 'E_He23S_out',
    'E_rot_H2', 'E_vib_H2', 'E_exc_H2_singlet', 'E_exc_H2_triplet_bound',
    'E_diss_H2_chem', 'E_exc_H_from_H2', 'E_attach_chem',
]
COUNT_KEYS = [
    'N_ion_HI', 'N_ion_HeI', 'N_ion_He23S', 'N_ion_H2', 'N_ion_H2_diss',
    'N_ion_HI_by_He_photons',
    'N_H_2s', 'N_H_2p', 'N_H_n3', 'N_H_n4', 'N_H_n5plus', 'N_Lya_caseB',
    'N_He_2_3S', 'N_He_2_3P', 'N_He_triplet_n3plus', 'N_He_2_1S',
    'N_He_2_1P', 'N_He_singlet_n3plus', 'N_He_singlet_P',
    'N_He23S_direct', 'N_He23S_superelastic', 'N_He23S_to_2_3P',
    'N_He23S_destroyed',
    'N_H2_vib1', 'N_H2_vib2', 'N_H2_B', 'N_H2_C', 'N_H2_singlet_u_other',
    'N_H2_singlet_g', 'N_H2_B_cascade', 'N_H2_C_cascade',
    'N_H2_triplet_bound', 'N_H2_b', 'N_H2_DE_excited', 'N_H2_attach',
    'N_H2_diss',
]
KEYS = ENERGY_KEYS + COUNT_KEYS
IDX = {k: i for i, k in enumerate(KEYS)}
NE = len(ENERGY_KEYS)
Q = len(KEYS)


def _vec(**kw):
    v = np.zeros(Q)
    for k, x in kw.items():
        v[IDX[k]] += x
    return v


@dataclass
class Channel:
    name: str
    sigma: np.ndarray            # on the grid, cm^2 (density applied later)
    density: float
    W: np.ndarray                # energy lost, on the grid (eV)
    fixed: np.ndarray            # deposit vector excluding the W remainder
    rem_key: str = None          # bucket of W - (energy in fixed) - spawn
    spawn: float = 0.0           # energy of an extra electron released
    absorbed: bool = False       # electron removed (attachment)


@dataclass
class IonChannel:
    name: str
    sigma: np.ndarray
    density: float
    I: float
    species: str                 # BED/Opal key
    deposit: np.ndarray


# --------------------------------------------------------------------------
# building the channels
# --------------------------------------------------------------------------
def _lsym(L):
    return 'SPDFGHIJKL'.index(L)


def build_channels(gas, opt, grid):
    ch, ion = [], []
    const = lambda v: np.full_like(grid, v, dtype=float)
    fs10 = (opt.level_set == 'fs10')
    nmax_H = 4 if fs10 else 10
    # ---------------- H ----------------
    if gas.n_HI > 0:
        dH = os.path.join(DATA, 'CCC_H')
        Plya = lya_probability_case_b(10)
        sig_n = {}
        for n in range(2, nmax_H + 1):
            nc = ':' if n == 10 else str(n)
            sig_n[n] = np.zeros_like(grid)
            for l in range(n):
                L = 'SPDFGHIJKL'[l]
                E, S = _read_table(os.path.join(dH, f'{nc}{L}.1S'))
                W = I_H * (1 - 1 / n**2)
                s = _on_grid(E, S, grid, W)
                sig_n[n] += s
                key = {2: 'N_H_2s' if l == 0 else 'N_H_2p', 3: 'N_H_n3',
                       4: 'N_H_n4'}.get(n, 'N_H_n5plus')
                if n == 2:
                    fixed = _vec(**{key: 1.0, 'N_Lya_caseB': float(l == 1)})
                else:
                    fixed = _vec(**{key: 1.0, 'N_Lya_caseB': Plya[(n, l)]})
                ch.append(Channel(f'H_{n}{L}', s, gas.n_HI, const(W), fixed,
                                  'E_exc_H'))
        # higher n
        if fs10:
            # FS10 section 4.1: sigma(n >= 5) = 1.43 sigma(n = 4), the n = 4
            # sublevel mixture and its Lya probabilities, at the n = 4 energy
            for l in range(4):
                L = 'SPDF'[l]
                E, S = _read_table(os.path.join(dH, f'4{L}.1S'))
                W = I_H * (1 - 1 / 16)
                s = 1.43 * _on_grid(E, S, grid, W)
                ch.append(Channel(f'H_n>4_{L}', s, gas.n_HI, const(W),
                                  _vec(N_H_n5plus=1.0,
                                       N_Lya_caseB=Plya[(4, l)]), 'E_exc_H'))
        else:
            n = np.arange(11, 2001)
            F = np.sum((10.0 / n)**3)
            Wm = np.sum(n**-3.0 * I_H * (1 - 1 / n**2)) / np.sum(n**-3.0)
            # Lya probability of the continuation: that of n = 10, weighted
            # by the n = 10 sublevel cross sections at 30 eV
            w10 = []
            for l in range(10):
                E, S = _read_table(os.path.join(dH, ':' + 'SPDFGHIJKL'[l]
                                                + '.1S'))
                w10.append(np.interp(30.0, E, S))
            p10 = sum(w * Plya[(10, l)] for l, w in enumerate(w10)) / sum(w10)
            s = F * sig_n[10]
            s[grid < Wm] = 0.0
            ch.append(Channel('H_n>10', s, gas.n_HI, const(Wm),
                              _vec(N_H_n5plus=1.0, N_Lya_caseB=p10),
                              'E_exc_H'))
        E, S = _read_table(os.path.join(dH, 'TICS.1S'))
        ion.append(IonChannel('H_ion', _on_grid(E, S, grid, I_H), gas.n_HI,
                              I_H, 'H', _vec(N_ion_HI=1.0, E_ion_HI=I_H)))
        if opt.neutral_recoil:
            E, S = _read_table(os.path.join(dH, '1S.1S'))
            ch.append(('recoil', _on_grid(E, S, grid), gas.n_HI, 1.00794))
    # ---------------- He ground ----------------
    if gas.n_HeI > 0:
        dHe = os.path.join(DATA, 'CCC_He')
        nmax = 4 if fs10 else 5
        n5 = {'s': np.zeros_like(grid), 't': np.zeros_like(grid)}
        for lev, W in HE_LEVELS.items():
            spin, n, L = lev[0], int(lev[1]), lev[2]
            if n > nmax:
                continue
            E, S = _read_table(os.path.join(dHe, f'{lev}.s1S'))
            s = _on_grid(E, S, grid, W)
            mult = 1.0
            if fs10 and n == 4:
                mult = 2.43          # n = 4 plus 1.43 n = 4 for n >= 5
            if n == nmax:
                n5[spin] += s
            ch.append(_he_exc_channel(lev, W, mult * s, gas, opt))
        if not fs10:
            n = np.arange(6, 2001)
            F = np.sum((5.0 / n)**3)
            Wm = I_HE - 13.6057 * np.sum(n**-5.0) / np.sum(n**-3.0)
            for spin in 'st':
                s = F * n5[spin]
                s[grid < Wm] = 0.0
                ch.append(_he_exc_channel(spin + '6+', Wm, s, gas, opt))
        E, S = _read_table(os.path.join(dHe, 'TICS.s1S'))
        ion.append(IonChannel('He_ion', _on_grid(E, S, grid, I_HE),
                              gas.n_HeI, I_HE, 'He',
                              _vec(N_ion_HeI=1.0, E_ion_HeI=I_HE)))
        if opt.neutral_recoil:
            E, S = _read_table(os.path.join(dHe, 's1S.s1S'))
            ch.append(('recoil', _on_grid(E, S, grid), gas.n_HeI, 4.002602))
    # ---------------- He 2^3S ----------------
    if gas.n_He23S > 0 and opt.he23s_targets:
        dHe = os.path.join(DATA, 'CCC_He')
        for lev in ['s1S'] + list(HE_LEVELS):
            if lev == 't2S':
                continue
            path = os.path.join(dHe, f'{lev}.t2S')
            if not os.path.exists(path):
                continue
            E, S = _read_table(path)
            if lev == 's1S':
                W = -E_HE23S
                s = _on_grid(E, S, grid)
                fixed = _vec(N_He23S_superelastic=1.0, N_He23S_destroyed=1.0)
                fixed[IDX['E_He23S_out']] += W
                ch.append(Channel('He23S->1S', s, gas.n_He23S, const(W),
                                  fixed, None))
                continue
            W = HE_LEVELS[lev] - E_HE23S
            s = _on_grid(E, S, grid, W)
            kw = {'N_He23S_to_2_3P': 1.0} if lev == 't2P' else \
                {'N_He23S_destroyed': 1.0}
            ch.append(Channel(f'He23S->{lev}', s, gas.n_He23S, const(W),
                              _vec(**kw), 'E_He23S_out'))
        E, S = _read_table(os.path.join(dHe, 'TICS.t2S'))
        I23 = I_HE - E_HE23S
        ion.append(IonChannel('He23S_ion', _on_grid(E, S, grid, I23),
                              gas.n_He23S, I23, 'He23S',
                              _vec(N_ion_He23S=1.0, E_ion_He23S=I23,
                                   N_He23S_destroyed=1.0)))
    # ---------------- H2 ----------------
    if gas.n_H2 > 0:
        dX = os.path.join(DATA, 'MCCC_H2', 'MCCC-el-H2-X1Sg')
        dexc = os.path.join(dX, 'excitation', 'vi=0')
        su = {'B1Su': 'N_H2_B', 'C1Pu': 'N_H2_C', 'Bp1Su':
              'N_H2_singlet_u_other', 'D1Pu': 'N_H2_singlet_u_other'}
        sg = ['EF1Sg', 'GK1Sg', 'H1Sg', 'I1Pg', 'J1Dg']
        Bcasc = {'B1Su', 'EF1Sg', 'GK1Sg', 'H1Sg'}
        Ccasc = {'C1Pu', 'I1Pg', 'J1Dg'}
        for fn in sorted(os.listdir(dexc)):
            m = re.match(r'MCCC-el-H2-(\w+?)_(vf=\d+|DE)\.X1Sg_vi=0\.txt', fn)
            if not m or m.group(1) == 'X1Sg':
                continue
            st, kind = m.group(1), m.group(2)
            path = os.path.join(dexc, fn)
            E, S = _read_table(path)
            thr = _header_value(path, 'Threshold')
            triplet = st[1] == '3'
            if kind == 'DE':
                if st == 'b3Su':
                    s = _on_grid(E, S, grid, thr)
                    W = np.minimum(D0_H2 + H2_B_FRAGMENT_KE,
                                   thr + 0.5 * (grid - thr))
                    W = np.maximum(W, thr)
                    fixed = _vec(N_H2_b=1.0, N_H2_diss=1.0,
                                 E_diss_H2_chem=D0_H2)
                    ch.append(Channel('H2_b_DE', s, gas.n_H2, W, fixed,
                                      'heat_fragments'))
                else:
                    s = _on_grid(E, S, grid, thr)
                    # products H(1s) + H(n = 2 or 3); the excited atom's
                    # energy (threshold - D0) is recorded as H excitation
                    eH = thr - D0_H2
                    fixed = _vec(N_H2_DE_excited=1.0, N_H2_diss=1.0,
                                 E_diss_H2_chem=D0_H2, E_exc_H_from_H2=eH)
                    ch.append(Channel(f'H2_{st}_DE', s, gas.n_H2, const(thr),
                                      fixed, 'heat_fragments'))
                continue
            s = _on_grid(E, S, grid, thr)
            if triplet:
                # bound triplets decay radiatively into b 3Sigma_u+ and
                # dissociate (Dalgarno et al. section 5.3): D0 chemical,
                # the rest continuum photon plus fragment kinetic energy
                fixed = _vec(N_H2_triplet_bound=1.0, N_H2_diss=1.0,
                             E_diss_H2_chem=D0_H2)
                ch.append(Channel(f'H2_{st}_{kind}', s, gas.n_H2, const(thr),
                                  fixed, 'E_exc_H2_triplet_bound'))
            else:
                fd = H2_B_RADDISS if st in Bcasc else \
                    (H2_C_RADDISS if st in Ccasc else 0.0)
                kw = {'N_H2_diss': fd}
                if st in su:
                    kw[su[st]] = 1.0
                else:
                    kw['N_H2_singlet_g'] = 1.0
                if st in Bcasc:
                    kw['N_H2_B_cascade'] = 1.0
                if st in Ccasc:
                    kw['N_H2_C_cascade'] = 1.0
                ch.append(Channel(f'H2_{st}_{kind}', s, gas.n_H2, const(thr),
                                  _vec(**kw), 'E_exc_H2_singlet'))
        # vibrational
        s = _on_grid(YOON_V01[:, 0], YOON_V01[:, 1] * 1e-16, grid, 0.516)
        m = (grid > 7.0) & (grid < 15.0)          # log interpolation gap
        s[m] = np.exp(np.interp(np.log(grid[m]), np.log([7.0, 15.0]),
                                np.log([0.21e-16, 0.0565e-16])))
        ch.append(Channel('H2_v1', s, gas.n_H2, const(0.516),
                          _vec(N_H2_vib1=1.0), 'E_vib_H2'))
        ch.append(Channel('H2_v2', janev_vib02(grid), gas.n_H2, const(1.003),
                          _vec(N_H2_vib2=1.0), 'E_vib_H2'))
        # dissociative attachment: electron removed
        # the table starts at 3.56 eV, below the v = 0 energetic threshold
        # D0 - EA = 3.72 eV (thermal spread of the measurement); the
        # electron cannot attach below it here
        s = _on_grid(YOON_DA[:, 0], YOON_DA[:, 1] * 1e-16, grid,
                     D0_H2 - EA_H, extrapolate='zero')
        ch.append(Channel('H2_DA', s, gas.n_H2, grid.copy(),
                          _vec(N_H2_attach=1.0, N_H2_diss=1.0,
                               E_attach_chem=D0_H2 - EA_H),
                          'heat_fragments', absorbed=True))
        # ionization, split into H2+ and H+ + H with Yoon Table 13
        E, S = _read_table(os.path.join(dX, 'ionisation', 'vi=0',
                                        'MCCC-el-H2-TICS.X1Sg_vi=0.txt'))
        stot = _on_grid(E, S, grid, I_H2)
        r = YOON_PARTIAL[:, 2] * 1e-18 / (YOON_PARTIAL[:, 1] * 1e-17 +
                                          YOON_PARTIAL[:, 2] * 1e-18)
        fdi = np.interp(grid, np.r_[I_H2_DI, YOON_PARTIAL[:, 0]],
                        np.r_[0.0, r])
        fdi[grid < I_H2_DI] = 0.0
        ion.append(IonChannel('H2_ion', stot * (1 - fdi), gas.n_H2, I_H2,
                              'H2', _vec(N_ion_H2=1.0, E_ion_H2=I_H2)))
        ion.append(IonChannel('H2_DI', stot * fdi, gas.n_H2, I_H2_DI, 'H2',
                              _vec(N_ion_H2_diss=1.0, N_H2_diss=1.0,
                                   E_ion_H2_diss=I_H2_DI)))
        if opt.neutral_recoil:
            s = np.interp(grid, YOON_MT[:, 0], YOON_MT[:, 1] * 1e-16)
            ch.append(('recoil', s, gas.n_H2, 2.01588))
    recoil = [c for c in ch if isinstance(c, tuple)]
    ch = [c for c in ch if not isinstance(c, tuple)]
    return ch, ion, recoil


def _he_exc_channel(lev, W, s, gas, opt):
    spin, n = lev[0], lev[1]
    triplet = spin == 't'
    if lev == 't2S':
        kw = {'N_He_2_3S': 1.0}
    elif lev == 't2P':
        kw = {'N_He_2_3P': 1.0}
    elif lev == 's2S':
        kw = {'N_He_2_1S': 1.0}
    elif lev == 's2P':
        kw = {'N_He_2_1P': 1.0, 'N_He_singlet_P': 1.0}
    elif triplet:
        kw = {'N_He_triplet_n3plus': 1.0}
    else:
        kw = {'N_He_singlet_n3plus': 1.0}
        if len(lev) == 3 and lev[2] == 'P':
            kw['N_He_singlet_P'] = 1.0
    if triplet:
        kw['N_He23S_direct'] = 1.0       # every triplet cascades to 2^3S
    key = 'E_exc_He_triplet' if triplet else 'E_exc_He_singlet'
    if opt.he_photon_recycling:
        # FS10 section 4: the excited atom decays at once to the ground
        # state and its photon (W < 24.6 eV) ionizes H I; the photoelectron
        # of energy W - I_H continues
        fixed = _vec(**kw)
        fixed += _vec(N_ion_HI_by_He_photons=1.0,
                      E_ion_HI_by_He_photons=I_H)
        return Channel('He_' + lev, s, gas.n_HeI,
                       np.full_like(s, W), fixed, key, spawn=W - I_H)
    return Channel('He_' + lev, s, gas.n_HeI, np.full_like(s, W),
                   _vec(**kw), key)


def rotational_loss(gas, opt, grid):
    """Net H2 rotational energy loss per unit length (eV/cm), MCCC N -> N+-2
    (v = 0) cross sections weighted by the rotational populations."""
    d = os.path.join(DATA, 'MCCC_H2', 'MCCC-el-H2-X1Sg',
                     'rovibrational_excitation', 'vi=0')
    T = gas.T_rot if gas.T_rot else gas.T_e
    lev = {}
    files = {}
    for N in range(10):
        for fn in os.listdir(os.path.join(d, f'Ni={N}')):
            m = re.search(r'Nf=(\d+)\.X1Sg_vi=0_Ni=(\d+)', fn)
            Nf = int(m.group(1))
            if abs(Nf - N) != 2:
                continue
            p = os.path.join(d, f'Ni={N}', fn)
            ei = _header_value(p, 'Initial-state energy')
            ef = _header_value(p, 'Final-state energy')
            lev[N] = ei
            lev[Nf] = ef
            files[(N, Nf)] = p
    if opt.h2_rot_pop is not None:
        pop = dict(opt.h2_rot_pop)
    else:
        e0 = lev[0]
        w = {N: (2 * N + 1) * (3 if N % 2 else 1) *
             math.exp(-(lev[N] - e0) * 27.211386 / (KB_EV * T))
             for N in range(10) if N in lev}
        z = sum(w.values())
        pop = {N: v / z for N, v in w.items()}
    L = np.zeros_like(grid)
    for (N, Nf), p in files.items():
        if pop.get(N, 0.0) <= 0:
            continue
        E, S = _read_table(p)
        dEr = (lev[Nf] - lev[N]) * 27.211386
        s = _on_grid(E, S, grid) if dEr > 0 else np.interp(grid, E, S)
        L += gas.n_H2 * pop[N] * s * dEr
    return np.maximum(L, 0.0)


def coulomb_loss(gas, opt, grid):
    """Coulomb drag on thermal electrons, eV/cm."""
    if gas.n_e <= 0:
        return np.zeros_like(grid)
    Ee = opt.E_e_swartz if opt.E_e_swartz is not None else KB_EV * gas.T_e
    with np.errstate(invalid='ignore', divide='ignore'):
        th = np.where(grid > Ee, ((grid - Ee) / (grid - 0.53 * Ee))**2.36,
                      0.0)
        if opt.coulomb == 'swartz':
            v = 5.9310e7 * np.sqrt(grid)
            L = 2.0e-4 * gas.n_e**0.97 / grid**0.44 * th / v
        elif opt.coulomb == 'spitzer':
            # FS10 eq. (11) (Spitzer 1962), binary collisions only,
            # Lambda = [(kT)^3/(pi n_e e^6 gamma_e)]^(1/2), ln gamma_e = 0.5772
            kT = KB_EV * gas.T_e
            lnL = 0.5 * np.log(kT**3 / (np.pi * gas.n_e * E4 * 1.43996448e-7
                                        * math.exp(0.5772157)))
            L = 2 * np.pi * E4 * gas.n_e * lnL / grid * th
        else:
            zeta = 7.40e-11 * gas.n_e
            lnL = np.log(np.maximum(4 * grid / zeta, 1.0))
            L = 2 * np.pi * E4 * gas.n_e * lnL / grid * th
    L[~np.isfinite(L)] = 0.0
    return L


def secondary_weights(species, E, I, dE, spectrum):
    """Normalized weights of the secondary energy eps_k = k dE, k = 0..K,
    eps <= (E - I)/2."""
    A = E - I
    K = int(math.floor(0.5 * A / dE + 1e-9))
    eps = np.arange(K + 1) * dE
    if spectrum == 'opal' or species == 'He23S':
        if species == 'He23S' and spectrum != 'opal':
            # BEB, Kim & Rudd (1994) eqs. (52)-(54), Q = 1, U = B
            B = I
            t, u, w = E / B, 1.0, eps / B
            F2 = 1.0 / (t + u + 1)
            F1 = -F2 / (t + 1)
            F3 = math.log(max(t, 1.0)) / (t + u + 1)
            p = sum(Fn * ((w + 1)**-n + (t - w)**-n)
                    for n, Fn in ((1, F1), (2, F2), (3, F3)))
        else:
            p = 1.0 / (1.0 + (eps / OPAL_EBAR[species])**2.1)
    else:
        b = BED[species]
        B, U, N = b['B'], b['U'], b['N']
        t, u, w = E / B, U / B, eps / B
        y = 1.0 / (w + 1.0)
        dfdw = sum(c * y**k for k, c in b['c'].items())
        Ni = sum(c / (k - 1) for k, c in b['c'].items())
        lnt = math.log(max(t, 1.0))
        p = ((Ni / N - 2) / (t + 1) * (1 / (w + 1) + 1 / (t - w))
             + (2 - Ni / N) * (1 / (w + 1)**2 + 1 / (t - w)**2)
             + lnt / (N * (w + 1)) * dfdw)
    p = np.maximum(p, 0.0)
    s = p.sum()
    if s <= 0:
        p = np.zeros(K + 1)
        p[0] = 1.0
        return p
    return p / s


# --------------------------------------------------------------------------
# solver
# --------------------------------------------------------------------------
class Solution:
    def __init__(self, grid, Y, gas, opt, timing):
        self.E, self.Y, self.gas, self.opt = grid, Y, gas, opt
        self.timing = timing

    def record(self, E0):
        """Deposition record at primary energy E0 (linear interpolation)."""
        y = np.array([np.interp(E0, self.E, self.Y[:, q]) for q in range(Q)])
        rec = {k: y[IDX[k]] for k in KEYS}
        rec['E0'] = E0
        rec['closure'] = y[:NE].sum() / E0 - 1.0
        return rec

    def fractions(self, E0):
        r = self.record(E0)
        heat = r['heat_coulomb'] + r['heat_recoil'] + r['heat_thermalized']
        return dict(
            E0=E0, f_heat_electrons=(r['heat_coulomb'] + r['heat_thermalized'])
            / E0, f_heat_heavy=(r['heat_recoil'] + r['heat_fragments']) / E0,
            f_heat=(heat + r['heat_fragments']) / E0,
            f_ion=(r['E_ion_HI'] + r['E_ion_HeI'] + r['E_ion_He23S'] +
                   r['E_ion_H2'] + r['E_ion_H2_diss'] +
                   r['E_ion_HI_by_He_photons']) / E0,
            f_exc_H=r['E_exc_H'] / E0, closure=r['closure'])


def solve(gas, opt=None, verbose=False):
    opt = opt or Options()
    t0 = time.time()
    dE = opt.dE
    N = int(round(opt.E_max / dE))
    grid = np.arange(N + 1) * dE
    ch, ion, recoil = build_channels(gas, opt, grid)
    t_build = time.time() - t0
    # continuous losses, eV/cm
    Lc = coulomb_loss(gas, opt, grid)
    Lr = np.zeros_like(grid)
    for (_, s, n, M) in recoil:
        Lr += n * s * (2 * ME_AMU / M) * grid
    Lrot = rotational_loss(gas, opt, grid) if gas.n_H2 > 0 else 0 * grid
    fC = opt.coulomb_discrete_fraction
    if fC > 0:
        disc = fC * grid >= 2 * dE
        sigC = np.where(disc, Lc / np.maximum(fC * grid, 1e-300), 0.0)
        ch.append(Channel('coulomb_discrete', sigC, 1.0, fC * grid,
                          np.zeros(Q), 'heat_coulomb'))
        Lc = np.where(disc, 0.0, Lc)
    Lt = Lc + Lr + Lrot
    # discrete rates
    nu = np.array([c.density * c.sigma for c in ch]) if ch else \
        np.zeros((0, N + 1))
    Wm = np.array([c.W for c in ch]) if ch else np.zeros((0, N + 1))
    nui = np.array([c.density * c.sigma for c in ion]) if ion else \
        np.zeros((0, N + 1))
    # fixed deposit including the W remainder (energy closure per event)
    fixed = np.array([c.fixed for c in ch]) if ch else np.zeros((0, Q))
    efix = fixed[:, :NE].sum(axis=1) if ch else np.zeros(0)
    remk = np.array([IDX[c.rem_key] if c.rem_key else -1 for c in ch], int)
    spawn = np.array([c.spawn for c in ch])
    absorbed = np.array([c.absorbed for c in ch], bool)
    if ch:
        pos = Wm[(nu > 0) & (Wm > 0)]
        assert pos.size == 0 or pos.min() >= dE, \
            'a discrete loss smaller than the bin width: reduce dE'
    floor = opt.floor if opt.floor is not None else \
        max(1.5 * KB_EV * gas.T_e, dE)
    i_floor = int(math.ceil(floor / dE - 1e-9))
    Y = np.zeros((N + 1, Q))
    Y[:, IDX['heat_thermalized']] = grid      # thermal below the floor;
    has_super = bool(ch) and np.any(Wm < 0)
    n_iter = opt.superelastic_iterations if has_super else 1
    # precompute ionization secondary weights lazily per bin
    sw_cache = {}
    for it in range(n_iter):
        Yold = Y.copy()
        for i in range(i_floor, N + 1):
            Ei = grid[i]
            rate_cont = Lt[i] / dE
            nu_i = nu[:, i] if ch else nu[:, 0]
            act = np.nonzero(nu_i > 0)[0]
            nui_i = nui[:, i]
            R = nu_i[act].sum() + nui_i.sum() + rate_cont
            if R <= 0:
                Y[i] = 0.0
                Y[i, IDX['heat_thermalized']] = Ei
                continue
            acc = np.zeros(Q)
            if rate_cont > 0:
                cont = Y[i - 1].copy()
                cont[IDX['heat_coulomb']] += dE * Lc[i] / Lt[i]
                cont[IDX['heat_recoil']] += dE * Lr[i] / Lt[i]
                cont[IDX['E_rot_H2']] += dE * Lrot[i] / Lt[i]
                acc += rate_cont * cont
            if act.size:
                W = Wm[act, i]
                dep = fixed[act].copy()
                rk = remk[act]
                rem = W - efix[act] - spawn[act]
                m = rk >= 0
                np.add.at(dep, (np.nonzero(m)[0], rk[m]), rem[m])
                # continuation of the scattered electron
                tgt = Ei - W
                cont = np.zeros((act.size, Q))
                ok = ~absorbed[act]
                x = tgt[ok] / dE
                j0 = np.floor(x + 1e-12).astype(int)
                f = x - j0
                rows_lo = _rows(Y, Yold, j0, i, N, grid)
                rows_hi = _rows(Y, Yold, j0 + 1, i, N, grid)
                cont[ok] = (1 - f)[:, None] * rows_lo + f[:, None] * rows_hi
                sp = spawn[act]
                if np.any(sp > 0):
                    ms = sp > 0
                    xs = sp[ms] / dE
                    k0 = np.floor(xs).astype(int)
                    fs = xs - k0
                    cont[ms] += (1 - fs)[:, None] * Y[k0] + \
                        fs[:, None] * Y[k0 + 1]
                acc += nu_i[act] @ (dep + cont)
            for t, c in enumerate(ion):
                if nui_i[t] <= 0:
                    continue
                key = (t, i)
                if key not in sw_cache:
                    sw_cache[key] = secondary_weights(c.species, Ei, c.I, dE,
                                                      opt.spectrum)
                w = sw_cache[key]
                K = w.size - 1
                a = (Ei - c.I) / dE
                mm = int(math.floor(a + 1e-9))
                fr = a - mm
                sec = w @ Y[0:K + 1]
                lo = w @ Y[mm - K:mm + 1][::-1]
                hi = w @ Y[mm - K + 1:mm + 2][::-1] if fr > 1e-12 else lo
                prim = (1 - fr) * lo + fr * hi
                acc += nui_i[t] * (c.deposit + sec + prim)
            Y[i] = acc / R
        if not has_super:
            break
        if np.max(np.abs(Y - Yold)) < 1e-12 * opt.E_max:
            break
    sol = Solution(grid, Y, gas, opt, dict(build=t_build,
                                           total=time.time() - t0,
                                           iterations=it + 1,
                                           n_channels=len(ch),
                                           n_ion=len(ion)))
    if verbose:
        print(sol.timing)
    return sol


def _rows(Y, Yold, j, i, N, grid):
    """Y rows at indices j: the current sweep below i, the previous
    iterate at or above i (superelastic landings), proportional
    extrapolation above the grid top."""
    out = np.empty((j.size, Y.shape[1]))
    below = j < i
    out[below] = Y[np.maximum(j[below], 0)]
    mid = (~below) & (j <= N)
    out[mid] = Yold[j[mid]]
    top = j > N
    if np.any(top):
        out[top] = Yold[N][None, :] * (j[top] / N)[:, None]
    return out


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--nH', type=float, default=1.0)
    ap.add_argument('--nHe', type=float, default=0.1)
    ap.add_argument('--nHe23S', type=float, default=0.0)
    ap.add_argument('--nH2', type=float, default=0.0)
    ap.add_argument('--xe', type=float, default=1e-2,
                    help='n_e / (n_H nuclei + n_He)')
    ap.add_argument('--ntot', type=float, default=1e10,
                    help='scale of all densities (cm^-3)')
    ap.add_argument('--Te', type=float, default=5000.0)
    ap.add_argument('--dE', type=float, default=0.02)
    ap.add_argument('--Emax', type=float, default=60.0)
    ap.add_argument('--E0', type=float, nargs='*',
                    default=[12, 15, 20, 25, 30, 40, 50])
    a = ap.parse_args()
    s = a.ntot / (a.nH + 2 * a.nH2 + a.nHe)
    nuc = (a.nH + 2 * a.nH2 + a.nHe) * s
    gas = Gas(n_HI=a.nH * s, n_HeI=a.nHe * s, n_He23S=a.nHe23S * s,
              n_H2=a.nH2 * s, n_e=a.xe * nuc, T_e=a.Te)
    sol = solve(gas, Options(dE=a.dE, E_max=a.Emax))
    print(f'# solve time {sol.timing["total"]:.2f} s, '
          f'{sol.timing["n_channels"]} discrete channels')
    for E0 in a.E0:
        r = sol.fractions(E0)
        print(' '.join(f'{k}={v:.4g}' for k, v in r.items()))


if __name__ == '__main__':
    main()
