"""He I effective collision strengths of Bray, Burgess, Fursa & Tully (2000).

Source: I. Bray, A. Burgess, D. V. Fursa & J. A. Tully (2000), A&AS 146,
481-498, Tables 2 (1s2 1S -> 1s nl 1,3L) and 3 (1s2s 3S -> 1s nl 1,3L),
n <= 5, READ from the published article (the text layer of the journal PDF,
checked against the printed page).  Convergent close-coupling (CCC) cross
sections, thermally averaged.  What the paper says of the range (p. 483):
"In Tables 2 - 5 we give effective collision strengths Upsilon over the
temperature range 3.75 <= log T <= 5.75.  Because of the uncertainty that
pseudo resonances introduce into our cross sections at energies close to
threshold, we limit the low temperature end of our tabulation to about
6 000 degrees."  And (Sect. 6): "For our tabulated range, we have
log T >= 3.75, where the error due to this threshold effect does not exceed
a few percent."

ONE MISPRINT is corrected: Table 2, 1^1S - 4^1P at log T = 5.25 is printed
1.713^-3 between 9.287^-3 (log T 5.00) and 3.040^-2 (5.50); the n = 3 and
n = 5 neighbors of the same row (4.340^-2, 8.605^-3) and the ratio 4^1P /
3^1P = 0.39-0.40 of every other row make it 1.713^-2, which is used.

THE ADOPTED STRENGTH OF ONE TERM-TO-TERM TRANSITION (adopted_upsilon), the
one construction used for every He I electron-impact excitation from 1^1S
and from 2^3S that the code carries:

  inside the table, 10^3.75 <= T <= 10^5.75 K: the Bray et al. values,
    interpolated by a cubic Hermite in (log10 T, log10 Upsilon) whose
    interior node slopes are the Fritsch-Carlson (PCHIP) ones, the harmonic
    mean of the two secants (zero at a local extremum), and whose two end
    node slopes are those of the CHIANTI shape below, so that the strength
    and its slope are continuous at the table ends (the end intervals can
    then overshoot the tabulated values slightly; overshoot_report gives
    the size);
  outside the table: the temperature dependence of the CHIANTI v11.0.2
    he_1 strength of the same term (the J levels summed; Burgess & Tully
    scaled knots descaled with the natural spline, chianti_cooling.upsilon)
    anchored to the Bray value at the table end, so value and slope are
    continuous there.  CHIANTI states it uses its he_1 strengths (Sawey &
    Berrington 1993, ADNDT 55, 81; Bray et al. 2000) over
    3.0 < log T < 5.75; below log T 3.75 only the Sawey & Berrington
    R-matrix data stand behind the shape, and that paper was not read.
    Where CHIANTI has no level of the term (2^3S -> 4^3S) the strength is
    held at the table end value.

Energies: the g-weighted observed level energies of CHIANTI he_1.elvlc (NIST
ASD).  The module is imported by fit_helium_triplet_cascade.py and
fit_hydrogen_helium_excitation.py; run as a script it prints the checks.
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chianti_cooling as cc

HC_K = 1.438776877            # h c / k [cm K]
BRAY_LOGT = 3.75 + 0.25 * np.arange(9)
T_LO, T_HI = 10.0 ** BRAY_LOGT[0], 10.0 ** BRAY_LOGT[-1]


def _block(text):
    return np.array([[float(v) for v in line.split()[1:]]
                     for line in text.strip().splitlines()])


# Table 2, 1^1S -> n^(2S+1)L.  Each block: log T, then the n = 2..5 columns of
# the singlet transition, then those of the triplet one, as printed.
_T2 = {
    "S": ("""
3.75 3.075e-2 8.431e-3 2.397e-3 1.000e-3 6.198e-2 1.642e-2 5.317e-3 2.593e-3
4.00 3.492e-2 8.308e-3 2.505e-3 1.084e-3 6.458e-2 1.584e-2 5.233e-3 2.542e-3
4.25 3.840e-2 8.364e-3 2.642e-3 1.189e-3 6.387e-2 1.530e-2 5.183e-3 2.477e-3
4.50 4.183e-2 8.696e-3 2.856e-3 1.334e-3 6.157e-2 1.486e-2 5.179e-3 2.448e-3
4.75 4.573e-2 9.313e-3 3.178e-3 1.526e-3 5.832e-2 1.434e-2 5.142e-3 2.421e-3
5.00 5.048e-2 1.022e-2 3.628e-3 1.771e-3 5.320e-2 1.346e-2 4.936e-3 2.325e-3
5.25 5.649e-2 1.147e-2 4.210e-3 2.068e-3 4.787e-2 1.198e-2 4.461e-3 2.108e-3
5.50 6.436e-2 1.318e-2 4.941e-3 2.426e-3 4.018e-2 9.968e-3 3.738e-3 1.774e-3
5.75 7.481e-2 1.553e-2 5.888e-3 2.868e-3 3.167e-2 7.741e-3 2.904e-3 1.384e-3
""", [2, 3, 4, 5], [2, 3, 4, 5]),
    "P": ("""
3.75 8.886e-3 2.783e-3 1.046e-3 6.597e-4 1.716e-2 6.041e-3 2.074e-3 1.193e-3
4.00 1.299e-2 3.678e-3 1.504e-3 9.203e-4 2.233e-2 6.906e-3 2.468e-3 1.466e-3
4.25 1.910e-2 5.128e-3 2.143e-3 1.223e-3 2.826e-2 7.931e-3 2.983e-3 1.762e-3
4.50 3.006e-2 7.760e-3 3.173e-3 1.701e-3 3.477e-2 9.178e-3 3.585e-3 2.070e-3
4.75 5.178e-2 1.295e-2 5.180e-3 2.673e-3 4.128e-2 1.059e-2 4.236e-3 2.378e-3
5.00 9.534e-2 2.343e-2 9.287e-3 4.706e-3 4.635e-2 1.185e-2 4.796e-3 2.625e-3
5.25 1.778e-1 4.340e-2 1.713e-2 8.605e-3 4.800e-2 1.240e-2 5.048e-3 2.704e-3
5.50 3.167e-1 7.716e-2 3.040e-2 1.521e-2 4.503e-2 1.184e-2 4.832e-3 2.545e-3
5.75 5.217e-1 1.274e-1 5.023e-2 2.515e-2 3.808e-2 1.020e-2 4.169e-3 2.166e-3
""", [2, 3, 4, 5], [2, 3, 4, 5]),
    "D": ("""
3.75 3.884e-3 1.416e-3 7.906e-4 1.863e-3 1.030e-3 8.023e-4
4.00 3.809e-3 1.442e-3 8.285e-4 2.060e-3 1.125e-3 8.630e-4
4.25 3.783e-3 1.548e-3 8.888e-4 2.202e-3 1.247e-3 9.306e-4
4.50 3.942e-3 1.764e-3 1.009e-3 2.322e-3 1.354e-3 9.867e-4
4.75 4.390e-3 2.152e-3 1.230e-3 2.396e-3 1.407e-3 9.957e-4
5.00 5.266e-3 2.790e-3 1.585e-3 2.372e-3 1.389e-3 9.480e-4
5.25 6.598e-3 3.642e-3 2.047e-3 2.195e-3 1.283e-3 8.485e-4
5.50 8.159e-3 4.521e-3 2.522e-3 1.860e-3 1.087e-3 7.044e-4
5.75 9.556e-3 5.222e-3 2.900e-3 1.434e-3 8.368e-4 5.371e-4
""", [3, 4, 5], [3, 4, 5]),
    "F": ("""
3.75 5.864e-4 4.221e-4 5.082e-4 5.039e-4
4.00 4.917e-4 3.813e-4 4.547e-4 4.501e-4
4.25 4.101e-4 3.488e-4 3.748e-4 3.664e-4
4.50 3.467e-4 3.149e-4 2.945e-4 2.801e-4
4.75 3.011e-4 2.797e-4 2.268e-4 2.078e-4
5.00 2.697e-4 2.472e-4 1.724e-4 1.525e-4
5.25 2.436e-4 2.179e-4 1.279e-4 1.111e-4
5.50 2.135e-4 1.881e-4 9.081e-5 7.878e-5
5.75 1.781e-4 1.564e-4 6.114e-5 5.335e-5
""", [4, 5], [4, 5]),
    "G": ("""
3.75 1.938e-4 3.037e-4
4.00 1.551e-4 2.389e-4
4.25 1.158e-4 1.717e-4
4.50 8.143e-5 1.157e-4
4.75 5.484e-5 7.485e-5
5.00 3.640e-5 4.734e-5
5.25 2.473e-5 2.968e-5
5.50 1.742e-5 1.860e-5
5.75 1.244e-5 1.197e-5
""", [5], [5]),
}

# Table 3, 2^3S -> n^(2S+1)L, same layout (no 2^3S - 2^3S column).
_T3 = {
    "S": ("""
3.75 2.389 3.544e-1 1.199e-1 5.624e-2 2.410 6.651e-1 2.321e-1
4.00 2.456 3.295e-1 1.062e-1 4.930e-2 2.286 5.911e-1 2.154e-1
4.25 2.275 2.832e-1 8.877e-2 4.129e-2 2.235 5.383e-1 2.001e-1
4.50 1.916 2.280e-1 7.059e-2 3.305e-2 2.370 5.328e-1 2.009e-1
4.75 1.496 1.730e-1 5.351e-2 2.517e-2 2.761 5.919e-1 2.265e-1
5.00 1.111 1.247e-1 3.870e-2 1.823e-2 3.397 7.092e-1 2.732e-1
5.25 8.003e-1 8.624e-2 2.682e-2 1.261e-2 4.187 8.591e-1 3.302e-1
5.50 5.660e-1 5.765e-2 1.793e-2 8.402e-3 5.013 1.016 3.889e-1
5.75 3.944e-1 3.747e-2 1.165e-2 5.419e-3 5.755 1.163 4.443e-1
""", [2, 3, 4, 5], [3, 4, 5]),
    "P": ("""
3.75 7.965e-1 1.306e-1 4.941e-2 3.272e-2 1.508e+1 1.606 4.476e-1 1.833e-1
4.00 9.579e-1 1.499e-1 5.791e-2 3.599e-2 2.580e+1 1.611 4.465e-1 1.936e-1
4.25 1.042 1.543e-1 6.089e-2 3.587e-2 4.185e+1 1.576 4.476e-1 2.015e-1
4.50 1.015 1.439e-1 5.740e-2 3.259e-2 6.455e+1 1.552 4.578e-1 2.105e-1
4.75 8.950e-1 1.234e-1 4.930e-2 2.730e-2 9.526e+1 1.615 4.938e-1 2.291e-1
5.00 7.265e-1 9.896e-2 3.939e-2 2.141e-2 1.351e+2 1.868 5.842e-1 2.707e-1
5.25 5.516e-1 7.494e-2 2.964e-2 1.586e-2 1.842e+2 2.410 7.603e-1 3.493e-1
5.50 3.948e-1 5.367e-2 2.106e-2 1.114e-2 2.401e+2 3.274 1.037 4.719e-1
5.75 2.677e-1 3.639e-2 1.419e-2 7.434e-3 2.976e+2 4.382 1.408 6.315e-1
""", [2, 3, 4, 5], [2, 3, 4, 5]),
    "D": ("""
3.75 2.479e-1 8.385e-2 4.285e-2 1.392 3.528e-1 1.678e-1
4.00 2.591e-1 8.375e-2 4.240e-2 1.954 4.605e-1 2.013e-1
4.25 2.558e-1 8.145e-2 4.053e-2 2.760 6.158e-1 2.496e-1
4.50 2.375e-1 7.635e-2 3.781e-2 3.900 8.407e-1 3.309e-1
4.75 2.071e-1 6.833e-2 3.398e-2 5.496 1.182 4.661e-1
5.00 1.694e-1 5.755e-2 2.874e-2 7.599 1.672 6.649e-1
5.25 1.292e-1 4.502e-2 2.253e-2 9.990 2.262 9.058e-1
5.50 9.188e-2 3.264e-2 1.635e-2 1.226e+1 2.841 1.143
5.75 6.137e-2 2.212e-2 1.108e-2 1.408e+1 3.319 1.339
""", [3, 4, 5], [3, 4, 5]),
    "F": ("""
3.75 4.106e-2 2.766e-2 3.325e-1 1.840e-1
4.00 4.019e-2 2.809e-2 3.889e-1 2.165e-1
4.25 3.838e-2 2.713e-2 4.620e-1 2.564e-1
4.50 3.475e-2 2.450e-2 5.534e-1 3.066e-1
4.75 2.954e-2 2.068e-2 6.671e-1 3.692e-1
5.00 2.341e-2 1.628e-2 7.814e-1 4.329e-1
5.25 1.719e-2 1.192e-2 8.556e-1 4.769e-1
5.50 1.175e-2 8.147e-3 8.707e-1 4.903e-1
5.75 7.577e-3 5.257e-3 8.385e-1 4.774e-1
""", [4, 5], [4, 5]),
    "G": ("""
3.75 1.409e-2 8.969e-2
4.00 1.175e-2 8.533e-2
4.25 9.224e-3 7.928e-2
4.50 6.900e-3 7.392e-2
4.75 4.932e-3 6.903e-2
5.00 3.359e-3 6.279e-2
5.25 2.181e-3 5.440e-2
5.50 1.360e-3 4.498e-2
5.75 8.284e-4 3.616e-2
""", [5], [5]),
}


def _unpack(blocks):
    out = {}
    for L, (text, n_singlet, n_triplet) in blocks.items():
        a = _block(text)
        assert a.shape == (9, len(n_singlet) + len(n_triplet))
        for j, n in enumerate(n_singlet):
            out[(n, 1, L)] = a[:, j]
        for j, n in enumerate(n_triplet):
            out[(n, 3, L)] = a[:, len(n_singlet) + j]
    return out


# BRAY[lower][(n, 2S+1, L)] = Upsilon at BRAY_LOGT; lower 1 = 1^1S, 2 = 2^3S
BRAY = {1: _unpack(_T2), 2: _unpack(_T3)}

# ---------------------------------------------------------------- CHIANTI
_d = cc.ion_dir("he", 1)
LEV = cc.read_elvlc(os.path.join(_d, "he_1.elvlc"))     # {i: (g, Eobs, Eth)}
TRS = cc.read_scups(os.path.join(_d, "he_1.scups"))
TERMS = {}                                               # (n, 2S+1, L): [levels]
with open(os.path.join(_d, "he_1.elvlc")) as fh:
    for line in fh:
        t = line.split()
        if not t or t[0] == "-1":
            break
        n = int(t[1].split(".")[-1][0]) if "." in t[1] else 1
        TERMS.setdefault((n, int(t[2]), t[3]), []).append(int(t[0]))
LOWER_LEVEL = {1: 1, 2: 2}                               # 1^1S, 2^3S levels
LOWER_TERM = {1: (1, 1, "S"), 2: (2, 3, "S")}


def term_energy_K(term):
    """g-weighted observed energy of a term above 1^1S [K]."""
    us = TERMS[term]
    g = sum(LEV[u][0] for u in us)
    return sum(LEV[u][0] * LEV[u][1] for u in us) / g * HC_K


def chianti_term_upsilon(lower, term, T):
    """CHIANTI v11 he_1 strength lower -> term, J levels summed; None if
    CHIANTI carries no level of the term for this lower level."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    ll = LOWER_LEVEL[lower]
    tot, found = np.zeros_like(T), False
    for u in TERMS[term]:
        t = [x for x in TRS if x["ll"] == ll and x["ul"] == u]
        if t:
            tot = tot + cc.upsilon(t[0], T)
            found = True
    return tot if found else None


def chianti_log_slope(lower, term, logT, h=1.0e-3):
    """d log10 Upsilon / d log10 T of the CHIANTI term strength."""
    up = chianti_term_upsilon(lower, term, 10.0 ** (logT + h))
    dn = chianti_term_upsilon(lower, term, 10.0 ** (logT - h))
    return float((np.log10(up) - np.log10(dn))[0] / (2.0 * h))


def adopted_upsilon(lower, term, T):
    """The adopted strength lower -> term (module docstring)."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    L = np.log10(BRAY[lower][term])
    ch = chianti_term_upsilon(lower, term, T) is not None
    # node slopes in d log10 Upsilon per table index (0.25 dex)
    m = np.zeros(9)
    for i in range(1, 8):
        sl, sr = L[i] - L[i - 1], L[i + 1] - L[i]
        m[i] = 0.0 if sl * sr <= 0.0 else 2.0 / (1.0 / sl + 1.0 / sr)
    if ch:
        m[0] = 0.25 * chianti_log_slope(lower, term, BRAY_LOGT[0])
        m[8] = 0.25 * chianti_log_slope(lower, term, BRAY_LOGT[8])
    lt = np.log10(T)
    out = np.empty_like(T)
    for j, x in enumerate((lt - BRAY_LOGT[0]) / 0.25):
        if x < 0.0 or x > 8.0:
            edge = 0 if x < 0.0 else 8
            if ch:
                r = (chianti_term_upsilon(lower, term, T[j])
                     / chianti_term_upsilon(lower, term, 10.0 ** BRAY_LOGT[edge]))
                out[j] = 10.0 ** L[edge] * r[0]
            else:
                out[j] = 10.0 ** L[edge]
            continue
        k = min(int(x), 7)
        f = x - k
        h00 = (1.0 + 2.0 * f) * (1.0 - f) ** 2
        h10 = f * (1.0 - f) ** 2
        h01 = f * f * (3.0 - 2.0 * f)
        h11 = f * f * (f - 1.0)
        out[j] = 10.0 ** (h00 * L[k] + h10 * m[k] + h01 * L[k + 1] + h11 * m[k + 1])
    return out


def upper_terms(lower):
    """Every term above the lower one that Table 2 (lower 1) or Table 3
    (lower 2) carries."""
    return sorted(BRAY[lower])


def smoothness_report():
    """Largest departure of each tabulated log10 Upsilon from the mean of its
    two neighbors: a misprint (or a transcription slip) stands out."""
    rows = []
    for lower in (1, 2):
        for term, v in BRAY[lower].items():
            L = np.log10(v)
            d = np.abs(L[1:-1] - 0.5 * (L[:-2] + L[2:]))
            rows.append((d.max(), lower, term, BRAY_LOGT[1 + d.argmax()]))
    return sorted(rows, reverse=True)


def overshoot_report():
    """Largest excursion of the adopted log10 Upsilon beyond the range of the
    two tabulated values bounding an interval [dex], over every term."""
    rows = []
    for lower in (1, 2):
        for term, v in BRAY[lower].items():
            L = np.log10(v)
            lt = BRAY_LOGT[0] + 0.25 * np.linspace(0.0, 8.0, 801)
            a = np.log10(adopted_upsilon(lower, term, 10.0 ** lt))
            k = np.minimum((np.arange(801) // 100), 7)
            ex = np.maximum(a - np.maximum(L[k], L[k + 1]),
                            np.minimum(L[k], L[k + 1]) - a)
            rows.append((ex.max(), lower, term))
    return sorted(rows, reverse=True)


if __name__ == "__main__":
    print("Largest second differences of log10 Upsilon (dex):")
    for d, lower, term, lt in smoothness_report()[:6]:
        print(f"  {d:.3f}  lower {lower} -> {term} at log T {lt:.2f}")
    print("Largest overshoot of the adopted strength inside the table (dex):")
    for d, lower, term in overshoot_report()[:4]:
        print(f"  {d:.4f}  lower {lower} -> {term}")
    print("CHIANTI / Bray at log T 3.75 ... 5.75, term by term:")
    Tb = 10.0 ** BRAY_LOGT
    for lower in (1, 2):
        for term in upper_terms(lower):
            c = chianti_term_upsilon(lower, term, Tb)
            if c is None:
                print(f"  {LOWER_TERM[lower]} -> {term}: not in CHIANTI he_1")
                continue
            r = c / BRAY[lower][term]
            print(f"  {LOWER_TERM[lower]} -> {term}: {r.min():.3f} - {r.max():.3f}")
