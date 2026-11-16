#!/usr/bin/env python3
"""Build the molecular infrared data module of EXHALE.

Writes ``src/modules/lower_atmosphere/molecular_infrared_data.f90``, the
tabulated input of ``molecular_infrared_cooling.f90``: the band-mean
absorption cross sections of H2O and CO, and the H2 infrared line list with
its rovibrational level ladder.

Physics the tables are meant to serve
-------------------------------------
The molecular layer between a 1 microbar base and the H2 -> H front exchanges
infrared photons with the (unresolved) atmosphere below it.  In the optically
thin, LTE limit the net radiative loss of one molecule of species X is

    L_net(T; T_rad, W) = 4 pi Int sigma_nu(T) [ B_nu(T) - W B_nu(T_rad) ] dnu

for a band absorber, and the equivalent line sum for H2.  Both vanish
identically at T = T_rad, W = 1, which is what makes the layer settle on a
radiative equilibrium temperature instead of radiating itself to zero.

H2O and CO --- band-mean cross sections
---------------------------------------
Source: the correlated-k coefficients distributed with Photochem
(Wogan et al. 2025, PSJ 6, 256), computed with a fork of HELIOS-K
(Grimm & Heng 2015, ApJ 808, 182; Grimm et al. 2021, ApJS 253, 30) from
HITEMP 2010 for H2O (Rothman et al. 2010, JQSRT 111, 2139) and HITEMP 2019
for CO (Li et al. 2015, ApJS 216, 15).  Pressures 1e-6 to 1e3 bar,
temperatures 50 to 2000 K.  The on-disk file is
``photochem_clima_data/photochem_clima_data/data/kdistributions/{H2O,CO}.h5``.

The k-distribution stores the sorted absorption cross section k(g) inside each
wavelength bin, so the bin-mean cross section is exactly Sum_g w_g k_g.  Line
strengths do not depend on the broadening, only the line shape does, so this
bin mean is nearly pressure independent: measured over 1e-6 to 1 bar it moves
by less than 2 percent at 1000 K, which is why only the 1e-6 bar slice is
carried.  That is also the lowest tabulated pressure, and the EXHALE domain
sits at and below it.

HITRAN/HITEMP line intensities carry the stimulated-emission factor
1 - exp(-h nu / k T), so Sum_g w_g k_g is the absorption coefficient corrected
for stimulated emission, which is the quantity that pairs with B_nu in the net
exchange above.

Cross-checks performed when this table was built (2026-08-31), all against
sources that share none of its line processing:

* The Planck-mean H2O cross section derived here agrees with the independent
  Planck-mean table ``aiolos/inputdata/1H2-16O_T400.aiopa`` to within 11
  percent over 600-2000 K (ratio 0.889-1.111 on the 14 shared grid
  temperatures).  BELOW 600 K the comparison is not meaningful: that table is
  not monotonic in T there (9.9e-21 at 100 K, 5.4e-20 at 200 K, 1.5e-20 at
  300 K, 7.9e-21 at 500 K, 1.2e-20 at 600 K), while the table built here is
  monotonic over its whole range.  The low-temperature end is covered by the
  Neufeld & Kaufman comparison below instead.
* Restricted to lambda > 9 um, i.e. to the pure rotational band, the H2O
  emission agrees with the optically thin LTE rate of Neufeld & Kaufman
  (1993), ApJ 418, 263, table 2 (Ntilde = 1e10) to a factor 0.70 (100 K),
  0.81 (200 K), 0.90 (400 K) and 1.02 (1000 K).  The cut is not a tuned
  number: moving it over 7-15 um leaves those ratios at 0.70 / 0.81 / 0.90
  unchanged and moves the 1000 K one only over 0.92-1.22.  The same cut on CO against
  their table 3 gives 0.81-0.93 over 300-2000 K.  At 100 K the CO ratio falls
  to 0.10 and the reason is NOT the coarse long-wavelength binning -- measured,
  97% of CO's 100 K emission lands in the eighteen 9-250 um bins and only 2.6%
  in the single 250-10000 um bin.  The likely cause is the g-quadrature: a
  k-distribution is built to reproduce TRANSMISSION, and its eight points
  resolve the linear bin mean poorly when the spectrum inside a bin is a few
  isolated narrow lines, which is what CO's rotational ladder is at 100 K and
  1e-6 bar.  The discrepancy shrinks as the ladder fills (0.81-0.93 above
  300 K), which is consistent with that reading.  It does not matter here: at
  100 K CO radiates 1.6e-20 erg/s per molecule against 1.8e-15 for H2O, five
  orders of magnitude less.
* The H2 line sum agrees with the LTE rates of Hollenbach & McKee (1979),
  ApJS 41, 555, eqs. (6.37) and (6.38) to 1-11 percent over 200-2000 K
  (rotational 1.05-1.11 above 400 K, vibrational 0.85-0.97 above 600 K).

Two numbers from those comparisons are the reason this module integrates the
whole band system instead of adopting a published rotational cooling function:
above 300 K, CO's total emission is four to five orders of magnitude above its
rotational part (the 4.7 um fundamental), and H2O's is 3.2 times its
rotational part at 1000 K (the 6.3 um bend and the 2.7 um stretches).  A
rotation-only coolant would have missed almost all of it.

H2 --- infrared line list
-------------------------
Source: Roueff, Abgrall, Czachorowski, Pachucki, Puchalski & Komasa (2019),
A&A 630, A58, table 2, retrieved from VizieR as catalog J/A+A/630/A58 and kept
byte-identical in ``roueff2019_h2_infrared_lines.dat`` (column description in
``roueff2019_h2_infrared_lines.ReadMe``).  It lists every infrared transition
within the X state of H2 with the electric quadrupole and magnetic dipole
Einstein coefficients, the transition wavenumber, the upper level term energy
measured from (v,J) = (0,0), and the upper level statistical weight
g = g_I (2J+1).

H2 has no permitted dipole spectrum; these quadrupole and magnetic dipole
lines are its whole infrared emission, and they are missing from every cooling
channel EXHALE carried before this table.

The (0,0) and (0,1) levels are never the upper level of an infrared
transition, so they do not appear in the file.  (0,0) is the origin.  (0,1) is
recovered exactly from the (0,3) -> (0,1) row as T(0,3) - hc sigma / k, which
is the ortho-para offset, and its weight is 3 x (2x1+1) = 9.

Usage
-----
    python3 cooling_data/molecular_infrared_bands.py

Run from the repository root.  Regenerating is only needed when the underlying
line lists or k-coefficients change.
"""

import os
import sys

import numpy as np

try:
    import h5py
except ImportError:  # pragma: no cover
    sys.exit("h5py is required to read the k-distribution files")

# ---------------------------------------------------------------- constants
H_PLANCK = 6.62607015e-27      # erg s        (CODATA 2018, exact)
C_LIGHT = 2.99792458e10        # cm s^-1      (exact)
K_BOLTZ = 1.380649e-16         # erg K^-1     (CODATA 2018, exact)
CM_TO_K = H_PLANCK * C_LIGHT / K_BOLTZ   # cm^-1 -> K

# numpy renamed trapz to trapezoid in 2.0; this machine carries 1.26.
TRAPZ = getattr(np, "trapezoid", None) or np.trapz

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
KDIR = os.path.join(ROOT, "photochem_clima_data", "photochem_clima_data",
                    "data", "kdistributions")
H2_LINES = os.path.join(HERE, "roueff2019_h2_infrared_lines.dat")
OUT = os.path.join(ROOT, "src", "modules", "lower_atmosphere",
                   "molecular_infrared_data.f90")

# Bands whose Planck-weighted share of the emission stays below this fraction
# at every gas temperature and every radiating temperature of the probe set
# are dropped.  The truncation error on the band sum is reported below.
BAND_KEEP = 1.0e-6
PROBE_T = np.array([100., 200., 400., 700., 1000., 1400., 2000., 3000.])

# H2 lines whose share of the LTE emission stays below this fraction at every
# probe temperature are dropped; the truncation error is reported below.
LINE_KEEP = 1.0e-5
H2_PROBE_T = np.array([100., 200., 400., 700., 1000., 1500., 2000.,
                       3000., 4000., 6000.])


# ------------------------------------------------------------------ helpers
def planck_band_integral(wl_um, T, nsub=128):
    """Int B_nu(T) dnu over each wavelength bin [erg cm^-2 s^-1 sr^-1]."""
    out = np.zeros(len(wl_um) - 1)
    for i in range(len(wl_um) - 1):
        lam = np.geomspace(wl_um[i], wl_um[i + 1], nsub) * 1.0e-4
        nu = C_LIGHT / lam
        x = H_PLANCK * nu / (K_BOLTZ * T)
        b = np.where(x < 600.0,
                     2.0 * H_PLANCK * nu ** 3 / C_LIGHT ** 2
                     / np.expm1(np.minimum(x, 600.0)), 0.0)
        out[i] = abs(TRAPZ(b, nu))
    return out


def band_mean_cross_section(species):
    """Bin-mean cross section [cm^2/molecule] at the lowest tabulated pressure."""
    with h5py.File(os.path.join(KDIR, species + ".h5"), "r") as f:
        temp = f["T"][:].astype(float)
        wl = f["wavelengths"][:].astype(float)
        wts = f["weights"][:].astype(float)
        k = 10.0 ** f["log10k"][:].astype(float)
        logp = f["log10P"][:].astype(float)
    ip = int(np.argmin(logp))           # 1e-6 bar
    sig = (k * wts[None, None, None, :]).sum(axis=3)[:, :, ip]   # (nband, nT)
    return temp, wl, sig, logp[ip]


def select_bands(temp, wl, sig):
    """Keep the bands that carry the emission; report the truncation error."""
    share = np.zeros(sig.shape[0])
    for tg in PROBE_T:
        s = np.array([np.interp(min(tg, temp[-1]), temp, sig[b])
                      for b in range(sig.shape[0])])
        for tr in PROBE_T:
            v = s * planck_band_integral(wl, tr)
            tot = v.sum()
            if tot > 0.0:
                share = np.maximum(share, v / tot)
    keep = np.where(share > BAND_KEEP)[0]
    err = 0.0
    for tg in PROBE_T:
        s = np.array([np.interp(min(tg, temp[-1]), temp, sig[b])
                      for b in range(sig.shape[0])])
        v = s * planck_band_integral(wl, tg)
        if v.sum() > 0.0:
            err = max(err, abs(1.0 - v[keep].sum() / v.sum()))
    return keep, err


def read_h2_lines():
    """Roueff et al. (2019) table 2 -> (lines, levels)."""
    lines = []
    levels = {}
    e01 = None
    with open(H2_LINES) as fh:
        for rec in fh:
            if len(rec) < 164:
                continue
            vu = int(rec[1:3])
            ju = int(rec[4:6])
            vl = int(rec[7:9])
            jl = int(rec[10:12])
            sigma = float(rec[13:29])        # cm^-1
            a_ul = float(rec[96:105])        # s^-1, quadrupole + magnetic dipole
            t_up = float(rec[148:159])       # K, from (0,0)
            g_up = float(rec[161:164])
            lines.append((sigma * CM_TO_K, a_ul, t_up, g_up))
            levels[(vu, ju)] = (t_up, g_up)
            if (vu, ju, vl, jl) == (0, 3, 0, 1):
                e01 = t_up - sigma * CM_TO_K
    if e01 is None:
        sys.exit("the (0,3) -> (0,1) row is missing; cannot place the "
                 "ortho-para offset")
    levels[(0, 0)] = (0.0, 1.0)
    levels[(0, 1)] = (e01, 9.0)
    return np.array(lines), levels


def select_h2_lines(lines, levels):
    lev_t = np.array([v[0] for v in levels.values()])
    lev_g = np.array([v[1] for v in levels.values()])
    share = np.zeros(len(lines))
    for tg in H2_PROBE_T:
        z = np.sum(lev_g * np.exp(-lev_t / tg))
        v = lines[:, 3] * np.exp(-lines[:, 2] / tg) / z * lines[:, 1] \
            * lines[:, 0] * K_BOLTZ
        share = np.maximum(share, v / v.sum())
    keep = np.where(share > LINE_KEEP)[0]
    err = 0.0
    for tg in H2_PROBE_T:
        z = np.sum(lev_g * np.exp(-lev_t / tg))
        v = lines[:, 3] * np.exp(-lines[:, 2] / tg) / z * lines[:, 1] \
            * lines[:, 0] * K_BOLTZ
        err = max(err, abs(1.0 - v[keep].sum() / v.sum()))
    return keep, err


# ------------------------------------------------------------- source emitter
def emit_data(fh, element, values, per_line=5, fmt="{:23.15e}"):
    """Write one `data` statement per chunk of at most CHUNK values.

    `element` is the array element designator with `i` as the running index,
    e.g. "h2o_sig(i,7)".  Chunking keeps every statement well inside the
    continuation-line limit of each of the three compilers the Makefile
    supports.
    """
    n = len(values)
    chunk = 150
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        fh.write("      data ({0}, i = {1}, {2}) /".format(
            element, start + 1, stop))
        for j in range(start, stop):
            if (j - start) % per_line == 0:
                fh.write("  &\n         ")
            fh.write(fmt.format(values[j]).replace("e", "d"))
            if j != stop - 1:
                fh.write(",")
        fh.write(" /\n")


def main():
    print("reading k-distributions from", KDIR)
    data = {}
    for sp in ("H2O", "CO"):
        temp, wl, sig, logp = band_mean_cross_section(sp)
        keep, err = select_bands(temp, wl, sig)
        data[sp] = dict(T=temp, wl_lo=wl[:-1][keep], wl_hi=wl[1:][keep],
                        sig=sig[keep, :], logp=logp)
        print("  {0}: {1} of {2} bands kept, {3:.3f}-{4:.1f} um, "
              "band-sum truncation error {5:.2e}".format(
                  sp, len(keep), sig.shape[0], wl[:-1][keep].min(),
                  wl[1:][keep].max(), err))

    print("reading the H2 line list from", H2_LINES)
    lines, levels = read_h2_lines()
    keep, lerr = select_h2_lines(lines, levels)
    lines_k = lines[keep, :]
    lev_t = np.array([v[0] for v in levels.values()])
    lev_g = np.array([v[1] for v in levels.values()])
    order = np.argsort(lev_t)
    lev_t = lev_t[order]
    lev_g = lev_g[order]
    print("  H2: {0} of {1} lines kept, {2} levels, line-sum truncation "
          "error {3:.2e}".format(len(keep), len(lines), len(lev_t), lerr))

    temp = data["H2O"]["T"]
    if not np.allclose(temp, data["CO"]["T"]):
        sys.exit("H2O and CO are tabulated on different temperature grids")

    with open(OUT, "w") as fh:
        fh.write(HEADER.format(
            n_t=len(temp), nb_h2o=data["H2O"]["sig"].shape[0],
            nb_co=data["CO"]["sig"].shape[0],
            n_lev=len(lev_t), n_line=lines_k.shape[0],
            t_lo=temp[0], t_hi=temp[-1],
            p_bar=10.0 ** data["H2O"]["logp"]))
        emit_data(fh, "mir_T(i)", temp)
        fh.write("\n")
        for sp, tag in (("H2O", "h2o"), ("CO", "co")):
            d = data[sp]
            emit_data(fh, tag + "_wl_lo(i)", d["wl_lo"], fmt="{:20.12e}")
            emit_data(fh, tag + "_wl_hi(i)", d["wl_hi"], fmt="{:20.12e}")
            fh.write("\n")
            for it in range(len(temp)):
                emit_data(fh, "{0}_sig(i,{1})".format(tag, it + 1),
                          d["sig"][:, it])
            fh.write("\n")
        emit_data(fh, "h2_lev_T(i)", lev_t)
        emit_data(fh, "h2_lev_g(i)", lev_g, per_line=10, fmt="{:8.1f}")
        fh.write("\n")
        emit_data(fh, "h2_line_dE(i)", lines_k[:, 0])
        emit_data(fh, "h2_line_A(i)", lines_k[:, 1])
        emit_data(fh, "h2_line_Tu(i)", lines_k[:, 2])
        emit_data(fh, "h2_line_gu(i)", lines_k[:, 3], per_line=10,
                  fmt="{:8.1f}")
        fh.write(FOOTER)
    print("wrote", OUT)


HEADER = '''\
!==============================================================================
! molecular_infrared_data -- tabulated input of molecular_infrared_cooling
!
! GENERATED FILE.  Do not edit by hand; rerun
!     python3 cooling_data/molecular_infrared_bands.py
! from the repository root.  That script carries the provenance, the validity
! ranges and the cross-checks; this file carries only the numbers.
!
! H2O, CO: band-mean absorption cross sections [cm^2/molecule], at the lowest
!   tabulated pressure ({p_bar:.1e} bar; the bin mean is pressure independent to
!   better than 2% up to 1 bar, so the whole EXHALE domain uses this slice).
!   Correlated-k coefficients distributed with Photochem (Wogan et al. 2025,
!   PSJ 6, 256), from HITEMP 2010 (H2O; Rothman et al. 2010, JQSRT 111, 2139)
!   and HITEMP 2019 (CO; Li et al. 2015, ApJS 216, 15) processed with HELIOS-K
!   (Grimm & Heng 2015, ApJ 808, 182).  Tabulated {t_lo:.0f}-{t_hi:.0f} K.
!   The cross sections carry the HITRAN stimulated-emission correction, so
!   they are the coefficients that pair with B_nu in the net exchange.
!
! H2: the infrared line list of Roueff et al. (2019), A&A 630, A58, table 2
!   (VizieR J/A+A/630/A58), electric quadrupole plus magnetic dipole.  Level
!   term energies and statistical weights g = g_I (2J+1) come from the same
!   table; (0,0) is the origin and (0,1) is placed by the (0,3) -> (0,1) row.
!==============================================================================
      module molecular_infrared_data

      implicit none
      public
      save

      integer, parameter :: n_mir_T   = {n_t}
      integer, parameter :: nb_h2o    = {nb_h2o}
      integer, parameter :: nb_co     = {nb_co}
      integer, parameter :: n_h2_lev  = {n_lev}
      integer, parameter :: n_h2_line = {n_line}

      ! temperature grid of the H2O and CO cross sections [K]
      real*8 :: mir_T(n_mir_T)

      ! band edges [micron] and band-mean cross sections [cm^2/molecule]
      real*8 :: h2o_wl_lo(nb_h2o), h2o_wl_hi(nb_h2o)
      real*8 :: co_wl_lo(nb_co),   co_wl_hi(nb_co)
      real*8 :: h2o_sig(nb_h2o, n_mir_T)
      real*8 :: co_sig(nb_co, n_mir_T)

      ! H2 rovibrational ladder: term energy [K] from (0,0), and weight
      real*8 :: h2_lev_T(n_h2_lev), h2_lev_g(n_h2_lev)

      ! H2 lines: transition energy [K], Einstein A [s^-1], upper level term
      ! energy [K] and upper level weight
      real*8 :: h2_line_dE(n_h2_line), h2_line_A(n_h2_line)
      real*8 :: h2_line_Tu(n_h2_line), h2_line_gu(n_h2_line)

      ! Implied-do index of the data statements below, and nothing else.
      ! PRIVATE: the module is otherwise public and molecular_infrared_cooling
      ! uses it without an only: list, so a public "i" would be in scope in
      ! every routine there -- a shared, saved integer that implicit none
      ! cannot catch as an undeclared loop index and that the OpenMP cell
      ! sweep would write from several threads.
      integer, private :: i

'''

FOOTER = '''
      end module molecular_infrared_data
'''


if __name__ == "__main__":
    main()
