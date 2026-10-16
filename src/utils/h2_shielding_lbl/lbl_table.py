"""Overlapping-line H2 Lyman-Werner table: sigma_pump, p_single, sigma_diss.

Same line data, populations and slab as lbl.py, but it returns the ABSOLUTE
quantities the code needs, not just a normalized shielding factor:

  sigma_pump(N)   pumps per unit incident band photon fluence [cm^2]
  p_single(N)     dissociations per pump with every fluorescent photon free to
                  escape, sum_i w_i p_i / sum_i w_i with w_i the pumping rate
                  of line i and p_i = D_u/(D_u + sum A_ul) from the level data
  sigma_diss(N)   = p_eff x sigma_pump

The trapping that turns p_single into p_eff is NOT computed here.  It is a
property of the fluorescent lines' escape probabilities and of the slab
geometry, and CLOUDY solves it (escape probabilities per transition, both
directions, iterated); reimplementing it would be a second, unvalidated
version of the same thing.  The ratio p_eff/p_single is therefore taken from
the CLOUDY runs node by node, and everything that line overlap touches -- the
pumping rate, hence sigma_pump and the weighting inside p_single -- is taken
from here.  See the module header of the generated Fortran.
"""
import numpy as np
from scipy.special import wofz
import lbl

CLIGHT = lbl.CLIGHT
HPL = lbl.HPL
PI_E2_MC = lbl.PI_E2_MC

# THE BAND, AND THE ONE NORMALIZATION IT FIXES.  The Lyman-Werner band of
# EXHALE is 912-1201 A: the H Lyman edge, below which atomic H absorbs
# everything, to the start of band B2 of oxygen_rates.f90.  The line list
# built below reaches the same 1201 A, so every line that absorbs a band
# photon is normalized per photon of the band that carries it.  This is one
# number: the column integral of sigma_pump then tends to the fraction of
# band photons the lines can take, which cannot pass 1, and the beam of
# lyman_werner.f90 loses exactly what this table rates.
#
# F_BAND_ERG is the flat band flux of the deck.  The deck of the earlier
# 912-1110 A normalization carried 343.0 erg cm^-2 s^-1 over 198 A, i.e.
# F_lambda = 343.0/198 erg cm^-2 s^-1 A^-1; over the 289 A of 912-1201 A the
# same F_lambda integrates to 343.0*289/198 = 500.6.  The value does not
# enter the result: sigma_pump below is pump/F_BAND_PHOT, pump is linear in
# F_lam_cm = F_BAND_ERG/(WL_HI - WL_LO), and F_BAND_PHOT = F_BAND_ERG/
# E_LW_PHOTON_ERG, so F_BAND_ERG cancels exactly and sigma_pump depends on
# the band only through its width and its mean photon energy.
#
# E_LW_PHOTON_ERG is the mean photon energy of a flat-F_lambda band,
# 2hc/(WL_LO + WL_HI) = 2hc/2113 A = 11.7354 eV, so that F_BAND_PHOT is the
# exact photon content of the band.  It is lyman_werner.f90's
# e_lw_photon_erg, which converts the user's band energy flux to the same
# photon flux; the two must be the same number.
F_BAND_ERG = 500.6             # flat 912-1201 A band flux of the deck
E_LW_PHOTON_ERG = 1.88021e-11  # lyman_werner.f90's mean band photon energy
F_BAND_PHOT = F_BAND_ERG / E_LW_PHOTON_ERG
WL_LO, WL_HI = 912.0, 1201.0   # the band the flux is normalized over


def absolute_curves(T, Ncols, nu_pts=400_000, wl_min=911.75, wl_max=1201.0,
                    core_hw=300.0, keep_frac=1.0 - 1e-7, verbose=False):
    """(sigma_pump, p_single) with overlapping lines, at each column."""
    ll = lbl.build_linelist(T, wl_min, wl_max)
    nu0, f, xl, Gam, pdis = (ll['nu0'], ll['f'], ll['xl'], ll['Gam'],
                             ll['pdis'])
    w0 = xl * f
    order = np.argsort(w0)[::-1]
    c = np.cumsum(w0[order]) / w0.sum()
    keep = order[:int(np.searchsorted(c, keep_frac) + 1)]
    nu0, f, xl, Gam, pdis = (a[keep] for a in (nu0, f, xl, Gam, pdis))
    nl = len(nu0)
    b = np.sqrt(2.0 * lbl.KB * T / lbl.M_H2)
    dnuD = nu0 * b / CLIGHT
    av = Gam / (4.0 * np.pi * dnuD)

    nu = np.linspace(CLIGHT / (wl_max * 1e-8), CLIGHT / (wl_min * 1e-8), nu_pts)
    dnu = nu[1] - nu[0]
    # flat F_lambda normalized to F_BAND_ERG over WL_LO-WL_HI:
    #   F_lambda = F_BAND_ERG/(WL_HI-WL_LO) [erg cm^-2 s^-1 A^-1]
    #   F_nu = F_lambda lambda^2/c ; photons per unit nu = F_nu/(h nu)
    F_lam_cm = F_BAND_ERG / (WL_HI - WL_LO) * 1.0e8      # per cm of wavelength
    lam = CLIGHT / nu
    phot = F_lam_cm * lam ** 2 / CLIGHT / (HPL * nu)      # photons cm^-2 s^-1 Hz^-1

    sig = np.zeros(nu_pts)
    for i in range(nl):
        sig += PI_E2_MC * f[i] * xl[i] * lbl._profiles(nu, nu0[i], dnuD[i],
                                                       av[i], Gam[i], core_hw)
    Ncols = np.atleast_1d(np.asarray(Ncols, dtype=float))
    att = np.exp(-np.outer(Ncols, sig))
    pump = np.zeros(len(Ncols))
    diss = np.zeros(len(Ncols))
    for i in range(nl):
        phi = lbl._profiles(nu, nu0[i], dnuD[i], av[i], Gam[i], core_hw)
        w = PI_E2_MC * f[i] * xl[i] * phi * phot * dnu
        r = att @ w
        pump += r
        diss += pdis[i] * r
        if verbose and i % 500 == 0:
            print('   line %d/%d' % (i, nl), flush=True)
    return pump / F_BAND_PHOT, diss / pump
