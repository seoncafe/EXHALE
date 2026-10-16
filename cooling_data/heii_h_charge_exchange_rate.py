#!/usr/bin/env python3
"""Rate coefficient of He+(1s) + H(1s) -> He(1s^2) + H+ (non-radiative) for
`nonradiative_electron_capture_Hep_from_H` in
src/modules/radiation/charge_exchange.f90.

Source: Loreau, Ryabchenko & Vaeck (2014), J. Phys. B 47, 135204, eq. (3)
and Table 2, process (1): the cross section of the REVERSE reaction
H+ + He(1s^2) -> H(1s) + He+(1s), log10(sigma / 1e-16 cm^2) as a Fourier
series in log10 E [eV/u], E the energy per nucleon of the proton (so that
E_cm = mu E with mu in u).  The forward cross section follows by
microscopic reversibility with statistical weights 4 (He+(1s) + H(1s)) and
1 (He(1^1S) + H+):
    sigma_f(E_f) = (1/4) (mu_r E_r) / (mu_f E_f) sigma_r(E_r),
    E_r = E_f + dE,  dE = I(He) - I(H) = 10.988953 eV.

Near threshold.  Loreau et al.'s calculation starts at 15 eV/u; the
reverse threshold lies at x_th = dE (m_p + m_He)/(m_p m_He) = 13.655 eV/u.
Between x_th and 15 eV/u the reverse cross section is taken as
    sigma_r(x) = sigma_r(15) sqrt((x - x_th)/(15 - x_th)),  0 below x_th,
which by reciprocity is sigma_f v -> constant at low forward energy
(capture on the ion-induced dipole potential with a transfer probability
independent of energy).  The two alternatives bracket it: sigma_r = 0
below 15 eV/u closes an open channel (0.46x at 1e4 K, 1.5e-5x at 1e3 K),
and the fit carried through the threshold, finite there, gives
sigma_f ~ 1/E_f (1.32x at 1e4 K, 3.6x at 1e3 K).  The assumption of a
constant transfer probability over forward energies 0 - 1.08 eV is not
tested by any calculation.  Adopted on 2026-09-26 in both EXHALE and MoCHII
so that the two codes carry the same rate.

The functions below follow the MoCHII generator
(/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tools/fitting/make_heii_h_charge_exchange.py,
same author), evaluated on the EXHALE temperature grid.  With low='zero'
they reproduce the EXHALE table carried before 2026-09-26 to 0.4 %.

Usage: python3 heii_h_charge_exchange_rate.py   (prints the Fortran table
and the two bracketing variants).
"""
import math
import numpy as np

# Loreau et al. (2014) Table 2, process (1)
A = [-5.167, -5.200, -0.3265, -0.3786, -0.04399, -0.05318, -0.06957, 0.01072]
B = [0.0, 0.1131, 1.279, 0.4867, 0.1944, 0.1305, 0.03501, -0.04530]
W = 0.7800

U = 1.66053906660e-24          # g
EV = 1.602176634e-12           # erg
KB = 8.617333262e-5            # eV/K
ME = 5.485799090e-4            # u
MP = 1.007276467               # u
MHE_ATOM = 4.001506179 + 2*ME  # u
MU_R = MP*MHE_ATOM/(MP + MHE_ATOM)                  # H+ + He
MU_F = (MP + ME)*(MHE_ATOM - ME)/(MP + MHE_ATOM)    # He+ + H
DE = 24.587387 - 13.598434     # eV
X_TH = DE*(MP + MHE_ATOM)/(MP*MHE_ATOM)             # reverse threshold [eV/u]
X_DATA = 15.0                                       # first quantal energy [eV/u]

# The temperature grid of the Fortran table [K]
T_GRID = [3.0e2, 5.0e2, 7.0e2, 1.0e3, 1.5e3, 2.0e3, 3.0e3, 4.0e3, 5.0e3,
          6.0e3, 7.0e3, 8.0e3, 1.0e4, 1.2e4, 1.5e4, 2.0e4, 2.5e4, 3.0e4,
          4.0e4, 5.0e4, 6.0e4, 8.0e4, 1.0e5]


def sigma_rev(x):
    """Loreau process (1) cross section [cm^2] at x eV/u."""
    L = np.log10(x)
    s = A[0] + sum(A[n]*np.cos(n*W*L) + B[n]*np.sin(n*W*L) for n in range(1, 8))
    return 10.0**s*1.0e-16


def k_forward(T, low="threshold", npts=200001):
    """Maxwellian <sigma_f v> [cm^3 s^-1]; below 15 eV/u, low = 'threshold'
    (the adopted form), 'zero' or 'fit' (the fit as given)."""
    kT = KB*T
    ef = np.exp(np.linspace(math.log(1e-7), math.log(60.0*kT + 60.0), npts))
    er = ef + DE
    x = er*(MP + MHE_ATOM)/(MP*MHE_ATOM)
    s = sigma_rev(x)
    if low == "zero":
        s = np.where(x < X_DATA, 0.0, s)
    elif low == "threshold":
        s15 = sigma_rev(np.array([X_DATA]))[0]
        s = np.where(x < X_DATA,
                     s15*np.sqrt(np.clip((x - X_TH)/(X_DATA - X_TH), 0.0, None)), s)
    sf = 0.25*(MU_R*er)/(MU_F*ef)*s
    f = ef*sf*np.exp(-ef/kT)*ef
    integ = np.sum(0.5*(f[1:] + f[:-1])*np.diff(np.log(ef)))
    return math.sqrt(8.0/(math.pi*MU_F*U))*(kT*EV)**-1.5*integ*EV**2


def main():
    k = [k_forward(t) for t in T_GRID]
    items = ['%.3fd%d' % (v/10**math.floor(math.log10(v)), math.floor(math.log10(v)))
             for v in k]
    print('      real*8, parameter  :: k_l(n_l) = [ ' + ', '.join(items) + ' ]')
    print('\n   T [K]    adopted     sigma_r=0 below 15   fit through threshold')
    for t in T_GRID:
        print('%9.0f  %.4e  %.4e  %.4e' % (t, k_forward(t), k_forward(t, low='zero'),
                                            k_forward(t, low='fit')))


if __name__ == '__main__':
    main()
