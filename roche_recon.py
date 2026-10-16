"""Roche-equipotential 3-D reconstruction for the transit transmission spectrum
(Phase 5b of the Huang+2023 WASP-121b reproduction).

The 1-D EXHALE Roche run gives the SUBSTELLAR profile: the state (rho, v, T, ion
densities) as a function of radius r along the planet->star axis (+x).  Assuming
the thermodynamic state is uniform on Roche-equipotential surfaces, the 3-D field
at any point (x,y,z) equals the substellar state at the "equivalent radius"
r_eff whose substellar potential equals the local 3-D potential:

    phi_sub(r_eff) = phi_3D(x,y,z)  ->  state(x,y,z) = state_sub(r_eff).

Positions are in units of the planet base radius R_p.  The dimensionless Roche
potential (units GM_p/R_p; the overall scale cancels in the mapping) mirrors
grav_field.f90 (the spherical_domain=.false. branch):

    phi(x,y,z) = -1/r - q/r_s - (1+q)/(2 a^3) [(x - x_cm)^2 + y^2],

with r = |(x,y,z)|, r_s = |(x,y,z)-(a,0,0)| the distance to the star, q = M*/M_p,
a = orbital distance [R_p], x_cm = q a/(1+q) the centre of mass.  Setting tidal
terms -> 0 (q->0) recovers phi = -1/r, i.e. the spherical geometry, so the
transit reconstruction reduces to the Phase-5a spherical case for Case A.

Run `python3 roche_recon.py` for the self-test (monotonicity, inversion round
trip, L1, triaxial-radii relations).
"""

import os
import sys

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq

# The Jupiter radius is defined once, in examples/exhale_io.py (RJ_CM), and
# imported rather than written down again, so that this module's self-test
# and the run it reconstructs mean one planet by `Planet radius [R_J]`.
# `__file__` is absent under exec(), hence the guard.
try:
    _HERE = os.path.dirname(os.path.realpath(__file__))
except NameError:
    _HERE = os.getcwd()
_EXAMPLES_DIR = os.path.join(_HERE, 'examples')
if _EXAMPLES_DIR not in sys.path:
    sys.path.append(_EXAMPLES_DIR)
from exhale_io import RJ_CM                                # noqa: E402


# --------------------------------------------------------------------------- #
#  Dimensionless Roche potential (R_p units, GM_p/R_p energy units)
# --------------------------------------------------------------------------- #
def roche_phi(x, y, z, q, a):
    """Full 3-D Roche potential at (x,y,z) [R_p].  Star at (a,0,0)."""
    r  = np.sqrt(x*x + y*y + z*z)
    rs = np.sqrt((x - a)**2 + y*y + z*z)
    xcm = q*a/(1.0 + q)
    return -1.0/r - q/rs - 0.5*(1.0 + q)/a**3*((x - xcm)**2 + y*y)


def roche_phi_sub(r, q, a):
    """Substellar potential phi(r,0,0) for r in (0, a)."""
    return roche_phi(r, 0.0, 0.0, q, a)


def roche_dphidr_sub(r, q, a):
    """d/dr of the substellar potential (= grav_field.f90 Dphi / b0)."""
    xcm = q*a/(1.0 + q)
    return 1.0/r**2 - q/(a - r)**2 + (1.0 + q)/a**3*(xcm - r)


def find_L1(q, a):
    """Inner Lagrange point along +x: root of dphi_sub/dr = 0 in (1, a)."""
    # dphi/dr > 0 near the planet, -> -inf approaching the star; bracket the root.
    lo, hi = 1.0 + 1.0e-6, a*(1.0 - 1.0e-6)
    f_lo = roche_dphidr_sub(lo, q, a)
    # scan inward from the star for the sign change closest to the planet
    rs = np.linspace(lo, hi, 4000)
    fs = roche_dphidr_sub(rs, q, a)
    sign_change = np.where(np.sign(fs[:-1]) != np.sign(fs[1:]))[0]
    if len(sign_change) == 0:
        return hi
    i = sign_change[0]
    return brentq(roche_dphidr_sub, rs[i], rs[i+1], args=(q, a))


# --------------------------------------------------------------------------- #
#  Reconstruction map built from a substellar radial grid
# --------------------------------------------------------------------------- #
class ReconMap:
    """Equipotential map: phi_3D(x,y,z) -> equivalent substellar radius r_eff.

    Built from the EXHALE substellar radial grid r_sub (in R_p).  phi_sub is
    monotone increasing on (1, L1), so it inverts cleanly.  Points deeper than
    the base map to r_eff = r_sub[0] (clamped); points above the outermost grid
    radius (or outside the Roche lobe, phi > phi(L1)) are flagged outside
    (r_eff = NaN) so the caller can zero the density there.
    """

    def __init__(self, r_sub, q, a):
        self.q, self.a = q, a
        self.r_sub = np.asarray(r_sub, dtype=float)
        self.phi_sub = roche_phi_sub(self.r_sub, q, a)
        # enforce strict monotonicity for the inverse interpolation
        if not np.all(np.diff(self.phi_sub) > 0):
            keep = np.concatenate(([True], np.diff(self.phi_sub) > 0))
            self.r_sub = self.r_sub[keep]
            self.phi_sub = self.phi_sub[keep]
        self.L1 = find_L1(q, a)
        self.phi_L1 = roche_phi_sub(self.L1, q, a)
        self._inv = interp1d(self.phi_sub, self.r_sub, kind='linear',
                             bounds_error=False, fill_value=np.nan)
        self.phi_base = self.phi_sub[0]
        self.phi_top = self.phi_sub[-1]

    def r_eff(self, x, y, z):
        """Equivalent substellar radius at (x,y,z); NaN where outside."""
        phi = roche_phi(x, y, z, self.q, self.a)
        r = self._inv(phi)
        # clamp points deeper than the base to the base radius
        r = np.where(phi <= self.phi_base, self.r_sub[0], r)
        # anything above the computed atmosphere (or outside the lobe) -> NaN
        r = np.where(phi > self.phi_top, np.nan, r)
        return r

    def triaxial_radii(self, r_ref):
        """x/y/z extents of the equipotential passing through (r_ref,0,0).

        R_px = r_ref by construction (substellar).  R_py, R_pz solve
        phi(0,Ry,0)=phi_ref and phi(0,0,Rz)=phi_ref.  Returns (R_px,R_py,R_pz).
        """
        phi_ref = roche_phi_sub(r_ref, self.q, self.a)

        def fy(yy):
            return roche_phi(0.0, yy, 0.0, self.q, self.a) - phi_ref

        def fz(zz):
            return roche_phi(0.0, 0.0, zz, self.q, self.a) - phi_ref

        # phi(0,s,0) and phi(0,0,s) -> -1/s as s->0 (very negative) and rise
        # with s; bracket between a tiny radius and a few L1.
        hi = 3.0*self.L1
        Rpy = brentq(fy, 1.0e-4, hi)
        Rpz = brentq(fz, 1.0e-4, hi)
        return r_ref, Rpy, Rpz

    def roche_lobe_radius(self):
        """Volume-equivalent Roche-lobe radius via the Eggleton (1983) fit."""
        q = self.q  # M*/M_p ; Eggleton uses the mass ratio of the lobe-filling body
        # lobe of the planet: mass ratio Mp/M* = 1/q
        Q = 1.0/q
        rL_over_a = 0.49*Q**(2.0/3.0)/(0.6*Q**(2.0/3.0) + np.log(1.0 + Q**(1.0/3.0)))
        return rL_over_a*self.a  # [R_p]


# --------------------------------------------------------------------------- #
#  Self-test
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    # WASP-121b Case D-ish geometry (q = M*/M_p, a = orbital distance in R_p)
    # SI (metres, kilograms); RJ_CM is in cm, hence the factor.
    RJ = RJ_CM*1.0e-2
    Rsun, AU, Msun, MJ = 6.957e8, 1.495978707e11, 1.989e30, 1.898e27
    Rp = 2.581*RJ
    q = (1.3521*Msun)/(1.1204*MJ)
    a = (0.02544*AU)/Rp
    r_sub = np.array([1.0**( (1.0)**0 )] )  # placeholder, replaced below
    r_sub = np.array([1.0* (1.35)**(j/499.0) for j in range(500)])  # 1 -> 1.35 R_p

    print('q = M*/M_p = %.1f   a = %.2f R_p' % (q, a))
    m = ReconMap(r_sub, q, a)
    print('L1 = %.4f R_p   phi(L1) = %.4f' % (m.L1, m.phi_L1))
    print('phi_sub monotone increasing :', bool(np.all(np.diff(m.phi_sub) > 0)))

    # inversion round-trip
    rt = m._inv(roche_phi_sub(np.array([1.05, 1.15, 1.25]), q, a))
    print('inversion round-trip r_eff(phi_sub(r)) =', np.round(rt, 4),
          '(expect 1.05,1.15,1.25)')

    # triaxial radii + relations at a level near the lobe
    for rref in (1.05, 1.15, 1.25):
        Rx, Ry, Rz = m.triaxial_radii(rref)
        print('  r_ref=%.2f -> R_px=%.3f R_py=%.3f R_pz=%.3f   '
              'R_py*R_pz=%.3f (R_px^2=%.3f)  R_py/R_pz=%.3f'
              % (rref, Rx, Ry, Rz, Ry*Rz, Rx*Rx, Ry/Rz))
    rL = m.roche_lobe_radius()
    print('Roche-lobe radius R_RL = %.3f R_p   R_RL/L1 = %.3f (expect ~2/3)'
          % (rL, rL/m.L1))
