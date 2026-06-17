#!/usr/bin/env python3
"""tpm_halpha_lart2d.py -- H-alpha transit transmission from a 2D (cylindrical)
Ly-alpha scattering rate computed by LaRT (spherical-illumination model).

Why this is separate from TPM.py
--------------------------------
TPM.py assumes a spherically symmetric atmosphere, n_2p = n_2p(r).  Under
stellar (spherical) illumination the Ly-alpha mean intensity -- and therefore the
H(2p) population that produces H-alpha -- is only *cylindrically* symmetric about
the star-planet axis, n_2p = n_2p(rho, z).  This script builds n_2p(rho,z) from the
LaRT scattering rate P_alpha(rho,z) and integrates the H-alpha optical depth along
the transit line of sight (= the star-planet axis = z), so the spherical assumption
is dropped.  He 10830 (set by the spherically symmetric He 2^3S ionization balance)
is unaffected and is still handled by TPM.py.

Pipeline
--------
  EXHALE run  -> n_HI(r), n_e(r), T(r), v(r)         (spherical profiles)
  LaRT h5     -> Pa_2D[iz, irho]  (scattering rate, *unit luminosity*)
  P_alpha(rho,z) = Pa_2D * L_lya            L_lya = F_lya*4pi a^2 / (h nu_lya)  [photons/s]
  J_lya,eff(rho,z) = P_alpha / B12_lya      (so B12*J = P_alpha = the 1s->2p pump)
  n_2p(rho,z) = Christie+2013 rate equilibrium (same as TPM.n2_populations)
  tau(b,nu)   = int n_2p(rho=b, z) sigma_Ha(nu, v_LOS, T) dz,   v_LOS = v(r) z/r
  disk-average over impact parameter b -> H-alpha transmission spectrum

Usage
-----
  python tpm_halpha_lart2d.py <exhale_run_dir> <lart_h5> [--obs FILE] [--out PNG]
"""

import os
import sys
import argparse
import numpy as np
from scipy.special import wofz
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exhale_io as aio   # noqa: E402

# ---- constants (SI, matching TPM.py) ----
kb   = 1.380649e-23
mp   = 1.672623e-27
me   = 9.109384e-31
c_l  = 2.99792458e8
E0   = 8.854188e-12
echg = -1.602176e-19
AU   = 1.495978707e11
hpl  = 6.62607015e-34
hp_eV = 4.135667696e-15

# H-alpha (n=2 -> n=3), air wavelength
l_Ha   = 6562.8e-10
nu_Ha  = c_l / l_Ha
f_Ha   = 0.6407
A12_Ha = 4.4101e7
Fadd_const = np.sqrt(np.pi) * echg**2 / (4.0 * np.pi * E0 * me * c_l)

# Ly-alpha 1s<->2p pumping (cgs B coefficients, as in TPM.py)
lA      = 1215.6701e-10
nu_Lya  = c_l / lA
A_2p1s  = 6.3e8
g1s, g2s, g2p = 2.0, 2.0, 6.0
c_cgs, h_cgs, kb_cgs = 2.99792458e10, 6.62607015e-27, 1.380649e-16
eV2Hz = 2.417989242e14
B21_lya = A_2p1s * c_cgs**2 / (2.0 * h_cgs * nu_Lya**3)
B12_lya = (g2p / g1s) * B21_lya
A_2s1s  = 8.26


def gamma_n2_balmer(T_star, R_over_a):
    """n=2 photoionization rate [s^-1] from a diluted stellar blackbody Balmer
    continuum (E > 3.4 eV).  Same as TPM.gamma_n2_balmer.  R_over_a = R_star/a."""
    if T_star <= 0.0:
        return 0.0
    E2   = 3.40
    nu2  = E2 * eV2Hz
    s2th = 1.4e-17
    Eg   = np.linspace(E2, 13.6, 400)
    nu   = Eg * eV2Hz
    Bnu  = (2.0 * h_cgs * nu**3 / c_cgs**2) / (np.exp(h_cgs * nu / (kb_cgs * T_star)) - 1.0)
    Fnu  = np.pi * Bnu * R_over_a**2
    sig2 = s2th * (nu2 / nu)**3
    return np.trapz(Fnu / (h_cgs * nu) * sig2, nu)


def n2_populations(T, n1s, ne, Jlya, G2s=0.0, G2p=0.0):
    """Christie+2013 2s/2p rate equilibrium -> (n2s, n2p, n2tot) [cm^-3].
    G2s, G2p are the n=2 Balmer-continuum photoionization sinks [s^-1]."""
    T  = np.maximum(T, 1.0)
    t4 = T / 1.0e4
    aB  = 2.54e-13 * t4**(-0.8163 - 0.0208 * np.log(t4))
    a2s = (0.282 + 0.047 * t4 - 0.006 * t4**2) * aB
    a2p = aB - a2s
    C1s2s = 1.21e-8 * (1.0 / t4)**0.455 * np.exp(-118400.0 / T)
    C1s2p = 1.71e-8 * (1.0 / t4)**0.077 * np.exp(-118400.0 / T)
    C2s2p = 6.21e-5 * (np.log(T / 1.02) - 0.57721) / np.sqrt(T)
    C2s1s = 1.21e-8 * (1.0 / t4)**0.455 * (g1s / g2s)
    C2p1s = 1.71e-8 * (1.0 / t4)**0.077 * (g1s / g2p)
    C2p2s = C2s2p * (g2s / g2p)
    Ppump = B12_lya * Jlya
    Pstim = B21_lya * Jlya
    L2p = A_2p1s + Pstim + (C2p1s + C2p2s) * ne + G2p
    L2s = (C2s1s + C2s2p) * ne + A_2s1s + G2s
    S2p = (Ppump + C1s2p * ne) * n1s + a2p * ne**2
    S2s = (C1s2s * ne) * n1s + a2s * ne**2
    M12 = C2s2p * ne
    M21 = C2p2s * ne
    det = L2p * L2s - M12 * M21
    det = np.where(np.abs(det) > 0.0, det, 1.0)
    n2p = np.maximum((S2p * L2s + M12 * S2s) / det, 0.0)
    n2s = np.maximum((L2p * S2s + M21 * S2p) / det, 0.0)
    return n2s, n2p, n2s + n2p


def read_lart_pa2d(h5file):
    """Return (Pa[iz,irho], z_grid[Rp], rho_grid[Rp]) for the cylindrical P_alpha
    (unit-luminosity).  Axes reconstructed from the LaRT cylindrical binning:
      z_center(iz)  = -zmax + (iz+0.5)*dz,    dz  = 2*zmax/nz
      rho_center(j) = j*dr,                   dr  = rmax/(nr-0.5)   (nr odd)."""
    with h5py.File(h5file, 'r') as f:
        Pa = f['Pa_2D/data'][:]                 # shape (nz, nr)  [Fortran (nr,nz) reversed]
        a = dict(f['Spectrum'].attrs)
    nz, nr = Pa.shape
    zmax = float(a['zmax']); rmax = float(a.get('xmax', zmax))
    dz = 2.0 * zmax / nz
    z = -zmax + (np.arange(nz) + 0.5) * dz
    dr = rmax / (nr - 0.5) if (nr % 2 == 1) else rmax / nr
    rho = np.arange(nr) * dr
    return Pa, z, rho, a


def build_n2p_2d(run, Pa, z, rho, L_lya, fluxfac=1.0, T_star=0.0, R_over_a=0.0):
    """n2p(rho,z) [cm^-3] from the physical scattering rate.

    P_alpha = Pa * L_lya * fluxfac.  The LaRT scattering-rate accumulation
    (add_to_Pnew) uses photon%wgt (=1 for stellar illumination), NOT
    photon%flux_factor, so the geometric solid-angle dilution Omega_star/4pi is
    *not* in Pa; it is written to the header as `fluxfac` and must be applied
    here (in addition to the physical Ly-alpha photon luminosity L_lya)."""
    r = np.asarray(run.r, float)
    nHI = np.asarray(run.ion['HI'], float)
    ne  = np.asarray(run.ion['HII'], float) + np.asarray(run.ion['HeII'], float) \
          + 2.0 * np.asarray(run.ion['HeIII'], float)
    T   = np.asarray(run.T, float)
    rmaxp = r.max()
    RR, ZZ = np.meshgrid(rho, z)              # (nz, nr), matches Pa[iz,irho]
    rad = np.sqrt(RR**2 + ZZ**2)              # spherical radius at each (rho,z)
    inside = rad <= rmaxp
    # interpolate spherical EXHALE profiles onto rad
    T_g   = np.interp(rad, r, T,   left=T[0],   right=T[-1])
    nHI_g = np.interp(rad, r, nHI, left=nHI[0], right=0.0)
    ne_g  = np.interp(rad, r, ne,  left=ne[0],  right=0.0)
    Palpha = Pa * L_lya * fluxfac              # [s^-1 atom^-1]  (fluxfac = Omega_star/4pi)
    Jlya_eff = Palpha / B12_lya
    G2 = gamma_n2_balmer(T_star, R_over_a)     # n=2 Balmer-continuum photoionization sink
    _, n2p, _ = n2_populations(T_g, nHI_g, ne_g, Jlya_eff, G2s=G2, G2p=G2)
    n2p = np.where(inside & (Pa > 0.0), n2p, 0.0)
    return n2p, rad


def halpha_transmission(run, Pa, z, rho, n2p, Rp_m, Rstar_m, lam_grid,
                        rib=None, nb=200, nz_los=801):
    """Disk-averaged H-alpha transmission spectrum.  Returns (lam[A], T_lambda)."""
    r = np.asarray(run.r, float)
    Tprof = np.asarray(run.T, float)
    vprof = np.asarray(run.v, float) / 100.0   # cm/s -> m/s
    rmaxp = r.max()
    rib = rib or rmaxp
    # 2D n2p interpolator (rho, z) -> n2p [cm^-3]
    from scipy.interpolate import RegularGridInterpolator
    n2p_interp = RegularGridInterpolator((z, rho), n2p, bounds_error=False, fill_value=0.0)

    b_grid = np.array([rib**(j / (nb - 1.0)) for j in range(nb)])   # 1..rib, log-spaced
    z_los  = np.linspace(-rib, rib, nz_los)                          # LOS coordinate [Rp]
    dz_m   = (z_los[1] - z_los[0]) * Rp_m
    nu_arr = c_l / (lam_grid * 1e-10)

    exp_tau = np.ones((nb, len(lam_grid)))
    for ip, b in enumerate(b_grid):
        rad = np.sqrt(b**2 + z_los**2)                  # spherical r along the LOS
        m = rad <= rmaxp
        if not m.any():
            continue
        zl = z_los[m]; radl = rad[m]
        T_l   = np.interp(radl, r, Tprof, left=Tprof[0], right=Tprof[-1])
        v_l   = np.interp(radl, r, vprof, left=vprof[0], right=vprof[-1])
        v_los = v_l * zl / radl                         # radial velocity projected on LOS
        n2_l  = n2p_interp(np.column_stack([zl, np.full_like(zl, b)])) * 1.0e6  # cm^-3 -> m^-3
        v_th  = np.sqrt(2.0 * kb * np.maximum(T_l, 1.0) / mp)
        Dnu   = nu_Ha * v_th / c_l
        a_v   = A12_Ha / (4.0 * np.pi * Dnu)
        for il, nu in enumerate(nu_arr):
            X    = (nu - nu_Ha) / Dnu
            arg  = X - v_los / v_th + 1j * a_v
            sig  = f_Ha * Fadd_const / Dnu * wofz(arg).real     # [m^2]
            integrand = n2_l * sig                              # [m^-1]
            tau = np.sum(0.5 * dz_m * (integrand[:-1] + integrand[1:]))
            exp_tau[ip, il] = np.exp(-tau)

    # disk-average (TPM.py convention)
    A_star   = np.pi * Rstar_m**2
    A_planet = np.pi * Rp_m**2
    A_atm    = np.pi * (rib * Rp_m)**2
    Tl = np.empty(len(lam_grid))
    for il in range(len(lam_grid)):
        prob = np.trapz(2.0 * exp_tau[:, il] * b_grid, b_grid) * A_planet / (A_atm - A_planet)
        avg  = ((A_star - A_atm) + (A_atm - A_planet) * prob) / A_star
        Tl[il] = avg * A_star / (A_star - A_planet)
    return lam_grid, Tl


def lya_photon_luminosity(run):
    """Stellar Ly-alpha photon luminosity [photons/s] from input.inp's
    'Stellar Lya flux [erg/cm2/s]' and the orbital distance: L = F*4pi a^2 / (h nu)."""
    raw   = run.inp.get('raw', {})
    F_lya = float(raw.get('Stellar Lya flux [erg/cm2/s]', '0').split()[0])  # erg/cm2/s
    a_cm  = run.inp['a_AU'] * AU * 100.0
    hnu   = 10.20 * 1.602176634e-12   # Ly-alpha photon energy [erg]
    return F_lya * 4.0 * np.pi * a_cm**2 / hnu


def compute_halpha(run_dir, lart_h5, lam_min=6561.0, lam_max=6564.6, nlam=201, adv=True,
                   lya_scale=1.0):
    """Full pipeline: EXHALE run + LaRT Pa_2D -> H-alpha transmission spectrum.
    lya_scale multiplies the incident Ly-alpha flux (hence P_alpha and the pump).
    Returns (lam[A], T_lambda, info-dict)."""
    run = aio.load_run(os.path.join(run_dir, 'output'),
                       os.path.join(run_dir, 'input.inp'), adv=adv)
    Rp_m    = run.inp['Rp_RJ'] * 6.9911e7
    raw     = run.inp.get('raw', {})
    Rstar_m = float(raw.get('Stellar radius [R_sun]', '1').split()[0]) * 6.96e8
    T_star  = float(raw.get('Stellar Teff [K]', '0').split()[0])
    a_m     = run.inp['a_AU'] * AU
    R_over_a = Rstar_m / a_m
    L_lya   = lya_photon_luminosity(run) * lya_scale
    Pa, z, rho, attrs = read_lart_pa2d(lart_h5)
    fluxfac = float(np.ravel(attrs.get('fluxfac', [1.0]))[0])
    n2p, _  = build_n2p_2d(run, Pa, z, rho, L_lya, fluxfac=fluxfac,
                           T_star=T_star, R_over_a=R_over_a)
    lam = np.linspace(lam_min, lam_max, nlam)
    lam, Tl = halpha_transmission(run, Pa, z, rho, n2p, Rp_m, Rstar_m, lam)
    info = dict(L_lya=L_lya, fluxfac=fluxfac, Palpha_max=(Pa * L_lya * fluxfac).max(),
                n2p_max=n2p.max(), depth_pct=(1.0 - Tl.min()) * 100.0)
    return lam, Tl, info


def main():
    ap = argparse.ArgumentParser(description='2D (cylindrical) H-alpha transit from a LaRT '
                                             'spherical-illumination scattering-rate output.')
    ap.add_argument('run_dir', help='EXHALE run directory (input.inp + output/)')
    ap.add_argument('lart_h5', help='LaRT output .h5 with Pa_2D')
    ap.add_argument('--lam_min', type=float, default=6561.0)
    ap.add_argument('--lam_max', type=float, default=6564.6)
    ap.add_argument('--nlam', type=int, default=201)
    ap.add_argument('--obs', default=None, help='observed Halpha file (lam, TS, e_TS, ...)')
    ap.add_argument('--out', default='halpha_lart2d.png')
    ap.add_argument('--eq', action='store_true')
    args = ap.parse_args()

    run = aio.load_run(os.path.join(args.run_dir, 'output'),
                       os.path.join(args.run_dir, 'input.inp'), adv=not args.eq)
    Rp_m    = run.inp['Rp_RJ'] * 6.9911e7
    a_m     = run.inp['a_AU'] * AU
    raw     = run.inp.get('raw', {})
    Rstar_m = float(raw.get('Stellar radius [R_sun]', '1').split()[0]) * 6.96e8
    F_lya   = float(raw.get('Stellar Lya flux [erg/cm2/s]', '0').split()[0])  # erg/cm2/s
    hnu_lya = hp_eV * 10.20 * 1.602176634e-12 / hp_eV  # = 10.2 eV in erg
    hnu_lya = 10.20 * 1.602176634e-12                  # [erg]
    L_lya   = F_lya * 4.0 * np.pi * (a_m * 100.0)**2 / hnu_lya   # photons/s (a in cm)

    Pa, z, rho, attrs = read_lart_pa2d(args.lart_h5)
    fluxfac = float(np.ravel(attrs.get('fluxfac', [1.0]))[0])
    n2p, rad = build_n2p_2d(run, Pa, z, rho, L_lya, fluxfac=fluxfac)
    print('L_lya = %.3e photons/s ;  fluxfac = %.4e ;  P_alpha max = %.3e s^-1 ;  n2p max = %.3e cm^-3'
          % (L_lya, fluxfac, (Pa * L_lya * fluxfac).max(), n2p.max()))

    lam = np.linspace(args.lam_min, args.lam_max, args.nlam)
    lam, Tl = halpha_transmission(run, Pa, z, rho, n2p, Rp_m, Rstar_m, lam)
    depth = (1.0 - Tl.min()) * 100.0
    print('H-alpha line-center excess absorption = %.3f %%' % depth)

    # plot
    import matplotlib
    matplotlib.use('Agg'); matplotlib.rcParams['text.usetex'] = False
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(lam, (1.0 - Tl) * 100.0, 'C3-', lw=2, label='LaRT 2D model (%.2f%%)' % depth)
    if args.obs and os.path.exists(args.obs):
        o = np.loadtxt(args.obs)
        ax.errorbar(o[:, 0], -o[:, 1] * 100.0, yerr=o[:, 2] * 100.0, fmt='o', ms=3,
                    color='0.4', ecolor='0.7', capsize=2, label='obs')
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set_xlabel('wavelength [A] (air)'); ax.set_ylabel('excess absorption [%]')
    ax.set_title('WASP-52b H-alpha: LaRT 2D (cylindrical Ly-a) coupled')
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(args.out, dpi=140)
    print('saved', args.out)


if __name__ == '__main__':
    main()
