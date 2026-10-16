#!/usr/bin/env python3
"""Run the k05 Photochem column (photochem_to_lower_profile.py, unchanged)
with the H2O condensate treatment set from the environment, and save the
final column state for the condensate diagnosis.

Environment:
  CD_RADIUS_CM     constant radius [cm] of every condensate particle
                   (unset: Photochem's default, 1e-5 cm)
  CD_KCOND_H2OAER  condensation rate coefficient A_cond of H2Oaer (unset:
                   Photochem's 100)
  CD_VDEP_H2OAER   deposition velocity [cm/s] added at the lower boundary of
                   H2Oaer, on top of its fall velocity (unset: 0, the
                   gas-giant extension's setting)
No file of src/ or of the Photochem installation is modified: the adapter
is imported as a module and the gas-giant class is wrapped at run time.
"""
import os
import sys

import numpy as np

SRC = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils'
sys.path.insert(0, SRC)
import photochem_to_lower_profile as ad                  # noqa: E402
from photochem.extensions import gasgiants               # noqa: E402

RADIUS = os.environ.get('CD_RADIUS_CM')
VDEP = os.environ.get('CD_VDEP_H2OAER')
KCOND = os.environ.get('CD_KCOND_H2OAER')

_orig_init_atm = gasgiants.EvoAtmosphereGasGiant._initialize_atmosphere


def _init_atm(self, *a, **k):
    _orig_init_atm(self, *a, **k)
    npart = int(self.dat.np)
    if RADIUS is not None:
        r = np.asarray(self.var.particle_radius).copy()
        r[:, :] = float(RADIUS)
        self.var.particle_radius = r
    if VDEP is not None:
        self.set_lower_bc('H2Oaer', bc_type='vdep', vdep=float(VDEP))
    if KCOND is not None:
        cp = self.var.cond_params[0]
        cp.k_cond = float(KCOND)
    print('  [condensate driver] H2Oaer k_cond %.3e' % self.var.cond_params[0].k_cond)
    print('  [condensate driver] particles %d, H2Oaer radius %.3e cm, '
          'H2Oaer extra lower vdep %s'
          % (npart, np.asarray(self.var.particle_radius)[0, 0], VDEP),
          flush=True)


gasgiants.EvoAtmosphereGasGiant._initialize_atmosphere = _init_atm

KB = 1.380649e-16
NAVO = 6.02214076e23


def dynamic_viscosity_air(T):
    # photochem_eqns.f90 dynamic_viscosity_air (Sutherland, modern Earth air)
    return 10.0*1.716e-5*(T/273.0)**1.5*(273.0 + 111.0)/(T + 111.0)


def slip(r, den):
    lam = 1.0/(den*6.0e-15)
    return 1.0 + (lam/r)*(1.257 + 0.4*np.exp(-1.1*r/lam))


def save_final(pc):
    names = list(pc.dat.species_names)
    npart = int(pc.dat.np)
    u = np.asarray(pc.wrk.usol)
    T = np.asarray(pc.var.temperature)
    z = np.asarray(pc.var.z)
    den = u[npart:].sum(0)
    mass = np.asarray(pc.dat.species_mass)
    mubar = (mass[npart:u.shape[0], None]*u[npart:]).sum(0)/den
    grav = np.asarray(pc.var.grav)
    r = np.asarray(pc.var.particle_radius)
    rho_air = den/NAVO*mubar
    w = (2.0/9.0)*grav*r[0]**2*(1.0 - rho_air)/dynamic_viscosity_air(T) \
        * slip(r[0], den)
    psat = np.array([pc.dat.particle_sat[0].sat_pressure(float(t))
                     for t in T])
    surf, top = pc.gas_fluxes()
    np.savez('final_state.npz', usol=u, names=np.array(names), np=npart,
             z=z, T=T, p=np.asarray(pc.wrk.pressure_hydro),
             den=den, mubar=mubar, grav=grav, edd=np.asarray(pc.var.edd),
             radius=r, wfall_H2Oaer=w, psat_H2O=psat,
             top_atmos=pc.var.top_atmos, tn=pc.wrk.tn,
             surf_names=np.array(list(surf.keys())),
             surf_flux=np.array(list(surf.values())),
             top_names=np.array(list(top.keys())),
             top_flux=np.array(list(top.values())),
             species_mass=mass)
    print('  [condensate driver] final state saved, tn = %.4e s'
          % pc.wrk.tn, flush=True)


_orig_measured = ad.measured_elemental_escape_flux


def _measured(pc, *a, **k):
    save_final(pc)
    return _orig_measured(pc, *a, **k)


ad.measured_elemental_escape_flux = _measured

if __name__ == '__main__':
    ad.main()
