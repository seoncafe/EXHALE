#!/usr/bin/env python3
"""Condensate diagnosis of a final Photochem column (final_state.npz written
by condensate_driver.py).  Usage: analyze.py <case_dir> [...]"""
import re
import sys

import numpy as np

NAVO = 6.02214076e23
KB = 1.380649e-16
MH2O = 18.01534
L_SUB_A, L_SUB_B = 27468840000.0, 4181527.0     # mechanism yaml, erg/g
P_FLOOR_BAR = 1.004                              # energy-domain floor (k05)


def load(case):
    d = np.load(case + '/final_state.npz', allow_pickle=True)
    return {k: d[k] for k in d.files}


def column_report(case):
    s = load(case)
    names = list(s['names'])
    u = s['usol']
    npart = int(s['np'])
    n = s['den']
    p = s['p']/1e6
    T = s['T']
    z = s['z']
    dz = np.gradient(z)
    iv = names.index('H2O')
    ia = names.index('H2Oaer')
    qv = u[iv]/n
    qa = u[ia]/n
    rh = u[iv]*KB*T/s['psat_H2O']
    rho_gas = n*s['mubar']/NAVO
    rho_c = u[ia]*MH2O/NAVO
    w = s['wfall_H2Oaer']
    out = dict(case=case, p=p, T=T, qv=qv, qa=qa, rh=rh,
               qNH3=u[names.index('NH3')]/n, qCH4=u[names.index('CH4')]/n,
               rho_gas=rho_gas, rho_c=rho_c, w=w, radius=s['radius'][0, 0])
    # elemental inventories [atoms cm^-2] with and without condensates
    comp = {}
    try:
        import yaml  # noqa: F401
    except Exception:
        pass
    # composition from species names via the mechanism yaml
    mech = case + '/photochem_work/zahnle_HHeNOC_rxns.yaml'
    import yaml
    y = yaml.safe_load(open(mech))
    for sp in y['species']:
        comp[sp['name']] = sp['composition']
    for sp in y.get('particles', []):
        comp[sp['name']] = sp['composition']
    inv = {}
    for el in ('H', 'O', 'C', 'N', 'He'):
        g = 0.0
        c = 0.0
        for i, sp in enumerate(names[:u.shape[0]]):
            k = comp.get(sp, {}).get(el, 0)
            if not k:
                continue
            col = (u[i]*dz).sum()*k
            if i < npart:
                c += col
            else:
                g += col
        inv[el] = (g, c)
    out['inventory'] = inv
    # water transport across faces: vapour eddy flux (mixing-ratio form,
    # Photochem's face averages) and particle fall flux (upwind)
    K = s['edd']
    Kf = np.sqrt(K[1:]*K[:-1])
    nf = np.sqrt(n[1:]*n[:-1])
    dzf = z[1:] - z[:-1]
    phi_v = -Kf*nf*(qv[1:] - qv[:-1])/dzf           # upward positive
    phi_a = -w[1:]*u[ia][1:]                          # face j+1/2, upwind
    pf = np.sqrt(p[1:]*p[:-1])
    out.update(pf=pf, phi_v=phi_v, phi_a=phi_a)
    # latent heat deposited above the energy-domain floor: the net
    # condensation above it is the condensate falling through it
    j = np.argmin(np.abs(pf - P_FLOOR_BAR))
    L = (L_SUB_A + L_SUB_B*T[j])*MH2O/NAVO        # erg / molecule
    out['floor'] = dict(p=pf[j], phi_v=phi_v[j], phi_a=phi_a[j],
                        latent=L*(-phi_a[j]), L=L)
    # bottom and top
    surf = dict(zip(list(s['surf_names']), s['surf_flux']))
    top = dict(zip(list(s['top_names']), s['top_flux']))
    out['bottom_H2O_flux'] = surf.get('H2O', np.nan)
    out['bottom_H2Oaer_loss'] = w[0]*u[ia][0]
    out['top'] = top
    out['tn'] = float(s['tn'])
    return out


def energy_record(case):
    rec = {}
    txt = open(case + '/conducted_heat_column.txt').read()
    m = re.search(r'# T_match_K (\S+)', txt)
    rec['T_match'] = float(m.group(1))
    m = re.search(r'emission excess (\S+), bottom flux (\S+), top flux (\S+),'
                  r' relative error (\S+)', txt)
    rec['excess'], rec['bottom'], rec['topflux'], rec['relerr'] = \
        [float(x.rstrip(',')) for x in m.groups()]
    rec['shares'] = re.findall(r'# band (\S+)\s+(.+?)\s+(lte|nlte)\s+share '
                               r'(\S+)', txt)
    m = re.search(r'thin_layer_max_absorbed_fraction (\S+)', txt)
    rec['absorbed_max'] = float(m.group(1)) if m else np.nan
    log = open(case + '/run.log').read()
    m = re.search(r'measured elemental escape at the model top.*?H (\S+) g/s,'
                  r' He (\S+) g/s', log)
    rec['top_H_gs'] = float(m.group(1).rstrip(',')) if m else np.nan
    rec['top_He_gs'] = float(m.group(2).rstrip(',')) if m else np.nan
    rec['steady'] = re.findall(r'(STEADY|GAVE UP)', log)
    return rec


if __name__ == '__main__':
    for case in sys.argv[1:]:
        o = column_report(case)
        print('=====', case, 'radius %.1e cm, tn %.3e s' % (o['radius'],
                                                          o['tn']))
        k = np.argmax(o['qa'])
        print('max condensate: p %.3f bar, T %.1f, n_c %.3e molec/cm3, '
              'q_c %.3e, rho_c %.3e g/cm3, rho_c/rho_gas %.3e, w %.3e cm/s'
              % (o['p'][k], o['T'][k], o['qa'][k]*o['rho_gas'][k]/o['rho_gas'][k]
                 * 0 + o['rho_c'][k]/(MH2O/NAVO), o['qa'][k], o['rho_c'][k],
                 o['rho_c'][k]/o['rho_gas'][k], o['w'][k]))
        for el, (g, c) in o['inventory'].items():
            print('inventory %-2s gas %.4e  condensed %.4e  total/gas %.4e'
                  % (el, g, c, (g + c)/g))
        f = o['floor']
        print('floor face %.3f bar: vapour flux %.3e, particle flux %.3e '
              'molec/cm2/s; latent heat above the floor %.3e erg/cm2/s'
              % (f['p'], f['phi_v'], f['phi_a'], f['latent']))
        print('bottom: H2O gas flux %.3e, H2Oaer loss %.3e molec/cm2/s'
              % (o['bottom_H2O_flux'], o['bottom_H2Oaer_loss']))
        e = energy_record(case)
        print('T_match %.4f K; budget excess %.6e bottom %.6e top %.6e; '
              'top H %.6e He %.6e g/s; %s'
              % (e['T_match'], e['excess'], e['bottom'], e['topflux'],
                 e['top_H_gs'], e['top_He_gs'], e['steady']))
        above = o['p'] < 0.6
        print('above 0.6 bar: q_H2O %.4e..%.4e, q_NH3 %.6e..%.6e, '
              'q_CH4 %.6e..%.6e'
              % (o['qv'][above].min(), o['qv'][above].max(),
                 o['qNH3'][above].min(), o['qNH3'][above].max(),
                 o['qCH4'][above].min(), o['qCH4'][above].max()))
