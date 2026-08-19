#!/usr/bin/env python3
"""Convert a Wind-AE windsoln.csv into EXHALE (EXHALE) IC files.

Writes output/Hydro_ioniz_IC.txt and output/Ion_species_IC.txt on the
exact EXHALE grid (taken from an IC_dump.txt of a fresh run of the same
input.inp), so 'Load IC? True' can ingest a Wind-AE steady solution.

Conventions:
- EXHALE rho = nH*mH*(1+4*HeH); load_IC rebuilds rho from the Ion file,
  v/p/T come from the Hydro file (cols 3,4,5), p is made EOS-consistent
  (p = (ntot+ne) kB T).
- Below Wind-AE's Rmin: hydrostatic-blend down to the EXHALE base anchor
  n_nuc(r=1)=10^lognbase, T->T0; v from mass-flux continuity.
- Above Wind-AE's Rmax: T,v frozen, rho ~ 1/(v r^2) (flux conservation).

Usage: windae_to_exhale_ic.py <windsoln.csv> <IC_dump.txt> <outdir>
          <lognbase> <T0[K]> <Rp[cm]> [HeH=0.083333333]
"""
import sys
import numpy as np

kB = 1.380649e-16
mH = 1.67353284e-24  # H atom, matching EXHALE's mu (was 1.6726e-24, the proton)

def read_windsoln(fn):
    scales = None
    rows = []
    cols = None
    for line in open(fn):
        if line.startswith('#scales:'):
            scales = [float(x) for x in line.split(':')[1].split(',')]
        elif line.startswith('#vars:'):
            cols = line.split(':')[1].strip().split(',')
        elif not line.startswith('#'):
            rows.append([float(x) for x in line.split(',')])
    d = np.array(rows)
    out = {c: d[:, i] for i, c in enumerate(cols)}
    # dimensionalize (scales align with the first len(scales) columns)
    for i, c in enumerate(cols[:len(scales)]):
        out[c] = out[c] * scales[i]
    return out

def main():
    wfile, gridfile, outdir, lognbase, T0, Rp = sys.argv[1:7]
    HeH = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0/12.0
    nbase = 10.0**float(lognbase)   # total H+He nuclei at r=1
    T0 = float(T0); Rp = float(Rp)

    w = read_windsoln(wfile)
    rw   = w['r'] / Rp                  # Wind-AE radius in EXHALE Rp units
    rhow = w['rho']                     # g/cm3
    vw   = w['v']                       # cm/s
    Tw   = w['T']                       # K
    fHIw  = np.clip(w['Ys_HI'],  0.0, 1.0)
    fHeIw = np.clip(w['Ys_HeI'], 0.0, 1.0)
    rmin_w, rmax_w = rw[0], rw[-1]

    # EXHALE grid (incl. ghosts) from IC_dump
    g = np.loadtxt(gridfile)
    r = g[:, 0]

    # --- nuclei density profile on the EXHALE grid ---
    # Wind-AE nuclei density rescaled to EXHALE composition:
    # rho = nH mH (1+4HeH)  =>  nH = rho/(mH (1+4HeH)); n_nuc = nH (1+HeH)
    nH_w = rhow / (mH * (1.0 + 4.0*HeH))
    nnuc_w = nH_w * (1.0 + HeH)

    lognnuc = np.interp(r, rw, np.log10(nnuc_w))
    T  = np.interp(r, rw, Tw)
    v  = np.interp(r, rw, vw)
    fHI  = np.interp(r, rw, fHIw)
    fHeI = np.interp(r, rw, fHeIw)

    # --- below Rmin: blend hydrostatic anchor (n(1)=nbase, T->T0) ---
    blo = r <= rmin_w
    if np.any(blo):
        # isothermal hydrostatic shape at T0 anchored at r=1 (uses the
        # local wind-ae-matched Jeans parameter for the exponent)
        # b0 = G Mp mu /(kB T0 Rp): recover from matching the windsoln
        # scale height just above Rmin instead of planet constants:
        # d ln n/d(1/r) ~ b0_eff
        i1 = np.searchsorted(rw, rmin_w + 0.1)
        b0_eff = ((np.log(nnuc_w[0]) - np.log(nnuc_w[i1]))
                  / (1.0/rw[0] - 1.0/rw[i1]))
        loghyd = np.log10(nbase) + b0_eff*(1.0/r[blo] - 1.0)/np.log(10.0)
        # blend exponent: 0 at r=1 -> 1 at rmin_w (linear in 1/r)
        s = (1.0/1.0 - 1.0/r[blo]) / (1.0/1.0 - 1.0/rmin_w)
        s = np.clip(s, 0.0, 1.0)
        logw_at = np.interp(rmin_w, rw, np.log10(nnuc_w))
        hyd_at  = np.log10(nbase) + b0_eff*(1.0/rmin_w - 1.0)/np.log(10.0)
        lognnuc[blo] = loghyd + (logw_at - hyd_at)*s
        T[blo] = T0 + (np.interp(rmin_w, rw, Tw) - T0)*s
        fHI[blo]  = 1.0 + (np.interp(rmin_w, rw, fHIw)  - 1.0)*s
        fHeI[blo] = 1.0 + (np.interp(rmin_w, rw, fHeIw) - 1.0)*s
        # velocity from mass-flux continuity (rho v r^2 = const at Rmin)
        flux0 = (np.interp(rmin_w, rw, np.log10(nnuc_w)), )
        nn_rmin = 10.0**np.interp(rmin_w, rw, np.log10(nnuc_w))
        v_rmin = np.interp(rmin_w, rw, vw)
        v[blo] = v_rmin * (nn_rmin * rmin_w**2) / (10.0**lognnuc[blo] * r[blo]**2)

    # --- above Rmax: frozen T,v; flux-conserving density ---
    bhi = r >= rmax_w
    if np.any(bhi):
        T[bhi] = Tw[-1]
        v[bhi] = vw[-1]
        lognnuc[bhi] = (np.log10(nnuc_w[-1])
                        + 2.0*np.log10(rmax_w) - 2.0*np.log10(r[bhi]))
        fHI[bhi]  = fHIw[-1]
        fHeI[bhi] = fHeIw[-1]

    nnuc = 10.0**lognnuc
    fHI  = np.clip(fHI, 0.0, 1.0)
    fHeI = np.clip(fHeI, 0.0, 1.0)

    nH  = nnuc / (1.0 + HeH)
    nHe = nH * HeH
    nHI    = fHI * nH
    nHII   = (1.0 - fHI) * nH
    nHeI   = fHeI * nHe
    nHeII  = (1.0 - fHeI) * nHe
    nHeIII = np.zeros_like(nHe)
    nHeTR  = np.zeros_like(nHe)

    ne   = nHII + nHeII + 2.0*nHeIII
    ntot = nH + nHe
    p = (ntot + ne) * kB * T
    nmass = nHI + nHII + 4.0*(nHeI + nHeII + nHeIII)   # writer col 2 analog

    with open(outdir + '/Hydro_ioniz_IC.txt', 'w') as f:
        f.write('# Wind-AE -> EXHALE IC (converted)\n')
        f.write('# columns r[Rp] n[cm-3] v[cm/s] p[cgs] T[K] '
                'heat[erg/cm3/s] cool[erg/cm3/s]\n')
        for j in range(len(r)):
            f.write(' %.10E %.10E %.10E %.10E %.10E %.10E %.10E\n'
                    % (r[j], nmass[j], v[j], p[j], T[j], 0.0, 0.0))

    with open(outdir + '/Ion_species_IC.txt', 'w') as f:
        f.write('# Wind-AE -> EXHALE IC (converted)\n')
        f.write('# columns r[Rp] HI HII HeI HeII HeIII HeITR\n')
        for j in range(len(r)):
            f.write(' %.10E %.10E %.10E %.10E %.10E %.10E %.10E\n'
                    % (r[j], nHI[j], nHII[j], nHeI[j], nHeII[j],
                       nHeIII[j], nHeTR[j]))

    print('grid points:', len(r), ' r:', r[0], '->', r[-1])
    print('Wind-AE coverage: [%.4f, %.4f] Rp; below=%d above=%d pts'
          % (rmin_w, rmax_w, blo.sum(), bhi.sum()))
    print('base anchor n_nuc(1)=%.3e  windae n_nuc(Rmin)=%.3e  b0_eff=%.1f'
          % (nbase, nnuc_w[0]*1.0, b0_eff if np.any(blo) else -1))

if __name__ == '__main__':
    main()
