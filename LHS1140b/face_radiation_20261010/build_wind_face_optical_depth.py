#!/usr/bin/env python3
"""STAND-IN for the wind's export of its optical depth at the matching face.

Builds `wind_face_optical_depth.txt` (format: md/matching_face_radiation_
python_20261010.md section 3) from the STORED outputs of a wind run, offline,
because the Fortran export does not exist yet.  The file it writes is labeled
`producer stand-in` in its header; the adapter reads it exactly as it will
read the wind's own file.

What it computes, and how it differs from what the wind will write:
  * radial columns N_X = sum_cells n_X dr from the face R0 to the last
    physical cell, from output/Ion_species.txt, cell faces at R0 and at the
    midpoints of the stored cell centres (the wind sums its own cell widths
    dr_j R0 and includes its outer ghost cells, calc_column_dens of
    utilities.f90; the difference is the stand-in's);
  * the continuum optical depth of every absorber the wind attenuates its
    XUV beam with (util_ion_eq.f90, photoionization_field: H I, He I, He II,
    He 2^3S, H2, and every photo-absorbing metal stage), with the wind's own
    cross sections, evaluated by the driver
    photoabsorption_cross_sections_on_grid.f90 (Cross_sections of
    src/modules/functions/cross_sec.f90, opacity model 'A');
  * the bins are the cells of the SED file's own wavelength samples (the
    wind exports on its photon grid, set_energy_vectors.f90);
  * the H2 Lyman-Werner lines are not a continuum and are not in tau_cont;
    the header carries the radial H2 column and its column-weighted
    temperature (output/Lyman_Werner.txt), from which the reader evaluates
    the lines at its own slant.

Usage:
  build_wind_face_optical_depth.py energies <run_dir> <energies.txt>
  build_wind_face_optical_depth.py build <run_dir> <xs_table> <out_file>
where xs_table is the driver's output for energies.txt.
"""
import hashlib
import os
import sys

import numpy as np

HC_EV_NM = 1239.841984            # h c [eV nm] (CODATA)
LAMBDA_MAX_NM = 120.1             # the face file's long-wavelength end
R_J = 7.1492e9                    # cm, as EXHALE's parameters.f90

# Cross_sections column order of the 17 photo-absorbing metal stages
# (cross_sec.f90 header, = species_table's mion_iphot) and their names in
# the Ion_species.txt header.
METAL_STAGES = ['CI', 'CII', 'OI', 'OII', 'NI', 'NII', 'MgI', 'MgII', 'SiI',
                'SiII', 'CaI', 'CaII', 'NaI', 'KI', 'SI', 'FeI', 'FeII']
ABSORBERS = ['HI', 'HeI', 'HeII', 'HeITR', 'H2'] + METAL_STAGES
H2_LAYERS = 5


def sed_cells(run_dir):
    """Wavelength samples of the SED file Photochem reads, and the cells
    (edges at the midpoints) that hold them, below LAMBDA_MAX_NM."""
    wl = np.loadtxt(os.path.join(run_dir, 'photochem_work',
                                 'flux_photochem.txt'), skiprows=1)[:, 0]
    edges = np.concatenate([[wl[0] - 0.5*(wl[1] - wl[0])],
                            0.5*(wl[1:] + wl[:-1]),
                            [wl[-1] + 0.5*(wl[-1] - wl[-2])]])
    keep = edges[1:] <= LAMBDA_MAX_NM + 1e-9
    lo, hi = edges[:-1][keep], edges[1:][keep]
    return wl[keep], lo, hi


def resolved(run_dir):
    out = {}
    for line in open(os.path.join(run_dir, 'EXHALE_resolved.out')):
        t = line.split()
        if len(t) >= 2 and not line.startswith('#'):
            out[t[0]] = t[1]
    return out


def radial_columns(run_dir):
    res = resolved(run_dir)
    Rp = float(res['planet_radius_RJ'])*R_J
    R0 = float(res['base_face_radius_cm'])
    path = os.path.join(run_dir, 'output', 'Ion_species.txt')
    hdr = [l for l in open(path) if l.startswith('# columns')][0].split()[2:]
    ion = np.loadtxt(path)[2:-2]                     # physical cells
    r = ion[:, 0]*Rp
    f = np.empty(r.size + 1)
    f[0] = R0
    f[1:-1] = 0.5*(r[1:] + r[:-1])
    f[-1] = 2.0*r[-1] - f[-2]
    dr = np.diff(f)
    N = {k: float(np.sum(ion[:, hdr.index(k)]*dr)) for k in hdr[1:]}
    # H2 Lyman-Werner: the radial H2 column split into H2_LAYERS layers of
    # equal column, outward from the face, each with its column-weighted
    # temperature (the LTE populations of a layer are those of its T)
    lw = np.loadtxt(os.path.join(run_dir, 'output', 'Lyman_Werner.txt'))[2:-2]
    w = lw[:, 3]*dr
    cum = np.cumsum(w)/np.sum(w)
    layers = []
    lo = 0
    for k in range(1, H2_LAYERS + 1):
        hi = (int(np.searchsorted(cum, k/H2_LAYERS - 1e-12)) + 1
              if k < H2_LAYERS else w.size)
        sl = slice(lo, max(hi, lo + 1))
        Nk = float(np.sum(w[sl]))
        Tk = float(np.sum(w[sl]*lw[sl, 1])/Nk)
        layers.append((Nk, Tk))
        lo = sl.stop
    if lo < w.size:
        Nk = float(np.sum(w[lo:]))
        Tk = float(np.sum(w[lo:]*lw[lo:, 1])/max(Nk, 1e-300))
        layers.append((Nk, Tk))
    return N, layers, R0, f[-1], res


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        h.update(fh.read())
    return h.hexdigest()


def main():
    mode, run_dir = sys.argv[1], sys.argv[2]
    wl, lo, hi = sed_cells(run_dir)
    E = HC_EV_NM/wl
    if mode == 'energies':
        with open(sys.argv[3], 'w') as fh:
            fh.write('%d\n' % E.size)
            for e in E:
                fh.write('%.10e\n' % e)
        return
    xs_path, out = sys.argv[3], sys.argv[4]
    xs = np.loadtxt(xs_path)
    if xs.shape[0] != E.size or np.max(np.abs(xs[:, 0]/E - 1.0)) > 1e-7:
        raise SystemExit('the cross-section table is not on the SED cells')
    sig = {k: xs[:, 1 + i]*1.0e-18 for i, k in enumerate(ABSORBERS)}
    N, layers, R0, r_out, res = radial_columns(run_dir)
    tau = np.zeros_like(wl)
    parts = {}
    for k in ABSORBERS:
        if k in N and N[k] > 0.0:
            parts[k] = sig[k]*N[k]
            tau += parts[k]
    sed = os.path.join(run_dir, 'photochem_work', 'flux_photochem.txt')
    with open(out, 'w') as fh:
        fh.write('# EXHALE wind_face_optical_depth 1\n')
        fh.write('# producer stand-in: build_wind_face_optical_depth.py from'
                 ' the stored outputs of %s (NOT the wind\'s export)\n'
                 % os.path.abspath(run_dir))
        fh.write('# lower_profile_solution_id %s\n'
                 % res.get('lower_profile_solution_id', 'unknown'))
        fh.write('# face_radius_cm %.10E\n' % R0)
        fh.write('# outer_radius_cm %.10E\n' % r_out)
        fh.write('# path radial, from the outer boundary to the face\n')
        fh.write('# dayside_dilution %.6f\n' % 1.0)
        fh.write('# spectrum_file %s sha256 %s\n' % (os.path.basename(sed),
                                                      sha256(sed)[:16]))
        fh.write('# cross_sections EXHALE Cross_sections (cross_sec.f90),'
                 ' opacity model A\n')
        fh.write('# absorbers %s\n' % ' '.join(sorted(parts)))
        for k in ABSORBERS:
            if k in N:
                fh.write('# N_%s_cm2 %.10E\n' % (k, N[k]))
        for Nk, Tk in layers:
            fh.write('# lyman_werner_lines H2 N_H2_cm2 %.10E T_H2_K %.4f\n'
                     % (Nk, Tk))
        fh.write('# columns: lambda_lo[nm] lambda_hi[nm] tau_cont\n')
        for a, b, t in zip(lo, hi, tau):
            fh.write('%.8E %.8E %.10E\n' % (a, b, t))
    print('wrote %s: %d bins %.3f-%.3f nm; N_HeI %.3e N_H2 %.3e; H2 layers %s'
          % (out, lo.size, lo[0], hi[-1], N['HeI'], N['H2'],
             ', '.join('%.3e at %.0f K' % l for l in layers)))
    for lam in (1.0, 2.0, 5.0, 10.0, 30.0, 50.0, 80.0, 90.0, 100.0, 110.0):
        i = int(np.argmin(np.abs(wl - lam)))
        top = sorted(parts, key=lambda k: -parts[k][i])[:3]
        print('  tau_cont(%6.2f nm) = %.4e  (%s)' % (
            wl[i], tau[i], ', '.join('%s %.2e' % (k, parts[k][i])
                                     for k in top)))


if __name__ == '__main__':
    main()
