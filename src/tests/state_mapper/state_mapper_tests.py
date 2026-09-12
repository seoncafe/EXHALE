#!/usr/bin/env python3
"""Assertions on map_state_to_grid.py on small manufactured states."""
import os, subprocess, sys
import numpy as np

TOOL, WORK = sys.argv[1], sys.argv[2]
n_fail = 0

def verdict(name, measured, reference, tol, ok):
    global n_fail
    print(f"{'PASS' if ok else 'FAIL'} {name} measured={measured} reference={reference} tol={tol}")
    if not ok: n_fail += 1

def check_rel(name, m, r, tol):
    verdict(name, m, r, tol, abs(m - r) <= tol*max(abs(r), 1e-300))

def check_eq(name, m, r):
    verdict(name, m, r, 0, m == r)

HCOLS = 'r[Rp] rho[mH/cm3] v[cm/s] p[cgs] T[K] heat[erg/cm3/s] cool[erg/cm3/s]'.split()
SCOLS = 'r[Rp] HI HII HeI HeII HeIII HeITR H2 H2p H3p HeHp'.split()

def write_state(d, r, rho, hehr=0.0793, coupling='mode=phys t_phys=123.0 certified=T', ic=False, r_species=None):
    os.makedirs(d, exist_ok=True)
    n = r.size
    v = 1e5*(r - 1.0); T = 1000.0*(1 + 0.5*(r - 1)); p = rho*T*1e-6
    heat = 1e-3*rho; cool = 2e-3*rho
    x_hp = 0.1*(r - 1)/(r.max() - 1 + 1e-30); nH = rho/1.1
    hi = nH*(1 - x_hp)*0.5; h2 = nH*(1 - x_hp)*0.25; hii = nH*x_hp
    hei = hehr*nH; he_rest = np.full(n, 1e-300)
    def hdr(cols):
        return ['# EXHALE schema 2', '# columns ' + ' '.join(cols),
                f'# rows {n}: 2 ghost cells at each end; physical cells are rows 3 to {n-2}',
                '# coupling: ' + coupling, '# provenance: git=test tree=clean run=test']
    rs = r if r_species is None else r_species
    with open(os.path.join(d, 'Hydro_ioniz' + ('_IC' if ic else '') + '.txt'), 'w') as f:
        f.write('\n'.join(hdr(HCOLS)) + '\n')
        for row in zip(r, rho, v, p, T, heat, cool):
            f.write(' '.join(f'{x:24.16E}' for x in row) + '\n')
    with open(os.path.join(d, 'Ion_species' + ('_IC' if ic else '') + '.txt'), 'w') as f:
        f.write('\n'.join(hdr(SCOLS)) + '\n')
        for row in zip(rs, hi, hii, hei, he_rest, he_rest, he_rest, h2, he_rest, he_rest, he_rest):
            f.write(' '.join(f'{x:24.16E}' for x in row) + '\n')

def read(path):
    hdr = [l.rstrip('\n') for l in open(path) if l.startswith('#')]
    a = np.loadtxt(path)
    return hdr, a

def run(src, tgt, out):
    p = subprocess.run([sys.executable, TOOL, src, tgt, out], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr

# ---- source: 24 rows (2 ghosts each end), geometric radii from 0.999 to 3
n = 400
r_src = 1.0*np.exp(np.linspace(np.log(0.9990), np.log(3.0), n))
rho_src = 1e12*np.exp(-6*(r_src - 1))
src = os.path.join(WORK, 'src'); write_state(src, r_src, rho_src)

# 1. identity mapping reproduces every value
out = os.path.join(WORK, 'identity'); rc, log = run(src, os.path.join(src, 'Hydro_ioniz.txt'), out)
check_eq('identity_mapping_exit_status', rc, 0)
if rc == 0:
    for fn in ('Hydro_ioniz', 'Ion_species'):
        a = read(os.path.join(src, fn + '.txt'))[1]; b = read(os.path.join(out, fn + '.txt'))[1]
        m = np.max(np.abs(b - a)/np.maximum(np.abs(a), 1e-300))
        check_rel(f'identity_mapping_{fn}_max_relative_change', m, 0.0, 1e-13) if False else verdict(f'identity_mapping_{fn}_max_relative_change', m, 0.0, 1e-12, m <= 1e-12)
    hdr = read(os.path.join(out, 'Hydro_ioniz.txt'))[0]
    coup = [h for h in hdr if h.startswith('# coupling:')][0]
    check_eq('seed_metadata_mode_init', 'mode=init' in coup, True)
    check_eq('seed_metadata_t_phys_zero', 't_phys=0' in coup, True)
    check_eq('seed_metadata_not_certified', 'certified=F' in coup, True)
    check_eq('source_coupling_kept_as_provenance', any(h.startswith('# mapped-from-coupling: mode=phys t_phys=123.0 certified=T') for h in hdr), True)
    check_eq('mapping_line_present', any(h.startswith('# mapped:') for h in hdr), True)

# 2. a grid shifted by 2e-4 in the centers: round trip and He/H invariance
r_tgt = r_src*(1 + 2e-4*np.sin(np.linspace(0, 3, n))); r_tgt[[0, 1, -2, -1]] = r_src[[0, 1, -2, -1]]
tgt = os.path.join(WORK, 'tgt'); write_state(tgt, r_tgt, rho_src)
out = os.path.join(WORK, 'shifted'); rc, log = run(src, os.path.join(tgt, 'Hydro_ioniz.txt'), out)
check_eq('shifted_grid_exit_status', rc, 0)
if rc == 0:
    b = read(os.path.join(out, 'Ion_species.txt'))[1]
    nH = b[:, 1] + b[:, 2] + 2*b[:, 7]; nHe = b[:, 3] + b[:, 4] + b[:, 5]
    ratio = nHe/nH
    verdict('he_to_h_ratio_invariant_after_mapping', float(np.max(np.abs(ratio - 0.0793))), 0.0793, 1e-12, np.max(np.abs(ratio - 0.0793)) <= 1e-12)
    h = read(os.path.join(out, 'Hydro_ioniz.txt'))[1]
    back = np.interp(np.log(r_src[2:-2]), np.log(h[2:-2, 0]), np.log(h[2:-2, 1]))
    m = np.max(np.abs(np.exp(back) - rho_src[2:-2])/rho_src[2:-2])
    verdict('shifted_grid_density_round_trip', m, 0.0, 2e-5, m <= 2e-5)

# 3. refusals
bad = os.path.join(WORK, 'bad_species_radii'); write_state(bad, r_src, rho_src, r_species=r_src*1.001)
rc, log = run(bad, os.path.join(src, 'Hydro_ioniz.txt'), os.path.join(WORK, 'o3')); check_eq('mismatched_species_radii_refused', rc, 2)
wide = os.path.join(WORK, 'wide'); r_w = 1.0*np.exp(np.linspace(np.log(0.9990), np.log(3.5), n)); write_state(wide, r_w, rho_src)
rc, log = run(src, os.path.join(wide, 'Hydro_ioniz.txt'), os.path.join(WORK, 'o4')); check_eq('physical_target_outside_support_refused', rc, 2)
nonm = os.path.join(WORK, 'nonmono'); r_n = r_src.copy(); r_n[10] = r_n[9]; write_state(nonm, r_n, rho_src)
rc, log = run(nonm, os.path.join(src, 'Hydro_ioniz.txt'), os.path.join(WORK, 'o5')); check_eq('non_monotone_source_radii_refused', rc, 2)
check_eq('refusals_write_nothing', os.path.exists(os.path.join(WORK, 'o4', 'Hydro_ioniz.txt')), False)

# 4. ghost rows outside the support are extrapolated, not clamped
r_g = r_src.copy(); r_g[0] = r_src[0]*0.9995; r_g[-1] = r_src[-1]*1.0005
gt = os.path.join(WORK, 'ghost'); write_state(gt, r_g, rho_src)
out = os.path.join(WORK, 'o6'); rc, log = run(src, os.path.join(gt, 'Hydro_ioniz.txt'), out); check_eq('ghost_outside_support_accepted', rc, 0)
if rc == 0:
    h = read(os.path.join(out, 'Hydro_ioniz.txt'))[1]
    check_eq('outer_ghost_not_clamped_to_the_last_cell', bool(h[-1, 1] < h[-2, 1]), True)

print(f'state_mapper: {n_fail} failure(s)')
sys.exit(1 if n_fail else 0)
