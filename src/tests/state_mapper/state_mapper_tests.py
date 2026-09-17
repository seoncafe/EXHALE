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

# ---- 5. --reservoir He/H: the helium carried onto another reservoir -------
# A state whose pressure is the particle count's, so that the temperature the
# tool writes can be checked against T = p/((n_tot + n_e) k_B), and whose
# He/H varies with radius the way an element-diffusion solution's does:
# pinned to the reservoir at the base rows, falling outward.
KB = 1.380649e-16
M_HE_OVER_M_H = 6.6464790722e-24/1.67353284e-24

def write_reservoir_state(d, r, reservoir=1.6, hehp=0.0, line=True, diffused=True):
    os.makedirs(d, exist_ok=True)
    n = r.size
    nH = 1e12*np.exp(-6*(r - 1.0))
    # He/H: the reservoir over the base rows (the two ghosts and cell 1),
    # then falling to a quarter of it at the top.
    q = np.full(n, float(reservoir))
    if diffused:
        f = np.clip((np.log(r/r[2]))/np.log(r[-1]/r[2]), 0.0, 1.0)
        q[3:] = reservoir*(1.0 - 0.75*f[3:])
    nHe = q*nH
    x = 0.1*(r - 1.0)/(r.max() - 1.0)
    HI, HII = nH*(1 - x), nH*x
    HeI, HeII, HeIII = 0.90*nHe, 0.08*nHe, 0.02*nHe
    HeTR = 1e-4*HeI
    HeHp = np.full(n, hehp*nH[0]) if hehp else np.full(n, 1e-300)
    zero = np.full(n, 1e-300)
    ntot_ne = HI + HII + HeI + HeII + HeIII + (HII + HeII + 2*HeIII)
    T = 1000.0*(1 + 0.5*(r - 1))
    p = ntot_ne*KB*T
    rho = HI + HII + M_HE_OVER_M_H*(HeI + HeII + HeIII)
    heat, cool = 1e-3*rho, 2e-3*rho
    def hdr(cols):
        h = ['# EXHALE schema 2', '# columns ' + ' '.join(cols),
             f'# rows {n}: 2 ghost cells at each end; physical cells are rows 3 to {n-2}']
        if line:
            h.append(f'# reservoir He/H {float(reservoir):.16E} C/H 2.7000000000000000E-04')
        h += [f'# constants set IAU2015+CODATA2018 RJ[cm] 7.1492000000000000E+09 kB[erg/K] {KB:.16E}',
              '# coupling: mode=phys t_phys=1.0 certified=T',
              '# provenance: git=test tree=clean run=test']
        return h
    with open(os.path.join(d, 'Hydro_ioniz.txt'), 'w') as f:
        f.write('\n'.join(hdr(HCOLS)) + '\n')
        for row in zip(r, rho, 1e5*(r - 1.0), p, T, heat, cool):
            f.write(' '.join(f'{x:24.16E}' for x in row) + '\n')
    with open(os.path.join(d, 'Ion_species.txt'), 'w') as f:
        f.write('\n'.join(hdr(SCOLS)) + '\n')
        for row in zip(r, HI, HII, HeI, HeII, HeIII, HeTR, zero, zero, zero, HeHp):
            f.write(' '.join(f'{x:24.16E}' for x in row) + '\n')
    return q

def run_opts(src, tgt, out, *opts):
    p = subprocess.run([sys.executable, TOOL, src, tgt, out, *opts],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr

res_src = os.path.join(WORK, 'res_src')
q_src = write_reservoir_state(res_src, r_src, reservoir=1.6)
res_tgt = os.path.join(res_src, 'Hydro_ioniz.txt')

# The control: the same state mapped without the option.  The columns the
# option does not touch have to come out of it unchanged to the last bit,
# which is a stronger statement than agreement with the source (the mapping
# itself round trips the density through log10 and moves it by 1e-15).
control = os.path.join(WORK, 'res_control')
rc_c, _ = run_opts(res_src, res_tgt, control)
check_eq('reservoir_He_control_exit_status', rc_c, 0)

# 5a. every helium column is multiplied by the ratio of the two reservoirs and
#     nothing else is, so He/H(r) is the source profile times that ratio
out = os.path.join(WORK, 'res_scaled')
rc, log = run_opts(res_src, res_tgt, out, '--reservoir', 'He/H', '2.13')
check_eq('reservoir_He_exit_status', rc, 0)
if rc == 0:
    a0 = read(os.path.join(control, 'Ion_species.txt'))[1]
    a1 = read(os.path.join(out, 'Ion_species.txt'))[1]
    fac = 2.13/1.6
    for j, name in ((3, 'HeI'), (4, 'HeII'), (5, 'HeIII'), (6, 'HeITR')):
        m = float(np.max(np.abs(a1[:, j]/(fac*a0[:, j]) - 1.0)))
        verdict(f'reservoir_He_{name}_scaled_by_the_ratio', m, 0.0, 1e-15, m <= 1e-15)
    for j, name in ((1, 'HI'), (2, 'HII')):
        m = float(np.max(np.abs(a1[:, j] - a0[:, j])))
        verdict(f'reservoir_He_{name}_untouched', m, 0.0, 0.0, m == 0.0)
    q1 = (a1[:, 3] + a1[:, 4] + a1[:, 5])/(a1[:, 1] + a1[:, 2])
    base = float(np.max(np.abs(q1[:3] - 2.13)))/2.13
    verdict('reservoir_He_base_rows_at_the_new_reservoir', base, 0.0, 1e-14, base <= 1e-14)
    shape = float(np.max(np.abs(q1/(fac*q_src) - 1.0)))
    verdict('reservoir_He_diffused_profile_shape_kept', shape, 0.0, 1e-12, shape <= 1e-12)
    verdict('reservoir_He_profile_is_not_uniform', float(q1.max()/q1.min()), 4.0, 0.1,
            abs(q1.max()/q1.min() - 4.0) <= 0.1)

# 5b. the header: both files state the new reservoir, the other elements of
#     the line are left alone, and the mapping line records the rescaling
if rc == 0:
    for fn in ('Hydro_ioniz', 'Ion_species'):
        hdr = read(os.path.join(out, fn + '.txt'))[0]
        res = [h for h in hdr if h.startswith('# reservoir')]
        check_eq(f'reservoir_He_line_rewritten_{fn}', len(res) == 1
                 and res[0].split()[3] == f'{2.13:.16E}', True)
        check_eq(f'reservoir_He_other_elements_kept_{fn}',
                 len(res) == 1 and res[0].split()[4:] == ['C/H', '2.7000000000000000E-04'], True)
    mapped = [h for h in read(os.path.join(out, 'Hydro_ioniz.txt'))[0] if h.startswith('# mapped:')]
    check_eq('reservoir_He_mapping_line_records_the_rescaling',
             len(mapped) == 1 and 'carried from the reservoir' in mapped[0]
             and 'initialization choice' in mapped[0], True)

# 5c. the pressure is kept, the density and the temperature follow
if rc == 0:
    h0 = read(os.path.join(control, 'Hydro_ioniz.txt'))[1]
    h1 = read(os.path.join(out, 'Hydro_ioniz.txt'))[1]
    mp = float(np.max(np.abs(h1[:, 3] - h0[:, 3])))
    verdict('reservoir_He_pressure_unchanged', mp, 0.0, 0.0, mp == 0.0)
    mv = float(np.max(np.abs(h1[:, 2] - h0[:, 2])))
    verdict('reservoir_He_velocity_unchanged', mv, 0.0, 0.0, mv == 0.0)
    a1 = read(os.path.join(out, 'Ion_species.txt'))[1]
    ntot_ne = (a1[:, 1] + a1[:, 2] + a1[:, 3] + a1[:, 4] + a1[:, 5]
               + a1[:, 2] + a1[:, 4] + 2*a1[:, 5])
    mt = float(np.max(np.abs(h1[:, 4]/(h1[:, 3]/(ntot_ne*KB)) - 1.0)))
    verdict('reservoir_He_temperature_is_p_over_the_particle_count', mt, 0.0, 1e-13, mt <= 1e-13)
    rho = a1[:, 1] + a1[:, 2] + M_HE_OVER_M_H*(a1[:, 3] + a1[:, 4] + a1[:, 5])
    mr = float(np.max(np.abs(h1[:, 1]/rho - 1.0)))
    verdict('reservoir_He_density_is_the_species_mass_sum', mr, 0.0, 1e-13, mr <= 1e-13)

# 5d. without the option nothing of this runs: the data rows are those of the
#     same state written without a reservoir line at all, and the reservoir
#     line is carried through as it stands
plain = os.path.join(WORK, 'res_plain')
write_reservoir_state(plain, r_src, reservoir=1.6, line=False)
o_res, o_plain = os.path.join(WORK, 'no_opt_res'), os.path.join(WORK, 'no_opt_plain')
rc1, _ = run_opts(res_src, res_tgt, o_res)
rc2, _ = run_opts(plain, res_tgt, o_plain)
check_eq('no_option_exit_status', (rc1, rc2), (0, 0))
if (rc1, rc2) == (0, 0):
    for fn in ('Hydro_ioniz', 'Ion_species'):
        rows_a = [l for l in open(os.path.join(o_res, fn + '.txt')) if not l.startswith('#')]
        rows_b = [l for l in open(os.path.join(o_plain, fn + '.txt')) if not l.startswith('#')]
        check_eq(f'no_option_data_rows_byte_identical_{fn}', rows_a == rows_b, True)
    hdr = read(os.path.join(o_res, 'Hydro_ioniz.txt'))[0]
    res = [h for h in hdr if h.startswith('# reservoir')]
    check_eq('no_option_reservoir_line_carried_through',
             len(res) == 1 and res[0].split()[3] == f'{1.6:.16E}', True)
    mapped = [h for h in hdr if h.startswith('# mapped:')]
    check_eq('no_option_mapping_line_states_no_rescaling',
             len(mapped) == 1 and 'carried from the reservoir' not in mapped[0], True)

# 5e. the refusals
rc, log = run_opts(plain, res_tgt, os.path.join(WORK, 'o7'), '--reservoir', 'He/H', '2.13')
check_eq('reservoir_He_without_a_reservoir_line_refused', rc, 2)
rc, log = run_opts(res_src, res_tgt, os.path.join(WORK, 'o8'), '--reservoir', 'He/H', '-1')
check_eq('reservoir_He_nonpositive_refused', rc, 2)
hehp = os.path.join(WORK, 'res_hehp')
write_reservoir_state(hehp, r_src, reservoir=1.6, hehp=0.05)
rc, log = run_opts(hehp, res_tgt, os.path.join(WORK, 'o9'), '--reservoir', 'He/H', '2.13')
check_eq('reservoir_He_with_HeHp_refused', rc, 2)
check_eq('reservoir_He_refusals_write_nothing',
         os.path.exists(os.path.join(WORK, 'o9', 'Hydro_ioniz.txt')), False)

# ---- 6. --reservoir <El>/H for a metal ------------------------------------
# The same rule carried to an element of the metal block: one factor over the
# ionization stages, the ratio taken from the '# reservoir' line, the mass at
# the weight calc_rho gives the element (amu_over_m_H x A_El) and the
# particles at 1 + the charge of the stage.
AMU_OVER_M_H = 1.66053906660e-24/1.67353284e-24
A_C, A_O = AMU_OVER_M_H*12.011, AMU_OVER_M_H*15.999
MCOLS = ('r[Rp] HI HII HeI HeII HeIII HeITR CI CII CIII OI OII OIII').split()
# Mg/H is stated but carries no column, N/H carries neither: the two refusals.
MRES = ('# reservoir He/H {heh:.16E} C/H {ch:.16E} O/H {oh:.16E}'
        ' Mg/H 3.8000000000000000E-05')

def write_metal_state(d, r, heh=1.6, ch=2.7e-4, oh=4.9e-4):
    os.makedirs(d, exist_ok=True)
    n = r.size
    nH = 1e12*np.exp(-6*(r - 1.0))
    x = 0.1*(r - 1.0)/(r.max() - 1.0)
    HI, HII = nH*(1 - x), nH*x
    nHe = heh*nH
    HeI, HeII, HeIII = 0.90*nHe, 0.08*nHe, 0.02*nHe
    HeTR = 1e-4*HeI
    nC, nO = ch*nH, oh*nH
    CI, CII, CIII = 0.5*nC, 0.3*nC, 0.2*nC
    OI, OII, OIII = 0.6*nO, 0.25*nO, 0.15*nO
    ntot_ne = (HI + HII + HeI + HeII + HeIII + CI + CII + CIII + OI + OII + OIII
               + HII + HeII + 2*HeIII + CII + 2*CIII + OII + 2*OIII)
    T = 1000.0*(1 + 0.5*(r - 1))
    p = ntot_ne*KB*T
    rho = (HI + HII + M_HE_OVER_M_H*(HeI + HeII + HeIII)
           + A_C*(CI + CII + CIII) + A_O*(OI + OII + OIII))
    def hdr(cols):
        return ['# EXHALE schema 2', '# columns ' + ' '.join(cols),
                f'# rows {n}: 2 ghost cells at each end; physical cells are rows 3 to {n-2}',
                MRES.format(heh=heh, ch=ch, oh=oh),
                f'# constants set IAU2015+CODATA2018 RJ[cm] 7.1492000000000000E+09 kB[erg/K] {KB:.16E}',
                '# coupling: mode=phys t_phys=1.0 certified=T',
                '# provenance: git=test tree=clean run=test']
    with open(os.path.join(d, 'Hydro_ioniz.txt'), 'w') as f:
        f.write('\n'.join(hdr(HCOLS)) + '\n')
        for row in zip(r, rho, 1e5*(r - 1.0), p, T, 1e-3*rho, 2e-3*rho):
            f.write(' '.join(f'{y:24.16E}' for y in row) + '\n')
    with open(os.path.join(d, 'Ion_species.txt'), 'w') as f:
        f.write('\n'.join(hdr(MCOLS)) + '\n')
        for row in zip(r, HI, HII, HeI, HeII, HeIII, HeTR,
                       CI, CII, CIII, OI, OII, OIII):
            f.write(' '.join(f'{y:24.16E}' for y in row) + '\n')

def metal_ratios(a):
    """(He/H, C/H, O/H) of every row of a species array in MCOLS order."""
    nH = a[:, 1] + a[:, 2]
    return ((a[:, 3] + a[:, 4] + a[:, 5])/nH,
            (a[:, 7] + a[:, 8] + a[:, 9])/nH,
            (a[:, 10] + a[:, 11] + a[:, 12])/nH)

def metal_state_is_consistent(tag, out):
    """The written state is ONE state: T = p/((n_tot + n_e) k_B) and rho the
    species mass sum at calc_rho's weights."""
    h = read(os.path.join(out, 'Hydro_ioniz.txt'))[1]
    a = read(os.path.join(out, 'Ion_species.txt'))[1]
    ntot_ne = (a[:, 1] + a[:, 2] + a[:, 3] + a[:, 4] + a[:, 5]
               + a[:, 7] + a[:, 8] + a[:, 9] + a[:, 10] + a[:, 11] + a[:, 12]
               + a[:, 2] + a[:, 4] + 2*a[:, 5]
               + a[:, 8] + 2*a[:, 9] + a[:, 11] + 2*a[:, 12])
    mt = float(np.max(np.abs(h[:, 4]/(h[:, 3]/(ntot_ne*KB)) - 1.0)))
    verdict(f'{tag}_temperature_is_p_over_the_particle_count', mt, 0.0, 1e-13, mt <= 1e-13)
    rho = (a[:, 1] + a[:, 2] + M_HE_OVER_M_H*(a[:, 3] + a[:, 4] + a[:, 5])
           + A_C*(a[:, 7] + a[:, 8] + a[:, 9])
           + A_O*(a[:, 10] + a[:, 11] + a[:, 12]))
    mr = float(np.max(np.abs(h[:, 1]/rho - 1.0)))
    verdict(f'{tag}_density_is_the_species_mass_sum', mr, 0.0, 1e-13, mr <= 1e-13)

met_src = os.path.join(WORK, 'met_src')
write_metal_state(met_src, r_src)
met_tgt = os.path.join(met_src, 'Hydro_ioniz.txt')
met_ctl = os.path.join(WORK, 'met_control')
rc_c, _ = run_opts(met_src, met_tgt, met_ctl)
check_eq('reservoir_metal_control_exit_status', rc_c, 0)

# 6a. one element: only its ionization stages move, by the ratio of the two
out = os.path.join(WORK, 'met_one')
rc, log = run_opts(met_src, met_tgt, out, '--reservoir', 'C/H', '3.0e-4')
check_eq('reservoir_metal_one_element_exit_status', rc, 0)
if rc == 0 and rc_c == 0:
    a0 = read(os.path.join(met_ctl, 'Ion_species.txt'))[1]
    a1 = read(os.path.join(out, 'Ion_species.txt'))[1]
    fac = 3.0e-4/2.7e-4
    for j, name in ((7, 'CI'), (8, 'CII'), (9, 'CIII')):
        m = float(np.max(np.abs(a1[:, j]/(fac*a0[:, j]) - 1.0)))
        verdict(f'reservoir_metal_{name}_scaled_by_the_ratio', m, 0.0, 1e-15, m <= 1e-15)
    for j, name in ((1, 'HI'), (3, 'HeI'), (10, 'OI'), (12, 'OIII')):
        m = float(np.max(np.abs(a1[:, j] - a0[:, j])))
        verdict(f'reservoir_metal_{name}_untouched', m, 0.0, 0.0, m == 0.0)
    q_he, q_c, q_o = metal_ratios(a1)
    dev = abs(q_c[2] - 3.0e-4)/3.0e-4
    verdict('reservoir_metal_first_physical_cell_at_the_new_ratio', dev, 0.0, 1e-14, dev <= 1e-14)
    flat = float(np.max(np.abs(q_c/q_c[2] - 1.0)))
    verdict('reservoir_metal_ratio_flat_over_the_column', flat, 0.0, 1e-14, flat <= 1e-14)
    h0 = read(os.path.join(met_ctl, 'Hydro_ioniz.txt'))[1]
    h1 = read(os.path.join(out, 'Hydro_ioniz.txt'))[1]
    mp = float(np.max(np.abs(h1[:, 3] - h0[:, 3])))
    verdict('reservoir_metal_pressure_unchanged', mp, 0.0, 0.0, mp == 0.0)
    mv = float(np.max(np.abs(h1[:, 2] - h0[:, 2])))
    verdict('reservoir_metal_velocity_unchanged', mv, 0.0, 0.0, mv == 0.0)
    metal_state_is_consistent('reservoir_metal_one_element', out)
    hdr = read(os.path.join(out, 'Hydro_ioniz.txt'))[0]
    res = [h for h in hdr if h.startswith('# reservoir')][0].split()[2:]
    check_eq('reservoir_metal_line_rewritten',
             res[2:4] == ['C/H', f'{3.0e-4:.16E}'], True)
    check_eq('reservoir_metal_other_elements_kept',
             res[0:2] == ['He/H', f'{1.6:.16E}']
             and res[4:6] == ['O/H', f'{4.9e-4:.16E}'], True)

# 6b. several elements in one call, helium among them
out = os.path.join(WORK, 'met_many')
rc, log = run_opts(met_src, met_tgt, out, '--reservoir', 'C/H', '3.0e-4',
                   '--reservoir', 'O/H', '5.5e-4', '--reservoir', 'He/H', '2.13')
check_eq('reservoir_metal_several_elements_exit_status', rc, 0)
if rc == 0:
    a1 = read(os.path.join(out, 'Ion_species.txt'))[1]
    q_he, q_c, q_o = metal_ratios(a1)
    for name, q, ref in (('He', q_he, 2.13), ('C', q_c, 3.0e-4), ('O', q_o, 5.5e-4)):
        dev = abs(q[2] - ref)/ref
        verdict(f'reservoir_metal_several_{name}_at_the_new_ratio', dev, 0.0, 1e-14, dev <= 1e-14)
    metal_state_is_consistent('reservoir_metal_several', out)
    hdr = read(os.path.join(out, 'Ion_species.txt'))[0]
    res = [h for h in hdr if h.startswith('# reservoir')][0].split()[2:]
    check_eq('reservoir_metal_several_line_rewritten',
             res == ['He/H', f'{2.13:.16E}', 'C/H', f'{3.0e-4:.16E}',
                     'O/H', f'{5.5e-4:.16E}', 'Mg/H', '3.8000000000000000E-05'], True)
    mapped = [h for h in hdr if h.startswith('# mapped:')]
    check_eq('reservoir_metal_mapping_line_records_every_element',
             len(mapped) == 1 and mapped[0].count('carried from the reservoir') == 3, True)

# 6c. the refusals
rc, log = run_opts(met_src, met_tgt, os.path.join(WORK, 'o10'),
                   '--reservoir', 'Mg/H', '4.0e-5')
check_eq('reservoir_metal_without_a_column_refused', rc, 2)
rc, log = run_opts(met_src, met_tgt, os.path.join(WORK, 'o11'),
                   '--reservoir', 'N/H', '8.0e-5')
check_eq('reservoir_metal_absent_from_the_reservoir_line_refused', rc, 2)
rc, log = run_opts(met_src, met_tgt, os.path.join(WORK, 'o12'),
                   '--reservoir', 'C/H', '-1')
check_eq('reservoir_metal_nonpositive_refused', rc, 2)
rc, log = run_opts(met_src, met_tgt, os.path.join(WORK, 'o13'),
                   '--reservoir', 'Xx/H', '1.0')
check_eq('reservoir_unknown_element_refused', rc, 2)
check_eq('reservoir_metal_refusals_write_nothing',
         os.path.exists(os.path.join(WORK, 'o10', 'Hydro_ioniz.txt')), False)

print(f'state_mapper: {n_fail} failure(s)')
sys.exit(1 if n_fail else 0)
