#!/usr/bin/env python3
"""Advective formation-energy inventory and flux diagnostic.

The thermal energy equation EXHALE integrates carries one material energy
flux, the hydrodynamic one that `Phys_flux` builds,

    F_E = v (E + p),   E = (1/2) rho v^2 + u_th        (Num_Fluxes.f90:434)

with `u_th` the caloric internal energy, so the rovibrational energy of bound
H2 is inside it (caloric_eos.f90, `energy_density_from_pressure`).  The
chemical (formation) energy u_form = sum_s n_s eps_s moves with the gas and
with the relative motion of the species, F_form = sum_s eps_s (n_s v +
Phi_s), and its balance follows from the species balances
d_t n_s + div(n_s v + Phi_s) = S_s:

    d_t u_form + div(F_form) = sum_s eps_s S_s .

THAT IS NOT A TERM MISSING FROM THE THERMAL EQUATION.  In the total energy
balance of the gas (kinetic, gravitational, sensible and formation energy,
with the work, the boundary fluxes and the radiation absorbed and emitted
all written in it) the chemical conversion is exchanged between u_form and
the rest; subtracting the balance above from it gives the sensible-energy
equation, in which the conversion appears ONCE, as -sum_s eps_s S_s, and
that is what the reaction heating and cooling rates of the code already are
(a photoionization heats by h nu - I, not by h nu; a recombination cools by
the kinetic energy of the captured electron, not by I).  Adding -div(F_form)
or -d_t u_form|chem to that equation would count the same conversion twice
(code audit of 2026-09-29, md/CODE_AUDIT_20260929.md F8; this header said
otherwise until then).  The sensible enthalpy carried by a species moving
relative to the mixture is a different quantity: binary_element_diffusion.f90
carries it for the elements, and the carriers are the open item of the same
audit (F1).

WHAT THIS SCRIPT MEASURES is the ADVECTIVE part only, F_form,adv = u_form v
from the cell-centred profiles; it reads no relative species flux Phi_s.  In
a stationary state

    div(F_form,adv) = sum_s eps_s S_s - div(sum_s eps_s Phi_s) ,

so its divergence equals the chemical conversion of each cell only where the
relative formation-energy transport is negligible, and only where the species
balances are stationary.  A difference between div(F_form,adv) and the
conversion is therefore a sum of distinct things -- the omitted relative
(diffusive) transport, the residual of the species balances of the saved
state, species with no eps_s entry, and the midpoint flux reconstruction of
this script -- and none of them can be singled out from the difference alone
(follow-up review md/CODE_AUDIT_20260929_review1.md, R2).  A complete ledger
test reads the operator's own fluxes and reaction sources and compares their
discrete balance; comparing div(u_form v) with the heating and cooling does
not establish it.

This script MEASURES `div(F_form)` on saved output directories and compares it
with the heating and cooling rates those same states carry.  It changes
nothing in the code.

`eps_s` is the formation and excitation energy of one particle of species `s`
above the single declared reference of B1 T1.2 (every element a neutral,
ground-state, free atom at rest; the electron carries no ionization energy,
the ion does).  Every value in EPS_EV below is READ from the tree, with the
site named in the table.

Discretization, stated because the output files carry cell centers and not
the grid's own faces `r_edg`:

  * the face between cells j and j+1 is placed at r_{j+1/2} = (r_j+r_{j+1})/2
    and carries the average of the two cell-centered flux values,
    F_{j+1/2} = (u_form,j v_j + u_form,j+1 v_{j+1})/2 ;
  * the divergence is the spherical conservative form
    div F = (A_+ F_+ - A_- F_-) / dV with A = r^2 and
    dV = (r_+^3 - r_-^3)/3 ;
  * the two ghost rows at each end supply the neighbours of the first and the
    last physical cell, so every physical cell has a two-sided stencil.

This is not the operator the code would use (that one rides on the Riemann
face mass flux), so the numbers below are the SIZE of the term, not its
discrete value in an implementation.  The script prints, beside every ratio,
the same operator applied to the mass flux, which is the same approximation
acting on a quantity whose exact divergence is zero in a converged state: it
bounds the part of the measured divergence that is the approximation rather
than the physics.

Usage
-----
    MPLBACKEND=Agg python3 src/utils/formation_energy_flux_diagnostic.py \
        [--plot-dir DIR] [label=outdir[,input=path] ...]

With no case arguments the default set is the regression goldens and the two
benchmark runs listed in DEFAULT_CASES.
"""

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, os.path.join(_ROOT, 'examples'))

import exhale_io as eio                                   # noqa: E402

EV_TO_ERG = 1.602176634e-12        # CODATA 2018, the eV_to_erg of
                                   # molecular_reaction_heat.f90
CM_TO_EV = 1.239841984e-4          # the cm_to_eV of the same module
KB = 1.380649e-23 * 1.0e7          # [erg/K], CODATA exact

# --------------------------------------------------------------------------
# The species formation/excitation energy table eps_s [eV].
#
# Reference state (B1 T1.2): every element a neutral, ground-state, free atom
# at rest; the free electron carries zero, the ionization energy is on the ion.
#
# H, He           the photoionization thresholds of parameters.f90:826, 842,
#                 843 (e_th_HI, e_th_HeI, e_th_HeII, NIST ASD).  Cumulative
#                 over the stages.
# He 2^3S         parameters.f90:844-849 states the metastable's excitation
#                 energy as 159855.9743 cm^-1 = 19.819614 eV above the ground
#                 singlet, and defines e_th_HeTR = 24.587389 - 19.819614.
#                 The HeI column of Ion_species.txt is the TOTAL He I density
#                 with the triplet inside it (exhale_io._NUC_HE note), so the
#                 triplet column adds only the excitation energy here and no
#                 nucleus is counted twice.
# H(n=2)          E_21 = e_th_HI (1 - 1/4), the Lyman-alpha energy; used only
#                 where Excited_H.txt exists.  The HI column likewise holds
#                 the total neutral hydrogen.
# metals          mion_ethr of species_table.f90:249-253, in the mion order,
#                 cumulative over the stages of each element.
# H2, H2+, H3+,   the species enthalpy table of
# HeH+            molecular_reaction_heat.f90 (species_enthalpies), which is
#                 the code's ONE definition of these energies and is written
#                 against exactly this reference ("the energy needed to build
#                 X out of ground-state H atoms, ground-state He atoms and
#                 free electrons at rest").  Its constants and their published
#                 sources: D0(H2) = 36118.11 cm^-1 (Huber & Herzberg 1979,
#                 through mol_rates::h2_dissociation_energy_eV),
#                 IP(H2) = 15.425927 eV (parameters.f90:828, NIST Chemistry
#                 WebBook), D0(H3+) = 35076 cm^-1 (Mizus et al. 2019,
#                 Mol. Phys. 117, 1663), D0(HeH+) = 16448.84 - 1566.6764
#                 cm^-1 (Coxon & Hajigeorgiou 1999, J. Mol. Spectrosc. 193,
#                 306, Tables 9 and 8).
# OH, H2O, CO     not entered: no case measured here runs the oxygen network,
#                 and B1 T1.9 is the increment that adds them.  A run carrying
#                 those columns is reported with the missing species named
#                 rather than silently given eps = 0.
# --------------------------------------------------------------------------
E_TH_HI = 13.598434599
E_TH_HEI = 24.587389
E_TH_HEII = 54.417765
E_TH_H2 = 15.425927
D0_H2_EV = 36118.11 * CM_TO_EV
D0_H3P_EV = 35076.0 * CM_TO_EV
D0_HEHP_EV = (16448.84 - 1566.6764) * CM_TO_EV
E_HE_TRIPLET = 19.819614
E_H_N2 = E_TH_HI * 0.75

# metal thresholds, mion order of species_table.f90
_METAL_THR = {
    'C': [11.26, 24.38], 'O': [13.62, 35.12], 'N': [14.53, 29.60],
    'Mg': [7.646, 15.035], 'Si': [8.152, 16.35], 'Ca': [6.113, 11.87],
    'Na': [5.139], 'K': [4.341], 'S': [10.36], 'Fe': [7.902, 16.199],
}
_ROMAN = ['I', 'II', 'III']

EPS_EV = {
    'HI': 0.0,
    'HII': E_TH_HI,
    'HeI': 0.0,
    'HeII': E_TH_HEI,
    'HeIII': E_TH_HEI + E_TH_HEII,
    'HeITR': E_HE_TRIPLET,
    'H2': -D0_H2_EV,
    'H2p': E_TH_H2 - D0_H2_EV,
    'H3p': -D0_H2_EV + E_TH_HI - D0_H3P_EV,
    'HeHp': E_TH_HI - D0_HEHP_EV,
}
for _sym, _thr in _METAL_THR.items():
    _cum = 0.0
    EPS_EV[_sym + 'I'] = 0.0
    for _k, _t in enumerate(_thr):
        _cum += _t
        EPS_EV[_sym + _ROMAN[_k + 1]] = _cum

# Species whose energy the table does not carry; a run with one of these
# columns is reported, not silently truncated.
_UNKNOWN_OK = set()


# --------------------------------------------------------------------------
# H2 rovibrational internal energy, for the enthalpy-flux comparison.
#
# The caloric equation of state sums the complete bound rovibrational ladder
# of H2 X^1 Sigma_g^+ (caloric_eos.f90 `h2_rovibrational_sum`) over the 302
# levels of `molecular_infrared_data.f90`, term energies in K and weights
# g_I (2J+1) (Roueff et al. 2019, A&A 630, A58).  This diagnostic READS that
# same table out of the Fortran source rather than writing a second ladder, so
# the u_rv it uses is the one the code's own energy flux carries.
# --------------------------------------------------------------------------
_H2_LEVEL_SOURCE = os.path.join(
    _ROOT, 'src', 'modules', 'lower_atmosphere', 'molecular_infrared_data.f90')


def _read_h2_levels(path):
    """(term energy [K], weight) of the H2 bound ladder, from the Fortran
    `data (h2_lev_T(i), ...)` and `data (h2_lev_g(i), ...)` statements."""
    text = open(path).read()
    out = {}
    for name in ('h2_lev_T', 'h2_lev_g'):
        vals = {}
        i = 0
        key = 'data (' + name + '(i)'
        while True:
            i = text.find(key, i)
            if i < 0:
                break
            head = text[i:text.index('/', i)]
            lo = int(head.split('i =')[1].split(',')[0])
            body = text[text.index('/', i) + 1:text.index('/', text.index('/', i) + 1)]
            nums = [float(t.strip().replace('d', 'e'))
                    for t in body.replace('&', ' ').split(',') if t.strip()]
            for k, val in enumerate(nums):
                vals[lo + k] = val
            i = text.index('/', text.index('/', i) + 1)
        out[name] = np.array([vals[k] for k in sorted(vals)])
    if len(out['h2_lev_T']) != len(out['h2_lev_g']):
        raise SystemExit('%s: h2_lev_T and h2_lev_g differ in length' % path)
    return out['h2_lev_T'], out['h2_lev_g']


_H2_ET, _H2_G = _read_h2_levels(_H2_LEVEL_SOURCE)


def h2_rovibrational_energy_K(T):
    """Mean rovibrational energy per bound H2 molecule [K], the u_rv of
    caloric_eos.f90, summed over the same ladder."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    ek = _H2_ET[None, :]
    x = -ek / np.maximum(T, 1.0e-30)[:, None]
    x = x - x.max(axis=1, keepdims=True)
    w = _H2_G[None, :] * np.exp(x)
    return (w * ek).sum(axis=1) / w.sum(axis=1)


# --------------------------------------------------------------------------
# geometry and the divergence operator
# --------------------------------------------------------------------------
def face_divergence(r, f):
    """Spherical divergence of a cell-centered flux `f` on cell centers `r`.

    Faces at the midpoints of the centers, face flux the average of the two
    neighbours, divergence (A_+ F_+ - A_- F_-)/dV with dV = (r_+^3-r_-^3)/3.
    Returns the divergence on cells 1..len(r)-2 (the cells with two
    neighbours) together with those cells' radii and volumes.
    """
    rf = 0.5 * (r[:-1] + r[1:])                      # faces, len n-1
    ff = 0.5 * (f[:-1] + f[1:])
    a = rf * rf
    rm, rp = rf[:-1], rf[1:]
    dv = (rp ** 3 - rm ** 3) / 3.0
    out = a[1:] * ff[1:]
    inn = a[:-1] * ff[:-1]
    div = (out - inn) / dv
    # Relative non-flatness of the flux across the cell: how much of the
    # difference survives the cancellation of the two face terms.  A quantity
    # whose exact spherical flux is constant (the mass flux of a steady state)
    # gives here only the error of the midpoint face reconstruction, so this
    # number on the mass flux is the operator's own floor and the same number
    # on F_form is what has to stand above it.
    rel = np.abs(out - inn) / np.maximum(np.abs(out) + np.abs(inn), 1e-300)
    return r[1:-1], dv, div, rel


# --------------------------------------------------------------------------
def u_form_density(ion, note):
    """u_form [erg/cm^3] = sum_s n_s eps_s, and the per-species contributions."""
    u = None
    parts = {}
    for name, n in ion.items():
        if name not in EPS_EV:
            note.append('species with no eps_s entry, contribution omitted: '
                        + name)
            continue
        c = n * (EPS_EV[name] * EV_TO_ERG)
        parts[name] = c
        u = c if u is None else u + c
    return u, parts


def analyse(label, outdir, inputfile, plot_dir=None):
    note = []
    hydro_path = os.path.join(outdir, 'Hydro_ioniz.txt')
    ion_path = os.path.join(outdir, 'Ion_species.txt')
    h = eio.load_hydro(hydro_path, ghost=True)
    r_all, ion = eio.load_ions(ion_path, ghost=True)
    inp = eio.read_input(inputfile)
    rp = inp['Rp_RJ'] * eio.RJ                       # [cm]

    if not np.allclose(r_all, h['r'], rtol=1e-12, atol=0.0):
        raise SystemExit(label + ': Hydro and Ion radius columns differ')

    # H(n=2) if the run wrote it: its excitation energy is in the reservoir
    # (B1 T1.4) and the HI column already holds the atom.
    exc_path = os.path.join(outdir, 'Excited_H.txt')
    n2 = None
    if os.path.exists(exc_path):
        e = eio.load_excited_H(exc_path, ghost=True)
        key2s = [k for k in e if k.lower() in ('n2s', 'n_2s')]
        key2p = [k for k in e if k.lower() in ('n2p', 'n_2p')]
        if key2s and key2p:
            n2 = e[key2s[0]] + e[key2p[0]]
            note.append('H(n=2) included from Excited_H.txt (%s + %s)'
                        % (key2s[0], key2p[0]))

    r = h['r'] * rp                                  # [cm]
    v = h['v']
    rho = h['n'] * eio.mu                            # [g/cm^3]
    p = h['p']
    T = h['T']
    heat, cool = h['heat'], h['cool']

    u_form, parts = u_form_density(ion, note)
    if n2 is not None:
        u_form = u_form + n2 * (E_H_N2 * EV_TO_ERG)
        parts['H(n=2)'] = n2 * (E_H_N2 * EV_TO_ERG)

    # --- the fluxes ---------------------------------------------------
    f_form = u_form * v
    f_mass = rho * v
    n_h2 = ion.get('H2')
    if n_h2 is not None:
        u_rv = n_h2 * KB * h2_rovibrational_energy_K(T)
    else:
        u_rv = np.zeros_like(p)
    u_th = 1.5 * p + u_rv
    e_tot = 0.5 * rho * v * v + u_th
    f_e_today = v * (e_tot + p)                      # Num_Fluxes.f90:434
    f_e_mono = 2.5 * p * v                           # monatomic enthalpy only
    f_kin = 0.5 * rho * v ** 3
    f_rv = v * u_rv

    rc, dv, div_form, rel_form = face_divergence(r, f_form)
    _, _, div_mass, rel_mass = face_divergence(r, f_mass)
    _, _, div_e_today, _ = face_divergence(r, f_e_today)
    _, _, div_e_mono, _ = face_divergence(r, f_e_mono)
    _, _, div_kin, _ = face_divergence(r, f_kin)
    _, _, div_rv, _ = face_divergence(r, f_rv)

    # interior slice of the cell-centered quantities
    s = slice(1, -1)
    heat_i, cool_i = heat[s], cool[s]
    net_i = heat_i - cool_i
    rho_i, u_form_i, v_i = rho[s], u_form[s], v[s]

    # the split of div(F_form) into the part the non-flat mass flux carries
    # and the part the composition gradient carries
    div_form_mass_part = (u_form_i / rho_i) * div_mass
    div_form_comp_part = div_form - div_form_mass_part

    # drop the ghost cells from the reported set: face_divergence already
    # consumed one row at each end, so one ghost row remains at each end.
    i0, i1 = eio.physical_cell_rows(hydro_path)
    keep = slice(i0 - 1, i1 - 1)                     # physical cells of rc

    out = dict(
        label=label, outdir=outdir, note=note, rp=rp,
        r=rc[keep] / rp, dv=dv[keep],
        div_form=div_form[keep], div_form_comp=div_form_comp_part[keep],
        div_form_mass=div_form_mass_part[keep],
        div_mass=div_mass[keep],
        rel_form=rel_form[keep], rel_mass=rel_mass[keep],
        heat=heat_i[keep], cool=cool_i[keep], net=net_i[keep],
        u_form=u_form_i[keep], rho=rho_i[keep], v=v_i[keep],
        div_e_today=div_e_today[keep], div_e_mono=div_e_mono[keep],
        div_kin=div_kin[keep], div_rv=div_rv[keep],
        u_rv=u_rv[s][keep], p=p[s][keep], T=T[s][keep],
        parts={k: val[s][keep] for k, val in parts.items()},
        x_hii=None,
    )
    if 'HI' in ion and 'HII' in ion:
        tot = ion['HI'] + ion['HII']
        out['x_hii'] = (ion['HII'] / np.where(tot > 0, tot, np.nan))[s][keep]
    return out


def _fmt(x):
    return '%.3e' % x


def report(a, fh):
    w = fh.write
    r, div, heat, cool, net = a['r'], a['div_form'], a['heat'], a['cool'], a['net']
    w('\n### %s  (%s)\n' % (a['label'], a['outdir']))
    for n in a['note']:
        w('  note: %s\n' % n)
    absnet = np.abs(net)
    ratio_net = np.abs(div) / np.where(absnet > 0, absnet, np.nan)
    ratio_heat = np.abs(div) / np.where(heat > 0, heat, np.nan)
    signed_heat = div / np.where(heat > 0, heat, np.nan)

    # A cell whose heating is a millionth of the peak carries no energy budget
    # to compare against, and its ratio is a division by a number the state
    # does not resolve.  Every statistic below is taken over the ENERGETIC
    # cells, those within 1e-4 of the peak heating rate; the count of the
    # excluded cells is printed so the restriction is visible.
    sig = heat >= 1.0e-4 * np.nanmax(heat)
    ncell = int(sig.sum())

    def frac(mask):
        m = mask & sig
        return '%d/%d (%.1f%%)' % (m.sum(), ncell, 100.0 * m.sum() / max(ncell, 1))

    def med(x):
        return np.nanmedian(x[sig])

    def mx(x):
        return np.nanmax(x[sig])

    w('  physical cells: %d, r = %.4f to %.4f R_p; energetic cells'
      ' (heat >= 1e-4 heat_max): %d, r = %.4f to %.4f R_p\n'
      % (len(r), r[0], r[-1], ncell, r[sig][0], r[sig][-1]))
    w('  operator floor, relative face-flux non-flatness |A+F+ - A-F-|/(|A+F+|+|A-F-|):\n')
    w('    mass flux  (exactly flat in a steady state): median %s  max %s\n'
      % (_fmt(med(a['rel_mass'])), _fmt(mx(a['rel_mass']))))
    w('    F_form                                     : median %s  max %s\n'
      % (_fmt(med(a['rel_form'])), _fmt(mx(a['rel_form']))))
    w('  max |u_form|/(3/2 p)      : %s\n'
      % _fmt(mx(np.abs(a['u_form']) / (1.5 * a['p']))))
    w('  |div F_form| / |heat-cool|: median %s  max %s\n'
      % (_fmt(med(ratio_net)), _fmt(mx(ratio_net))))
    w('     > 0.1 in %s ; > 1 in %s\n'
      % (frac(ratio_net > 0.1), frac(ratio_net > 1.0)))
    w('  |div F_form| / heat       : median %s  max %s\n'
      % (_fmt(med(ratio_heat)), _fmt(mx(ratio_heat))))
    w('     > 0.1 in %s ; > 1 in %s\n'
      % (frac(ratio_heat > 0.1), frac(ratio_heat > 1.0)))
    w('  div F_form / heat (signed): min %s  max %s\n'
      % (_fmt(np.nanmin(signed_heat[sig])), _fmt(np.nanmax(signed_heat[sig]))))
    k = int(np.nanargmax(np.where(sig, ratio_net, -np.inf)))
    w('  worst cell r = %.4f R_p: div=%s heat=%s cool=%s'
      % (r[k], _fmt(div[k]), _fmt(heat[k]), _fmt(cool[k])))
    if a['x_hii'] is not None:
        w('  x_HII=%.4f' % a['x_hii'][k])
    w('\n')

    # where the ionization front is, for the "where does it matter" statement
    if a['x_hii'] is not None:
        x = a['x_hii']
        g = np.abs(np.gradient(x, r))
        jf = int(np.nanargmax(g))
        w('  steepest x_HII gradient at r = %.4f R_p (x_HII = %.3f):'
          ' |div|/|heat-cool| = %s, |div|/heat = %s\n'
          % (r[jf], x[jf], _fmt(ratio_net[jf]), _fmt(ratio_heat[jf])))

    # split of the divergence
    den = np.abs(div)
    den = np.where(den > 0, den, np.nan)
    w('  composition-gradient share of div F_form: median %.3f\n'
      % med(np.abs(a['div_form_comp']) / den))
    w('  the part the non-flat mass flux carries, (u_form/rho) div(rho v),'
      ' as a share: median %.3f\n'
      % med(np.abs(a['div_form_mass']) / den))

    # radial breakdown
    edges = [1.0, 1.1, 1.3, 2.0, 1e9]
    w('  by radial band, median |div F_form|/heat (energetic cells only):\n')
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = sig & (r >= lo) & (r < hi)
        if m.sum() == 0:
            continue
        w('    %5.2f <= r/R_p < %-6s n=%3d  |div|/heat %s  |div|/|heat-cool| %s'
          '  rel_form %s  rel_mass %s\n'
          % (lo, ('%.2f' % hi) if hi < 1e8 else 'top', m.sum(),
             _fmt(np.nanmedian(ratio_heat[m])), _fmt(np.nanmedian(ratio_net[m])),
             _fmt(np.nanmedian(a['rel_form'][m])),
             _fmt(np.nanmedian(a['rel_mass'][m]))))

    # integrated budget
    dv = a['dv']
    iform = np.sum(div * dv)
    iheat = np.sum(heat * dv)
    icool = np.sum(cool * dv)
    w('  integrated over the physical cells [erg/s per steradian x 4pi omitted,'
      ' the same dV in all three]:\n')
    w('    int div(F_form) dV = %s\n' % _fmt(iform))
    w('    int heat dV        = %s   ratio %s\n'
      % (_fmt(iheat), _fmt(iform / iheat if iheat else np.nan)))
    w('    int cool dV        = %s   ratio %s\n'
      % (_fmt(icool), _fmt(iform / icool if icool else np.nan)))
    w('    int (heat-cool) dV = %s   ratio %s\n'
      % (_fmt(iheat - icool),
         _fmt(iform / (iheat - icool) if iheat != icool else np.nan)))

    # enthalpy flux comparison
    dth, dmo = a['div_e_today'], a['div_e_mono']
    dnet = np.where(absnet > 0, absnet, np.nan)
    w('  thermal flux, v(E+p) as built by Num_Fluxes.f90 against (5/2) p v:\n')
    w('    |div(v(E+p)) - div((5/2)pv)| / |heat-cool|: median %s  max %s\n'
      % (_fmt(med(np.abs(dth - dmo) / dnet)), _fmt(mx(np.abs(dth - dmo) / dnet))))
    w('    of which the kinetic part (1/2 rho v^3): median %s  max %s\n'
      % (_fmt(med(np.abs(a['div_kin']) / dnet)),
         _fmt(mx(np.abs(a['div_kin']) / dnet))))
    w('    of which the H2 rovibrational part (v u_rv): median %s  max %s,'
      ' max u_rv/(3/2 p) = %s\n'
      % (_fmt(med(np.abs(a['div_rv']) / dnet)),
         _fmt(mx(np.abs(a['div_rv']) / dnet)),
         _fmt(mx(a['u_rv'] / (1.5 * a['p'])))))

    # the largest species contributions to u_form
    tot = np.abs(a['u_form'])
    peak = int(np.nanargmax(np.where(sig, tot, -np.inf)))
    contrib = sorted(((np.abs(c[peak]), k) for k, c in a['parts'].items()),
                     reverse=True)[:5]
    w('  u_form at its peak (r = %.4f R_p, %s erg/cm^3), top species:'
      % (a['r'][peak], _fmt(a['u_form'][peak])))
    for val, k in contrib:
        w(' %s(%.0f%%)' % (k, 100.0 * val / max(tot[peak], 1e-300)))
    w('\n')


def plot(cases, plot_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams['text.usetex'] = False
    n = len(cases)
    fig, axes = plt.subplots(1, n, figsize=(6.0 * n, 4.6))
    if n == 1:
        axes = [axes]
    paths = []
    for ax, a in zip(axes, cases):
        r, d, h, c = a['r'], a['div_form'], a['heat'], a['cool']
        ax.plot(r, h, color='C3', lw=1.4, label='heat')
        ax.plot(r, c, color='C0', lw=1.4, label='cool')
        ax.plot(r, np.abs(d), color='k', lw=1.2,
                label=r'$|\nabla\cdot F_{\rm form}|$')
        ax.plot(r, np.where(d > 0, d, np.nan), 'k.', ms=2.5,
                label=r'$\nabla\cdot F_{\rm form} > 0$ (sink of $u_{\rm th}$)')
        ax.set_yscale('log')
        ax.set_xlabel(r'$r\ [R_p]$')
        ax.set_ylabel(r'erg cm$^{-3}$ s$^{-1}$')
        ax.set_title(a['label'])
        ax.legend(fontsize=8, loc='best')
        ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(plot_dir, 'formation_energy_flux.png')
    fig.savefig(out, dpi=130)
    paths.append(out)
    plt.close(fig)
    return paths


DEFAULT_CASES = [
    ('wasp_full', 'backup/regression/golden/wasp_full',
     'backup/regression/wasp_full/input.inp'),
    ('wasp_full_newton', 'backup/regression/golden/wasp_full_newton',
     'backup/regression/wasp_full_newton/input.inp'),
    ('mol_base_handoff', 'backup/regression/golden/mol_base_handoff',
     'backup/regression/mol_base_handoff/input.inp'),
    ('mol_carrier', 'backup/regression/golden/mol_carrier',
     'backup/regression/mol_carrier/input.inp'),
    ('lower_profile', 'backup/regression/golden/lower_profile',
     'backup/regression/lower_profile/input.inp'),
    ('hd209', 'benchmarks/hd209/output', 'benchmarks/hd209/input.inp'),
    ('wasp52', 'benchmarks/wasp52/output', 'benchmarks/wasp52/input.inp'),
]


def main(argv):
    plot_dir = '.'
    cases = []
    args = []
    i = 0
    while i < len(argv):
        if argv[i] == '--plot-dir':
            plot_dir = argv[i + 1]
            i += 2
            continue
        args.append(argv[i])
        i += 1
    if args:
        for a in args:
            label, rest = a.split('=', 1)
            parts = rest.split(',')
            outdir = parts[0]
            inputfile = None
            for q in parts[1:]:
                if q.startswith('input='):
                    inputfile = q[6:]
            if inputfile is None:
                inputfile = os.path.join(os.path.dirname(outdir), 'input.inp')
            cases.append((label, outdir, inputfile))
    else:
        cases = [(l, os.path.join(_ROOT, o), os.path.join(_ROOT, f))
                 for l, o, f in DEFAULT_CASES]

    print('# eps_s table [eV], reference: neutral ground-state free atoms')
    for k in sorted(EPS_EV):
        print('#   %-6s %+12.6f' % (k, EPS_EV[k]))
    print('#   %-6s %+12.6f  (only where Excited_H.txt exists)'
          % ('H(n=2)', E_H_N2))

    done = []
    for label, outdir, inputfile in cases:
        if not os.path.exists(os.path.join(outdir, 'Hydro_ioniz.txt')):
            print('\n### %s: no Hydro_ioniz.txt in %s, skipped'
                  % (label, outdir))
            continue
        a = analyse(label, outdir, inputfile)
        report(a, sys.stdout)
        done.append(a)

    want = [a for a in done if a['label'] in ('wasp_full', 'mol_base_handoff')]
    if want:
        for p in plot(want, plot_dir):
            print('\nplot: %s' % p)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
