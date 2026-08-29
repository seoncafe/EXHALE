#!/usr/bin/env python3
"""VULCAN -> EXHALE lower-atmosphere profile (the cross-check arm of E2).

Same schema, same command line and same fingerprint rule as
`photochem_to_lower_profile.py`; the only difference is where the chemistry
comes from.  VULCAN has no climate model, so this arm can never be the path
away from a prescribed T(p) (`docs/vulcan_photochem_comparison.md` P1.7);
what it is for is stating, on the same planet and the same matching level,
how much of the handoff is the chemistry code rather than the chemistry
(P1 measured 1.70x in q_H on an 864 K base and 1.08x on a 2331 K base).

A VULCAN run is not launched from here: VULCAN is configured by editing
`vulcan_cfg.py` and regenerating `chem_funs.py`, and repeating that from a
command line would be a second, divergent way to run it.  This adapter reads
a finished `.vul` output, which is what `VULCAN_run_*/` and
`EXHALE_v1.00/VULCAN/output/` hold.

The scalar generator `vulcan_to_base.py` stays: it writes the single-level
`base.inp` special case, which four regression cases pin.  This one writes
the profile, and the two must not both be read into the same run -- with a
profile in use EXHALE refuses every scalar physics key beside it.

Example (HD 209458 b, the P1 configuration):

  python3 src/utils/vulcan_to_lower_profile.py <run_dir> \\
      --vulfile vulcan_work/hd209_vulcan/output/HD209.vul \\
      --mp 0.720 --r-ref 1.36 --p-match 1e-6
"""
import argparse
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lower_profile_schema as sch                       # noqa: E402


def load_vul(path):
    with open(path, 'rb') as fh:
        d = pickle.load(fh)
    species = list(d['variable']['species'])
    ymix = np.asarray(d['variable']['ymix'], dtype=float)   # (nz, nsp)
    p_dyn = np.asarray(d['atm']['pco'], dtype=float)
    T = np.asarray(d['atm']['Tco'], dtype=float)
    kzz = np.asarray(d['atm']['Kzz'], dtype=float)
    pico = np.asarray(d['atm']['pico'], dtype=float)
    # The condensed carriers VULCAN carries as their own species ('H2O_l_s',
    # 'NH3_l_s').  They are nuclei that have left the gas, so they belong to
    # the cold trap and not to the gas EXHALE's base inherits; the label
    # alone cannot say so, because `formula_elements` strips the phase
    # suffix and would count 'H2O_l_s' as a second H2O.
    # `gas_indx` is the index list VULCAN itself sums the gas over
    # (`store.py:149`); its complement is the condensed set.
    gas_indx = d['atm'].get('gas_indx')
    if gas_indx is None:
        non_gas = []
    else:
        gi = set(int(i) for i in gas_indx)
        non_gas = [sp for i, sp in enumerate(species) if i not in gi]
    # VULCAN carries K_zz on the cell INTERFACES, excluding the two domain
    # boundaries; put it on the cell centres the composition is stated on.
    if len(kzz) == len(p_dyn):
        p_kzz = p_dyn
    elif len(kzz) == len(pico) - 2:
        p_kzz = pico[1:-1]
    elif len(kzz) == len(pico):
        p_kzz = pico
    else:
        sch.refuse('the K_zz array has %d entries, which matches neither the'
                   ' %d cell centres nor the %d interfaces'
                   % (len(kzz), len(p_dyn), len(pico)))
    o = np.argsort(np.log(p_kzz))
    Kzz = np.interp(np.log(p_dyn), np.log(p_kzz)[o], kzz[o])
    return species, ymix, p_dyn, T, Kzz, non_gas


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sch.add_common_arguments(ap)
    ap.add_argument('--vulfile', required=True,
                    help='a finished VULCAN output (.vul pickle)')
    ap.add_argument('--network', default=None,
                    help='the reaction network the run used, for the header'
                         ' (a path is hashed)')
    args = ap.parse_args()

    species, ymix, p_dyn, T, Kzz_src, non_gas = load_vul(args.vulfile)
    p_bar = p_dyn/sch.BAR

    counts = {sp: sch.formula_elements(sp) for sp in species}
    counts = {sp: c for sp, c in counts.items() if c}
    mixing = {sp: ymix[:, species.index(sp)] for sp in counts}
    # The gas phase is what sets the mean molecular weight and hence the
    # scale height; VULCAN excludes the condensed species from mu the same
    # way (`build_atm.py:267`).
    gas = {sp: y for sp, y in mixing.items() if sp not in non_gas}
    mu = sch.mean_molecular_weight(gas, counts)
    ratios, _ = sch.element_ratios(
        mixing if args.count_condensates else gas, counts)

    # A stated eddy coefficient (--kzz-const, or --kzz-power with --kzz-ref)
    # replaces the one the VULCAN run was mixed at; with neither stated the
    # file keeps VULCAN's own K_zz(p).  This arm cannot re-solve the
    # chemistry, so a stated profile here changes only what is handed over.
    Kzz = sch.eddy_diffusion_coefficient(p_bar, args)
    if Kzz is None:
        Kzz = Kzz_src
    n_tot = p_dyn/(sch.KB*T)
    rho = n_tot*mu*sch.MAMU
    r_RJ = sch.hydrostatic_radius(p_bar, T, mu, args.mp, args.r_ref,
                                  args.p_ref)

    def y(sp):
        return (ymix[:, species.index(sp)] if sp in species
                else np.zeros_like(p_bar))

    cols = {'p': p_bar, 'r': r_RJ, 'T': T, 'n_tot': n_tot, 'rho': rho,
            'Kzz': Kzz, 'q_H2': y('H2'), 'q_H': y('H')}
    order = ['p', 'r', 'T', 'n_tot', 'rho', 'Kzz', 'q_H2', 'q_H']
    for el in sch.PROFILE_ELEMENTS:
        if el in ratios:
            cols['X_' + el] = ratios[el]
            order.append('X_' + el)
    cols['F_H'] = np.full_like(p_bar, args.trial_flux_H)
    cols['F_He'] = np.full_like(p_bar, args.trial_flux_He)
    order += ['F_H', 'F_He']
    for sp in sch.DIAGNOSTIC_MOLECULES:
        if sp in species:
            cols['q_' + sp] = y(sp)
            order.append('q_' + sp)

    p_deep = args.p_deep if args.p_deep else p_bar.max()
    p_top = args.p_top if args.p_top else p_bar.min()
    cols = sch.clip_table(sch.sort_deep_to_shallow(cols), p_deep, p_top)

    net = args.network or 'stated by the VULCAN configuration'
    net_sha = sch.sha256_file(args.network) if args.network else 'absent'
    fp = dict(
        source_code='vulcan',
        vulfile=os.path.basename(args.vulfile),
        vulfile_sha=sch.sha256_file(args.vulfile),
        network=os.path.basename(net), network_sha=net_sha,
        mp_MJ=args.mp, r_ref_RJ=args.r_ref, p_ref_bar=args.p_ref,
        kzz_const=args.kzz_const, count_condensates=args.count_condensates,
        p_match_bar=args.p_match,
        trial_flux_H=args.trial_flux_H, trial_flux_He=args.trial_flux_He,
    )
    if args.kzz_power is not None:
        # Only when it is in use: a fingerprint entry that is always present
        # would move the solution_id of every constant-K_zz run.
        fp.update(kzz_power=args.kzz_power, kzz_ref=args.kzz_ref,
                  kzz_ref_bar=args.kzz_ref_bar)
    version = 'vulcan (output %s)' % os.path.basename(args.vulfile)
    conden = ('; condensed carriers %s counted in the El/H ratios'
              % (', '.join(non_gas) if non_gas else 'none')
              if args.count_condensates else
              ('; condensed carriers %s excluded from the El/H ratios (the'
               ' cold trap)' % ', '.join(non_gas)) if non_gas else '')
    header = dict(
        solution_id=sch.solution_fingerprint(fp),
        source_code='vulcan',
        source_version=version + ', vulcan_to_lower_profile.py',
        mechanism='%s sha256 %s' % (os.path.basename(net), net_sha[:16]),
        stellar_flux='stated by the VULCAN configuration of ' +
                     os.path.basename(args.vulfile),
        p_match_bar=args.p_match, p_top_bar=0.0, p_deep_bar=0.0,
        trial_flux_H=args.trial_flux_H, trial_flux_He=args.trial_flux_He,
        iteration=args.iteration, reached_steady_state=True,
        notes=('cross-check arm: a finished VULCAN solution, T(p) prescribed'
               ' and no climate step; F_H and F_He are the STATED trial'
               ' fluxes, not measured') + conden)

    # The El/H ratios VULCAN was given are a property of its own
    # configuration, not of this file, so there is no input vector to check
    # the deepest level against here; what the two arms are compared on is
    # the matching level (docs/phase_e_flux_closure_design.md test T-E9).
    sch.write_handoff(args, cols, header, order)


if __name__ == '__main__':
    main()
