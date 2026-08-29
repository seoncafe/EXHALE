#!/usr/bin/env python3
"""Where the deepest-level El/H departure is born.

Runs the SAME adapter call the closure ladder makes, but records the
elemental ratios at three points of the chain:

  (a) `mix1`, the equilibrium composition on the CLIMATE grid, as
      `initialize_to_climate_equilibrium_PT` hands it to
      `_initialize_atmosphere` (after the quench overwrite and the
      renormalization Photochem applies to CH4/CO/NH3/HCN/H2/CO2);
  (b) the photochemical model's own initial state, i.e. after that column
      has been interpolated onto the photochemical altitude grid;
  (c) the converged steady state -- the state the adapter's conservation
      check is formed on.

Nothing in EXHALE is modified: the adapter is imported and its main() run,
with Photochem's gas-giant class subclassed to take the two snapshots.

usage: probe_initial_vs_steady.py <run_dir> <He/H>
"""
import os
import sys

import numpy as np

EX = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00'
sys.path.insert(0, os.path.join(EX, 'src/utils'))

REC = {}


def install_probe():
    import photochem.extensions.gasgiants as gg
    base = gg.EvoAtmosphereGasGiant

    class ProbedGasGiant(base):
        def _initialize_atmosphere(self, P1, T1, Kzz1, z1, mix1):
            REC['clima_P'] = np.array(P1)
            REC['clima_mix'] = {k: np.array(v) for k, v in mix1.items()}
            out = base._initialize_atmosphere(self, P1, T1, Kzz1, z1, mix1)
            REC['photo_init'] = self.return_atmosphere()
            import photochem_to_lower_profile as _ad
            # Element counts from Photochem's own composition matrix, never
            # parsed from the label: O1D and N2D are excited states.
            REC['counts'] = _ad.species_element_counts(self)
            return out

    gg.EvoAtmosphereGasGiant = ProbedGasGiant


def ratios(mix, counts, idx):
    """El/H nuclei ratios at level `idx` of a {species: profile} dict."""
    nuc = {}
    for sp, y in mix.items():
        for el, k in counts.get(sp, {}).items():
            nuc[el] = nuc.get(el, 0.0) + k*float(np.asarray(y)[idx])
    nH = nuc['H']
    return {el: v/nH for el, v in nuc.items() if el != 'H'}, nH, nuc


def report(tag, mix, counts, idx, want):
    r, nH, nuc = ratios(mix, counts, idx)
    print('  %s (level %d):  n_H = %.9e' % (tag, idx, nH))
    for el in ('He', 'C', 'N', 'O'):
        if el in r and el in want:
            print('     %-3s %.12e   dev %+.4e'
                  % (el, r[el], r[el]/want[el] - 1.0))
    return r, nuc


def main():
    run_dir, heh = sys.argv[1], sys.argv[2]
    os.makedirs(run_dir, exist_ok=True)
    install_probe()
    import photochem_to_lower_profile as ad

    sed = os.path.join(EX, 'LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt')
    sys.argv = ['photochem_to_lower_profile.py', run_dir,
                '--mp', '0.0176220', '--r-ref', '0.157692', '--p-ref', '1.0',
                '--p-match', '1.0e-6', '--climate', '--climate-p-deep', '20.0',
                '--boa-pressure-factor', '1.0', '--stellar-flux', sed,
                '--flux-at-planet', '--wavelength-unit', 'A', '--toa', '1.0e-2',
                '--atoms', 'H,He,N,O,C', '--abundances', 'He=' + heh,
                '--kzz-const', '1.0e9',
                '--trial-flux-H', '4768206.9882',
                '--trial-flux-He', '20433648.55',
                '--iteration', '0', '--abundance-tol', '1.0']

    # main() itself needs the converged model; re-run its body far enough to
    # keep the object.  Simplest: let main() write the files, then read the
    # deepest row of the profile it wrote for point (c).
    ad.main()

    want = dict(ad.SOLAR_ABUNDANCES)
    want['He'] = float(heh)
    want = {el: want[el]/want['H'] for el in want if el != 'H'}

    print()
    print('=' * 70)
    print('DEEPEST-LEVEL ELEMENTAL RATIOS, He/H = %s' % heh)
    print('=' * 70)

    # Element counts as the adapter forms them: Photochem's composition
    # matrix, recorded by the probe.
    counts = REC['counts']

    cm = {k: v for k, v in REC['clima_mix'].items() if k in counts}
    r_a, nuc_a = report('(a) equilibrium+quench on the climate grid',
                        cm, counts, 0, want)
    pi = {k: np.asarray(v) for k, v in REC['photo_init'].items()
          if k in counts}
    r_b, nuc_b = report('(b) photochemical initial state', pi, counts, 0,
                        want)

    prof = os.path.join(run_dir, 'lower_atmosphere_profile.dat')
    cols = None
    rows = []
    with open(prof) as fh:
        for line in fh:
            if line.startswith('#'):
                if line.startswith('# columns:'):
                    cols = line.split(':', 1)[1].split()
                continue
            rows.append([float(x) for x in line.split()])
    a = np.array(rows)
    print('  (c) converged steady state, from the written profile:')
    for el in ('He', 'C', 'N', 'O'):
        k = 'X_' + el
        if cols and k in cols:
            got = a[0, cols.index(k)]
            print('     %-3s %.12e   dev %+.4e'
                  % (el, got, got/want[el] - 1.0))

    print()
    print('  nitrogen carriers at the deepest level (mixing ratio):')
    for tag, mix in (('(a) clima', cm), ('(b) photo init', pi)):
        tot = 0.0
        parts = []
        for sp in sorted(mix, key=lambda s: -counts[s].get('N', 0)
                         * float(np.asarray(mix[s])[0])):
            k = counts[sp].get('N', 0)
            if k:
                v = k*float(np.asarray(mix[sp])[0])
                tot += v
                if v > 1e-30:
                    parts.append('%s %.6e' % (sp, v))
        print('     %-14s sum(N nuclei) %.9e  ::  %s'
              % (tag, tot, ', '.join(parts[:8])))

    print()
    print('  stoichiometry of the (a)->(b) and (b)->(c) change:')
    for el in ('H', 'He', 'C', 'N', 'O'):
        if el in nuc_a and el in nuc_b:
            print('     d%-3s (a->b) %+.6e  (relative %+.4e)'
                  % (el, nuc_b[el] - nuc_a[el],
                     nuc_b[el]/nuc_a[el] - 1.0))


if __name__ == '__main__':
    main()
