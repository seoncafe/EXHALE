#!/usr/bin/env python3
"""T2: hold the lower-atmosphere handoff fixed and solve only the wind.

The two converged states of `kzz_profile_scan/kzz1e7` (He/H = 3.4, cool, and
He/H = 3.399, hot) carry handoffs that state the same base to within 1.5e-4
in He/H and 4e-7 in T0.  If the two states survive when the SAME handoff is
imposed on both seeds, the two states are a property of the steady wind
equations at a fixed boundary condition, and not an attractor of the
composition Picard loop that wraps them.

Runs the 2x2: profile in {cool, hot} times seed in {cool, hot}.  Each cell is
one `solve_escape_wind` call -- the same JFNK-under-PTC solve the closure
uses, with no outer composition pass.

Usage:  python3 T2_fixed_composition.py
"""
import os
import sys
import json

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.abspath(os.path.join(HERE, '../../..'))
sys.path.insert(0, os.path.join(EX, 'src/utils'))
import shutil                                               # noqa: E402
import element_flux_closure as efc                          # noqa: E402

SCAN = os.path.join(EX, 'LHS1140b/exhale/kzz_profile_scan/kzz1e7')
COOL = os.path.join(SCAN, 'heh3p4/k02')       # log10 Mdot 7.54, EW 0.4287
HOT = os.path.join(SCAN, 'heh3p399/k03')      # log10 Mdot 7.87, EW 1.3248


def main():
    root = os.path.join(HERE, 'T2_fixedBC')
    os.makedirs(root, exist_ok=True)
    results = []
    for pname, pdir in (('coolprof', COOL), ('hotprof', HOT)):
        cfgsrc = os.path.join(os.path.dirname(pdir), 'closure.json')
        for sname, sdir in (('coolseed', COOL), ('hotseed', HOT)):
            tag = '%s_%s' % (pname, sname)
            d = os.path.join(root, tag)
            if os.path.isfile(os.path.join(d, 'DONE')):
                print('%s already done, skipping' % tag)
                results.append(json.load(
                    open(os.path.join(d, 'DONE'))))
                continue
            os.makedirs(d, exist_ok=True)
            cfg = efc.read_configuration(cfgsrc)
            # the handoff: profile + its base.inp, taken as one pair
            for f in ('lower_atmosphere_profile.dat', 'base.inp'):
                shutil.copyfile(os.path.join(pdir, f), os.path.join(d, f))
            print('=== %s : profile from %s, seed from %s' % (tag, pdir, sdir))
            try:
                outd, info, mdot = efc.solve_escape_wind(
                    cfg, d, os.path.join(sdir, 'output'), print)
            except Exception as exc:                        # noqa: BLE001
                print('   FAILED: %s' % exc)
                results.append({'tag': tag, 'info': None, 'log10_Mdot': None,
                                'error': str(exc)})
                continue
            rec = {'tag': tag, 'profile': pdir, 'seed': sdir,
                   'info': info, 'log10_Mdot': mdot}
            json.dump(rec, open(os.path.join(d, 'DONE'), 'w'), indent=2)
            results.append(rec)
            print('   info=%s  log10 Mdot = %.3f' % (info, mdot))

    print()
    print('%-24s %6s %12s' % ('cell', 'info', 'log10 Mdot'))
    for r in results:
        print('%-24s %6s %12s'
              % (r['tag'], r['info'],
                 'n/a' if r.get('log10_Mdot') is None
                 else '%.3f' % r['log10_Mdot']))


if __name__ == '__main__':
    main()
