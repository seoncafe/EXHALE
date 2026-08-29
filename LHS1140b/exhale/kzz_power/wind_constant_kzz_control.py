#!/usr/bin/env python3
"""Branch A of the K_zz(p) test: chemistry on K_zz(p), wind on K_zz = const.

The lower-atmosphere profile is solved on the pressure-dependent eddy
coefficient, and then its `Kzz` column alone is overwritten with the
constant the reference arm used before the wind reads it.  EXHALE takes
K_zz from that column (`eddy_diffusion_on_grid`), and above the profile's
top level every cell inherits the top value, so with the column made
constant the wind sees exactly the eddy coefficient it saw in the reference
run.  What is left different is the composition and the temperature the
chemistry produced -- which is the quantity the test is about.

Usage:
  python3 wind_constant_kzz_control.py <source_kdir> <dest_kdir> <K_zz>
      <seed_output_dir> <closure.json>
"""
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../../../src/utils'))
import element_flux_closure as efc                          # noqa: E402


def constant_eddy_column(src, dst, kzz):
    """Copy the profile, replacing the `Kzz` column with one value."""
    icol, out = None, []
    for line in open(src):
        if line.startswith('#'):
            if line.startswith('# columns:'):
                icol = line.split(':', 1)[1].split().index('Kzz')
            out.append(line)
            continue
        f = line.split()
        f[icol] = '%25.17E' % kzz
        out.append(' '.join(f) + '\n')
    out.insert(1, '# NOTE the Kzz column was overwritten with the constant'
                  ' %.6E cm^2/s after the chemistry was solved on the'
                  ' pressure-dependent profile: this file states the'
                  ' composition of a K_zz(p) column and the eddy coefficient'
                  ' of a constant-K_zz one, on purpose, so that the wind is'
                  ' held fixed while the chemistry changes.\n' % kzz)
    open(dst, 'w').writelines(out)
    return icol


def main():
    srck, dstk, kzz, seed, cfgpath = sys.argv[1:6]
    kzz = float(kzz)
    cfg = efc.read_configuration(cfgpath)
    os.makedirs(os.path.join(dstk, 'output'), exist_ok=True)
    constant_eddy_column(os.path.join(srck, cfg['profile_name']),
                         os.path.join(dstk, cfg['profile_name']), kzz)
    shutil.copyfile(os.path.join(srck, 'base.inp'),
                    os.path.join(dstk, 'base.inp'))
    out, info, mdot = efc.solve_escape_wind(cfg, dstk, seed, print)
    print('info = %s, log10 Mdot = %.3f' % (info, mdot))
    return 0 if info == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
