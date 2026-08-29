#!/usr/bin/env python3
"""exhale_to_lart.py -- build a LaRT spherical-illumination Ly-alpha input from an
EXHALE run.

Reads an EXHALE run directory (``input.inp`` + ``output/*_adv.txt``) and writes the
LaRT input files for the ``geometry='spherical_atmosphere'`` +
``source_geometry='stellar_illumination'`` model:

  dens_profile.txt    r [R_p],  n_HI  [cm^-3]   (Ly-alpha scatterer = ground-state H)
  temp_profile.txt    r [R_p],  T     [K]
  velo_profile.txt    r [R_p],  v_r   [km/s]
  line_profile.txt    lambda [A], normalized incident stellar Ly-alpha profile
  <name>_lya.in       LaRT namelist (geometry, star-planet geometry, P_alpha output)
  run.sh              optional multi-host mpirun launcher

Length unit: ``distance2cm = R_p`` so that all LaRT radii/lengths are in R_p, matching
the EXHALE radial grid.  The star-planet geometry (a/R_p, R_star/R_p) is taken from
``input.inp``.

Incident stellar Ly-alpha line profile: a double-Gaussian with peaks at +/-m km/s and
Gaussian width s km/s, following Huang et al. (2017) / Yan et al. (2022)
(defaults m=74, s=49 km/s; total integrated flux ~37,300 erg cm^-2 s^-1 for WASP-52b
is only a normalization and is recorded as a comment).

Typical use
-----------
  python exhale_to_lart.py WASP-52b/fxuv0p25_he98 -o WASP-52b/LaRT_lya --name WASP-52b
  # then:  cd WASP-52b/LaRT_lya && ./run.sh

The scattering rate is requested as a cylindrical 2D grid (``geometry_JPa=2``): the
stellar illumination breaks spherical symmetry, so P_alpha(rho,z) is axisymmetric about
the star-planet axis -- the grid needed to build the 2D H(2p) population for the
H-alpha transit.  Build LaRT with ``make CALCPnew=1`` (-> LaRT_calcJPP.x) to populate it.
"""

import os
import sys
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exhale_io as aio   # noqa: E402

C_KMS = 299792.458
LYA0  = 1215.6701   # Ly-alpha rest wavelength [A], vacuum


def _raw_num(run, label, default=None):
    """Pull a numeric value from input.inp by its label (e.g. 'Stellar radius [R_sun]')."""
    val = run.inp.get('raw', {}).get(label)
    if val is None:
        return default
    try:
        return float(val.split()[0])
    except (ValueError, IndexError):
        return default


def double_gaussian_lya(m_kms=74.0, s_kms=49.0, lam_min=1214.0, lam_max=1217.4,
                        nwav=341, width_is_fwhm=False):
    """Incident stellar Ly-alpha double-Gaussian (Huang 2017 / Yan 2022).

    Two Gaussians centered at +/-m_kms with Gaussian width s_kms.  If
    width_is_fwhm, s_kms is interpreted as FWHM (sigma = FWHM/2.3548); otherwise
    s_kms is the Gaussian standard deviation (sigma), the literal Huang/Yan
    convention.  Returns (lambda [A], normalized profile, peak to 1)."""
    sig = s_kms / 2.35482 if width_is_fwhm else s_kms
    lam = np.linspace(lam_min, lam_max, nwav)
    vel = (lam - LYA0) / LYA0 * C_KMS
    prof = np.exp(-(vel - m_kms)**2 / (2.0 * sig**2)) + \
           np.exp(-(vel + m_kms)**2 / (2.0 * sig**2))
    prof /= prof.max()
    return lam, prof


def write_profiles(run, outdir, rmin=1.0, rmax=None):
    """Write dens/temp/velo profiles (n_HI, T, v_kms) over [rmin, rmax]."""
    r = np.asarray(run.r, float)
    if rmax is None:
        rmax = r.max()
    m = (r >= rmin) & (r <= rmax)
    r = r[m].copy()
    nHI = np.asarray(run.ion['HI'], float)[m]
    T   = np.asarray(run.T, float)[m]
    vkm = np.asarray(run.v, float)[m] / 1.0e5
    r[0] = rmin   # pin the base exactly
    np.savetxt(os.path.join(outdir, 'dens_profile.txt'), np.c_[r, nHI], fmt='%.6e  %.6e')
    np.savetxt(os.path.join(outdir, 'temp_profile.txt'), np.c_[r, T],   fmt='%.6e  %.6e')
    np.savetxt(os.path.join(outdir, 'velo_profile.txt'), np.c_[r, vkm], fmt='%.6e  %.6e')
    return r, nHI, T, vkm, rmax


IN_TEMPLATE = """&parameters
 !=====================================================================
 ! LaRT spherical-illumination Ly-alpha RT for {name}.
 ! Structure (n_HI, T, v_radial) from the EXHALE run:
 !   {src}
 ! Incident stellar Ly-alpha: double-Gaussian, peaks +/-{m:g} km/s,
 ! width(sigma) {s:g} km/s (Huang 2017 / Yan 2022).
 ! distance2cm = R_p, so all radii/lengths are in R_p.
 ! Build LaRT with CALCPnew=1 (-> LaRT_calcJPP.x) for the scattering rate P_alpha.
 !=====================================================================

 par%no_photons  = {nphotons:g}

 par%DGR         = 0.0
 par%use_stokes  = .false.
 par%comoving_source = .false.
 par%save_Jin        = .true.
 par%save_direc0     = .true.
 par%recoil          = .true.
 par%save_all        = .false.

 !-- scattering rate P_alpha output geometry: 2 = cylindrical (axisymmetric about
 !   the star-planet axis), needed for the 2D H(2p) population (H-alpha transit).
 par%geometry_JPa = {geometry_jpa:d}

 !-- incident stellar Ly-alpha spectrum (double-Gaussian; see line_profile.txt) --
 par%spectral_type       = 'line_prof_file'
 par%line_prof_file      = './line_profile.txt'
 par%line_prof_file_type = 1

 !-- spherical planetary atmosphere illuminated by the spherical star --
 par%geometry        = 'spherical_atmosphere'
 par%source_geometry = 'stellar_illumination'
 par%dens_file       = './dens_profile.txt'   ! n_HI  [cm^-3]
 par%temp_file       = './temp_profile.txt'   ! T     [K]
 par%velo_file       = './velo_profile.txt'   ! v_r   [km/s]

 !-- length unit: 1 code unit = R_p = {Rp_cm:.4e} cm --
 par%distance2cm  = {Rp_cm:.4e}

 !-- spectral sampling around Ly-alpha 1215.67 A --
 par%wavelength_min = {lam_min:g}
 par%wavelength_max = {lam_max:g}
 par%nwavelength    = {nwav:d}

 !-- Cartesian grid covering r = {rmin:g}..{rmax:g} R_p --
 par%nx     = {ngrid:d}
 par%ny     = {ngrid:d}
 par%nz     = {ngrid:d}
 par%xmax   = {rmax:g}
 par%ymax   = {rmax:g}
 par%zmax   = {rmax:g}
 par%rmin   = {rmin:g}
 par%rmax   = {rmax:g}

 !-- star-planet geometry (in R_p units) --
 par%distance_star_to_planet = {a_Rp:.4f}
 par%stellar_radius          = {Rstar_Rp:.4f}

 !-- distance to the system (only scales output images, not internal P_alpha) --
 par%distance = {distance_cm:.3e}

 !-- observer & detector (peel-off; transit view along the star-planet axis) --
 par%save_peeloff_2D = .true.
 par%save_peeloff_3D = .false.
 par%nxim = {ngrid:d}
 par%nyim = {ngrid:d}
 par%out_bitpix = -64
 par%alpha = 0.0
 par%beta  = 90.0
 par%gamma = 90.0

 par%out_merge = .true.
 par%nprint = 1e6
/
"""

RUN_SH = """#!/bin/bash
# Run LaRT spherical-illumination Ly-alpha RT for {name}.
exec < /dev/null 2>&1
trap "" HUP

EXEC={exe}

# Hosts to run on, and the MPI slots each provides.  Both are site-specific:
# override with LART_HOSTS / LART_SLOTS, or per host as "name:slots".
HOSTS=${{LART_HOSTS:-{hosts}}}
SLOTS=${{LART_SLOTS:-{slots}}}
host_file=/tmp/host_file_$RANDOM
for host in $(echo $HOSTS | tr "," "\\n"); do
   case $host in
      *:*) echo $host >> $host_file ;;
      *)   echo $host:$SLOTS >> $host_file ;;
   esac
done

echo "Running $EXEC on $HOSTS"
echo "   machinefile $host_file:"; cat $host_file
mpirun -machinefile $host_file $EXEC {infile}
"""


def build(run_dir, outdir, name=None, m_kms=74.0, s_kms=49.0, width_is_fwhm=False,
          rmax=None, ngrid=201, nwav=161, nphotons=1e7, geometry_jpa=2,
          exe=None, hosts=None, slots=None,
          distance_pc=174.0, adv=True, write_run_sh=True):
    """Build a complete LaRT spherical-illumination input set from an EXHALE run.

    LaRT lives outside this repository, so the executable and the MPI hosts are
    site-specific: they come from --exe/--hosts/--slots, or from the
    environment variables LART_EXE, LART_HOSTS and LART_SLOTS.
    """
    exe = exe or os.environ.get('LART_EXE', 'LaRT_calcJPP.x')
    hosts = hosts or os.environ.get('LART_HOSTS', 'localhost')
    slots = slots or os.environ.get('LART_SLOTS', '1')
    run_dir = os.path.abspath(run_dir)
    name = name or os.path.basename(run_dir.rstrip('/'))
    run = aio.load_run(os.path.join(run_dir, 'output'),
                       os.path.join(run_dir, 'input.inp'), adv=adv)
    os.makedirs(outdir, exist_ok=True)

    Rp_cm = run.inp['Rp_RJ'] * aio.RJ
    a_cm  = run.inp['a_AU'] * aio.AU
    Rstar = _raw_num(run, 'Stellar radius [R_sun]', 1.0) * 6.957e10  # R_sun [cm]
    a_Rp     = a_cm / Rp_cm
    Rstar_Rp = Rstar / Rp_cm

    r, nHI, T, vkm, rmax = write_profiles(run, outdir, rmax=rmax)

    lam, prof = double_gaussian_lya(m_kms, s_kms, nwav=341, width_is_fwhm=width_is_fwhm)
    np.savetxt(os.path.join(outdir, 'line_profile.txt'), np.c_[lam, prof], fmt='%8.3f  %.6e')

    infile = '%s_lya.in' % name
    with open(os.path.join(outdir, infile), 'w') as f:
        f.write(IN_TEMPLATE.format(
            name=name, src=run_dir, m=m_kms, s=s_kms, nphotons=nphotons,
            geometry_jpa=geometry_jpa, Rp_cm=Rp_cm, lam_min=lam.min(), lam_max=lam.max(),
            nwav=nwav, ngrid=ngrid, rmin=1.0, rmax=rmax, a_Rp=a_Rp, Rstar_Rp=Rstar_Rp,
            distance_cm=distance_pc * 3.0856775814913673e18))

    if write_run_sh:
        runp = os.path.join(outdir, 'run.sh')
        with open(runp, 'w') as f:
            f.write(RUN_SH.format(name=name, exe=exe, hosts=hosts, slots=slots,
                                  infile=infile))
        os.chmod(runp, 0o755)

    print('EXHALE -> LaRT input written to %s/' % outdir)
    print('  source run : %s' % run_dir)
    print('  n_HI: %.2e (base) -> %.2e (r=%.1f);  T: %.0f..%.0f K;  v: %.1f..%.1f km/s'
          % (nHI[0], nHI[-1], rmax, T.min(), T.max(), vkm.min(), vkm.max()))
    print('  distance2cm = R_p = %.4e cm;  a/R_p = %.3f;  R_star/R_p = %.3f'
          % (Rp_cm, a_Rp, Rstar_Rp))
    print('  Ly-alpha line: double-Gaussian +/-%g km/s, %s=%g km/s'
          % (m_kms, 'FWHM' if width_is_fwhm else 'sigma', s_kms))
    print('  namelist: %s   (geometry_JPa=%d, %g photons, grid %d^3)'
          % (infile, geometry_jpa, nphotons, ngrid))
    return outdir


def main():
    p = argparse.ArgumentParser(description='Build a LaRT spherical-illumination Ly-alpha '
                                            'input from an EXHALE run.')
    p.add_argument('run_dir', help='EXHALE run directory (contains input.inp and output/)')
    p.add_argument('-o', '--outdir', required=True, help='output directory for LaRT input')
    p.add_argument('--name', default=None, help='model name (default: run-dir basename)')
    p.add_argument('--m', type=float, default=74.0, help='Ly-alpha double-Gaussian peak |v| [km/s]')
    p.add_argument('--s', type=float, default=49.0, help='Ly-alpha Gaussian width [km/s]')
    p.add_argument('--fwhm', action='store_true', help='interpret --s as FWHM (else sigma)')
    p.add_argument('--rmax', type=float, default=None, help='outer radius [R_p] (default: grid max)')
    p.add_argument('--ngrid', type=int, default=201, help='Cartesian grid size per axis')
    p.add_argument('--nwav', type=int, default=161, help='number of wavelength bins')
    p.add_argument('--nphotons', type=float, default=1e7, help='number of MC photons')
    p.add_argument('--geometry_jpa', type=int, default=2, choices=(1, 2, 3),
                   help='P_alpha output geometry: 1 spherical, 2 cylindrical, 3 full 3D')
    p.add_argument('--exe', default=None,
                   help='LaRT executable (default: $LART_EXE, else LaRT_calcJPP.x on PATH)')
    p.add_argument('--hosts', default=None,
                   help='comma-separated MPI hosts (default: $LART_HOSTS, else localhost)')
    p.add_argument('--slots', default=None,
                   help='MPI slots per host (default: $LART_SLOTS, else 1)')
    p.add_argument('--distance_pc', type=float, default=174.0, help='system distance [pc]')
    p.add_argument('--eq', action='store_true', help='use eq (non-_adv) EXHALE files')
    p.add_argument('--no-run-sh', action='store_true', help='do not write run.sh')
    a = p.parse_args()
    build(a.run_dir, a.outdir, name=a.name, m_kms=a.m, s_kms=a.s, width_is_fwhm=a.fwhm,
          rmax=a.rmax, ngrid=a.ngrid, nwav=a.nwav, nphotons=a.nphotons,
          geometry_jpa=a.geometry_jpa, exe=a.exe, hosts=a.hosts, slots=a.slots,
          distance_pc=a.distance_pc,
          adv=not a.eq, write_run_sh=not a.no_run_sh)


if __name__ == '__main__':
    main()
