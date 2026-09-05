#!/usr/bin/env python3
"""Build src/modules/lower_atmosphere/h2_self_shielding_table.f90 from a
LINE-BY-LINE calculation in which the H2 Lyman-Werner lines absorb each
other's beam.

Usage
-----
    python3 h2_shielding_table_line_by_line.py <abs_dir> <cloudy_run_dir> \\
            [<out.f90>]

<abs_dir> holds abs_T####.npz, one per temperature of the grid, each carrying
sigma_pump(N) and p_single(N) from lbl_table.absolute_curves.
<cloudy_run_dir> holds the same 27 "save h2 rates" files the CLOUDY table was
built from; they are read ONLY for the trapping ratio p_eff/p_single.

WHY THIS REPLACED THE CLOUDY TABLE.  CLOUDY shields one line at a time
(Federman, Glassgold & Kwan 1979, applied one transition at a time), and at
the columns
a planetary molecular base reaches that is not a correction but the dominant
term: at T = 1300 K, n_H = 1e13 cm^-3, letting the lines absorb each other's
beam lowers the surviving pumping by a factor 24 at N_H2 = 1e21 cm^-2 and 630
at 4.4e21.  docs/p38_line_overlap_shielding.md is the measurement, and
Update_EXHALE.md section 135 the changelog entry.

WHAT IS COMPUTED HERE AND WHAT IS NOT.  Everything line overlap touches is
computed here: the pumping rate of every line in the attenuated beam, hence
sigma_pump, and the pump-weighted branching p_single.  The fluorescent
trapping that turns p_single into p_eff is NOT: it is a property of the
escape probabilities of the decay lines and of the slab geometry, both of
which CLOUDY solves (one transition at a time, both directions, iterated), and
a second unvalidated implementation of it would be worse than reusing the
first.  The ratio p_eff/p_single is read node by node from the CLOUDY runs and
multiplied in.  Its geometry is plane-parallel, illuminated on one face and
closed on the other, while the layer is spherical and open outward, so the
trapping is the one ingredient of this table computed in a geometry the layer
does not have; its provenance name is slab_surrogate.  Measured, from a matched
slab/sphere calculation on the same line data (docs/e2_lw_geometry.md): the
plane-parallel trapping over-predicts p_eff by 0.2 per cent at r/H = 41 and by
up to 20 per cent at r/H = 16, one-signed, 6.5 per cent weighted by the
dissociation rate over a hot-Uranus molecular layer -- smaller than the 20-40
per cent between two escape-probability methods in the same slab, which is why
it is not rebuilt.  The density response of the same quantity is 13 per cent
between n_H = 1e12 and 1e14 at fixed T = 1300 K.
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import importlib.util
_spec = importlib.util.spec_from_file_location(
    'h2_shielding_table_from_cloudy',
    os.path.join(_HERE, 'h2_shielding_table_from_cloudy.py'))
_cl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_cl)

T_GRID, LOGN_GRID, LOGCOL_GRID, NCOL = (_cl.T_GRID, _cl.LOGN_GRID,
                                        _cl.LOGCOL_GRID, _cl.NCOL)


def trapping_ratio(run_dir, T, ln):
    """p_eff/p_single of the CLOUDY run at this node, on the column axis."""
    name = 't%04d_n%d' % (int(round(T)), int(round(ln)))
    N, sig, ps, pe = _cl.curves(os.path.join(run_dir, name + '.rat'))
    return 10.0 ** np.interp(LOGCOL_GRID, np.log10(N), np.log10(pe / ps))


def build(abs_dir, run_dir):
    shape = (NCOL, len(T_GRID), len(LOGN_GRID))
    lsig = np.zeros(shape); lps = np.zeros(shape); lpe = np.zeros(shape)
    for it, T in enumerate(T_GRID):
        d = np.load(os.path.join(abs_dir, 'abs_T%04d.npz' % int(round(T))))
        if not np.allclose(d['logN'], LOGCOL_GRID):
            raise SystemExit('abs_T%04d.npz is on a different column axis'
                             % int(round(T)))
        sp, ps = d['sigma_pump'], d['p_single']
        for jn, ln in enumerate(LOGN_GRID):
            pe = np.minimum(ps * trapping_ratio(run_dir, T, ln), 1.0)
            lsig[:, it, jn] = np.log10(sp * pe)
            lps[:, it, jn] = np.log10(ps)
            lpe[:, it, jn] = np.log10(pe)
    return lsig, lps, lpe


def emit(lsig, lps, lpe, out_path):
    L = [HEADER.format(nt=len(T_GRID), nn=len(LOGN_GRID), ncol=NCOL), '']
    L.append('      ! Grid axes, log10 of T [K], n_H [cm^-3] and N_H2 [cm^-2].')
    L.append(_cl.fortran_array('h2_shield_log_temp', np.log10(T_GRID)))
    L.append(_cl.fortran_array('h2_shield_log_dens', LOGN_GRID))
    L.append(_cl.fortran_array('h2_shield_log_col', LOGCOL_GRID))
    L.append('')
    L.append('      ! log10 sigma_diss [cm^2] on (N_H2, T, n_H): dissociations per unit')
    L.append('      ! INCIDENT band photon fluence, with the lines absorbing each other')
    L.append('      ! and with the fluorescent trapping inside it.')
    L.append(_cl.cube('h2_shield_log_sigma_diss', lsig))
    L.append('')
    L.append('      ! log10 p_single: dissociations per pump with every fluorescent')
    L.append('      ! photon free to escape.  The vibrational heat return needs this one.')
    L.append(_cl.cube('h2_shield_log_p_single', lps))
    L.append('')
    L.append('      ! log10 p_eff: dissociations per band photon removed from the beam.')
    L.append('      ! The band photon ledger needs this one.')
    L.append(_cl.cube('h2_shield_log_p_eff', lpe))
    L.append(FOOTER)
    open(out_path, 'w').write('\n'.join(L))
    print('wrote %s: 3 x %d x %d x %d = %d values'
          % (out_path, NCOL, len(T_GRID), len(LOGN_GRID),
             3 * NCOL * len(T_GRID) * len(LOGN_GRID)))




HEADER = """      module h2_self_shielding_table
      ! H2 Lyman-Werner photodissociation with OVERLAPPING LINES: the absolute
      ! dissociation cross section per incident band photon, and the two
      ! branching ratios that go with it.
      !
      ! GENERATED FILE.  Rebuild with
      !     python3 src/utils/h2_shielding_table_line_by_line.py <abs_dir> <run_dir>
      ! Edit the script, not this file.
      !
      ! ---------------------------------------------------------------
      ! WHAT IS TABULATED, AND WHY THREE THINGS AND NOT ONE
      !
      ! A Lyman-Werner band photon absorbed by H2 lifts it into B, C, B' or D.
      ! From there the molecule either dissociates into the X continuum or
      ! fluoresces back into a bound level of X.  Once the layer is optically
      ! thick the fluorescent photon can be re-absorbed, which starts the cycle
      ! again, so THREE different counts exist and the code needs all three:
      !
      !   sigma_diss  dissociations per unit INCIDENT band photon fluence
      !               [cm^2].  Line self-shielding, LINE OVERLAP and the
      !               trapping of the fluorescent photons are all inside it.
      !               This sets the rate, and with it the 0.4 eV fragment
      !               heating.
      !   p_single    dissociations per PUMP with every fluorescent photon
      !               free to escape, D/(D + sum A_ul).  The number of
      !               fluorescent decays per dissociation is
      !               (1 - p_single)/p_single WHETHER OR NOT the photons are
      !               trapped -- a trapped decay plus its re-absorption moves
      !               one molecule down and another up and deposits nothing,
      !               so the trapping cancels out of this count.  This is the
      !               vibrational heat return.
      !   p_eff       dissociations per photon REMOVED FROM THE BEAM, i.e.
      !               per fresh pump, D/(D + sum A_ul Ploss).  Only fresh
      !               pumps take photons out of the stellar beam, so this is
      !               the band photon ledger.
      !
      ! sigma_diss = p_eff x sigma_pump, and sigma_pump -- the pump cross
      ! section -- carries no trapping at all.
      !
      ! ---------------------------------------------------------------
      ! WHY A LINE-BY-LINE CALCULATION AND NOT CLOUDY, AND NOT A FIT
      !
      ! Both published fits and the level-resolved CLOUDY calculation that
      ! stood here until 2026-09-03 shield ONE LINE AT A TIME.  At the columns
      ! a planetary molecular base reaches that is not a correction but the
      ! dominant term.  Measured, at n_H = 1e13 cm^-3, as the ratio of the
      ! surviving pumping with the lines absorbing each other's beam to the
      ! same calculation with each line alone:
      !
      !     N_H2 [cm^-2]     900 K    1300 K    1800 K    2700 K
      !     1e19             0.73      0.62      0.47      0.30
      !     1e20             0.45      0.33      0.22      0.096
      !     1e21             0.080     0.042     0.019     0.0040
      !     4.4e21           0.0047    0.0016    0.00056   0.000048
      !
      ! and the base of the molecular layer sits at N_H2 = 4.4e21.  Two
      ! independent codes agree on this: the calculation tabulated here, and
      ! the Meudon PDR code (Le Petit et al. 2006) with its exact UV line
      ! transfer, which solves the transfer equation on one wavelength grid
      ! carrying every H2 line at once.  They agree to 5-22 per cent over
      ! 1e19 <= N_H2 <= 1e21 at 900, 1300 and 1800 K.
      !
      ! ---------------------------------------------------------------
      ! WHAT WAS CALCULATED
      !
      ! Plane-parallel slab, normally incident beam, flat F_lambda
      ! (f_nu ~ nu^-2, the Draine & Bertoldi 1996 eq. 24 spectrum that
      ! sigma_lw is calibrated to) normalized to F = 343 erg cm^-2 s^-1 over
      ! 912-1110 A, the same band and the same mean band photon energy
      ! e_lw_photon_erg that lyman_werner.f90 converts F_LW with.
      !
      !   tau(nu, N) = N sum_i x_i sigma_i(nu)
      !   pump_i(N)  = x_i (pi e^2/m c) f_i int phi_i(nu) F_nu e^-tau /(h nu) dnu
      !
      ! summed over every Lyman and Werner transition with 911.75 A <= lambda
      ! <= 1200 A, on a 4e5-point frequency grid, with a Voigt profile
      ! (scipy wofz inside +-300 Doppler widths, the exact Lorentzian outside)
      ! and b = sqrt(2kT/m_H2), the same Doppler parameter as
      ! h2_doppler_parameter.  Level populations are LTE at the imposed T,
      ! which is what CLOUDY's own runs show at these densities (dep coef =
      ! 1.000 at every level) and what the Meudon runs show as well.
      !
      ! THE LINE DATA.  Abgrall, Roueff & Drira (2000), A&AS 141, 297,
      ! "Total transition probability and spontaneous radiative dissociation
      ! of B, C, B' and D states of molecular hydrogen": the level energies,
      ! total decay rates and dissociation probabilities of the six excited
      ! electronic states, and the individual line probabilities of their
      ! electronic companion.  These are the SAME data CLOUDY's dissprob_*
      ! and transprob_* files carry and the same the Meudon PDR code uses,
      ! so no line-list difference enters the comparison.
      !
      ! ---------------------------------------------------------------
      ! THE ONE PIECE THAT IS NOT LINE-BY-LINE: THE TRAPPING
      !
      ! p_eff/p_single -- the suppression of the branching by re-absorption of
      ! the fluorescent decay photons -- is taken node by node from the same
      ! CLOUDY runs the previous table was built from, and multiplied into
      ! p_single.  Reimplementing CLOUDY's escape probabilities would be a
      ! second, unvalidated version of the same thing.  ITS GEOMETRY IS
      ! PLANE-PARALLEL, ILLUMINATED ON ONE FACE AND CLOSED ON THE OTHER, and
      ! it is the one ingredient of this table computed in a geometry the
      ! layer does not have: the layer is spherical and open outward.  The
      ! provenance name for it is SLAB_SURROGATE.
      !
      ! HOW MUCH THAT IS WORTH, MEASURED.  Slab and sphere were computed from
      ! this same line list, these same LTE populations and this same
      ! frequency integral, with a single-flight escape probability under
      ! complete redistribution, so that only the set of ray columns over
      ! solid angle differs (docs/e2_lw_geometry.md).  The plane-parallel
      ! trapping OVER-PREDICTS p_eff, by 0.2 per cent at r/H = 41 and by up
      ! to 20 per cent at r/H = 16, one-signed everywhere, and by 6.5 per
      ! cent weighted by the dissociation rate over a hot-Uranus molecular
      ! layer.  Escape is dominated by the near-radial directions, where the
      ! two geometries have the SAME column, which is why the effect is
      ! small; the difference is the near-tangential rays, which a sphere
      ! caps at sqrt(pi r / 2H) times the radial column where a slab lets
      ! them diverge, plus the lit sliver past the terminator.  The
      ! controlling parameter is therefore r/H, not the column.
      !
      ! WHY IT IS NOT REBUILT.  Two escape-probability methods in the SAME
      ! slab -- CLOUDY's, taken one transition at a time and iterated, and
      ! the single-flight one above -- differ by 20 to 40 per cent, and not
      ! with one sign (900 K: this work 20-40 per cent higher; 2700 K: up to
      ! 35 per cent lower at depth).  The method is the leading term and the
      ! geometry is not, so a spherical rebuild of the trapping would move
      ! sigma_diss by less than the uncertainty already inside it.  The
      ! density axis remains the only axis the trapping ratio has; an r/H
      ! axis would span 0.90-1.00 over the H2-bearing layer and is not
      ! carried.  The earlier statement of this sensitivity -- 13 per cent
      ! between n_H = 1e12 and 1e14 at fixed T = 1300 K -- is the density
      ! response of the same quantity and still holds.
      !
      ! ---------------------------------------------------------------
      ! CHECKS
      !
      !  - At 100 K the overlapping-line calculation reproduces the Draine &
      !    Bertoldi (1996) self-shielding function to 2-9 per cent over
      !    1e16 <= N_H2 <= 1e21.  That fit was built from a calculation that
      !    includes overlap (their sec. 5.2), and it is the one published
      !    overlap-including curve there is to check against.
      !  - sigma_pump at the BOTTOM of the column axis (N_H2 = 1e12 cm^-2,
      !    the optically thin end) is 2.59e-17 (700 K) to 2.88e-17 cm^2
      !    (3200 K), against DB96's own 3.4520e-18/0.135 =
      !    2.557e-17 -- 1 to 13 per cent.
      !  - p_single at the same place is 0.1425 (700 K) to 0.1592 (3200 K),
      !    against 0.1465 from the CLOUDY runs at 1300 K, where this
      !    calculation gives 0.1475: 0.7 per cent.
      !  - sigma_diss at the same place agrees with the CLOUDY-based table it
      !    replaces to 1.3-3.6 per cent, which is the statement that the two
      !    differ in the shielding and not in the normalization.
      !
      ! ---------------------------------------------------------------
      ! VALIDITY, AND WHERE IT STOPS
      !
      !  - TEMPERATURE 700-3200 K and DENSITY 1e12-1e14 cm^-3 are the grid.
      !    Outside it the edge value is used (clamped), which is a statement
      !    that nothing was calculated there, not that the value holds.  The
      !    density axis enters ONLY through the trapping ratio; the
      !    line-by-line pumping has no density in it at all, because the
      !    populations are LTE.
      !  - COLUMN 1e12-5e21 cm^-2 is the grid.  h2_shield_max_column() is the
      !    top of it; above that the edge value is returned.  THERE IS NO
      !    LONGER A COLUMN ABOVE WHICH THE TABLE IS AN UPPER BOUND: line
      !    overlap is inside it.  The table that stood here until 2026-09-03
      !    carried h2_shielding_overlap_column for that purpose; it is gone.
      !  - RADIATION FIELD.  One G_0 only.  Dropping F_LW by four decades
      !    moves the shielding by at most 9% anywhere on the grid.
      !  - The ortho/para ratio is the LTE one at the imposed temperature.
      !
      ! ---------------------------------------------------------------
      ! Reference: Draine & Bertoldi (1996) ApJ 468, 269 (the band, the
      ! normalization, the equivalent width, and the line-overlap statement of
      ! their sec. 5.2); Abgrall, Roueff & Drira (2000) A&AS 141, 297 (the
      ! line data); Le Petit et al. (2006) ApJS 164, 506 (the Meudon PDR code
      ! this was checked against); Ferland et al. (2017) RMxAA 53, 385
      ! (CLOUDY, the trapping).  docs/p38_line_overlap_shielding.md is the
      ! measurement.

      implicit none
      private

      public :: h2_lw_dissociation_cross_section,                        &
                h2_lw_dissociation_per_pump,                             &
                h2_lw_dissociation_per_absorbed_photon,                  &
                h2_self_shielding_level_resolved,                        &
                h2_shield_max_column

      integer, parameter :: n_temp_sh = {nt}, n_dens_sh = {nn}, n_col_sh = {ncol}
"""

FOOTER = """
      contains

      ! Trilinear interpolation of a log10 table on (N_H2, T, n_H), clamped at
      ! every edge rather than extrapolated.
      double precision function h2_shield_interp(tab, N_H2, T, n_H)         &
                                result(v)
      real*8, intent(in) :: tab(n_col_sh,n_temp_sh,n_dens_sh)
      real*8, intent(in) :: N_H2, T, n_H
      real*8  :: wt, wn, wc, v00, v01, v10, v11, v0, v1
      integer :: it, in, ic

      call h2_shield_locate(h2_shield_log_temp, n_temp_sh,                  &
                            log10(max(T, 1.0d0)), it, wt)
      call h2_shield_locate(h2_shield_log_dens, n_dens_sh,                  &
                            log10(max(n_H, 1.0d-30)), in, wn)
      call h2_shield_locate(h2_shield_log_col, n_col_sh,                    &
                            log10(max(N_H2, 1.0d-30)), ic, wc)

      v00 = (1.0d0 - wc)*tab(ic,it,in)     + wc*tab(ic+1,it,in)
      v01 = (1.0d0 - wc)*tab(ic,it,in+1)   + wc*tab(ic+1,it,in+1)
      v10 = (1.0d0 - wc)*tab(ic,it+1,in)   + wc*tab(ic+1,it+1,in)
      v11 = (1.0d0 - wc)*tab(ic,it+1,in+1) + wc*tab(ic+1,it+1,in+1)
      v0  = (1.0d0 - wn)*v00 + wn*v01
      v1  = (1.0d0 - wn)*v10 + wn*v11
      v   = 10.0d0**((1.0d0 - wt)*v0 + wt*v1)
      end function h2_shield_interp

      ! Dissociations per unit incident band photon fluence [cm^2] behind a
      ! star-ward H2 column N_H2 [cm^-2] at gas temperature T [K] and hydrogen
      ! nucleus density n_H [cm^-3].  Line self-shielding and the trapping of
      ! the fluorescent decay photons are both in it; the continuum of the
      ! same interval is not (that is tau_cont in lyman_werner.f90).
      double precision function h2_lw_dissociation_cross_section(N_H2, T,   &
                                n_H) result(sigma)
      real*8, intent(in) :: N_H2, T, n_H
      sigma = h2_shield_interp(h2_shield_log_sigma_diss, N_H2, T, n_H)
      end function h2_lw_dissociation_cross_section

      ! Dissociations per PUMP with every fluorescent photon free to escape,
      ! D/(D + sum A_ul).  The count of fluorescent decays per dissociation is
      ! (1 - p)/p with THIS p and not with the effective one: a trapped decay
      ! and the re-absorption that follows it move one molecule down and
      ! another up, depositing nothing, so trapping cancels out of the count
      ! (module header).  Used by the vibrational heat return.
      double precision function h2_lw_dissociation_per_pump(N_H2, T, n_H)   &
                                result(p)
      real*8, intent(in) :: N_H2, T, n_H
      p = h2_shield_interp(h2_shield_log_p_single, N_H2, T, n_H)
      end function h2_lw_dissociation_per_pump

      ! Dissociations per band photon REMOVED FROM THE BEAM, i.e. per fresh
      ! pump, D/(D + sum A_ul Ploss).  Only fresh pumps take photons out of
      ! the stellar beam, so the band photon ledger divides the dissociation
      ! rate by this and not by the single-pump branching.
      double precision function h2_lw_dissociation_per_absorbed_photon(     &
                                N_H2, T, n_H) result(p)
      real*8, intent(in) :: N_H2, T, n_H
      p = h2_shield_interp(h2_shield_log_p_eff, N_H2, T, n_H)
      end function h2_lw_dissociation_per_absorbed_photon

      ! The suppression of the dissociation rate relative to the top of the
      ! same column, sigma_diss(N)/sigma_diss(N_min).  It is a DIAGNOSTIC --
      ! the rate itself now comes from the cross section above -- and it is
      ! what output/Lyman_Werner.txt reports.  At the bottom of the column
      ! axis the tabulated factor is 1 to 0.1 per cent, so the ratio is the
      ! optically thin limit there.
      double precision function h2_self_shielding_level_resolved(N_H2, T,   &
                                n_H) result(f_shield)
      real*8, intent(in) :: N_H2, T, n_H
      real*8 :: N_bottom
      N_bottom = 10.0d0**h2_shield_log_col(1)
      f_shield = h2_shield_interp(h2_shield_log_sigma_diss, N_H2, T, n_H)   &
               / h2_shield_interp(h2_shield_log_sigma_diss, N_bottom, T,    &
                                  n_H)
      end function h2_self_shielding_level_resolved

      ! Top of the tabulated column axis [cm^-2].  Above it the table is
      ! CLAMPED -- the edge value is returned, which is a statement that
      ! nothing was calculated further in, not that the value holds.  There is
      ! no longer a column above which the table is an upper bound: line
      ! overlap is inside it (module header).
      double precision function h2_shield_max_column() result(N_max)
      N_max = 10.0d0**h2_shield_log_col(n_col_sh)
      end function h2_shield_max_column


      ! Bracket x in the ascending grid g(1:n): returns the lower index i and
      ! the fractional position w in [0,1] within [g(i), g(i+1)].  Outside the
      ! grid w is clamped to 0 or 1, so the edge value is returned rather than
      ! an extrapolation.
      subroutine h2_shield_locate(g, n, x, i, w)
      integer, intent(in)  :: n
      real*8,  intent(in)  :: g(n), x
      integer, intent(out) :: i
      real*8,  intent(out) :: w
      integer :: k
      if (x .le. g(1)) then
         i = 1
         w = 0.0d0
         return
      endif
      if (x .ge. g(n)) then
         i = n - 1
         w = 1.0d0
         return
      endif
      i = 1
      do k = 1,n-1
         if (x .ge. g(k) .and. x .le. g(k+1)) then
            i = k
            exit
         endif
      enddo
      w = (x - g(i))/(g(i+1) - g(i))
      end subroutine h2_shield_locate

      ! End of module
      end module h2_self_shielding_table
"""


if __name__ == '__main__':
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    out = (sys.argv[3] if len(sys.argv) > 3 else
           os.path.join(_HERE, '..', 'modules', 'lower_atmosphere',
                        'h2_self_shielding_table.f90'))
    emit(*build(sys.argv[1], sys.argv[2]), out)
