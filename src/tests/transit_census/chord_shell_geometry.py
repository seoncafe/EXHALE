"""Chord geometry of the transit post-processor, and what the census of it says.

A ray at impact parameter b crosses every shell of radius r >= b, at
line-of-sight coordinate |x| = sqrt(r^2 - b^2).  The selection has no upper
radial limit, so a shell lying outside the stellar disk still absorbs on the
rays that fall on the disk.  These rows pin that geometry on the ONE selection
the optical-depth integrals and the refusal census share, and pin the
line-center optical-depth share the census reports for the rows whose steady
advective correction was refused.

Every row prints one
    PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
line; lines beginning with two spaces or with DIAGNOSTIC are context, not
verdicts.  Exit status is nonzero if any row fails.
"""

import os
import sys

import numpy as np

# The rows that need the library.  The script census row stands on its own, so
# it is not repeated here when the import fails.
ROWS = ('outer_shell_crosses_the_inner_ray',
        'chord_crossing_of_the_outer_shell',
        'shell_inside_the_ray_is_not_crossed',
        'census_counts_the_shell_above_the_disk_cap',
        'integral_absorbs_in_shells_above_the_disk_cap',
        'refused_share_is_one_when_every_row_is_refused',
        'refused_share_is_zero_when_no_row_is_refused',
        'refused_and_retained_shares_sum_to_one',
        'refused_share_sees_a_refusal_above_the_disk_cap',
        'saved_spectral_files_still_read')

FAIL = 0


def say(ok, name, measured, reference, tol):
	global FAIL
	if not ok:
		FAIL = 1
	print('%s %s measured=%s reference=%s tol=%s'
	      % ('PASS' if ok else 'FAIL', name, measured, reference, tol))


# ---------------------------------------------------------------------------
# Row 1: the census in the production script selects with the shared function
# and imposes no upper radius.  A radial mask (r >= 1) & (r <= Rib) describes
# less material than the integrals sample, since it drops every shell above the
# stellar-disk cap while the chords still cross them, so this row measures the
# source text of the script it is run against.
# ---------------------------------------------------------------------------
MODDIR = os.environ.get('TRANSIT_CENSUS_MODULE_DIR', '.')
SCRIPT = os.path.join(MODDIR, 'EXHALE_transit.py')
if os.path.exists(SCRIPT):
	src = open(SCRIPT).read()
	uses_shared = 'chord_shell_indices(r, r_grid[0])' in src
	has_cap     = '(r <= Rib)' in src
	say(uses_shared and not has_cap, 'script_census_uses_the_shared_selection',
	    'shared selection %s, upper radius mask %s'
	    % ('present' if uses_shared else 'absent',
	       'present' if has_cap else 'absent'),
	    'shared selection present, upper radius mask absent', '0')
else:
	say(False, 'script_census_uses_the_shared_selection',
	    'no %s' % SCRIPT, 'the script the census lives in', '0')


try:
	from exhale_transit_lib import (chord_shell_indices,
	                                refused_line_center_tau_share,
	                                resonance_depth,
	                                mp, f_la, A12_HI, lA)
except ImportError as _exc:
	# The census and the integrals have no shared selection to test.
	for _row in ROWS:
		say(False, _row, 'no shared selection (%s)' % _exc,
		    'exhale_transit_lib.chord_shell_indices', '0')
	sys.exit(1)


# ---------------------------------------------------------------------------
# Rows 2 to 4: the crossing condition itself, on the geometry of the review
# (shell r = 12, ray b = 2, ray grid capped at Rib = 10).
# ---------------------------------------------------------------------------
SHELLS = np.array([0.5, 1.5, 2.0, 5.0, 12.0])
B_RAY  = 2.0
RIB    = 10.0

sel = chord_shell_indices(SHELLS, B_RAY)
say(12.0 in SHELLS[sel], 'outer_shell_crosses_the_inner_ray',
    'r = 12 %s' % ('counted' if 12.0 in SHELLS[sel] else 'missing'),
    'r = 12 counted at b = 2', '0')

x_cross = float(np.sqrt(12.0**2.0 - B_RAY**2.0))
x_ref   = 11.832159566199232          # sqrt(140)
say(abs(x_cross - x_ref) <= 1.0e-12*x_ref, 'chord_crossing_of_the_outer_shell',
    '%.15f' % x_cross, '%.15f' % x_ref, '1.0e-12 relative')

say(1.5 not in SHELLS[sel], 'shell_inside_the_ray_is_not_crossed',
    'r = 1.5 %s' % ('counted' if 1.5 in SHELLS[sel] else 'not counted'),
    'r = 1.5 not counted at b = 2', '0')

# ---------------------------------------------------------------------------
# Row 5: the census of a ray grid 1 <= b <= Rib is the union of the crossings,
# which is the crossing condition of the innermost ray.  A radial mask
# r <= Rib, which the census carried until this suite existed, drops the shell
# above the cap although the integrals keep it.
# ---------------------------------------------------------------------------
n_census = int(chord_shell_indices(SHELLS, 1.0).size)
n_capped = int(np.count_nonzero((SHELLS >= 1.0) & (SHELLS <= RIB)))
say(n_census == 4, 'census_counts_the_shell_above_the_disk_cap',
    '%d shells' % n_census, '4 shells (r = 1.5, 2, 5, 12)', '0')
print('  DIAGNOSTIC a radial mask (r >= 1) & (r <= %g) counts %d of them and '
      'drops r = 12' % (RIB, n_capped))

# ---------------------------------------------------------------------------
# A synthetic atmosphere, used by the rows below.  Uniform temperature, no
# wind, absorbers ONLY in the shells above the disk cap: whatever depth comes
# out is depth those shells produced.
# ---------------------------------------------------------------------------
Rp      = 1.0e8                                    # planet radius [m]
r_shell = np.geomspace(1.0, 12.0, 60)              # shell radii [Rp]
T_shell = np.full(r_shell.size, 1.0e4)             # [K]
v_shell = np.zeros(r_shell.size)                   # [m/s]
n_above = np.where(r_shell > RIB, 1.0e8, 0.0)      # [m^-3], unsaturated

# Mirrored (night + day) chord arrays, as the integrals receive them.
data_r = np.append(-np.flip(r_shell), r_shell)*Rp
data_T = np.append(np.flip(T_shell), T_shell)
data_v = np.append(-np.flip(v_shell), v_shell)
data_n_above = np.append(np.flip(n_above), n_above)

Grid_Number = 40
r_grid    = np.array([RIB**(j/(Grid_Number - 1)) for j in range(Grid_Number)])
A_planet  = np.pi*Rp**2.0
A_atm     = np.pi*(RIB*Rp)**2.0
A_star    = np.pi*(RIB*Rp)**2.0                    # the cap is the disk edge

# ---------------------------------------------------------------------------
# Row 6: the optical-depth integral does absorb in those shells.  This is the
# control the census has to match: the review found the integral correct and
# the census narrower than the integral.
# ---------------------------------------------------------------------------
d_lc, _d_band, _l, _t = resonance_depth(lA*1e10, f_la, A12_HI, mp,
                                        data_n_above, 5.0e4,
                                        Grid_Number, r_grid, Rp,
                                        data_r, data_v, data_T,
                                        A_star, A_atm, A_planet)
say(d_lc > 0.0, 'integral_absorbs_in_shells_above_the_disk_cap',
    'line-center depth %.6g percent' % d_lc, '> 0', '0')

# ---------------------------------------------------------------------------
# Rows 7 to 10: the line-center optical-depth share the census reports.
# The absorbers now fill the column, so every ray carries optical depth.
# ---------------------------------------------------------------------------
n_all  = np.full(r_shell.size, 1.0e10)
data_n = np.append(np.flip(n_all), n_all)
LINE   = [('Ly-alpha', lA*1e10, [(lA*1e10, f_la, A12_HI, mp, data_n)])]

refused_all = np.ones(r_shell.size, dtype=bool)
data_ref    = np.append(np.flip(refused_all), refused_all)
tau, share  = refused_line_center_tau_share(LINE, r_grid, Rp, data_r, data_v,
                                            data_T, data_ref)['Ly-alpha']
lit  = tau > 0.0
worst = float(np.max(np.abs(share[lit] - 1.0)))
say(worst == 0.0 and np.count_nonzero(lit) > 1,
    'refused_share_is_one_when_every_row_is_refused',
    'largest departure from 1 over %d lit rays: %.3e'
    % (int(np.count_nonzero(lit)), worst), '0 over more than one ray', '0')

refused_none = np.zeros(r_shell.size, dtype=bool)
data_ref     = np.append(np.flip(refused_none), refused_none)
tau0, share0 = refused_line_center_tau_share(LINE, r_grid, Rp, data_r, data_v,
                                             data_T, data_ref)['Ly-alpha']
worst0 = float(np.max(np.abs(share0)))
say(worst0 == 0.0, 'refused_share_is_zero_when_no_row_is_refused',
    'largest share %.3e' % worst0, '0', '0')

# Half the column refused: the shares of the two disjoint subsets add up to
# one ray by ray, because the quadrature is a sum over rows.
refused_half = r_shell < 3.0
data_ref     = np.append(np.flip(refused_half), refused_half)
_t1, s_ref   = refused_line_center_tau_share(LINE, r_grid, Rp, data_r, data_v,
                                             data_T, data_ref)['Ly-alpha']
data_keep    = np.append(np.flip(~refused_half), ~refused_half)
_t2, s_keep  = refused_line_center_tau_share(LINE, r_grid, Rp, data_r, data_v,
                                             data_T, data_keep)['Ly-alpha']
lit   = _t1 > 0.0
worst = float(np.max(np.abs(s_ref[lit] + s_keep[lit] - 1.0)))
partial = bool(np.any(s_ref[lit] > 0.0) and np.any(s_keep[lit] > 0.0))
say(worst <= 1.0e-13 and partial, 'refused_and_retained_shares_sum_to_one',
    'largest departure %.3e, both subsets nonzero: %s' % (worst, partial),
    '<= 1.0e-13 with both subsets nonzero', '1.0e-13')

# Refusal ONLY above the disk cap: an inner ray still takes optical depth
# from it, which is the whole content of the geometry above.
refused_cap = r_shell > RIB
data_ref    = np.append(np.flip(refused_cap), refused_cap)
_t3, s_cap  = refused_line_center_tau_share(LINE, r_grid, Rp, data_r, data_v,
                                            data_T, data_ref)['Ly-alpha']
say(s_cap[0] > 0.0, 'refused_share_sees_a_refusal_above_the_disk_cap',
    'share at b = 1 Rp is %.3e' % s_cap[0], '> 0', '0')

# ---------------------------------------------------------------------------
# Row 11: the saved spectral files.  Four numerical columns, and the column
# line among the comment header, wherever in it: a curve saved with the
# validity metadata carries that block above its column line, and one saved
# before the metadata existed carries the column line alone.  Every reader of
# these files skips comments, so both read the same.
# ---------------------------------------------------------------------------
root  = os.environ.get('TRANSIT_CENSUS_ROOT', '.')
names = ['WASP-121b/tpm_He10830.txt', 'WASP-121b/tpm_MgII.txt',
         'WASP-52b/tpm_Lya.txt', 'HD209458b/tpm_He10830.txt',
         'HD189733b/tpm_Halpha.txt']
found, bad = 0, []
for nm in names:
	path = os.path.join(root, nm)
	if not os.path.exists(path):
		continue
	found += 1
	try:
		col = np.loadtxt(path)
		head = [ln for ln in open(path) if ln.startswith('#')]
	except Exception as exc:                       # noqa: BLE001
		bad.append('%s: %s' % (nm, exc))
		continue
	named = [ln for ln in head
	         if ln.lstrip('#').strip().startswith(
	             'lambda[A]  T_theo  T_instr  T_rot+instr')]
	if col.ndim != 2 or col.shape[1] != 4:
		bad.append('%s: shape %s' % (nm, col.shape))
	elif len(named) != 1:
		bad.append('%s: %d column lines in the header' % (nm, len(named)))
if found == 0:
	print('  DIAGNOSTIC no saved spectral file found under %s; nothing to read'
	      % root)
	say(True, 'saved_spectral_files_still_read', 'no file present',
	    'skipped where the planet directories are absent', '0')
else:
	say(not bad, 'saved_spectral_files_still_read',
	    '%d files, %d unreadable' % (found, len(bad)),
	    '4 columns and one column line in the header of every file', '0')
	for msg in bad:
		print('  ' + msg)

sys.exit(FAIL)
