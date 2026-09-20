"""WHICH STATE OF A RUN A TRANSIT SPECTRUM STANDS ON, AND HOW IT SAYS SO.

A run can hold two descriptions of one column and they are not the same
numbers.  The solved pair, `Hydro_ioniz.txt` and `Ion_species.txt`, is the
state the wind solver produced and the certification judged; where the
ionization stages are transported its composition is the partition the stage
transport itself produced.  The `_adv` pair is a post-process of a marching
state: a steady advective correction applied cell by cell, first-order upwind
at the bulk velocity, with no eddy term and no element drift, refused where an
assumption of the steady equations fails.  A stationary solve writes no `_adv`
products at all, so a solved run has no `adv` state to read, and one is not to
be manufactured for it.

These rows pin that:

* each selection returns the pair it names, and the two files of a pair are
  read together, since the temperature of one state with the composition of
  the other solves neither set of equations;
* a missing `_adv` pair under the default selection is REFUSED by name, and
  the solved pair is never read in its place;
* the solution selection works on a run that has no `_adv` products;
* the saved curve names the files actually read, by path and by md5, since a
  run directory is written again by the next run and a path alone identifies
  nothing, together with the statements the state makes about itself;
* the difference the tool reports between the two states is the difference a
  direct comparison of the same columns gives, and is None, not zero, when the
  other state is not there.

Every row prints one
    PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
line; lines beginning with two spaces or with DIAGNOSTIC are context, not
verdicts.  Exit status is nonzero if any row fails.
"""

import os
import subprocess
import sys
import tempfile

import numpy as np

ROWS = ('adv_selection_names_the_post_process_pair',
        'solution_selection_names_the_solved_pair',
        'unknown_selection_refused',
        'missing_adv_pair_refused_by_name',
        'missing_adv_pair_not_substituted_by_the_solved_pair',
        'solution_selection_works_without_adv_products',
        'oi_levels_follow_the_selection',
        'file_identity_is_the_digest_of_the_file',
        'provenance_header_names_the_files_read',
        'provenance_header_carries_the_state_statements',
        'state_pair_difference_agrees_with_a_direct_comparison',
        'state_pair_difference_absent_is_unknown_not_zero')

FAIL = 0


def say(ok, name, measured, reference, tol):
	global FAIL
	if not ok:
		FAIL = 1
	print('%s %s measured=%s reference=%s tol=%s'
	      % ('PASS' if ok else 'FAIL', name, measured, reference, tol))


try:
	import exhale_transit_lib as lib
except ImportError as _exc:
	for _row in ROWS:
		say(False, _row, 'no module (%s)' % _exc, 'exhale_transit_lib', '0')
	sys.exit(1)


def entry(name):
	"""The module's function of this name, or None if the text has none."""
	return getattr(lib, name, None)


# ---------------------------------------------------------------------------
# Two runs on disk.  `marched` carries both states, `solved` carries only the
# state the solver wrote, which is what a stationary route leaves behind.
# Six rows of which two at each end are ghosts, the layout the writer
# produces, so the readers return the two physical rows.
# ---------------------------------------------------------------------------
TMP = tempfile.mkdtemp(prefix='transit_state_')

HYD_COLS = ('columns r[Rp] rho[mH/cm3] v[cm/s] p[cgs] T[K] '
            'heat[erg/cm3/s] cool[erg/cm3/s]')
ION_COLS = 'columns r HI HII HeI HeII HeIII HeITR'
ROWS_HDR = 'rows 6: 2 ghost cells at each end; physical cells are rows 3 to 4'
PROV = 'provenance: git=0123456789ab tree=clean run=2026-09-19T00:00:00'
PROV2 = 'provenance: recon=WENO3 base_bc=characteristic carrier=none N=2'
BMODEL = 'boundary_model characteristic_face_ps_reservoir_C_minus_contact_upwind_v2'
BRESV = ('boundary_reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp] '
         '1.0E+000 1.0E+000 8.0E-001 1.0E+000')
COUPL = ('coupling: sec_ion=T recon=WENO3 iontrans=T certified=T '
         'cert_reason=certified_in_wind mode=init')

# T of the six rows, solved state and post-process.
T_SOL = [900.0, 950.0, 6000.0, 8000.0, 9000.0, 9500.0]
T_ADV = [900.0, 950.0, 6300.0, 8000.0, 9000.0, 9500.0]
# n(H I) and n(He 2^3S) of the six rows in each state.
HI_SOL = [1.0e10, 9.0e9, 1.0e9, 5.0e8, 1.0e8, 1.0e7]
HI_ADV = [1.0e10, 9.0e9, 1.0e9, 6.0e8, 1.0e8, 1.0e7]
TR_SOL = [1.0e2, 1.0e2, 2.0e3, 1.0e3, 5.0e2, 1.0e2]
TR_ADV = [1.0e2, 1.0e2, 2.0e3, 1.0e3, 5.0e2, 1.0e2]
RG = [1.0, 1.05, 1.1, 1.2, 1.3, 1.4]


def write_hydro(path, header, T):
	with open(path, 'w') as fh:
		for line in header:
			fh.write('# %s\n' % line)
		for k in range(6):
			fh.write(' '.join('%16.8e' % x
			                  for x in (RG[k], 1.0e9, 1.0e4, 1.0e-4, T[k],
			                            1.0e-10, 1.0e-11)) + '\n')


def write_ioniz(path, header, HI, TR):
	with open(path, 'w') as fh:
		for line in header:
			fh.write('# %s\n' % line)
		for k in range(6):
			fh.write(' '.join('%16.8e' % x
			                  for x in (RG[k], HI[k], 1.0e9, 1.0e8, 1.0e7,
			                            1.0e6, TR[k])) + '\n')


def write_oi(path):
	with open(path, 'w') as fh:
		fh.write('# columns r OI\n')
		for k in range(6):
			fh.write('%16.8e %16.8e\n' % (RG[k], 1.0e3))


SOLVED_HDR_H = [HYD_COLS, ROWS_HDR, COUPL, PROV, PROV2, BMODEL, BRESV]
SOLVED_HDR_I = [ION_COLS, ROWS_HDR, PROV]
ADV_HDR_H = [HYD_COLS, ROWS_HDR, 'adv_schema 0', PROV]
ADV_HDR_I = [ION_COLS, ROWS_HDR, PROV]

marched = os.path.join(TMP, 'marched')
solved = os.path.join(TMP, 'solved')
for d in (marched, solved):
	os.makedirs(os.path.join(d, 'output'))
write_hydro(os.path.join(marched, 'output', 'Hydro_ioniz.txt'),
            SOLVED_HDR_H, T_SOL)
write_ioniz(os.path.join(marched, 'output', 'Ion_species.txt'),
            SOLVED_HDR_I, HI_SOL, TR_SOL)
write_hydro(os.path.join(marched, 'output', 'Hydro_ioniz_adv.txt'),
            ADV_HDR_H, T_ADV)
write_ioniz(os.path.join(marched, 'output', 'Ion_species_adv.txt'),
            ADV_HDR_I, HI_ADV, TR_ADV)
write_oi(os.path.join(marched, 'output', 'OI_levels.txt'))
write_oi(os.path.join(marched, 'output', 'OI_levels_adv.txt'))
write_hydro(os.path.join(solved, 'output', 'Hydro_ioniz.txt'),
            SOLVED_HDR_H, T_SOL)
write_ioniz(os.path.join(solved, 'output', 'Ion_species.txt'),
            SOLVED_HDR_I, HI_SOL, TR_SOL)


def select(path, selection, **kw):
	"""transit_state_files, returning the exception instead of raising."""
	try:
		return lib.transit_state_files(path, selection, **kw)
	except Exception as exc:                                    # noqa: BLE001
		return exc


# --- the pair each selection names ----------------------------------------
got = select(marched, 'adv')
want = (os.path.join(marched, 'output', 'Hydro_ioniz_adv.txt'),
        os.path.join(marched, 'output', 'Ion_species_adv.txt'), 'adv')
say(got == want, ROWS[0],
    'None' if isinstance(got, Exception) else '/'.join(
        os.path.basename(q) for q in got[:2]),
    'Hydro_ioniz_adv.txt/Ion_species_adv.txt', 'exact')

got = select(marched, 'solution')
want = (os.path.join(marched, 'output', 'Hydro_ioniz.txt'),
        os.path.join(marched, 'output', 'Ion_species.txt'), 'solution')
say(got == want, ROWS[1],
    'None' if isinstance(got, Exception) else '/'.join(
        os.path.basename(q) for q in got[:2]),
    'Hydro_ioniz.txt/Ion_species.txt', 'exact')

got = select(marched, 'the corrected one')
say(isinstance(got, ValueError), ROWS[2],
    type(got).__name__ if isinstance(got, Exception) else 'returned a pair',
    'ValueError', 'exact')

# --- a missing post-process pair is refused, and not replaced -------------
got = select(solved, 'adv')
msg = str(got) if isinstance(got, Exception) else ''
named = (isinstance(got, Exception)
         and 'Hydro_ioniz_adv.txt' in msg and 'Ion_species_adv.txt' in msg)
say(named, ROWS[3],
    ('returned a pair' if not isinstance(got, Exception)
     else 'raised %s, names both missing files: %s'
          % (type(got).__name__, named)),
    'a refusal naming Hydro_ioniz_adv.txt and Ion_species_adv.txt', 'exact')

substituted = (not isinstance(got, Exception)
               and os.path.basename(got[0]) == 'Hydro_ioniz.txt')
say(not substituted, ROWS[4],
    'solved pair returned under the adv selection' if substituted
    else 'refused', 'never the other state', 'exact')

# --- the solution selection on a run that has no post-process -------------
got = select(solved, 'solution')
ok = (not isinstance(got, Exception)
      and os.path.exists(got[0]) and os.path.exists(got[1])
      and got[2] == 'solution')
say(ok, ROWS[5],
    type(got).__name__ if isinstance(got, Exception) else 'the solved pair',
    'the solved pair, read on a run with no _adv products', 'exact')

# --- the O I level populations come from the selected state ---------------
oi = entry('transit_state_oi_levels')
if oi is None:
	say(False, ROWS[6], 'no transit_state_oi_levels',
	    'the level file follows the selection', 'exact')
else:
	a = os.path.basename(oi(marched, 'adv'))
	b = os.path.basename(oi(marched, 'solution'))
	say(a == 'OI_levels_adv.txt' and b == 'OI_levels.txt', ROWS[6],
	    '%s / %s' % (a, b), 'OI_levels_adv.txt / OI_levels.txt', 'exact')

# --- the identity of the files a curve was built from ---------------------
ident = entry('file_identity')
hpath = os.path.join(solved, 'output', 'Hydro_ioniz.txt')
if ident is None:
	say(False, ROWS[7], 'no file_identity', 'path, md5 and size', 'exact')
	ident_line = ''
else:
	ident_line = ident(hpath)
	ref = subprocess.run(['md5sum', hpath], capture_output=True,
	                     text=True).stdout.split()[0]
	say(('md5=%s' % ref) in ident_line and hpath in ident_line, ROWS[7],
	    ident_line, 'path and md5=%s' % ref, 'exact')

# --- what the saved curve says about the state it stands on ---------------
prov = entry('state_provenance_statements')
ipath = os.path.join(solved, 'output', 'Ion_species.txt')
if ident is None or prov is None:
	for _row in (ROWS[8], ROWS[9]):
		say(False, _row, 'no state identity or provenance reader',
		    'the header names the files read and the state statements',
		    'exact')
else:
	record = {'selection': 'solution', 'hydro': hpath, 'ioniz': ipath,
	          'identity': [ident(hpath), ident(ipath)],
	          'provenance': prov(hpath), 'difference': None}
	block = lib.transit_metadata_block(None, 'git=none', [], None,
	                                   state=record)
	idl = [q for q in block if q.startswith('state_file_identity')]
	ok = (len(idl) == 2 and hpath in idl[0] and ipath in idl[1]
	      and all('md5=' in q for q in idl))
	say(ok, ROWS[8], '%d identity lines' % len(idl),
	    'one naming each file read, with its digest', 'exact')
	want = ('state_provenance: git=0123456789ab tree=clean '
	        'run=2026-09-19T00:00:00')
	has_bm = any(q.startswith('state_boundary_model') for q in block)
	has_cp = any(q.startswith('state_coupling') for q in block)
	say(want in block and has_bm and has_cp, ROWS[9],
	    'provenance=%s boundary_model=%s coupling=%s'
	    % (want in block, has_bm, has_cp),
	    'the state provenance, its boundary model and its coupling', 'exact')

# --- the reported distance between the two states -------------------------
spd = entry('state_pair_difference')
if spd is None:
	def spd(*a, **k):
		return None
d = spd(
    os.path.join(marched, 'output', 'Hydro_ioniz.txt'),
    os.path.join(marched, 'output', 'Ion_species.txt'),
    os.path.join(marched, 'output', 'Hydro_ioniz_adv.txt'),
    os.path.join(marched, 'output', 'Ion_species_adv.txt'))


def rel(x, y):
	x = np.asarray(x, dtype=float)
	y = np.asarray(y, dtype=float)
	return float(np.max(np.abs(x - y)
	                    /np.maximum(np.maximum(np.abs(x), np.abs(y)), 1e-99)))


# The physical rows are the middle pair, rows 3 and 4 of the six.
eT = rel(T_SOL[2:4], T_ADV[2:4])
eH = rel(HI_SOL[2:4], HI_ADV[2:4])
eR = rel(TR_SOL[2:4], TR_ADV[2:4])
if d is None:
	say(False, ROWS[10], 'None', 'T %.6e HI %.6e HeI_2_3S %.6e' % (eT, eH, eR),
	    '0')
else:
	err = max(abs(d['T'] - eT), abs(d['HI'] - eH), abs(d['HeI_2_3S'] - eR))
	say(err == 0.0 and d['rows'] == 2, ROWS[10],
	    'T %.6e HI %.6e HeI_2_3S %.6e rows %d'
	    % (d['T'], d['HI'], d['HeI_2_3S'], d['rows']),
	    'T %.6e HI %.6e HeI_2_3S %.6e rows 2' % (eT, eH, eR), '0')

d0 = spd(
    os.path.join(solved, 'output', 'Hydro_ioniz.txt'),
    os.path.join(solved, 'output', 'Ion_species.txt'),
    os.path.join(solved, 'output', 'Hydro_ioniz_adv.txt'),
    os.path.join(solved, 'output', 'Ion_species_adv.txt'))
say(d0 is None, ROWS[11], repr(d0), 'None, the distance is UNKNOWN', 'exact')

sys.exit(FAIL)
