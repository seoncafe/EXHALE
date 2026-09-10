"""What a transit spectrum knows about the validity of the rows it stands on.

The `_adv` profile states, row by row, whether the steady advective correction
was adopted, and the schema of that statement is versioned in the file's own
header.  These rows pin the three cases the transit post-processor has to
distinguish, and pin the comment metadata every saved spectrum carries:

schema 2  two fields, `adv_T_status` and `adv_comp_status`, plus the
          `adv_mass_row` column and the `adv_conditional_tol` line: a
          corrected row is a CONDITIONAL correction accurate to that fraction
          of itself in the mass flux, and both the fraction and each row's own
          measure have to reach the spectrum
schema 1  the single `adv_status` field
schema 0  neither: the validity of every row is UNKNOWN, which is not the same
          statement as no row refused

The metadata rows check that the block is written, that it parses, and that it
leaves the numerical columns of the file untouched, since a spectrum whose
numbers move when a comment is added is not the same spectrum.

Every row prints one
    PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
line; lines beginning with two spaces or with DIAGNOSTIC are context, not
verdicts.  Exit status is nonzero if any row fails.
"""

import os
import sys
import tempfile

import numpy as np

ROWS = ('adv_schema_2_read_by_name',
        'adv_schema_2_header_travels',
        'adv_mass_row_and_tolerance_read_by_name',
        'conditional_tolerance_travels_into_the_metadata',
        'adv_schema_1_read_as_one_field',
        'adv_legacy_read_as_unknown',
        'legacy_metadata_says_unknown_not_zero_refused',
        'metadata_block_parses_in_a_saved_file',
        'metadata_leaves_the_columns_bitwise',
        'script_reads_the_validity_by_header')

FAIL = 0


def say(ok, name, measured, reference, tol):
	global FAIL
	if not ok:
		FAIL = 1
	print('%s %s measured=%s reference=%s tol=%s'
	      % ('PASS' if ok else 'FAIL', name, measured, reference, tol))


try:
	from exhale_transit_lib import (read_adv_validity, transit_metadata_block,
	                                transit_tool_identity,
	                                transit_environment_overrides,
	                                TRANSIT_SCHEMA)
except ImportError as _exc:
	for _row in ROWS:
		say(False, _row, 'no validity reader (%s)' % _exc,
		    'exhale_transit_lib.read_adv_validity', '0')
	sys.exit(1)


# ---------------------------------------------------------------------------
# Synthetic profiles.  Six rows, of which two at each end are ghosts, exactly
# the layout the writer produces, so the readers return the two physical rows
# and the status columns line up with them.
# ---------------------------------------------------------------------------
PHYS  = 'r[Rp] rho[mH/cm3] v[cm/s] p[cgs] T[K] heat[erg/cm3/s] cool[erg/cm3/s]'
BODY  = [(1.0, 1.0e10, 0.0, 1.0e-3, 1000.0, 1.0e-10, 1.0e-11),
         (1.1, 5.0e9, 1.0e4, 5.0e-4, 3000.0, 2.0e-10, 2.0e-11),
         (1.2, 1.0e9, 5.0e4, 1.0e-4, 6000.0, 3.0e-10, 3.0e-11),
         (1.3, 5.0e8, 1.0e5, 5.0e-5, 8000.0, 4.0e-10, 4.0e-11),
         (1.4, 1.0e8, 2.0e5, 1.0e-5, 9000.0, 5.0e-10, 5.0e-11),
         (1.5, 1.0e7, 3.0e5, 1.0e-6, 9500.0, 6.0e-10, 6.0e-11)]
ROWS_HDR = ('rows %d: 2 ghost cells at each end; physical cells are rows 3 to %d'
            % (len(BODY), len(BODY) - 2))
PROV = ('provenance: git=0123456789ab tree=clean '
        'run=2026-09-09T00:00:00')

# The status of the six rows; the two physical rows are the middle pair.
T_ALL    = [0, 0, 1, 0, 0, 0]
COMP_ALL = [0, 0, 2, 3, 0, 0]
ONE_ALL  = [0, 0, 2, 3, 0, 0]
# The measure each row was decided by, in the same row order: the retained
# row is the one above the fraction below.
MROW_ALL = [1.0e-9, 1.0e-9, 4.0e-2, 3.0e-4, 1.0e-9, 1.0e-9]
COND_TOL = 1.0e-2
COND_TEXT = ('%.1E the correction of a corrected row is accurate to this '
             'fraction of itself in the mass flux' % COND_TOL)

CERT_TEXT = ('F the certification found failing entries in the mass balance '
             'of the input state')
OPER_TEXT = 'face_flux_mass'
REST_TEXT = ('molecular and oxygen columns are not advection corrected; the '
             'product is a one-way correction on a fixed density and velocity '
             'field')

TMP = tempfile.mkdtemp(prefix='transit_census_')


def write_profile(name, header, status_cols, real_cols=()):
	path = os.path.join(TMP, name)
	with open(path, 'w') as fh:
		for line in header:
			fh.write('# %s\n' % line)
		for k, row in enumerate(BODY):
			fh.write(' '.join('%16.8e' % x for x in row))
			for col in status_cols:
				fh.write(' %3d' % col[k])
			for col in real_cols:
				fh.write(' %16.8e' % col[k])
			fh.write('\n')
	return path


p2 = write_profile(
	'schema2_adv.txt',
	['EXHALE schema 2',
	 'columns %s adv_T_status adv_comp_status adv_mass_row' % PHYS,
	 'adv_schema 2',
	 'adv_T_status 0 corrected. 1 retained. 2 failed. 3 unsupported. '
	 '4 not_evaluated.',
	 'adv_comp_status 0 corrected. 1 retained. 2 failed. 3 unsupported. '
	 '4 not_evaluated.',
	 'adv_status_counts T 5 1 0 0 0 comp 4 0 1 1 0',
	 'adv_input_certified %s' % CERT_TEXT,
	 'adv_stationarity_operator %s' % OPER_TEXT,
	 'adv_conditional_tol %s' % COND_TEXT,
	 'adv_model_restrictions %s' % REST_TEXT,
	 ROWS_HDR, PROV],
	[T_ALL, COMP_ALL], [MROW_ALL])

p1 = write_profile(
	'schema1_adv.txt',
	['EXHALE schema 2',
	 'columns %s adv_status' % PHYS,
	 'adv_status why the row is not the full steady advective correction; '
	 '0 corrected. 1 T kept. 2 T and composition kept. 3 composition kept. '
	 '4 T kept.',
	 'adv_status_counts 4 0 1 1 0',
	 ROWS_HDR, PROV],
	[ONE_ALL])

p0 = write_profile(
	'legacy_adv.txt',
	['EXHALE schema 2', 'columns %s' % PHYS, ROWS_HDR, PROV],
	[])

# ---------------------------------------------------------------------------
# Row 1: schema 2 is read by column name, both fields, and a row is refused on
# the temperature field while its composition field is reported beside it.
# ---------------------------------------------------------------------------
a2 = read_adv_validity(p2)
want_T    = np.array(T_ALL[2:-2], dtype=float)
want_comp = np.array(COMP_ALL[2:-2], dtype=float)
ok = (a2['schema'] == 2
      and a2['T_status'] is not None and a2['comp_status'] is not None
      and np.array_equal(np.asarray(a2['T_status']).ravel(), want_T)
      and np.array_equal(np.asarray(a2['comp_status']).ravel(), want_comp)
      and np.array_equal(np.asarray(a2['refused']).ravel(), want_T != 0))
say(ok, 'adv_schema_2_read_by_name',
    'schema %s T %s comp %s refused %s'
    % (a2['schema'],
       None if a2['T_status'] is None
       else list(np.asarray(a2['T_status']).ravel().astype(int)),
       None if a2['comp_status'] is None
       else list(np.asarray(a2['comp_status']).ravel().astype(int)),
       None if a2['refused'] is None
       else list(np.asarray(a2['refused']).ravel())),
    'schema 2 T %s comp %s refused %s'
    % (list(want_T.astype(int)), list(want_comp.astype(int)),
       list(want_T != 0)), '0')

# ---------------------------------------------------------------------------
# Row 2: the certification, the stationarity operator and the model
# restrictions of the input travel out of its header verbatim.
# ---------------------------------------------------------------------------
ok = (a2['input_certified'] == CERT_TEXT
      and a2['stationarity_operator'] == OPER_TEXT
      and a2['model_restrictions'] == REST_TEXT
      and PROV in a2['provenance'])
say(ok, 'adv_schema_2_header_travels',
    'certified %r operator %r restrictions %r provenance %s'
    % (a2['input_certified'][:24], a2['stationarity_operator'],
       a2['model_restrictions'][:24], 'kept' if PROV in a2['provenance']
       else 'lost'),
    'the four header statements of the input, verbatim', '0')

# ---------------------------------------------------------------------------
# Row 2b: the measure of every row and the fraction a corrected row is
# accurate to are read by name, both of them, so a caller has the CONDITION of
# each row and not only its verdict.
# ---------------------------------------------------------------------------
want_m = np.array(MROW_ALL[2:-2], dtype=float)
# .get, so a reader that does not carry these keys at all fails the row
# with its own line instead of raising.
got_m = (None if a2.get('mass_row') is None
         else np.asarray(a2['mass_row']).ravel())
ok = (got_m is not None and np.array_equal(got_m, want_m)
      and a2.get('conditional_tol') == COND_TOL
      and a2.get('conditional_tol_text') == COND_TEXT)
say(ok, 'adv_mass_row_and_tolerance_read_by_name',
    'measure %s, tol %s'
    % (None if got_m is None else list(got_m), a2.get('conditional_tol')),
    'measure %s, tol %s' % (list(want_m), COND_TOL), '0')

# ---------------------------------------------------------------------------
# Row 2c: that fraction travels into the metadata of a saved curve, verbatim.
# A spectrum built on corrected rows is accurate to no better than the
# fraction those rows were corrected to, and a reader of the curve alone
# cannot know it otherwise.
# ---------------------------------------------------------------------------
block_c = transit_metadata_block(a2, 'git=unavailable', [], None)
tol_lines = [ln for ln in block_c if ln.startswith('adv_conditional_tol')]
mrow_lines = [ln for ln in block_c if ln.startswith('adv_mass_row')]
ok = (tol_lines == ['adv_conditional_tol %s' % COND_TEXT]
      and len(mrow_lines) == 1
      and '%.3e' % float(want_m.max()) in mrow_lines[0])
say(ok, 'conditional_tolerance_travels_into_the_metadata',
    'tol line %r, measure line %r'
    % (tol_lines[0][:40] if tol_lines else None,
       mrow_lines[0][:40] if mrow_lines else None),
    'the tolerance line verbatim and the largest measure', '0')

# ---------------------------------------------------------------------------
# Row 3: schema 1 is one field.  There is no composition field to report, and
# the refusal is the nonzero value of that one field.
# ---------------------------------------------------------------------------
a1 = read_adv_validity(p1)
want_1 = np.array(ONE_ALL[2:-2], dtype=float)
ok = (a1['schema'] == 1 and a1['comp_status'] is None
      and np.array_equal(np.asarray(a1['T_status']).ravel(), want_1)
      and np.array_equal(np.asarray(a1['refused']).ravel(), want_1 != 0))
say(ok, 'adv_schema_1_read_as_one_field',
    'schema %s comp %s status %s'
    % (a1['schema'], a1['comp_status'],
       None if a1['T_status'] is None
       else list(np.asarray(a1['T_status']).ravel().astype(int))),
    'schema 1, no composition field, status %s'
    % list(want_1.astype(int)), '0')

# ---------------------------------------------------------------------------
# Row 4: a profile stating no row validity leaves it UNKNOWN.  `refused` is
# None, and nothing in the reader turns silence into 'corrected'.
# ---------------------------------------------------------------------------
a0 = read_adv_validity(p0)
ok = (a0['schema'] == 0 and a0['refused'] is None
      and a0['T_status'] is None and a0['input_certified'] == 'unknown')
say(ok, 'adv_legacy_read_as_unknown',
    'schema %s refused %s certified %r'
    % (a0['schema'], a0['refused'], a0['input_certified']),
    'schema 0, refused None, certified unknown', '0')

# ---------------------------------------------------------------------------
# Row 5: the metadata of a spectrum built on such a profile says UNKNOWN, and
# does not report zero refused rows.
# ---------------------------------------------------------------------------
block0 = transit_metadata_block(a0, 'git=unavailable', [], None)
text0  = '\n'.join(block0)
ok = ('UNKNOWN' in text0 and 'refused_rows' not in text0
      and 'transit_schema %d' % TRANSIT_SCHEMA in text0)
say(ok, 'legacy_metadata_says_unknown_not_zero_refused',
    'UNKNOWN %s, refused_rows %s'
    % ('present' if 'UNKNOWN' in text0 else 'absent',
       'present' if 'refused_rows' in text0 else 'absent'),
    'UNKNOWN present, refused_rows absent', '0')

# ---------------------------------------------------------------------------
# Rows 6 and 7: the block in a saved file.  The spectrum here is a stand-in
# for a curve: two columns of what the writer would save.
# ---------------------------------------------------------------------------
lam   = np.linspace(10828.2, 10831.2, 41)
curve = np.c_[lam, 1.0 - 0.01*np.exp(-((lam - 10830.3)/0.2)**2.0),
              1.0 - 0.008*np.exp(-((lam - 10830.3)/0.3)**2.0),
              1.0 - 0.006*np.exp(-((lam - 10830.3)/0.4)**2.0)]
census = {'sampled': 4, 'refused': 1, 'above_cap': 1,
          'reasons': 'retained: 1 rows', 'comp': 'corrected: 3 rows'}
line_census = ([1.0, 1.25, 1.5], [0.5, 0.25, 0.0], 0.5, 1.0)
block = transit_metadata_block(a2, transit_tool_identity([__file__]),
                               transit_environment_overrides(), census,
                               line_census)
COLS  = 'lambda[A]  T_theo  T_instr  T_rot+instr  (He I 10830)'
f_meta = os.path.join(TMP, 'tpm_He10830.txt')
f_bare = os.path.join(TMP, 'tpm_He10830_nometa.txt')
np.savetxt(f_meta, curve, header='\n'.join(block + [COLS]))
np.savetxt(f_bare, curve, header=COLS)

head = [ln.lstrip('#').strip() for ln in open(f_meta) if ln.startswith('#')]
keys = [ln.split()[0] for ln in head if ln]
need = ('transit_schema', 'adv_schema', 'adv_input_certified',
        'adv_stationarity_operator', 'adv_conditional_tol', 'adv_mass_row',
        'adv_model_restrictions',
        'transit_tool', 'census', 'census_share_b[Rp]', 'census_share',
        'census_share_max', 'census_note')
missing = [k for k in need if k not in keys]
has_prov = any(k.startswith('input_provenance') for k in keys)
schema_val = [ln.split()[1] for ln in head if ln.split()[0] == 'transit_schema']
contribution = any('CONTRIBUTION' in ln for ln in head)
data = np.loadtxt(f_meta)
ok = (not missing and has_prov and contribution
      and schema_val == [str(TRANSIT_SCHEMA)]
      and head[-1].startswith('lambda[A]')
      and data.shape == curve.shape)
say(ok, 'metadata_block_parses_in_a_saved_file',
    'missing %s, provenance %s, CONTRIBUTION %s, columns line last %s, '
    'shape %s' % (missing or 'none', 'present' if has_prov else 'absent',
                  'present' if contribution else 'absent',
                  head[-1].startswith('lambda[A]'), data.shape),
    'every key present, columns line last, %s columns of data'
    % curve.shape[1], '0')

d_meta = [ln for ln in open(f_meta) if not ln.startswith('#')]
d_bare = [ln for ln in open(f_bare) if not ln.startswith('#')]
say(d_meta == d_bare, 'metadata_leaves_the_columns_bitwise',
    '%d data lines, %s' % (len(d_meta),
                           'identical' if d_meta == d_bare else 'differ'),
    'the same bytes as a save without the metadata', '0')

# ---------------------------------------------------------------------------
# Row 8: the production script decides validity from the profile's header, and
# a legacy profile reaches the UNKNOWN branch there rather than a count of
# zero refused rows.
# ---------------------------------------------------------------------------
MODDIR = os.environ.get('TRANSIT_CENSUS_MODULE_DIR', '.')
SCRIPT = os.path.join(MODDIR, 'EXHALE_transit.py')
if os.path.exists(SCRIPT):
	src = open(SCRIPT).read()
	by_header = 'read_adv_validity(Hydro_file)' in src
	unknown   = ('if _adv_refused is None:' in src and 'UNKNOWN' in src)
	writes    = 'transit_metadata_block(' in src
	say(by_header and unknown and writes,
	    'script_reads_the_validity_by_header',
	    'header reader %s, UNKNOWN branch %s, metadata written %s'
	    % ('present' if by_header else 'absent',
	       'present' if unknown else 'absent',
	       'present' if writes else 'absent'),
	    'all three present', '0')
else:
	say(False, 'script_reads_the_validity_by_header',
	    'no %s' % SCRIPT, 'the script the census lives in', '0')

sys.exit(FAIL)
