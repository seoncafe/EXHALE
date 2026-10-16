#!/usr/bin/env python3
"""Write src/modules/radiation/low_energy_electron_tables.f90 from the
electronic tables of Furlanetto & Stoever (2010, MNRAS 404, 1869).

The tables are the authors' own: the directory x_int_tables/ of the archive
elec_interp.tar.gz that the paper's footnote 9 points to
(http://www.astro.ucla.edu/~sfurlane/xray.htm). That archive is no longer
served there (404); the copy used is the Wayback Machine capture of
2017-11-02,

    http://web.archive.org/web/20171102020104id_/http://www.astro.ucla.edu/~sfurlane/docs/elec_interp.tar.gz

whose SHA-1 equals the CDX digest of the capture (base32
BXITPPSGCYDO7IYBTMM7UXYCK5WOU4FY). The script checks both before it writes
anything: the SHA-1 of the archive in the table directory, and that the .dat
files it reads are byte for byte the ones inside that archive.

One file per ionized fraction x_i (14 of them, 1e-4 to 0.999), 258 rows of
electron energy E (10 eV to 9.9 keV, logarithmic), nine columns:

    E  f_ion  f_heat  f_exc  n_Lya  n_{ion,HI}  n_{ion,HeI}  n_{ion,HeII}  Shull

The rows written are those with E_MIN <= E <= E_MAX (default 10 to 32.2 eV,
which brackets the 30 eV at which the code's partition of fast photoelectrons
takes over, E_sec_ion in parameters.f90), every column but the last (the
Shull & van Steenberg heating fit the authors print for comparison). The
values are copied as the files print them, as double-precision literals;
the authors' interpolator elec_interp.c reads the same text into C floats,
so the two differ by the float rounding of the printed decimals and nothing
else.

Usage:
    python3 src/utils/furlanetto2010_tables_to_fortran.py [--tables DIR]
            [--emin 10.0] [--emax 32.2] [--out PATH]
"""
import argparse, hashlib, io, os, sys, tarfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_TABLES = os.path.normpath(os.path.join(ROOT, '..', 'references',
                                               'Furlanetto_2010_tables'))
DEFAULT_OUT = os.path.join(ROOT, 'src', 'modules', 'radiation',
                           'low_energy_electron_tables.f90')

ARCHIVE_SHA1 = '0dd137be461606efa3019b19fa5f02576cea70b8'
ARCHIVE_CDX = 'BXITPPSGCYDO7IYBTMM7UXYCK5WOU4FY'
ARCHIVE_URL = ('http://web.archive.org/web/20171102020104id_/'
               'http://www.astro.ucla.edu/~sfurlane/docs/elec_interp.tar.gz')

# The ionized fractions and the file names, in the order and with the values
# of initialize_interp_arrays() in elec_interp.c.
XI = [1.0e-4, 2.318e-4, 4.677e-4, 1.0e-3, 2.318e-3, 4.677e-3, 1.0e-2,
      2.318e-2, 4.677e-2, 1.0e-1, 0.5, 0.9, 0.99, 0.999]
FILES = ['log_xi_-4.0.dat', 'log_xi_-3.6.dat', 'log_xi_-3.3.dat',
         'log_xi_-3.0.dat', 'log_xi_-2.6.dat', 'log_xi_-2.3.dat',
         'log_xi_-2.0.dat', 'log_xi_-1.6.dat', 'log_xi_-1.3.dat',
         'log_xi_-1.0.dat', 'xi_0.500.dat', 'xi_0.900.dat', 'xi_0.990.dat',
         'xi_0.999.dat']
COLS = ['E', 'f_ion', 'f_heat', 'f_exc', 'n_Lya', 'n_HI', 'n_HeI', 'n_HeII']
NAMES = {'E': 'lee_E', 'f_ion': 'lee_f_ion', 'f_heat': 'lee_f_heat',
         'f_exc': 'lee_f_exc', 'n_Lya': 'lee_n_Lya', 'n_HI': 'lee_n_HI',
         'n_HeI': 'lee_n_HeI', 'n_HeII': 'lee_n_HeII'}


def sha1(data):
    return hashlib.sha1(data).hexdigest()


def read_table(text):
    lines = text.splitlines()
    head = lines[1].split()
    rows = [l.split() for l in lines[3:] if l.strip()]
    if len(rows) != 258 or any(len(r) != 9 for r in rows):
        raise SystemExit('unexpected table shape')
    return head, rows


def flit(s):
    """A printed decimal as a Fortran double literal, digits unchanged."""
    s = s.strip()
    if 'e' in s or 'E' in s:
        m, e = s.lower().split('e')
        if '.' not in m:
            m += '.0'
        return m + 'd' + e
    if '.' not in s:
        s += '.0'
    return s + 'd0'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tables', default=DEFAULT_TABLES)
    ap.add_argument('--emin', type=float, default=10.0)
    ap.add_argument('--emax', type=float, default=32.2)
    ap.add_argument('--out', default=DEFAULT_OUT)
    a = ap.parse_args()

    arch = os.path.join(a.tables, 'elec_interp.tar.gz')
    raw = open(arch, 'rb').read()
    if sha1(raw) != ARCHIVE_SHA1:
        raise SystemExit('archive SHA-1 %s is not %s' % (sha1(raw), ARCHIVE_SHA1))
    with tarfile.open(fileobj=io.BytesIO(raw)) as tf:
        packed = {os.path.basename(m.name): tf.extractfile(m).read()
                  for m in tf.getmembers()
                  if m.isfile() and m.name.endswith('.dat')
                  and not os.path.basename(m.name).startswith('._')}

    tabs, heads, digests = [], [], []
    for f in FILES:
        data = open(os.path.join(a.tables, 'x_int_tables', f), 'rb').read()
        if packed.get(f) != data:
            raise SystemExit('%s differs from the copy inside the archive' % f)
        head, rows = read_table(data.decode())
        tabs.append(rows)
        heads.append(head)
        digests.append(sha1(data))

    energies = [r[0] for r in tabs[0]]
    for t in tabs[1:]:
        if [r[0] for r in t] != energies:
            raise SystemExit('the files do not share one energy column')
    sel = [i for i, e in enumerate(energies) if a.emin <= float(e) <= a.emax]
    i0, i1 = sel[0], sel[-1]
    if sel != list(range(i0, i1 + 1)):
        raise SystemExit('selected rows are not contiguous')
    ne, nx = len(sel), len(XI)

    o = []
    w = o.append
    w('      module low_energy_electron_tables')
    w('      ! GENERATED by src/utils/furlanetto2010_tables_to_fortran.py; do not')
    w('      ! edit by hand, rerun the script.')
    w('      !')
    w('      ! Energy deposition of a fast electron in a primordial H/He gas,')
    w('      ! Furlanetto & Stoever (2010, MNRAS 404, 1869): the authors\' own')
    w('      ! electronic tables (their footnote 9), rows %d to %d (0-based) of'
      % (i0, i1))
    w('      ! the 258 energies, E = %s to %s eV, at all %d ionized fractions.'
      % (energies[i0], energies[i1], nx))
    w('      !')
    w('      ! Source archive: elec_interp.tar.gz, retrieved from')
    cut = ARCHIVE_URL.index('id_/') + 4
    w('      !   %s' % ARCHIVE_URL[:cut])
    w('      !   %s' % ARCHIVE_URL[cut:])
    w('      ! (the paper\'s http://www.astro.ucla.edu/~sfurlane/xray.htm points')
    w('      ! to it; the original address answers 404).')
    w('      !   SHA-1 %s' % ARCHIVE_SHA1)
    w('      !   = CDX digest %s (base32) of the capture' % ARCHIVE_CDX)
    w('      ! Local copy read: <workspace>/references/Furlanetto_2010_tables/')
    w('      ! x_int_tables/*.dat, each checked byte for byte against the')
    w('      ! archive. SHA-1 of each file, and its header line 2 (neutral')
    w('      ! fractions f_HI, f_HeI, the He II fraction, z, T[K]):')
    for f, d, h in zip(FILES, digests, heads):
        w('      !   %-16s %s  %s' % (f, d, ' '.join(h)))
    w('      !')
    w('      ! Columns kept (the ninth, the Shull & van Steenberg heating fit the')
    w('      ! authors print for comparison, is dropped): electron energy E [eV],')
    w('      ! the energy fractions f_ion, f_heat, f_exc, and the numbers per')
    w('      ! primary of H I Ly-alpha photons and of H I, He I and He II')
    w('      ! ionizations. Values are the printed decimals, digit for digit.')
    w('      ! Array index (energy row, ionized fraction).')
    w('')
    w('      implicit none')
    w('      public')
    w('')
    w('      integer, parameter :: n_lee_E = %d' % ne)
    w('      integer, parameter :: n_lee_x = %d' % nx)
    w('      ! First and last row of the authors\' 258, 0-based as in elec_interp.c.')
    w('      integer, parameter :: lee_row_first = %d' % i0)
    w('      integer, parameter :: lee_row_last  = %d' % i1)
    w('      ! The ionized fractions, as initialize_interp_arrays() sets them.')
    w('      real*8, parameter :: lee_x(n_lee_x) = (/                          &')
    xs = ['%s' % flit(repr(v)) for v in XI]
    for k in range(0, nx, 5):
        chunk = ', '.join(xs[k:k + 5])
        tail = ' /)' if k + 5 >= nx else ','
        w('           ' + chunk + tail + ('' if k + 5 >= nx else '  &'))

    def array(name, col):
        vals = []
        for jx in range(nx):
            for i in sel:
                vals.append(flit(tabs[jx][i][col]))
        w('      real*8, parameter :: %s(n_lee_E,n_lee_x) = reshape( (/   &' % name)
        per = 4
        for k in range(0, len(vals), per):
            chunk = ', '.join(vals[k:k + per])
            last = k + per >= len(vals)
            w('           ' + chunk + (' /), (/n_lee_E,n_lee_x/) )' if last else ', &'))

    w('      real*8, parameter :: lee_E(n_lee_E) = (/                          &')
    vals = [flit(energies[i]) for i in sel]
    for k in range(0, ne, 5):
        chunk = ', '.join(vals[k:k + 5])
        last = k + 5 >= ne
        w('           ' + chunk + (' /)' if last else ', &'))
    for c in COLS[1:]:
        array(NAMES[c], COLS.index(c))
    w('')
    w('      end module low_energy_electron_tables')

    with open(a.out, 'w') as fh:
        fh.write('\n'.join(o) + '\n')
    print('wrote %s: rows %d-%d (E = %s-%s eV), %d fractions'
          % (a.out, i0, i1, energies[i0], energies[i1], nx))


if __name__ == '__main__':
    main()
