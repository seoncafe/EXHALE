#!/usr/bin/env python3
"""Write an XUV-scaled copy of a spectrum: the whole flux column times a
factor, on the same wavelength grid.

    scale_sed.py <factor> <tag> ["<note>"]     e.g.  scale_sed.py 0.06 0p06

The scaling is ACHROMATIC and redefines no band: every wavelength keeps its
value and every band integral is the reference integral times the factor, so
a model series built on these files is a one-parameter family in the
irradiation strength alone.  The flux the run sees is the file's: with
`Spectrum type: Load from file..` the code overwrites LX and LEUV with this
file's own band integrals (sed_read.f90, "Calculate LEUV and LX
luminosities"), so the luminosity lines of input.inp do not enter.

The rule is the one the existing grid was written by, and this script
reproduces every one of those files byte for byte (`--check`).
"""
import sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
REF  = 'lhs1140_sed_gj1132_at_b.txt'
NHEAD = 9                      # the reference's own header lines


def scaled_text(factor, note):
    lines = open(os.path.join(HERE, REF)).read().splitlines()
    head, body = lines[:NHEAD], lines[NHEAD:]
    out = list(head) + [l for l in note.splitlines()]
    for l in body:
        w, f = l.split()
        out.append('%s %.6e' % (w, float(f)*factor))
    return '\n'.join(out) + '\n'


def note_of(factor, extra):
    return ('# flux column scaled by %g relative to %s (XUV grid,\n# %s)'
            % (factor, REF, extra))


if sys.argv[1] == '--check':
    # Every scaled file on the tree, with the note it actually carries, so
    # that only the numbers are compared.
    bad = 0
    for tag in ('0p33', '0p30', '0p28', '0p26', '0p25', '0p20', '0p18',
                '0p16', '0p15', '0p10', '0p07', '0p065', '0p06', '0p05',
                '0p03', '0p02', '0p015', '0p01'):
        path = os.path.join(HERE, 'lhs1140_sed_gj1132_at_b_xuv%s.txt' % tag)
        if not os.path.exists(path):
            continue
        have = open(path).read().splitlines()
        note = '\n'.join(have[NHEAD:NHEAD+2])
        factor = float(tag.replace('p', '.'))
        want = scaled_text(factor, note).splitlines()
        n = sum(1 for a, b in zip(have, want) if a != b)
        bad += (n > 0) or (len(have) != len(want))
        print('%-6s rows %5d  differing lines %d' % (tag, len(have), n
                                                     + abs(len(have)-len(want))))
    sys.exit(1 if bad else 0)

factor = float(sys.argv[1])
tag    = sys.argv[2]
extra  = sys.argv[3] if len(sys.argv) > 3 else 'no reason recorded'
out = os.path.join(HERE, 'lhs1140_sed_gj1132_at_b_xuv%s.txt' % tag)
if os.path.exists(out):
    sys.exit('%s exists; a spectrum of this tree is never rewritten' % out)
open(out, 'w').write(scaled_text(factor, note_of(factor, extra)))
print('wrote %s' % out)
