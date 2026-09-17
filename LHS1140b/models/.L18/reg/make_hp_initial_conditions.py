#!/usr/bin/env python3
"""Build the IC pairs of the three transported-proton cases from a
mol_base_handoff output.

The three cases hp_front, hp_zero_seed and hp_trace_seed restart the
mol_base_handoff state with `Ionization transport: True` (their READMEs).
Their IC pairs must sit on the grid the binary builds, so whenever the grid
moves (cell width, R_J, ...) they are rebuilt from a fresh mol_base_handoff
run with this script:

    ./make_hp_initial_conditions.py [mol_base_handoff/output]

For every case `IC/Hydro_ioniz_IC.txt` is the source `Hydro_ioniz.txt`
unchanged and `IC/Ion_species_IC.txt` is the source `Ion_species.txt` with
the H II column edited row by row, the removed protons moved into H I so
that the hydrogen nuclei of every row are unchanged (load_IC refuses a
changed He/H, see the case READMEs):

    hp_front       H II as written
    hp_zero_seed   H II = 0
    hp_trace_seed  H II = 1e-12 of the hydrogen nuclei of the row
                   (H I + H II + 2 H2 + 2 H2+ + 3 H3+ + HeH+)

Comment lines (schema, columns, rows, coupling, provenance) are copied
verbatim; numbers are rewritten with repr(), which round-trips a double.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "mol_base_handoff", "output")

CASES = {"hp_front": None, "hp_zero_seed": 0.0, "hp_trace_seed": 1.0e-12}


def read_table(path):
    header, rows = [], []
    with open(path) as f:
        for line in f:
            if line.lstrip().startswith("#"):
                header.append(line.rstrip("\n"))
            elif line.strip():
                rows.append([float(x) for x in line.split()])
    return header, rows


def column_index(header, name):
    for h in header:
        if h.startswith("# columns"):
            cols = h.split()[2:]
            return cols.index(name) + 0  # cols[0] is r
    raise SystemExit("no '# columns' header in the source file")


def main():
    hyd = os.path.join(SRC, "Hydro_ioniz.txt")
    ion = os.path.join(SRC, "Ion_species.txt")
    header, rows = read_table(ion)
    iHI, iHII = column_index(header, "HI"), column_index(header, "HII")
    cols = [h for h in header if h.startswith("# columns")][0].split()[2:]
    def idx(name):
        return cols.index(name) if name in cols else None
    iH2, iH2p, iH3p, iHeHp = idx("H2"), idx("H2p"), idx("H3p"), idx("HeHp")
    for case, seed in CASES.items():
        icdir = os.path.join(HERE, case, "IC")
        os.makedirs(icdir, exist_ok=True)
        with open(hyd) as f, open(os.path.join(icdir, "Hydro_ioniz_IC.txt"), "w") as g:
            g.write(f.read())
        with open(os.path.join(icdir, "Ion_species_IC.txt"), "w") as g:
            for h in header:
                g.write(h + "\n")
            for row in rows:
                r = list(row)
                if seed is not None:
                    nuc = r[iHI] + r[iHII]
                    for k, mult in ((iH2, 2.0), (iH2p, 2.0), (iH3p, 3.0), (iHeHp, 1.0)):
                        if k is not None:
                            nuc += mult * r[k]
                    new_hii = seed * nuc
                    r[iHI] = r[iHI] + r[iHII] - new_hii
                    r[iHII] = new_hii
                g.write("  " + "  ".join(repr(x) for x in r) + "\n")
        print(f"{case}: IC pair written from {SRC} (H II seed: {'as written' if seed is None else seed})")


if __name__ == "__main__":
    main()
