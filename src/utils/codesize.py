#!/usr/bin/env python3
"""Code-size and diff measurement behind the Update_EXHALE_stage1.tex appendices
("Code-size summary" and "Fortran source inherited unchanged from ATES").

Methodology (as stated in the appendix):
- SLOC excludes blank lines and comments; the comment stripper is
  string-aware, so a '!' inside a quoted string is not treated as a comment.
- The line-by-line diff additionally normalizes whitespace and case, so
  reformatting and comment-only edits are not counted as changes.
- Byte-identical files are compared with no normalization at all.

Usage:
    python3 src/utils/codesize.py [ATES_ROOT]

ATES_ROOT defaults to ../ATES/ATES-Code-main relative to the EXHALE
repository root (the pristine upstream tree this appendix diffs against).
The EXHALE source root is derived from this script's own location.
"""
import os, sys, difflib, filecmp

HERE = os.path.dirname(os.path.abspath(__file__))
EX_ROOT = os.path.dirname(os.path.dirname(HERE))          # .../EXHALE
EX = os.path.join(EX_ROOT, 'src')
ATES = sys.argv[1] if len(sys.argv) > 1 else \
    os.path.join(os.path.dirname(EX_ROOT), 'ATES', 'ATES-Code-main')

# Known renames between the trees (ATES basename -> EXHALE basename).
RENAME = {'ATES_main.f90': 'EXHALE_main.f90'}


def strip_comment(line):
    out, q = [], None
    for ch in line:
        if q:
            out.append(ch)
            if ch == q:
                q = None
        else:
            if ch in ('"', "'"):
                q = ch
                out.append(ch)
            elif ch == '!':
                break
            else:
                out.append(ch)
    return ''.join(out)


def code_lines(path):
    lines = []
    with open(path, errors='replace') as f:
        for raw in f:
            s = strip_comment(raw).strip()
            if s:
                lines.append(s)
    return lines


def norm_lines(path):
    return [' '.join(l.lower().split()) for l in code_lines(path)]


def ffiles(root):
    out = []
    for dp, _, fns in os.walk(root):
        for fn in fns:
            if fn.endswith('.f90'):
                out.append(os.path.join(dp, fn))
    return sorted(out)


def sloc(files):
    return sum(len(code_lines(f)) for f in files)


def rel(path):
    return path.replace(EX + '/', '').replace('modules/', '')


def main():
    ates = ffiles(ATES)
    ex_all = ffiles(EX)
    ex_wae = [f for f in ex_all if '/wind_ae/' in f]
    ex_core = [f for f in ex_all if '/wind_ae/' not in f]

    print(f"ATES ({ATES}): {len(ates)} files, {sloc(ates)} SLOC")
    print(f"EXHALE core: {len(ex_core)} files, {sloc(ex_core)} SLOC")
    print(f"wind_ae: {len(ex_wae)} files, {sloc(ex_wae)} SLOC")
    print(f"EXHALE total: {len(ex_all)} files, {sloc(ex_all)} SLOC")

    amap = {os.path.basename(f): f for f in ates}
    emap = {os.path.basename(f): f for f in ex_core}

    added = removed = 0
    modified, unchanged, byteident, missing = [], [], [], []
    for aname, apath in sorted(amap.items()):
        ename = RENAME.get(aname, aname)
        if ename not in emap:
            missing.append(aname)
            removed += len(code_lines(apath))
            continue
        epath = emap[ename]
        if filecmp.cmp(apath, epath, shallow=False):
            byteident.append(ename)
            unchanged.append(ename)
            continue
        a, e = norm_lines(apath), norm_lines(epath)
        if a == e:
            unchanged.append(ename)
            continue
        modified.append(ename)
        sm = difflib.SequenceMatcher(a=a, b=e, autojunk=False)
        for tag, i1, i2, j1, j2 in sm.get_opcodes():
            if tag in ('replace', 'delete'):
                removed += i2 - i1
            if tag in ('replace', 'insert'):
                added += j2 - j1

    new_core = sorted(set(emap) - {RENAME.get(a, a) for a in amap})
    new_sloc = {n: len(code_lines(emap[n])) for n in new_core}
    added += sum(new_sloc.values())

    print(f"\noriginal files: modified {len(modified)}, "
          f"unchanged(normalized) {len(unchanged)}, removed {len(missing)}")
    print(f"byte-identical: {len(byteident)}")
    print(f"diff core: added {added}, removed {removed}, "
          f"total {added + removed}")

    print(f"\nnew core modules ({len(new_core)}), SLOC (desc):")
    for n in sorted(new_core, key=lambda x: -new_sloc[x]):
        print(f"  {rel(emap[n]):55s} {new_sloc[n]}")
    print(f"new-module total SLOC: {sum(new_sloc.values())}")

    print("\nbyte-identical files:")
    for n in sorted(byteident):
        print("  ", rel(emap[n]))
    only_norm = sorted(set(unchanged) - set(byteident))
    if only_norm:
        print("\nnormalized-unchanged but not byte-identical:")
        for n in only_norm:
            print("  ", rel(emap[n]))
    if missing:
        print(f"\nmissing (in ATES, no EXHALE match): {missing}")


if __name__ == '__main__':
    main()
