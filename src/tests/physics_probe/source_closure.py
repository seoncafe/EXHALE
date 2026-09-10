#!/usr/bin/env python3
"""Dependency-ordered list of the production sources a test driver needs.

Given one or more driver files, scan their `use` statements, follow them
through the modules under src/modules (and build/build_stamp.f90 when the
Makefile has generated it), and print the closure in an order in which each
file compiles after every module it uses.  The drivers under
src/tests/physics_probe link the production sources themselves, never a
copy of them, and this is what tells the build script which ones and in
which order.

Usage:  source_closure.py <repo root> <driver.f90> [driver.f90 ...]
"""
import os
import re
import sys

module_definition = re.compile(r'^\s*module\s+([a-zA-Z]\w*)', re.I)
module_procedure = re.compile(r'^\s*module\s+procedure\b', re.I)
module_use = re.compile(r'^\s*use\s*(?:,\s*intrinsic\s*)?(?:::)?\s*([a-zA-Z]\w*)', re.I)


def scan(path):
    provides, uses = set(), set()
    with open(path, encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if module_procedure.match(line):
                continue
            hit = module_definition.match(line)
            if hit:
                provides.add(hit.group(1).lower())
            hit = module_use.match(line)
            if hit:
                uses.add(hit.group(1).lower())
    return provides, uses


def main():
    root, drivers = sys.argv[1], sys.argv[2:]
    bases = [os.path.join(root, 'src', 'modules'), os.path.join(root, 'build')]
    bases += [os.path.dirname(os.path.abspath(d)) for d in drivers]
    pool = []
    for base in bases:
        if not os.path.isdir(base):
            continue
        for here, _, names in os.walk(base):
            if '.ipynb_checkpoints' in here.split(os.sep):
                continue   # editor auto-save copies are not source
            for name in sorted(names):
                if name.endswith('.f90'):
                    pool.append(os.path.join(here, name))
    provides, uses, owner = {}, {}, {}
    pool = [p for p in pool if os.path.abspath(p) not in
            [os.path.abspath(d) for d in drivers]]
    for path in pool + drivers:
        p, u = scan(path)
        provides[path], uses[path] = p, u
        for m in p:
            owner.setdefault(m, path)

    needed, stack = set(), list(drivers)
    while stack:
        path = stack.pop()
        if path in needed:
            continue
        needed.add(path)
        for m in uses[path]:
            src = owner.get(m)
            if src is not None and src not in needed:
                stack.append(src)

    ordered, placed = [], set()
    remaining = sorted(needed)
    while remaining:
        ready = [p for p in remaining
                 if all(owner.get(m) in placed or owner.get(m) is None
                        or owner.get(m) == p for m in uses[p])]
        if not ready:
            sys.exit('circular or unresolved module dependency in: %s' % remaining)
        for p in ready:
            ordered.append(p)
            placed.add(p)
        remaining = [p for p in remaining if p not in placed]
    # the drivers are programs, not modules: they always compile last
    ordered = [p for p in ordered if p not in drivers] + list(drivers)
    print('\n'.join(ordered))


if __name__ == '__main__':
    main()
