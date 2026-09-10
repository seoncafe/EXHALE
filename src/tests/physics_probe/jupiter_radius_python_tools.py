#!/usr/bin/env python3
"""One Jupiter radius and one Jupiter mass for the Fortran and for the Python.

`Planet radius [R_J]` is read by the code and, separately, by the tools that
prepare its input (the lower-atmosphere and VULCAN drivers, which turn the key
into a base radius and a surface gravity) and by the tools that analyse its
output (the plotter, the transit library, the Roche reconstruction, the
heating-efficiency estimator).  If those hold their own copy of R_J, the run
and its analysis are about two different planets, and a regenerated `base.inp`
is built for a third.

So this test asserts three things:

  1. the Fortran definition is the IAU 2015 nominal equatorial radius
     R_J^N(eq) = 7.1492e9 cm (Prsa et al. 2016, AJ 152, 41, Table 1), read out
     of src/modules/init/parameters.f90 as text;
  2. the one Python definition, RJ_CM of examples/exhale_io.py, is that same
     number, to round-off;
  3. every listed tool resolves its own R_J to that number, and none of them
     still carries the volumetric-mean literal 6.9911e9 cm (or its value in
     metres, 6.9911e7) anywhere in its source;
  4. the Jupiter and solar MASSES of the Python definition are the IAU 2015
     nominal values of parameters.f90 as well, since the same key pair
     `Planet mass [M_J]` / `Parent star mass [M_sun]` is read on both sides.

Assertion 3 is made on the VALUE the tool ends up with, not on the presence of
an import line: a tool that imported the constant and then overwrote it would
pass a grep and fail here.  The tools that work in SI are checked against
RJ_CM/100.

Verdict lines follow the convention of the directory:
    PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
"""

import ast
import importlib
import importlib.util
import os
import re
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))

RJ_IAU_2015_NOMINAL_EQUATORIAL_CM = 7.1492e9
MJ_IAU_2015_NOMINAL_G = 1.8982e30
MSUN_IAU_2015_NOMINAL_G = 1.98842e33
# The volumetric-mean Jupiter radius, 6.9911e9 cm, which is a different
# radius of the same planet and is what the tools used to carry.
RJ_VOLUMETRIC_MEAN_CM = 6.9911e9

failures = 0


def check_relative(name, measured, reference, tol):
    global failures
    if reference == 0.0:
        ok = measured == 0.0
        rel = abs(measured)
    else:
        rel = abs(measured / reference - 1.0)
        ok = rel <= tol
    if not ok:
        failures += 1
    print("%s %s measured=%.15e reference=%.15e tol=%.3e relative"
          % ("PASS" if ok else "FAIL", name, measured, reference, tol))


def check_true(name, ok, detail):
    global failures
    if not ok:
        failures += 1
    print("%s %s measured=%s reference=%s tol=exact"
          % ("PASS" if ok else "FAIL", name, detail, "none"))


# ---- 1. the Fortran definition, read as text ------------------------------ #
params = open(os.path.join(ROOT, 'src', 'modules', 'init',
                           'parameters.f90')).read()
m = re.search(r'^\s*real\*8\s*,\s*parameter\s*::\s*RJ\s*=\s*([0-9.dDeE+-]+)',
              params, re.MULTILINE)
if m is None:
    print("FAIL jupiter_radius_declared_in_parameters_f90 measured=absent "
          "reference=one_declaration tol=exact")
    sys.exit(1)
rj_fortran = float(m.group(1).lower().replace('d', 'e'))
check_relative('jupiter_radius_of_parameters_f90', rj_fortran,
               RJ_IAU_2015_NOMINAL_EQUATORIAL_CM, 0.0)


def fortran_constant(name):
    """The value of a `real*8, parameter :: <name> = ...` of parameters.f90."""
    hit = re.search(r'^\s*real\*8\s*,?\s*parameter\s*::\s*%s\s*=\s*'
                    r'([0-9.dDeE+-]+)' % name, params, re.MULTILINE)
    if hit is None:
        return float('nan')
    return float(hit.group(1).lower().replace('d', 'e'))


mj_fortran = fortran_constant('MJ')
msun_fortran = fortran_constant('Msun')
check_relative('jupiter_mass_of_parameters_f90', mj_fortran,
               MJ_IAU_2015_NOMINAL_G, 0.0)
check_relative('solar_mass_of_parameters_f90', msun_fortran,
               MSUN_IAU_2015_NOMINAL_G, 0.0)

# ---- 2. the one Python definition ---------------------------------------- #
sys.path.append(os.path.join(ROOT, 'examples'))
sys.path.append(os.path.join(ROOT, 'src', 'utils'))
sys.path.append(ROOT)
exhale_io = importlib.import_module('exhale_io')
# Absent means there is no single Python definition yet; that is the failure
# this assertion exists to report, so it is reported and not raised, and the
# tool comparisons below then run against the Fortran value.
rj_python = getattr(exhale_io, 'RJ_CM', None)
check_true('jupiter_radius_of_the_python_definition_exists', rj_python
           is not None, 'examples/exhale_io.py defines RJ_CM: %s'
           % (rj_python is not None))
if rj_python is None:
    rj_python = rj_fortran
else:
    check_relative('jupiter_radius_of_the_python_definition', rj_python,
                   rj_fortran, 0.0)
check_relative('jupiter_mass_of_the_python_definition', exhale_io.MJ,
               mj_fortran, 0.0)
check_relative('solar_mass_of_the_python_definition', exhale_io.Msun,
               msun_fortran, 0.0)

# ---- 3. every tool resolves to it ---------------------------------------- #
# (module, attribute) pairs whose value is read by importing the module.  The
# SI lists hold the tools that work in metres, whose radius must be RJ_CM/100.
IMPORTED_CGS = [
    ('run_lower', 'RJ'),
    ('vulcan_to_base', 'RJ'),
    ('lower_profile_schema', 'RJ'),
]
IMPORTED_SI = [
    ('exhale_transit_lib', 'RJ'),
]

# EXHALE_plots reads ./input.inp at import, eta_approx asks the terminal for a
# radius, and roche_recon binds its R_J only inside its self-test, so none of
# the three can simply be imported here.  Each takes its value from a single
# module-level assignment, which is read as text and evaluated with RJ_CM
# bound: that is the value the tool runs with.
TEXT_CGS = [
    ('EXHALE_plots', 'RJ'),
    ('eta_approx', 'RJ'),
]
TEXT_SI = [
    ('roche_recon', 'RJ'),
]

PATHS = {
    'run_lower': 'src/utils/run_lower.py',
    'vulcan_to_base': 'src/utils/vulcan_to_base.py',
    'lower_profile_schema': 'src/utils/lower_profile_schema.py',
    'interface_state': 'src/utils/EXHALE_interface_state.py',
    'vulcan_driver': 'src/utils/vulcan_driver.py',
    'EXHALE_plots': 'EXHALE_plots.py',
    'eta_approx': 'eta_approx.py',
    'exhale_transit_lib': 'exhale_transit_lib.py',
    'roche_recon': 'roche_recon.py',
}


def value_from_text(mod, attr):
    """Evaluate the tool's own assignment of `attr`, with RJ_CM bound."""
    txt = open(os.path.join(ROOT, PATHS[mod])).read()
    pat = re.compile(r'^\s*%s\s*=\s*([^#\n]+)' % attr, re.MULTILINE)
    hits = pat.findall(txt)
    if len(hits) != 1:
        return float('nan')
    return eval(hits[0].strip(), {'RJ_CM': rj_python})


for mod, attr in IMPORTED_CGS:
    value = getattr(importlib.import_module(mod), attr)
    check_relative('jupiter_radius_of_' + mod, value, rj_python, 0.0)

for mod, attr in TEXT_CGS:
    check_relative('jupiter_radius_of_' + mod, value_from_text(mod, attr),
                   rj_python, 0.0)

for mod, attr in IMPORTED_SI:
    value = getattr(importlib.import_module(mod), attr)
    check_relative('jupiter_radius_of_' + mod, value,
                   rj_python*1.0e-2, 1.0e-15)

for mod, attr in TEXT_SI:
    check_relative('jupiter_radius_of_' + mod, value_from_text(mod, attr),
                   rj_python*1.0e-2, 1.0e-15)

# EXHALE_interface_state.py fills its constants inside init(), so the module
# has to be executed and init() called before RJ exists.  It lives in
# src/utils/, which is not on this driver's import path, so it is loaded by
# path like the other tools of PATHS.
_spec = importlib.util.spec_from_file_location(
    'exhale_interface_state', os.path.join(ROOT, PATHS['interface_state']))
interface_state = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(interface_state)
interface_state.init()
check_relative('jupiter_radius_of_interface_state', interface_state.RJ,
               rj_python, 0.0)

vulcan_driver_txt = open(os.path.join(ROOT, PATHS['vulcan_driver'])).read()
check_true('jupiter_radius_of_vulcan_driver',
           re.search(r'from\s+run_lower\s+import[^\n]*\bRJ\b',
                     vulcan_driver_txt) is not None,
           'imports RJ from run_lower')

# ---- the tools that take the constants from exhale_io -------------------- #
# These are analysis and input-preparation tools whose radius is now a single
# module-level assignment fed by the import.  None of them can be imported in
# a test process without side effects (figures, HDF5 files, a profile write),
# so the assignment is read from the source with the ast module and evaluated
# with the names the import binds.  It is still the VALUE that is asserted:
# an assignment that imported the constant and then overwrote it fails here.
_aio = types.SimpleNamespace(RJ_CM=rj_python, MJ=exhale_io.MJ,
                             Msun=exhale_io.Msun)
_NAMESPACE = {'RJ_CM': rj_python, 'MJ_G': exhale_io.MJ,
              'aio': _aio, 'exhale_io': _aio}

IMPORTING_TOOLS = [
    # (path, name, expected value)
    ('examples/lya_insitu_emissivity.py', 'RJ', rj_python),
    ('examples/tpm_halpha_lart2d.py', 'RJ_M', rj_python*1.0e-2),
    ('examples/17_lower_profile/make_example_profile.py', 'RJ', rj_python),
    ('examples/17_lower_profile/make_example_profile.py', 'MJ',
     exhale_io.MJ),
    ('docs/compare_vulcan_photochem.py', 'RJ', rj_python),
    ('docs/compare_vulcan_photochem.py', 'MJ', exhale_io.MJ),
    ('LHS1140b/compare_exhale_pwinds.py', 'RJ', rj_python),
    ('LHS1140b/make_memo_figures.py', 'R_JUP_CM', rj_python),
    ('LHS1140b/make_memo_figures.py', 'R_JUP', rj_python),
]


def module_level_value(rel, name):
    """Value of the single module-level assignment of `name` in `rel`.

    Handles both `name = expr` and the tuple form `a, name, b = x, y, z`.
    Returns NaN when the name is assigned zero times or more than once, so
    an ambiguous tool fails rather than being silently skipped.
    """
    tree = ast.parse(open(os.path.join(ROOT, rel)).read())
    found = []
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id == name:
                found.append(node.value)
            elif isinstance(target, ast.Tuple) and isinstance(node.value,
                                                             ast.Tuple):
                for elt, val in zip(target.elts, node.value.elts):
                    if isinstance(elt, ast.Name) and elt.id == name:
                        found.append(val)
    if len(found) != 1:
        return float('nan')
    return eval(compile(ast.Expression(found[0]), '<tool>', 'eval'),
                dict(_NAMESPACE))


for rel, name, expected in IMPORTING_TOOLS:
    label = os.path.basename(rel)[:-3]
    check_relative('%s_of_%s' % (name.lower(), label),
                   module_level_value(rel, name), expected, 1.0e-15)

for rel, _name, _expected in IMPORTING_TOOLS:
    PATHS.setdefault(os.path.basename(rel)[:-3] + '_' + rel.split('/')[0],
                     rel)

# ---- the volumetric-mean literal is gone from all of them ---------------- #
stale = []
for mod, rel in sorted(PATHS.items()):
    txt = open(os.path.join(ROOT, rel)).read()
    if re.search(r'%s[eE][+]?[79]' % re.escape('%.4f'
                 % (RJ_VOLUMETRIC_MEAN_CM/1.0e9)), txt):
        stale.append(rel)
check_true('volumetric_mean_jupiter_radius_absent_from_the_tools',
           not stale, 'files still carrying 6.9911e9/6.9911e7: %s'
           % (', '.join(stale) if stale else 'none'))

if failures:
    print('jupiter_radius_python_tools: %d assertion(s) failed' % failures)
    sys.exit(1)
print('jupiter_radius_python_tools: all assertions passed')
