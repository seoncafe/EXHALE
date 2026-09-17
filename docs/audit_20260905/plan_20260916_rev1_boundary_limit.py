"""Analytical probe of the proposed normalized-RMS boundary gate.

This is not an execution of EXHALE. The cubic weight follows
base_boundary.f90; the RMS replacement is the proposal in PLAN rev1.
Flux is normalized so that its mean is the wind Mach number.
The fixed local branch weight is zero. No files are written.
"""

from math import sqrt


def weight(x):
    if x <= -1.0:
        return 1.0
    if x >= 1.0:
        return 0.0
    return 0.5 - 0.75 * x + 0.25 * x**3


def boundary_weight(values):
    mean = sum(values) / len(values)
    rms = sqrt(sum((x - mean)**2 for x in values) / len(values))
    relative_spread = rms / abs(mean)
    window_weight = weight((2 * relative_spread - 0.03) / 0.01)
    return window_weight * weight(mean / 1e-8)


print("epsilon       uniform_path       nonuniform_path")
for eps in (1e-10, 1e-12, 1e-14, 1e-16):
    uniform = boundary_weight([eps, eps])
    nonuniform = boundary_weight([eps, 2 * eps])
    print(f"{eps:.1e}       {uniform:.12f}     {nonuniform:.12f}")
    assert nonuniform == 0.0
assert abs(boundary_weight([1e-20, 1e-20]) - 0.5) < 1e-10
print("Distinct directional limits: changing only the value at zero cannot repair continuity.")
