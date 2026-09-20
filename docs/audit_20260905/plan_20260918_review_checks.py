"""Arithmetic counterexamples for the September 18 plan review.

These checks demonstrate mathematical implications of inspected expressions.
They do not execute the EXHALE solver, grid generator, or chemistry network.
Run with Python 3; no third-party packages are required.
"""
import math
import struct


def branch_weight(x):
    if x >= 1.0:
        return 0.0
    if x <= -1.0:
        return 1.0
    return 0.5 - 0.75 * x + 0.25 * x**3


legacy_width = struct.unpack("f", struct.pack("f", 2.0e-4))[0]
explicit_width = float("2.0e-4")
print("D1 binary32 literal widened to binary64:", format(legacy_width, ".17g"))
print("D1 binary64 parsed decimal:", format(explicit_width, ".17g"))
print("D1 relative width difference:", abs(legacy_width / explicit_width - 1.0))
print("D1 legacy explicit spelling round-trips:",
      float("0.00019999999494757503") == legacy_width)

print("D2 original zero weight:", branch_weight(0.0))
for h in (1.0e-2, 1.0e-4, 1.0e-6):
    print("D2 h, original weight, slope if only zero changes to 0:",
          h, branch_weight(h), branch_weight(h) / h)
print("D2 Euler characteristic speeds at rest, sound speed 1:", (-1.0, 0.0, 1.0))

# A fast reaction with an exactly balanced, nonzero transport divergence.
production, loss, divergence, face_magnitudes = 1.0e8, 1.0e8 - 1.0, 1.0, 3.0
net_source = production - loss
residual = divergence - net_source
print("D3 exact balanced row / current-style scale:",
      abs(residual) / (face_magnitudes + abs(net_source)))

# A turnover scale alone can hide a large transport-relative defect.
wrong_divergence = 2.0
wrong_residual = wrong_divergence - net_source
print("D3 wrong row / turnover scale:",
      abs(wrong_residual) / (face_magnitudes + production + loss))
print("D3 wrong row / transport-source scale:",
      abs(wrong_residual) / (face_magnitudes + abs(net_source)))

# The inspected helium projection leaves an ion sum below one unchanged.
n_he_total, n_hehp, x_heii, x_heiii = 1.0, 0.1, 0.6, 0.35
ion_sum = x_heii + x_heiii
neutral = max(n_he_total * (1.0 - ion_sum) - n_hehp, 0.0)
returned_total = n_he_total * ion_sum + n_hehp + neutral
print("Helium stage sum admitted by unit simplex:", ion_sum)
print("Helium nuclei after inspected clipping/write-back expressions:", returned_total)
print("Helium physically available ion-fraction bound:", 1.0 - n_hehp / n_he_total)

# Convergence does not require monotone changes in a particular norm.
state = (10.0, 1.1)
changes = []
for _ in range(4):
    next_state = (10.0 * state[1], 0.01 * state[0])
    changes.append(max(abs(a - b) for a, b in zip(next_state, state)))
    state = next_state
print("D4 convergent linear map spectral radius:", math.sqrt(0.1))
print("D4 nonmonotone successive infinity-norm changes:", changes)
