#!/usr/bin/env python3
"""Analytical counterexamples supporting the September 6 plan review.

These are mathematical checks, not executions of EXHALE's Fortran routines.
They write no files and require only the Python standard library.
"""

from math import exp


def main():
    # An exact implicit solve need not accurately resolve the physical time.
    rate = 1.0
    duration = 10.0
    initial = 1.0
    implicit = initial / (1.0 + rate * duration)
    residual = implicit - initial + duration * rate * implicit
    exact = initial * exp(-rate * duration)
    assert abs(residual) < 1.0e-14
    assert abs(implicit - exact) > 0.09
    print("DECAY: dy/dt=-y, y(0)=1, final time=10")
    print(f"  BE residual={residual:.9e}; BE={implicit:.9e}; exact={exact:.9e}")
    previous_error = abs(implicit - exact)
    for steps in (2, 4, 8, 16, 32, 64):
        value = initial / (1.0 + rate * duration / steps) ** steps
        error = abs(value - exact)
        assert error < previous_error
        previous_error = error
        print(f"  steps={steps:2d}; value={value:.9e}; absolute error={error:.9e}")

    # A sum over species can hide an unbalanced species with smaller rates.
    row_residuals = (0.0, 1.0)
    row_scales = (1.0e8, 1.0)
    aggregate = sum(map(abs, row_residuals)) / sum(row_scales)
    worst_row = max(abs(r) / s for r, s in zip(row_residuals, row_scales))
    assert aggregate < 1.0e-6 and worst_row > 1.0e-6
    print(f"AGGREGATION: combined={aggregate:.9e}; worst row={worst_row:.9e}")

    # Arithmetic admitted by the current relative-drop stagnation criterion.
    iteration, start, previous, current = 3, 1.0e8, 100.0, 60.0
    floor, drop = 1.0e-8, 1.0e-6
    accepted_by_drop = (iteration >= 3 and current > 0.5 * previous
                        and (current <= floor or current <= drop * start))
    assert accepted_by_drop and current > floor
    print(f"RELATIVE DROP: accepted={accepted_by_drop}; residual={current:g}; "
          f"absolute floor={floor:g}")
    print("All mathematical assertions passed; no EXHALE simulation was run.")


if __name__ == "__main__":
    main()
