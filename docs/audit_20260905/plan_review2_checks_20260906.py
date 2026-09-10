#!/usr/bin/env python3
"""Mathematical examples for PLAN_20260906_rev1, not Fortran regression tests.

Run with Python 3. This script writes no files and uses only the standard library.
"""


def steady_residual(x):
    return x * x - 2.0


def main():
    tolerance = 1.0e-8
    initial = 1.0
    direction = -steady_residual(initial) / (2.0 * initial)
    trial = initial + direction
    residual = steady_residual(trial)
    merit_ratio = residual ** 2 / steady_residual(initial) ** 2
    assert trial > 0.0 and merit_ratio < 1.0
    assert abs(residual) > tolerance
    trials = [initial + direction * 0.5 ** k for k in range(20)]
    assert all(abs(steady_residual(x)) > tolerance for x in trials)
    print(f"NEWTON: trial={trial:g}, residual={residual:g}, "
          f"squared-merit ratio={merit_ratio:g}")
    print("  All 20 first-iteration backtracks fail final convergence,")
    print("  although the full Newton step is admissible and decreases merit.")
    x = initial
    for iteration in range(1, 11):
        x -= steady_residual(x) / (2.0 * x)
        if abs(steady_residual(x)) <= tolerance:
            break
    assert abs(steady_residual(x)) <= tolerance
    print(f"  Ordinary Newton reaches the tolerance in {iteration} iterations.")

    # A valid transient solves a time-discrete balance, not the steady equation.
    old, dt, source = 1.0, 0.1, -1.0
    new = old + dt * source
    time_residual = (new - old) / dt - source
    assert abs(time_residual) < 1.0e-14 and abs(source) > tolerance
    print(f"TRANSIENT: time-discrete residual={time_residual:.9e}, "
          f"steady source={source:g}")

    # Unit volumes, a single shared internal face, no external boundary flux.
    left, right, face_flux = 1.0, 1.0, 1.0
    dt_left, dt_right = 0.1, 0.2
    unequal_sum = left - dt_left * face_flux + right + dt_right * face_flux
    common_sum = left - dt_left * face_flux + right + dt_left * face_flux
    assert abs(common_sum - 2.0) < 1.0e-14
    assert abs(unequal_sum - 2.1) < 1.0e-14
    print(f"LOCAL CLOCKS: initial sum=2, unequal-step sum={unequal_sum:g}, "
          f"common-step sum={common_sum:g}")
    print("All mathematical assertions passed; no EXHALE simulation was run.")


if __name__ == "__main__":
    main()
