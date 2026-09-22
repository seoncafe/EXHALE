"""Small algebra checks for the published CIP formulas, not a hydro solver.

Source: Frelikh & Murray-Clay, ApJ 996, 96, published Appendix A,
equations A9, A10, and A26-A27. No EXHALE source is imported or changed.
All numerical examples use normalized units and exactly representable inputs.
"""

import json


def main():
    # A linear profile must be reproduced by a cubic Hermite interpolant.
    h, f0, f1, g0, g1 = 1.0, 2.0, 3.0, 1.0, 1.0
    a = (g0 + g1) / h**2 - 2.0 * (f1 - f0) / h**3
    b_printed = 3.0 * (f1 - f0) / h**2 - (g1 - 2.0 * g0) / h
    b_consistent = 3.0 * (f1 - f0) / h**2 - (g1 + 2.0 * g0) / h

    def endpoint_errors(b):
        return {
            "value": a * h**3 + b * h**2 + g0 * h + f0 - f1,
            "gradient": 3.0 * a * h**2 + 2.0 * b * h + g0 - g1,
        }

    printed = endpoint_errors(b_printed)
    consistent = endpoint_errors(b_consistent)
    assert printed == {"value": 4.0, "gradient": 8.0}
    assert consistent == {"value": 0.0, "gradient": 0.0}

    # Constant radial velocity still has spherical divergence 2 v / r.
    radius, dr, velocity, density = 2.0, 1.0, 0.5, 1.0
    rp, rm = radius + dr / 2.0, radius - dr / 2.0
    density_rate_printed = -density / radius**2 * (velocity - velocity) / dr
    density_rate_consistent = -density / radius**2 * (
        rp**2 * velocity - rm**2 * velocity
    ) / dr
    assert density_rate_printed == 0.0
    assert density_rate_consistent == -2.0 * density * velocity / radius

    # With the stellar term and pressure gradient removed, gravity is inward.
    gm = 1.0
    acceleration_printed = gm / radius**2
    acceleration_consistent = -gm / radius**2
    assert acceleration_printed == 0.25
    assert acceleration_consistent == -0.25

    print(json.dumps({
        "scope": "Published-formula algebra only; no atmospheric simulation",
        "A27_linear_profile_endpoint_errors_printed": printed,
        "A27_linear_profile_endpoint_errors_consistent": consistent,
        "A9_constant_radial_velocity_density_rate_printed": density_rate_printed,
        "A9_constant_radial_velocity_density_rate_consistent": density_rate_consistent,
        "A10_planetary_gravity_printed": acceleration_printed,
        "A10_planetary_gravity_consistent": acceleration_consistent,
        "assertions": "passed",
    }, indent=2))


if __name__ == "__main__":
    main()
