#!/usr/bin/env python3
"""Read-only SED diagnostics for the September 6 decision review.

Run with Python 3 from any directory. No spectra or simulation outputs are
written. Integration uses piecewise-linear F_lambda with interpolated band
boundaries, not the production reader's F_E quadrature. Published comparison
values are from Huang et al. (2023), ApJ 951, 123, Section 2.2.
"""

from bisect import bisect_left, bisect_right
from hashlib import sha256
from pathlib import Path


SED_DIR = Path(__file__).resolve().parents[2] / "inputdata" / "sed"


def read_spectrum(name):
    path = SED_DIR / name
    payload = path.read_bytes()
    rows = [tuple(map(float, line.split())) for line in payload.decode().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    w, f = zip(*rows)
    assert all(b > a for a, b in zip(w, w[1:]))
    print(name, "sha256=" + sha256(payload).hexdigest())
    return w, f


def flux_at(w, f, x):
    i = min(max(bisect_right(w, x) - 1, 0), len(w) - 2)
    return f[i] + (f[i + 1] - f[i]) * (x - w[i]) / (w[i + 1] - w[i])


def band(w, f, lo, hi):
    assert w[0] <= lo < hi <= w[-1], "Requested band is not covered"
    points = [lo] + [x for x in w if lo < x < hi] + [hi]
    values = [flux_at(w, f, x) for x in points]
    return sum((b - a) * (fa + fb) / 2
               for a, b, fa, fb in zip(points, points[1:], values, values[1:]))


def main():
    print("Band flux units: erg cm^-2 s^-1; wavelength units: A")
    stored = {}
    for name in ("wasp121b_solar_huang2023.txt",
                 "wasp121b_wasp17_huang2023.txt",
                 "hd189733b_epseri_salz2016.txt",
                 "hd189733b_bourrier2020.txt"):
        w, f = read_spectrum(name)
        stored[name] = w, f
        for lo, hi in ((10, 912), (100, 912), (912, 1110), (1700, 2583)):
            print(f"  F({lo}-{hi}) = {band(w, f, lo, hi):.8e}")
        if name.startswith("wasp121"):
            ratios = []
            for lo, hi, reference in ((700, 800, 1.4e4),
                                      (800, 912, 3.6e4),
                                      (912, 1170, 6.3e4)):
                value = band(w, f, lo, hi)
                ratios.append(value / reference)
                print(f"  Huang F({lo}-{hi}): stored={value:.8e} "
                      f"published={reference:.8e} ratio={value / reference:.6f}")
            print("  band-shape ratios relative to the middle band:",
                  " ".join(f"{x / ratios[1]:.6f}" for x in ratios))
            i = bisect_left(w, 1700)
            print(f"  adjacent-node join ratio = {f[i] / f[i - 1]:.6f}")
    a = band(*stored["wasp121b_solar_huang2023.txt"], 1700, 2583)
    b = band(*stored["wasp121b_wasp17_huang2023.txt"], 1700, 2583)
    print(f"WASP-121 NUV integrated-flux A/B = {a / b:.6f}")
    print("No He 2^3S photoionization rate or wind solution was evaluated.")


if __name__ == "__main__":
    main()
