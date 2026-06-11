"""Fetch NIST ASD Fe I atomic data for the Fe I line-cooling builder.

CHIANTI v11 has NO Fe I model atom (only fe_1.diparams), so - exactly as in
Huang et al. (2023) Section 2.5 - Fe I line cooling is built from NIST
oscillator strengths + the Van Regemorter approximation. This script downloads
the two auditable raw inputs:

  fe1_nist_lines.tsv   permitted (E1) Fe I lines: Aki, g_i, g_k, E_i, E_k.
                       Fetched over a wide wavelength range, chunked to dodge
                       any per-query row cap. Forbidden (M1/E2) lines and lines
                       with no transition probability are kept in the raw dump
                       but filtered out downstream.
  fe1_nist_levels.tsv  every Fe I level (g, E) for the partition function U(T)
                       and the Boltzmann lower-level populations.

Provenance: NIST Atomic Spectra Database (ASD), https://physics.nist.gov/asd .
Re-run to refresh; the cooling builder reads only these local files.
"""
import io
import os
import sys
import time
import urllib.parse
import urllib.request

OUTDIR = os.path.dirname(os.path.abspath(__file__))
LINES_URL = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl"
LEVELS_URL = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl"

# Wide enough to capture every Fe I coolant: the strong near-UV resonance
# multiplets (~2300-4500 A) dominate, but cool down to the IR is included.
WL_LO, WL_HI, WL_STEP = 1500.0, 13000.0, 1000.0


def _get(url, params):
    q = urllib.parse.urlencode(params)
    req = urllib.request.Request(url + "?" + q,
                                 headers={"User-Agent": "ATES-FeI-cooling/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


def fetch_lines():
    rows = []
    header = None
    for w0 in _frange(WL_LO, WL_HI, WL_STEP):
        w1 = w0 + WL_STEP
        params = [
            ("spectra", "Fe I"), ("limits_type", "0"),
            ("low_w", f"{w0:.0f}"), ("upp_w", f"{w1:.0f}"), ("unit", "0"),
            ("de", "0"), ("format", "3"), ("line_out", "0"),
            ("remove_js", "on"), ("en_unit", "0"), ("output", "0"),
            ("page_size", "15"),
            ("show_obs_wl", "1"), ("show_calc_wl", "1"),
            ("order_out", "0"), ("show_av", "2"), ("A_out", "0"),
            ("allowed_out", "1"), ("forbid_out", "1"),
            ("conf_out", "on"), ("term_out", "on"),
            ("enrg_out", "on"), ("J_out", "on"), ("g_out", "on"),
            ("submit", "Retrieve Data"),
        ]
        txt = _get(LINES_URL, params)
        lines = [ln for ln in txt.splitlines() if ln.strip()]
        if not lines or "\t" not in lines[0]:
            print(f"  [{w0:.0f}-{w1:.0f} A] no tabular data "
                  f"({'Input Error' if 'Input Error' in txt else 'empty'})")
            time.sleep(0.5)
            continue
        if header is None:
            header = lines[0]
        body = [ln for ln in lines[1:] if ln != header]
        rows.extend(body)
        print(f"  [{w0:.0f}-{w1:.0f} A] {len(body)} rows")
        time.sleep(0.5)
    path = os.path.join(OUTDIR, "fe1_nist_lines.tsv")
    with open(path, "w") as fh:
        fh.write("# NIST ASD Fe I lines (E1+forbidden, with energies/g). "
                 f"Fetched {WL_LO:.0f}-{WL_HI:.0f} A.\n")
        fh.write("# Type blank = E1 (permitted); downstream keeps blank-Type "
                 "rows that carry an Aki.\n")
        fh.write(header + "\n")
        fh.write("\n".join(rows) + "\n")
    print(f"wrote {path}  ({len(rows)} line rows)")


def fetch_levels():
    params = [
        ("de", "0"), ("spectrum", "Fe I"), ("units", "0"),
        ("format", "3"), ("output", "0"), ("page_size", "15"),
        ("multiplet_ordered", "0"),
        ("conf_out", "on"), ("term_out", "on"), ("level_out", "on"),
        ("j_out", "on"), ("g_out", "on"), ("temp", ""),
        ("submit", "Retrieve Data"),
    ]
    txt = _get(LEVELS_URL, params)
    lines = [ln for ln in txt.splitlines() if ln.strip()]
    path = os.path.join(OUTDIR, "fe1_nist_levels.tsv")
    with open(path, "w") as fh:
        fh.write("# NIST ASD Fe I energy levels (Configuration, Term, J, g, "
                 "Level cm-1).\n")
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {path}  ({len(lines) - 1} level rows)")


def _frange(a, b, step):
    x = a
    while x < b - 1e-9:
        yield x
        x += step


if __name__ == "__main__":
    print("Fetching Fe I lines from NIST ASD ...")
    fetch_lines()
    print("Fetching Fe I levels from NIST ASD ...")
    fetch_levels()
