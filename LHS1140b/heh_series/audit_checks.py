#!/usr/bin/env python3
"""Phase C acceptance checks on a set of converged He/H runs.

For each case directory with a finished run (output/ + run.log), report:
  1. finiteness of every column of Hydro_ioniz, Ion_species, and the two
     breakdown files;
  2. charge neutrality: n_e written by the code (breakdown col 3) against
     n_e rebuilt from the ion densities (HII + HeII + 2 HeIII + metal ions);
  3. energy closure: heating and cooling channel sums against the totals
     (the files promise this; the audit verifies it in the He-rich limit);
  4. the dominant heating and cooling channel at the base, at 1.5 R_p and
     at the top, and the He-carried share of each budget;
  5. the steady-state Mdot recomputed from the profile (4 pi rho v r^2 near
     the outer boundary, with the 2D-approximation factor), the fractional
     radial spread of rho v r^2, and the JFNK exit state from run.log;
  6. how far the smallest densities of the solution sit above the numerical
     guard constants that divide by them.

Usage:
    audit_checks.py                       # the default tutorial-planet cases,
                                          # summary written next to this script
    audit_checks.py CASE [CASE ...]       # summary written to the parent
                                          # directory of the first case
    audit_checks.py -o FILE CASE [...]    # summary written to FILE
"""
import os, re, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'examples'))
import exhale_io

CASES = ['heh_1', 'heh_10', 'heh_100', 'heh_1000']
HEAT_COLS = ['HI', 'HeI', 'HeII', 'He23S', 'H2', 'metals', 'Hpe', 'Hdx',
             'He_recomb', 'He23S_Penning', 'He23S_H2_Penning', 'H2_LW']
COOL_COLS = ['reco', 'coio', 'coex_HI', 'coex_HeI', 'coex_HeII', 'brem',
             'H3p_IR']
ION_CHARGE = {'HII': 1, 'HeII': 1, 'HeIII': 2,
              'CII': 1, 'CIII': 2, 'OII': 1, 'OIII': 2, 'NII': 1, 'NIII': 2,
              'MgII': 1, 'MgIII': 2, 'SiII': 1, 'SiIII': 2, 'CaII': 1,
              'CaIII': 2, 'NaII': 1, 'KII': 1, 'SII': 1, 'FeII': 1, 'FeIII': 2}

# Numerical guards that divide by a density on the H/He (metals-off,
# molecules-off) path. Each entry is (label, guard value, source location).
# The audit reports how far the solution's smallest value of that density
# sits above the guard.
GUARDS = [
    ('n_H + n_He', 1.0e-99, 'ionization_equilibrium.f90:246, '
                            'util_ion_eq.f90:1081, post_process_adv.f90:284'),
    ('n_H',        1.0e-30, 'post_process_adv.f90:445 (metal electrons per H)'),
    ('n_HI',       1.0e-30, 'lya_rt.f90:159,205 (Lya source, off by default)'),
    ('n_e',        1.0e-30, 'Cool_coeff.f90:1359-1448 (metal fine-structure, '
                            'off without metals)'),
]


def finite_report(name, arr):
    bad = ~np.isfinite(arr)
    return f"{name}: {'all finite' if not bad.any() else f'{bad.sum()} non-finite entries'}"


def load_breakdown(path):
    d = np.loadtxt(path)
    return d


def dominant(cols, vals):
    i = int(np.argmax(vals))
    tot = vals.sum()
    return cols[i], (vals[i]/tot if tot > 0 else np.nan)


def jfnk_state(logpath):
    """Exit info and final ||R|| of the JFNK finish, from run.log."""
    info, res = None, None
    if not os.path.isfile(logpath):
        return info, res
    with open(logpath, errors='replace') as fh:
        for ln in fh:
            m = re.search(r'done info=\s*(-?\d+)\s+\|\|R\|\|=\s*([0-9.eE+-]+)', ln)
            if m:
                info, res = int(m.group(1)), float(m.group(2))
    return info, res


def audit(case):
    out = os.path.join(case, 'output')
    name = os.path.basename(os.path.normpath(case))
    lines = [f"## {name}"]
    if not os.path.isfile(os.path.join(out, 'Hydro_ioniz.txt')):
        lines.append("no output yet"); return lines, None
    _, ions = exhale_io.load_ions(os.path.join(out, 'Ion_species.txt'))
    hyd = np.loadtxt(os.path.join(out, 'Hydro_ioniz.txt'))
    ion = np.loadtxt(os.path.join(out, 'Ion_species.txt'))
    heat = load_breakdown(os.path.join(out, 'Heating_breakdown.txt'))
    cool = load_breakdown(os.path.join(out, 'Cooling_breakdown.txt'))

    # 1. finiteness
    nbad = 0
    for nm, a in (('Hydro_ioniz', hyd), ('Ion_species', ion),
                  ('Heating_breakdown', heat), ('Cooling_breakdown', cool)):
        nbad += int((~np.isfinite(a)).sum())
        lines.append("- " + finite_report(nm, a))

    # 2. charge neutrality
    ne_code = heat[:, 2]
    ne_ions = sum(ions[k]*q for k, q in ION_CHARGE.items() if k in ions)
    ok = ne_code > 0
    rel = np.abs(ne_ions[ok] - ne_code[ok])/ne_code[ok]
    lines.append(f"- charge neutrality: max |n_e(ions) - n_e(code)|/n_e = {rel.max():.2e}"
                 f" (median {np.median(rel):.2e})")

    # 3. energy closure
    hs = heat[:, 4:4+len(HEAT_COLS)].sum(axis=1); ht = heat[:, 3]
    cs = cool[:, 4:4+len(COOL_COLS)].sum(axis=1)
    if cool.shape[1] > 4+len(COOL_COLS):            # metal-ion columns
        cs = cs + cool[:, 4+len(COOL_COLS):].sum(axis=1)
    ct = cool[:, 3]
    hrel = np.abs(hs - ht)/np.maximum(np.abs(ht), 1e-300)
    crel = np.abs(cs - ct)/np.maximum(np.abs(ct), 1e-300)
    hmax, cmax = hrel[ht > 0].max(), crel[ct > 0].max()
    lines.append(f"- heating channel closure: max rel. residual {hmax:.2e}")
    lines.append(f"- cooling channel closure: max rel. residual {cmax:.2e}")

    # 4. dominant channels and He share
    r = heat[:, 0]
    doms = {}
    for lab, j in (('base', 1), ('1.5 R_p', int(np.argmin(np.abs(r-1.5)))),
                   ('top', len(r)-1)):
        hch, hfr = dominant(HEAT_COLS, heat[j, 4:4+len(HEAT_COLS)])
        cvals = cool[j, 4:]
        ccols = COOL_COLS + [f"metal{k}" for k in range(cool.shape[1]-4-len(COOL_COLS))]
        cch, cfr = dominant(ccols, cvals)
        he_heat = heat[j, [5, 6, 7, 12, 13, 14]].sum()/max(heat[j, 3], 1e-300)
        doms[lab] = (hch, hfr, he_heat, cch, cfr)
        lines.append(f"- {lab} (r = {r[j]:.3f}, T = {heat[j,1]:.0f} K): heating led by {hch}"
                     f" ({100*hfr:.0f}%), He-channel share {100*he_heat:.0f}%;"
                     f" cooling led by {cch} ({100*cfr:.0f}%)")

    # 5. Mdot recomputed from the profile, mass-flux spread, JFNK exit state
    run = exhale_io.load_run(out, os.path.join(case, 'input.inp'), adv=False)
    mdot = exhale_io.mdot_log10(run)
    flux = run.n*run.v*run.r**2
    win = (run.r > 2.0) & (run.v > 0)
    du = (flux[win].max() - flux[win].min())/np.abs(np.median(flux[win])) \
        if win.sum() > 2 else np.nan
    info, res = jfnk_state(os.path.join(case, 'run.log'))
    lines.append(f"- log10 Mdot [g/s] = {mdot:.3f}; mass-flux spread above 2 R_p"
                 f" = {du:.2e}; JFNK exit info={info} ||R||={res}")

    # 6. guard margins
    nh = ions['HI'] + ions['HII']
    nhe = ions['HeI'] + ions['HeII'] + ions['HeIII']
    vals = {'n_H + n_He': nh + nhe, 'n_H': nh, 'n_HI': ions['HI'],
            'n_e': ne_code}
    nzero = int((ions['HI'] <= 0.0).sum())
    if nzero:
        rz = ions['HI'] <= 0.0
        lines.append(f"- n_HI reaches exactly 0 in {nzero} of {len(rz)} cells"
                     f" (r >= {r[rz].min():.3f}): hydrogen is fully ionized there,"
                     " an exact zero rather than a clamped floor")
    gl = []
    for lab, g, _src in GUARDS:
        mn = float(np.min(vals[lab]))
        gl.append(f"{lab} min {mn:.3e} ({mn/g:.1e}x guard {g:.0e})")
    lines.append("- guard margins: " + "; ".join(gl))

    summ = dict(case=name, nbad=nbad, neut=rel.max(), heat=hmax, cool=cmax,
                mdot=mdot, du=du, info=info, res=res, doms=doms,
                guards={lab: float(np.min(vals[lab])) for lab, _g, _s in GUARDS})
    return lines, summ


def main(argv):
    out_path = None
    args = list(argv)
    if args and args[0] in ('-o', '--out'):
        out_path = args[1]; args = args[2:]
    if args:
        cases = args
        if out_path is None:
            out_path = os.path.join(
                os.path.dirname(os.path.abspath(os.path.normpath(cases[0]))),
                'audit_summary.md')
    else:
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        cases = CASES
        if out_path is None:
            out_path = 'audit_summary.md'

    report = ["# He/H audit summary (Phase C)", "",
              "Generated by audit_checks.py on the cases: "
              + ", ".join(os.path.basename(os.path.normpath(c)) for c in cases)
              + ".", ""]
    rows = []
    for c in cases:
        lines, summ = audit(c)
        report += lines + [""]
        if summ:
            rows.append(summ)

    if rows:
        report += ["## Summary table", "",
                   "| case | non-finite | max charge-neutrality rel. err. |"
                   " max heating closure | max cooling closure | log10 Mdot |"
                   " mass-flux spread (r > 2) | JFNK info |",
                   "|---|---|---|---|---|---|---|---|"]
        for s in rows:
            report.append(f"| {s['case']} | {s['nbad']} | {s['neut']:.2e} |"
                          f" {s['heat']:.2e} | {s['cool']:.2e} |"
                          f" {s['mdot']:.3f} | {s['du']:.2e} | {s['info']} |")
        report += ["", "## Guard margins", "",
                   "Numerical guards on the H/He path (metals and molecules "
                   "off) are division denominators, not density floors: the "
                   "solution is never clamped onto them. Minimum of each "
                   "density over the domain, and the guard it is compared "
                   "with:", "",
                   "| case | " + " | ".join(g[0] + " min" for g in GUARDS) + " |",
                   "|---|" + "---|"*len(GUARDS)]
        for s in rows:
            report.append("| " + s['case'] + " | " + " | ".join(
                f"{s['guards'][g[0]]:.3e} ({s['guards'][g[0]]/g[1]:.1e}x)"
                for g in GUARDS) + " |")
        report += ["", "Guard sources:", ""]
        for lab, g, src in GUARDS:
            report.append(f"- `{lab}` vs {g:.0e}: {src}")
        report.append("")

    with open(out_path, 'w') as fh:
        fh.write("\n".join(report))
    print("\n".join(report))
    print(f"\n[written to {out_path}]")


if __name__ == '__main__':
    main(sys.argv[1:])
