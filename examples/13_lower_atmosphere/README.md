# 13_lower_atmosphere — Tier-1/3 lower-atmosphere connection examples

Four planets (HD 209458 b, HD 189733 b, WASP-121 b, WASP-52 b), each subfolder
holding the planet's `input.inp` plus two driver-generated handoff files:

- `base_iso.inp`     — isothermal-Teq analytic column (Koskinen+2022 style)
- `base_guillot.inp` — Guillot (2010) semi-grey T(p) column

Regenerate (from the EXHALE root; `--r1bar` = broadband transit radius in R_J):

```bash
python3 src/utils/run_lower.py examples/13_lower_atmosphere/HD209458b --r1bar 1.36
python3 src/utils/run_lower.py examples/13_lower_atmosphere/HD189733b --r1bar 1.138
python3 src/utils/run_lower.py examples/13_lower_atmosphere/WASP-121b --r1bar 1.753
python3 src/utils/run_lower.py examples/13_lower_atmosphere/WASP-52b  --r1bar 1.27
# add --guillot for the semi-grey T(p) variant
```

To USE a handoff: copy the chosen file to `base.inp` in the run directory —
EXHALE reads it on startup and overrides T0 / base radius / He/H / K_zz
(echoed to stdout; absent file = no-op). The in-code consistency report is
the one-line alternative: add `Lower column: <R_1bar>` to `input.inp`.

Results, interpretation, and the full description:
`docs/lower_atmosphere_coupling.pdf` (§Examples).

Headline numbers (isothermal column; r in R_J):

| planet | input R0 | r0(1 μbar) eq | atomic bracket | q_H2(base) | verdict |
|---|---|---|---|---|---|
| HD 209458 b | 1.401 | 1.472 | 1.582 | 0.83 | input ~5% below bracket |
| HD 189733 b | 1.193 | 1.174 | 1.205 | 0.84 | input INSIDE bracket ✓ |
| WASP-121 b  | 2.2075 | 2.010 | 2.133 | 0.005 | atomic base VERIFIED (UHJ) |
| WASP-52 b   | 1.270 | 1.437 | 1.615 | 0.84 | base too deep; recomputed: Mdot ×1.5, He 10830 ~unchanged (see docs) |
