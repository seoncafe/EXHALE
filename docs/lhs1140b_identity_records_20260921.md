# LHS 1140 b: Stage A0 identity records of five checkpoints

2026-09-21. Written from `models/identity_record.py`, which is new with
this document and reads only what the tree already stores.

## 1. What this is, and what it is not

Stage A0 of `PLAN_20260920_rev9.md` section 4 asks that the IDENTITY of a
checkpoint be recorded before anything is measured on it. This document is
that record for the five checkpoints the plan works on.

**Classification, rev9 section 3: every record below is a HISTORICAL
AUDIT.** It reads stored manifests, certificates, state files and logs, so
it carries the HISTORICAL executable identity and says nothing about the
current operator. No binary was run, no residual was evaluated and no state
was loaded, so no part of this document is Mode R. The only numbers the
tool produces itself are md5 sums of files on disk, labeled MEASURED;
everything a stored record states is labeled READ.

The command behind every section below, run in `LHS1140b/models/`:

```bash
models/identity_record.py <case dir> [--generation <id>] [--json <path>]
```

It writes nothing into a case directory and refuses a `--json` path inside
the case it is reading. It does not need the binary to exist: where the
executable is gone, the md5 the manifest recorded is READ and stands as the
identity, and the record says that it was not measured here.

## 2. The two tables of rev9 section 4.1, rebuilt from these records

A stored state and an iterate that exists only in a log never share a row.
The two tables below are kept apart for that reason and for no other.

### 2.1 The stored states the indexes select

Each row is one generation: the pair of files whose md5 the tool
recomputed and matched against the manifest. Generation, binary, boundary,
phase, ending and refusing rows are READ from
`states/<generation>/manifest.json` and `states/<generation>/certification.txt`.

| case | generation | binary | boundary | phase, ending | refusing rows (READ) |
|---|---|---|---|---|---|
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | `g0002_20260919T004843Z_f7485b14` | `7670f310` | `_v2` | evaluate, `evaluate` | hydrodynamic mass row: row measure 9.806E-01 above 1.9E-08 at cell 2; hydrodynamic momentum row: row measure 8.724E-06 above 1.0E-08 at cell 1; hydrodynamic energy row: row measure 1.071E+00 above 1.0E-06 at cell 1 |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | `g0004_20260919T094145Z_259fe9c3` | `2c3b0acc` | `_v2` | stationary alternation, `hydrodynamic_refusal` | hydrodynamic mass row: row measure 5.880E-01 above 1.2E-07 at cell 11; hydrodynamic momentum row: row measure 7.774E-07 above 1.0E-08 at cell 495; hydrodynamic energy row: row measure 1.949E+00 above 1.0E-06 at cell 7; elemental transport He/H partition: gated row measure 3.538E-02 above 1.0E-05 at cell 389 (a wind cell) |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | `g0002_20260919T004834Z_d59ec9ff` | `7670f310` | `_v2` | evaluate, `evaluate` | hydrodynamic mass row: row measure 9.589E-01 above 3.7E-08 at cell 2; hydrodynamic momentum row: row measure 5.539E-06 above 1.0E-08 at cell 1; hydrodynamic energy row: row measure 1.112E+00 above 1.0E-06 at cell 1 |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | `g0004_20260920T053734Z_79b42a03` | `75d55d9d` | `_v3` | evaluate, `evaluate` | carrier balance H2: gated row measure 8.150E-04 above 1.0E-05 at cell 499 (a wind cell); elemental transport He/H partition: gated row measure 2.503E-04 above 1.0E-05 at cell 280 (a wind cell) |
| `molecular_photochem_gj1132_kzzprofile/HeH9` | `g0005_20260920T053700Z_f59de41d` | `75d55d9d` | `_v3` | evaluate, `evaluate` | hydrodynamic mass row: row measure 1.041E+00 above 3.4E-10 at cell 2; hydrodynamic momentum row: row measure 1.450E-03 above 1.0E-08 at cell 1; hydrodynamic energy row: row measure 1.021E+00 above 1.0E-06 at cell 1; carrier balance H2: gated row measure 1.000E+00 above 1.0E-05 at cell 432 (a wind cell); elemental transport He/H partition: gated row measure 1.494E-04 above 1.0E-05 at cell 217 (a wind cell) |

Every state pair in this table MEASURED equal to the md5 the manifest
recorded: there is no mismatch anywhere in the five.

### 2.2 The iterates that exist only in a log

An iterate belongs here when the last iterate of the case `run.log` stands
after the last certificate of a published generation in that same log, so
no state was written from it. The tool decides that by locating each
published generation's certificate in the log BY ITS TEXT, never by a file
name or a modification time.

| case | binary | boundary | last completed solve (READ from `run.log`) | the log's last line | ending class |
|---|---|---|---|---|---|
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | `75d55d9d` | `_v3` | (JFNK) done info=2 \|\|R\|\|= 9.978E-01 flux spread= 4.424E-12 non-monotone accepts=16 solve=4 | line 1835: (EXHALE_main) outer pass 4: hydro info=2, worst gated species row 9.77E-04 of 1.0E-05 at cell 354 (elemental transport He/H partition), mass 7.33E-08, momentum 8.25E-14, energy 9.98E-01, omega 0.500, trust 1.0E-02, 282.01 s | stopped_by_wall_ceiling |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | `75d55d9d` | `_v3` | (JFNK) done info=2 \|\|R\|\|= 1.163E+00 flux spread= 9.356E-03 non-monotone accepts=24 solve=4 | line 1740: (JFNK) it 12 \|\|R\|\|= 1.171E+00 \|\|Fs\|\|2= 5.63E-04 lam= 1.00E+00 dtau= 1.18E-03 gm= 1 worst r= 1.000 worst row: mass of cell 1 solve=5 | stopped_by_wall_ceiling |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | `75d55d9d` | `_v3` | (JFNK) done info=2 \|\|R\|\|= 2.458E-01 flux spread= 1.956E-02 non-monotone accepts=12 solve=5 | line 2267: (JFNK) it 26 \|\|R\|\|= 2.542E-01 \|\|Fs\|\|2= 1.97E-05 lam= 1.00E+00 dtau= 6.47E-04 gm= 1 worst r= 1.261 worst row: mass of cell 232 solve=6 | stopped_by_wall_ceiling |

The two molecular cases have NO row here. The last solve of each of their
case `run.log` files ended in a published generation:

- `molecular_scalar_gj1132_kzz1e9/HeH0.083`: the last certificate in `run.log` is that of `g0001_20260916T023329Z_32b34e57`, at lines 8769 to 8821, and the last iterate of the log stands at line 8682, before it.

- `molecular_photochem_gj1132_kzzprofile/HeH9`: the last certificate in `run.log` is that of `g0003_20260917T222518Z_c55563e3`, at lines 439 to 495, and the last iterate of the log stands at line 412, before it.

## 3. The records

### 3.1 `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7`

```
==============================================================================
STAGE A0 IDENTITY RECORD
==============================================================================
case            atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7
case directory  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7
generation      g0002_20260919T004843Z_f7485b14
selected by     the reference `latest_complete` of state_index.json
latest_complete g0002_20260919T004843Z_f7485b14   (state_index.json)
latest_certified (none)  (state_index.json)
                the index names NO certified generation: no state
                of this case carries a certificate that covers its
                own bytes.
assembled       2026-09-21T08:14:57 by models/identity_record.py
classification  HISTORICAL AUDIT (PLAN_20260920_rev9 section 3). Every
                stored value below is READ; every md5 this tool
                recomputed is MEASURED. No binary was run and no
                state was evaluated, so nothing here is Mode R.

1. THE STATE PAIR
-----------------
directory       /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/states/g0002_20260919T004843Z_f7485b14

  Hydro_ioniz.txt
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/states/g0002_20260919T004843Z_f7485b14/Hydro_ioniz.txt
    bytes             93749
    md5 MEASURED      25e11562ec3cd1f6d204bf36f8e503eb
    md5 READ          25e11562ec3cd1f6d204bf36f8e503eb   [states/g0002_20260919T004843Z_f7485b14/manifest.json `components`]
    verdict           matches the recorded md5

  Ion_species.txt
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/states/g0002_20260919T004843Z_f7485b14/Ion_species.txt
    bytes             447309
    md5 MEASURED      78ac49d7e1ce90a3a4f271785400efda
    md5 READ          78ac49d7e1ce90a3a4f271785400efda   [states/g0002_20260919T004843Z_f7485b14/manifest.json `components`]
    verdict           matches the recorded md5

  the header of Hydro_ioniz.txt, which the md5 covers:
    # EXHALE schema 2
    # columns r[Rp] rho[mH/cm3] v[cm/s] p[cgs] T[K] heat[erg/cm3/s] cool[erg/cm3/s]
    # rows 504: 2 ghost cells at each end; physical cells are rows 3 to 502
    # coupling: sec_ion=T sec_ion_step=0 recon=WENO3 certified=F cert_reason=failing_entries mode=init
    # provenance: git=3c73905ca8a7 tree=dirty run=2026-09-19T09:48:43
    # provenance: ck_input=4027990851733187225 ck_base=-1 ck_metals=-1
    # provenance: recon=WENO3 base_bc=characteristic carrier=none restart_schema=3 resid_def=145 N=500
    # boundary_model characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
    # boundary_reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp] 1.0000000001000000E+000 1.0000000000000000E+000 2.7072294627320281E-001 1.0000000000000000E+000
    # restart_schema 1
    # reservoir He/H 9.6999999999999993E+00
    # species_columns 34 r HI HII HeI HeII HeIII HeITR CI CII CIII OI OII OIII NI NII NIII MgI MgII MgIII SiI SiII SiIII CaI CaII CaIII NaI NaII KI KII SI SII FeI FeII FeIII
    # grid N 500 R0[cm] 1.1273716464000001E+09 r_min[Rp] 1.0001933187778382E+00 r_max[Rp] 2.9031224061195065E+01 mode Mixed
    # constants set IAU2015+CODATA2018 RJ[cm] 7.1492000000000000E+09 kB[erg/K] 1.3806490000000000E-16
    # options He23S=T metals=F eos_metals=T mol=F molbase=F oxychem=F carrier=F carrier_newton=F iontrans=F he_diff=T he_metal_diff=F sec_ion=T caloric_mono=F excH=F base_ir=F mol_ir=F mol_heat=T visc=F cond=F jlya=0 wellbal=T
    # t_phys[s] 0.0000000000000000E+00
    # source git=3c73905ca8a7 tree=dirty restart_input=provenance_unknown

2. ORIGIN
---------
generation id       g0002_20260919T004843Z_f7485b14   [states/<gen>/manifest.json]
iteration_phase     evaluate   [manifest.json]
published_at        2026-09-19T09:48:44   [manifest.json]
state_written       2026-09-19T09:48:43   [manifest.json]
published_by        models/publish_state.py on lart4   [manifest.json]
source_directory    runs/r20260919T004843Z_15406/eval/output   [manifest.json]
run_id              (none)   [manifest.json]
certification       NOT CERTIFIED   [manifest.json `certification.status`]
state claim         certified=False cert_reason=failing_entries mode=init   [manifest]

parent              g0001_20260917T172820Z_3a57917b
  source            (none)
  compatibility     stated by the run that took it

seed                [provenance/g0002_20260919T004843Z_f7485b14/provenance_recovered*.json (RECOVERED, not recorded by the run)]
  established       True
  kind              a published generation
  case              atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7
  generation        g0001_20260917T172820Z_3a57917b
  directory         (none)
  chosen by         the manifest of this generation names a parent generation of this same case, which is the state the pass was started from (MODELS.md section 9.2)
  transformation    (none)
  confidence        established
  evidence          atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/states/g0002_20260919T004843Z_f7485b14/manifest.json
  note              generation g0001_20260917T172820Z_3a57917b of this same case

ending class        evaluate   [manifest.json `ending.class`]
ending reason       the certified generation g0001_20260917T172820Z_3a57917b handed back and measured as it stands under the D9 step 3 catalog binary; an evaluation has no solve ending
index ending_class  evaluate   [state_index.json]
solver verdict      info = (none), ||R|| = (none)   [manifest.json `ending.solver`]

the case `ENDING` file (/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/ENDING, written 2026-09-20T07:32:35):
    stopped_by_wall_ceiling: stopped by the wall ceiling of 30m that models/run_campaign.sh applies (SIGTERM 1800 s after the campaign started the case), in the wind pass after outer pass 4 had completed; that pass had written no complete state of its own, so nothing new is published
    NOTE: this file is the ending of the LAST run of the case
    directory. It is not in general the ending of the run that
    published this generation, which is the `ending` block above.

3. THE OPERATOR
---------------
binary path         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_7670f310.x   [manifest.json `source_identity`]
  on disk now       yes
  md5 READ          7670f31031fb4db91d27b44cb0da6f70   [manifest.json]
  md5 MEASURED      7670f31031fb4db91d27b44cb0da6f70
  source manifest   BINARY_MANIFEST_7670f31031fb.txt
  its md5 READ      0376f5510808d64d32af864a9cf55742
  its md5 MEASURED  0376f5510808d64d32af864a9cf55742

model options       He23S=T metals=F eos_metals=T mol=F molbase=F oxychem=F carrier=F carrier_newton=F iontrans=F he_diff=T he_metal_diff=F sec_ion=T caloric_mono=F excH=F base_ir=F mol_ir=F mol_heat=T visc=F cond=F jlya=0 wellbal=T   [manifest.json `model_identity.options`]
reservoir           He/H 9.6999999999999993E+00

boundary model of the STORED state:
  state header      characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
                    [# boundary_model of Hydro_ioniz.txt, part of the hashed bytes]
  manifest          characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
  boundary reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp] 1.0000000001000000E+000 1.0000000000000000E+000 2.7072294627320281E-001 1.0000000000000000E+000

boundary model a log of this case names for its own restart:
  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/run.log:20
    (load_IC) boundary model of the restart: characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3 (this run's).

  *** THE TWO ARE DIFFERENT EXPERIMENTS.
  *** the stored state was written under
  ***   characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
  *** a log of this case ran under
  ***   characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3
  *** A run that rebuilds the boundary from the physical column
  *** and its own reservoir is the CURRENT operator on an OLD
  *** state (rev9 section 3); its residual is not the residual
  *** the stored state was written with, and the two may not
  *** share a row of any table.

4. THE CONFIGURATION
--------------------
  input.inp                      matches the recorded md5
    md5 MEASURED      e43c6697786711fb916d3b8eb3c113aa
    md5 READ          e43c6697786711fb916d3b8eb3c113aa   [manifest.json `configuration_identity`]
  base.inp                       absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  metals.inp                     absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  opacity.inp                    absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  lower_atmosphere_profile.dat   absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]

  spectrum
    named in input    ../../../sed/lhs1140_sed_gj1132_at_b_xuv0p01.txt   [input.inp `Spectrum file:`]
    resolved to       /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/sed/lhs1140_sed_gj1132_at_b_xuv0p01.txt
    resolution rule   relative to the case directory, which is the working directory of a run: the binary reads ./input.inp and opens the string as written
    md5 MEASURED      7da5fa41a0706feacac74b7d245fbec3
    md5 READ          7da5fa41a0706feacac74b7d245fbec3   [manifest.json]
    verdict           matches the recorded md5

  EXHALE_resolved.out
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/EXHALE_resolved.out
    md5 MEASURED      a24d52aa7971b8d1366a9a1591d716e9
    written           2026-09-20T07:02:40
    keys              38
    route keys it DOES carry:
      well_balanced              T
      carrier_transport          F
      carrier_in_newton          F
      carrier_newton_on_stall    F
      ionization_transport       F
      oxygen_chemistry           F
    route keys it does NOT carry, taken from elsewhere:
      - the boundary model identity (the state file header and the manifest carry it; the resolved writer never emits it, write_setup_report.f90:906-1120)
      - the numerical flux and the reconstruction method (stated in prose by EXHALE_setup.out)
      - the residual assembly selector EXHALE_RESID_QUAD, which is an environment variable read at hydrodynamic_rows.f90:177 and is in no resolved record

  EXHALE_setup.out (the resolved route, in prose)
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/EXHALE_setup.out
    md5 MEASURED      24cad0f8914a95509252ab99b6057327
    written           2026-09-20T07:02:40
    numerical flux    ROE   [line 46]
    reconstruction method PLM   [line 47]
    well balanced     on   [line 48]
    base boundary     characteristic condition at the face r_edg(0)   [line 53]

    NOTE: EXHALE_resolved.out and EXHALE_setup.out of this case
    were written AFTER this generation was published, so they are
    the record of a LATER run in the same directory. The input
    md5 comparison above is what ties the configuration to this
    generation; the route lines are read from the later report
    and are the route of that run.

5. THE ROW REGISTRY
-------------------
certificate         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/states/g0002_20260919T004843Z_f7485b14/certification.txt
  md5 MEASURED      f800a42985aee3831be0a9b399df31c1
  active equations 6 of 30 in the inventory   [line 2]

  rows evaluated (6):
    hydrodynamic mass row                         ABOVE  max= 9.806E-01 at cell 2 (r = 1.00039), tol= 1.9E-08   [line 14]
    hydrodynamic momentum row                     ABOVE  max= 8.724E-06 at cell 1 (r = 1.00019), tol= 1.0E-08   [line 17]
    hydrodynamic energy row                       ABOVE  max= 1.071E+00 at cell 1 (r = 1.00019), tol= 1.0E-06   [line 19]
    elemental transport He/H partition            within max= 9.267E-01 at cell 2 (r = 1.00039), tol= 1.0E-05   [line 21]
    level balance He 2^3S                         within max= 7.201E-21 at cell 482 (r = 21.48368), tol= 1.0E-06   [line 24]
    eliminated-species closure System_HeH_TR      within max= 5.623E-17 at cell 482 (r = 21.48368), tol= 1.0E-06   [line 26]

  verdict           NOT CERTIFIED: 3 entry/entries of the inventory refuse it   [line 38]

  refusing rows (3):

    hydrodynamic mass row: row measure  9.806E-01 above  1.9E-08 at cell 2
      row max         9.806E-01 at cell 2 (r = 1.00039, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.9E-08
      volume measure  1.846E-01
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      verdict cell    2 (r = 1.00039E+00 READ, 1.00039 MEASURED): 9.806E-01 against 1.9E-08, distance 5.049E+07
                      [line 16]

    hydrodynamic momentum row: row measure  8.724E-06 above  1.0E-08 at cell 1
      row max         8.724E-06 at cell 1 (r = 1.00019, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-08
      volume measure  8.185E-07
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      the certificate prints no tolerance-normalized verdict
      cell for this row (see part 7)

    hydrodynamic energy row: row measure  1.071E+00 above  1.0E-06 at cell 1
      row max         1.071E+00 at cell 1 (r = 1.00019, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-06
      volume measure  3.199E-01
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      the certificate prints no tolerance-normalized verdict
      cell for this row (see part 7)

6. THE LOG BLOCK
----------------
the log the manifest names for this generation
  path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/runs/r20260919T004843Z_15406/pp.log
  named in          states/g0002_20260919T004843Z_f7485b14/manifest.json `certification.source`
  on disk now       yes
  md5 MEASURED      63d5fcc87727714f0d39bb3d7ab156a5
  md5 READ          (none)
  written           2026-09-19T09:48:44

where the certificate of this generation actually stands, located by
its own text and not by a file name:
  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/runs/r20260919T004843Z_15406/pp.log
  lines 169 to 213 (45 lines, the whole of certification.txt)

  the block:
     169 |  (certification) FINAL STATE AS WRITTEN -- this run is NOT a certified stationary solution
     170 |    active equations 6 of 30 in the inventory
     171 |    tolerances (convergence study, contract sections 9 and 10; the hydrodynamic values were
     172 |    re-anchored on the self-consistent JFNK residual after B5a):
     173 |      hydrodynamic mass  3.00E-12, momentum  1.00E-08, energy  1.00E-06
     174 |      the mass value is a FLOOR: where the cell's own rounding of the flux difference stands above it, 10.0 times that floor is the tolerance
     175 |      of that cell (cert_tol_mass_at), and the verdict on the mass row is taken cell by cell
     176 |      closure  1.00E-06, level  1.00E-06; and IN THE WIND, at r >= 1.20E+00, carrier  1.00E-05
     177 |      and elemental transport  1.00E-05; the species rows of the cells below that radius are
     178 |      REPORTED AND DO NOT GATE: the layer r < 1.10E+00, whose element fluxes
     179 |      are not conserved (spread 9.8), and the band above it, whose discretization (4.6e-4 at the
     180 |      binding radius) stands above any tolerance one could set there
     181 |      the run's own "Resid tol" is  5.00E-05, which is the SOLVER's target and not a certification tolerance
     182 |    hydrodynamic mass row                                evaluated       max= 9.806E-01  vol= 1.846E-01  cell=2  tol= 1.9E-08  ABOVE   
     183 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     184 |         verdict at cell 2 (r= 1.00039E+00):  9.806E-01 against  1.9E-08, distance  5.049E+07 -- the rounding anchor of that cell
     185 |    hydrodynamic momentum row                            evaluated       max= 8.724E-06  vol= 8.185E-07  cell=1  tol= 1.0E-08  ABOVE   
     186 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     187 |    hydrodynamic energy row                              evaluated       max= 1.071E+00  vol= 3.199E-01  cell=1  tol= 1.0E-06  ABOVE   
     188 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     189 |    elemental transport He/H partition                   evaluated       max= 9.267E-01  vol= 4.892E-02  cell=2  tol= 1.0E-05  within  
     190 |         scale: sum of the row's own terms [g cm^-3 s^-1]; floor 1e-20 rho X_base/(R0/v0)
     191 |         gated at r >= 1.200E+00 against  1.00E-05:  6.723E-06 at cell 395; the cells below that radius are reported only
     192 |    level balance He 2^3S                                evaluated       max= 7.201E-21  vol= 2.046E-21  cell=482  tol= 1.0E-06  within  
     193 |         scale: the row's own turnover rate; dimensionless, scale 1
     194 |    eliminated-species closure System_HeH_TR             evaluated       max= 5.623E-17  vol= 1.813E-17  cell=482  tol= 1.0E-06  within  
     195 |         scale: the row's own turnover rate; dimensionless, scale 1
     196 |    cells without a chemical root: 0 (acceptance class 4 or 6)
     197 |    mass closure of the composition, max |sum_i f_i A_i - 1| =  5.578E-16 at cell 32, r =   1.00619
     198 |      reported only: the sum is one by the definition of the mass fractions, so this is the arithmetic of the state and no equation of it
     199 |    validity states (B1a section 4):
     200 |      active unvalidated physics: not produced
     201 |      out-of-domain closure activations: 0
     202 |      rejected trials with no adopted contribution: not produced
     203 |      unbudgeted accepted corrections: 0
     204 |      specified external reservoirs: not produced
     205 |    attempts (not validity states): 0 energy-update and 0 conduction cell(s) reached the lower bracket end
     206 |    NOT CERTIFIED: 3 entry/entries of the inventory refuse it
     207 |      hydrodynamic mass row: row measure  9.806E-01 above  1.9E-08 at cell 2
     208 |      hydrodynamic momentum row: row measure  8.724E-06 above  1.0E-08 at cell 1
     209 |      hydrodynamic energy row: row measure  1.071E+00 above  1.0E-06 at cell 1
     210 |  
     211 |  (certification) the state was written in full; the run exits with status 2 because it is not certified.
     212 | Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
     213 | STOP 2

run.log of this case
  path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/run.log
  lines             1835
  md5 MEASURED      1a0b0f2037caa769673b29bda5073f3a
  written           2026-09-20T07:32:30
  runs in the file  1 (a run begins at line 2)
  is the log the manifest names     no
  holds this generation's block     NO
  boundary model of its restart     characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3   [line 20]

  published generations whose certificate stands in this log,
  located by the certificate text (a pass certification of a
  state the run then moved away from is not one):
    NONE. No pass recorded in this log ended in a state that
    was published as a generation of this case.

  the binary that wrote this log:
    `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_75d55d9d.x`, md5 `75d55d9d4fd0e748cd01d6e35713e34c`
    [/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/not_solved.md, written by models/run_case.sh for the run that produced this log]

  the last iterate this log records:
    the last outer pass   [line 1835]
      (EXHALE_main) outer pass 4: hydro info=2, worst gated species row  9.77E-04 of  1.0E-05 at cell 354 (elemental transport He/H partition), mass  7.33E-08, momentum  8.25E-14, energy  9.98E-01, omega 0.500, trust 1.0E-02,  282.01 s
    the last inner solve verdict   [line 1807]
      (JFNK) done info=2 ||R||=  9.978E-01  flux spread= 4.424E-12  non-monotone accepts=16  solve=4
    the last residual gate   [line 1834]
      (JFNK) gate NOT met: residual ||R|| = 9.978E-01 >= 5.000E-05

  the last lines of the file:
    1833 |  (JFNK) composition elimination: 7285 sweeps over 1597 residual evaluations,   4.56 per evaluation; the Newton model carried 6 Picard term(s)
    1834 |  (JFNK) gate NOT met: residual ||R|| = 9.978E-01 >= 5.000E-05
    1835 |  (EXHALE_main) outer pass 4: hydro info=2, worst gated species row  9.77E-04 of  1.0E-05 at cell 354 (elemental transport He/H partition), mass  7.33E-08, momentum  8.25E-14, energy  9.98E-01, omega 0.500, trust 1.0E-02,  282.01 s

the identity of the run that wrote the case `run.log`
  [/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/not_solved.md, written 2026-09-20T07:32:35]
  A log never names its own executable, so this is where the
  identity of a run that published nothing is kept. It belongs to
  the LOG-ONLY table of rev9 section 4.1 and never to a row that
  also carries a stored state.
    case               `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7`
    ending class       `stopped_by_wall_ceiling`
    solver verdict     info = none printed
    outer passes       4 of the 40 this run allowed
    pseudo-time start  1.0
    binary             `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_75d55d9d.x`, md5 `75d55d9d4fd0e748cd01d6e35713e34c`
    host, threads      `lart4`, OMP_NUM_THREADS=8
    seed               /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7/states/g0004_20260919T212852Z_de543034  tier3:models/atomic_scalar_gj1132x0.10_kzz1e9  HeH=9.7  target=9.7  dlog10=0.0000  candidates=101  (a certified case at another XUV normalization)  compat=same system (scalar, read from the state)

the publication test, MEASURED here
  Hydro_ioniz.txt
    output/ md5       25e11562ec3cd1f6d204bf36f8e503eb   (written 2026-09-19T09:48:43)
    generation md5    25e11562ec3cd1f6d204bf36f8e503eb
  Ion_species.txt
    output/ md5       78ac49d7e1ce90a3a4f271785400efda   (written 2026-09-19T09:48:43)
    generation md5    78ac49d7e1ce90a3a4f271785400efda
  verdict             the pair in `output/` is BITWISE the stored generation

  *** The last iterate of the case `run.log` (line 1835) stands AFTER every
  *** line of that log, none of which is the certificate of a published
  *** generation of this case, so no state was written from it. The pair in
  *** `output/` is bitwise the stored generation
  *** g0002_20260919T004843Z_f7485b14, which confirms that nothing was
  *** published after it. It exists ONLY in that log: it cannot be reloaded,
  *** it has no generation, no manifest and no certificate of its own, and it
  *** may not share a row of any table with a stored state (rev9 section 4.1).

7. FIELDS THIS RECORD COULD NOT FILL
------------------------------------
 1. the tolerance-normalized verdict of the row 'hydrodynamic momentum row'
    the certificate states the row measure and its tolerance but prints
    neither a `verdict at cell` line nor a `gated at r >=` line for it
 2. the tolerance-normalized verdict of the row 'hydrodynamic energy row'
    the certificate states the row measure and its tolerance but prints
    neither a `verdict at cell` line nor a `gated at r >=` line for it
 3. a recorded md5 for the log /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/runs/r20260919T004843Z_15406/pp.log
    the manifest records no `source_md5`, so the file on disk cannot be
    shown to be the file the certificate was read from
```

### 3.2 `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7`

```
==============================================================================
STAGE A0 IDENTITY RECORD
==============================================================================
case            atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7
case directory  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7
generation      g0004_20260919T094145Z_259fe9c3
selected by     the reference `latest_complete` of state_index.json
latest_complete g0004_20260919T094145Z_259fe9c3   (state_index.json)
latest_certified (none)  (state_index.json)
                the index names NO certified generation: no state
                of this case carries a certificate that covers its
                own bytes.
assembled       2026-09-21T08:14:57 by models/identity_record.py
classification  HISTORICAL AUDIT (PLAN_20260920_rev9 section 3). Every
                stored value below is READ; every md5 this tool
                recomputed is MEASURED. No binary was run and no
                state was evaluated, so nothing here is Mode R.

1. THE STATE PAIR
-----------------
directory       /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/states/g0004_20260919T094145Z_259fe9c3

  Hydro_ioniz.txt
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/states/g0004_20260919T094145Z_259fe9c3/Hydro_ioniz.txt
    bytes             93834
    md5 MEASURED      c19b3df722fb19c91662c633bc5f1bc1
    md5 READ          c19b3df722fb19c91662c633bc5f1bc1   [states/g0004_20260919T094145Z_259fe9c3/manifest.json `components`]
    verdict           matches the recorded md5

  Ion_species.txt
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/states/g0004_20260919T094145Z_259fe9c3/Ion_species.txt
    bytes             447390
    md5 MEASURED      9f04c20d6c365f87067c081ec2054b53
    md5 READ          9f04c20d6c365f87067c081ec2054b53   [states/g0004_20260919T094145Z_259fe9c3/manifest.json `components`]
    verdict           matches the recorded md5

  the header of Hydro_ioniz.txt, which the md5 covers:
    # EXHALE schema 2
    # columns r[Rp] rho[mH/cm3] v[cm/s] p[cgs] T[K] heat[erg/cm3/s] cool[erg/cm3/s]
    # rows 504: 2 ghost cells at each end; physical cells are rows 3 to 502
    # coupling: sec_ion=T sec_ion_step=0 recon=WENO3 certified=F cert_reason=no_stationary_claim mode=init
    # provenance: git=3c73905ca8a7 tree=dirty run=2026-09-19T18:41:45
    # provenance: ck_input=2729954298256964254 ck_base=-1 ck_metals=-1
    # provenance: recon=WENO3 base_bc=characteristic carrier=none restart_schema=3 resid_def=145 N=500
    # boundary_model characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
    # boundary_reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp] 1.0000336674112269E+000 1.0000000000000000E+000 2.7068048521511734E-001 1.0000000000000000E+000
    # restart_schema 1
    # reservoir He/H 9.7110732820766632E+00 C/H 2.7780271379679837E-04 O/H 9.0277252501755680E-07 N/H 8.1907551439822828E-05
    # species_columns 34 r HI HII HeI HeII HeIII HeITR CI CII CIII OI OII OIII NI NII NIII MgI MgII MgIII SiI SiII SiIII CaI CaII CaIII NaI NaII KI KII SI SII FeI FeII FeIII
    # grid N 500 R0[cm] 1.1592997810038679E+09 r_min[Rp] 1.0001933187778382E+00 r_max[Rp] 2.9031224061195065E+01 mode Mixed
    # constants set IAU2015+CODATA2018 RJ[cm] 7.1492000000000000E+09 kB[erg/K] 1.3806490000000000E-16
    # options He23S=T metals=T eos_metals=T mol=F molbase=F oxychem=F carrier=F carrier_newton=F iontrans=F he_diff=T he_metal_diff=F sec_ion=T caloric_mono=F excH=F base_ir=F mol_ir=F mol_heat=T visc=F cond=F jlya=0 wellbal=T
    # t_phys[s] 0.0000000000000000E+00
    # source git=3c73905ca8a7 tree=dirty restart_input=provenance_unknown

2. ORIGIN
---------
generation id       g0004_20260919T094145Z_259fe9c3   [states/<gen>/manifest.json]
iteration_phase     stationary alternation   [manifest.json]
published_at        2026-09-19T18:41:50   [manifest.json]
state_written       2026-09-19T18:41:45   [manifest.json]
published_by        models/publish_state.py on lart4   [manifest.json]
source_directory    output   [manifest.json]
run_id              (none)   [manifest.json]
certification       NOT CERTIFIED   [manifest.json `certification.status`]
state claim         certified=False cert_reason=no_stationary_claim mode=init   [manifest]

parent              g0003_20260919T091147Z_279c76bd
  source            (none)
  compatibility     stated by the run that took it

seed                [provenance/g0004_20260919T094145Z_259fe9c3/provenance_recovered*.json (RECOVERED, not recorded by the run)]
  established       True
  kind              a published generation
  case              atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7
  generation        g0003_20260919T091147Z_279c76bd
  directory         (none)
  chosen by         the manifest of this generation names a parent generation of this same case, which is the state the pass was started from (MODELS.md section 9.2)
  transformation    (none)
  confidence        established
  evidence          atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/states/g0004_20260919T094145Z_259fe9c3/manifest.json
  note              generation g0003_20260919T091147Z_279c76bd of this same case

ending class        hydrodynamic_refusal   [manifest.json `ending.class`]
ending reason       the solve returned info=1 with at least one hydrodynamic row above its tolerance: hydrodynamic mass row: row measure  5.880E-01 above  1.2E-07 at cell 11
index ending_class  hydrodynamic_refusal   [state_index.json]
solver verdict      info = 1, ||R|| = 1.949   [manifest.json `ending.solver`]

the case `ENDING` file (/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/ENDING, written 2026-09-20T08:02:36):
    stopped_by_wall_ceiling: stopped by the wall ceiling of 30m that models/run_campaign.sh applies (SIGTERM 1800 s after the campaign started the case), in the wind pass after outer pass 4 had completed; that pass had written no complete state of its own, so nothing new is published
    NOTE: this file is the ending of the LAST run of the case
    directory. It is not in general the ending of the run that
    published this generation, which is the `ending` block above.

3. THE OPERATOR
---------------
binary path         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_2c3b0acc.x   [manifest.json `source_identity`]
  on disk now       yes
  md5 READ          2c3b0acc9983aec03bb4f844294fed18   [manifest.json]
  md5 MEASURED      2c3b0acc9983aec03bb4f844294fed18
  source manifest   BINARY_MANIFEST_2c3b0acc9983.txt
  its md5 READ      54809f71fa208097c540c334c111f3e6
  its md5 MEASURED  54809f71fa208097c540c334c111f3e6

model options       He23S=T metals=T eos_metals=T mol=F molbase=F oxychem=F carrier=F carrier_newton=F iontrans=F he_diff=T he_metal_diff=F sec_ion=T caloric_mono=F excH=F base_ir=F mol_ir=F mol_heat=T visc=F cond=F jlya=0 wellbal=T   [manifest.json `model_identity.options`]
reservoir           He/H 9.7110732820766632E+00 C/H 2.7780271379679837E-04 O/H 9.0277252501755680E-07 N/H 8.1907551439822828E-05

boundary model of the STORED state:
  state header      characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
                    [# boundary_model of Hydro_ioniz.txt, part of the hashed bytes]
  manifest          characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
  boundary reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp] 1.0000336674112269E+000 1.0000000000000000E+000 2.7068048521511734E-001 1.0000000000000000E+000

boundary model a log of this case names for its own restart:
  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/run.log:32
    (load_IC) boundary model of the restart: characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3 (this run's).
  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/attempt_1/run.log:32
    (load_IC) boundary model of the restart: characteristic_face_ps_reservoir_C_minus_contact_upwind_v2 (this run's).
  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/attempt_1/solve_1/run.log:32
    (load_IC) boundary model of the restart: characteristic_face_ps_reservoir_C_minus_contact_upwind_v2 (this run's).

  *** THE TWO ARE DIFFERENT EXPERIMENTS.
  *** the stored state was written under
  ***   characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
  *** a log of this case ran under
  ***   characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3
  *** A run that rebuilds the boundary from the physical column
  *** and its own reservoir is the CURRENT operator on an OLD
  *** state (rev9 section 3); its residual is not the residual
  *** the stored state was written with, and the two may not
  *** share a row of any table.

4. THE CONFIGURATION
--------------------
  input.inp                      matches the recorded md5
    md5 MEASURED      6fc953139ef22e7e88397b465d997b64
    md5 READ          6fc953139ef22e7e88397b465d997b64   [manifest.json `configuration_identity`]
  base.inp                       absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  metals.inp                     absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  opacity.inp                    absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  lower_atmosphere_profile.dat   matches the recorded md5
    md5 MEASURED      3f96da22717e79667a5c93e932ded7f3
    md5 READ          3f96da22717e79667a5c93e932ded7f3   [manifest.json `configuration_identity`]

  spectrum
    named in input    ../../../sed/lhs1140_sed_gj1132_at_b_xuv0p01.txt   [input.inp `Spectrum file:`]
    resolved to       /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/sed/lhs1140_sed_gj1132_at_b_xuv0p01.txt
    resolution rule   relative to the case directory, which is the working directory of a run: the binary reads ./input.inp and opens the string as written
    md5 MEASURED      7da5fa41a0706feacac74b7d245fbec3
    md5 READ          7da5fa41a0706feacac74b7d245fbec3   [manifest.json]
    verdict           matches the recorded md5

  EXHALE_resolved.out
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/EXHALE_resolved.out
    md5 MEASURED      6d60d3febd20180d0d18c8b0f2811d95
    written           2026-09-20T07:32:40
    keys              55
    route keys it DOES carry:
      well_balanced              T
      carrier_transport          F
      carrier_in_newton          F
      carrier_newton_on_stall    F
      ionization_transport       F
      oxygen_chemistry           F
    route keys it does NOT carry, taken from elsewhere:
      - the boundary model identity (the state file header and the manifest carry it; the resolved writer never emits it, write_setup_report.f90:906-1120)
      - the numerical flux and the reconstruction method (stated in prose by EXHALE_setup.out)
      - the residual assembly selector EXHALE_RESID_QUAD, which is an environment variable read at hydrodynamic_rows.f90:177 and is in no resolved record

  EXHALE_setup.out (the resolved route, in prose)
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/EXHALE_setup.out
    md5 MEASURED      86b2a08eafa922734c82c05b78501a6b
    written           2026-09-20T07:32:40
    numerical flux    ROE   [line 47]
    reconstruction method PLM   [line 48]
    well balanced     on   [line 49]
    base boundary     characteristic condition at the face r_edg(0)   [line 54]

    NOTE: EXHALE_resolved.out and EXHALE_setup.out of this case
    were written AFTER this generation was published, so they are
    the record of a LATER run in the same directory. The input
    md5 comparison above is what ties the configuration to this
    generation; the route lines are read from the later report
    and are the route of that run.

5. THE ROW REGISTRY
-------------------
certificate         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/states/g0004_20260919T094145Z_259fe9c3/certification.txt
  md5 MEASURED      c5ed3d071f9be48ffccc3d15510f6589
  active equations 6 of 30 in the inventory   [line 2]

  rows evaluated (6):
    hydrodynamic mass row                         ABOVE  max= 6.925E-01 at cell 9 (r = 1.00174), tol= 1.2E-07   [line 14]
    hydrodynamic momentum row                     ABOVE  max= 7.774E-07 at cell 495 (r = 26.69350), tol= 1.0E-08   [line 18]
    hydrodynamic energy row                       ABOVE  max= 1.949E+00 at cell 7 (r = 1.00135), tol= 1.0E-06   [line 20]
    elemental transport He/H partition            ABOVE  max= 3.538E-02 at cell 389 (r = 5.04842), tol= 1.0E-05   [line 22]
    level balance He 2^3S                         within max= 9.338E-21 at cell 491 (r = 24.96303), tol= 1.0E-06   [line 25]
    eliminated-species closure System_HeH_TR_metals within max= 4.512E-16 at cell 480 (r = 20.78185), tol= 1.0E-06   [line 27]

  verdict           NOT CERTIFIED: 4 entry/entries of the inventory refuse it   [line 39]

  refusing rows (4):

    hydrodynamic mass row: row measure  5.880E-01 above  1.2E-07 at cell 11
      row max         6.925E-01 at cell 9 (r = 1.00174, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.2E-07
      volume measure  5.745E-02
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      verdict cell    11 (r = 1.00213E+00 READ, 1.00213 MEASURED): 5.880E-01 against 1.2E-07, distance 5.075E+06
                      [line 16]

    hydrodynamic momentum row: row measure  7.774E-07 above  1.0E-08 at cell 495
      row max         7.774E-07 at cell 495 (r = 26.69350, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-08
      volume measure  2.295E-08
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      the certificate prints no tolerance-normalized verdict
      cell for this row (see part 7)

    hydrodynamic energy row: row measure  1.949E+00 above  1.0E-06 at cell 7
      row max         1.949E+00 at cell 7 (r = 1.00135, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-06
      volume measure  3.179E-01
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      the certificate prints no tolerance-normalized verdict
      cell for this row (see part 7)

    elemental transport He/H partition: gated row measure  3.538E-02 above  1.0E-05 at cell 389 (a wind cell)
      row max         3.538E-02 at cell 389 (r = 5.04842, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-05
      volume measure  1.140E-02
      scale           sum of the row's own terms [g cm^-3 s^-1]; floor 1e-20 rho X_base/(R0/v0)
      gated at r >= 1.200E+00 against 1.00E-05: 3.538E-02 at cell 389 (r = 5.04842 MEASURED)
                      [line 24]

6. THE LOG BLOCK
----------------
the log the manifest names for this generation
  path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/run.log
  named in          states/g0004_20260919T094145Z_259fe9c3/manifest.json `certification.source`
  on disk now       yes
  md5 MEASURED      4a82d4fcc667ced98feb59fe0d74553c
  md5 READ          (none)
  written           2026-09-20T08:01:52

where the certificate of this generation actually stands, located by
its own text and not by a file name:
  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/attempt_1/run.log
  lines 2318 to 2368 (51 lines, the whole of certification.txt)

  *** this is NOT the log the manifest names. The named file
  *** has been written by a later run of this directory; the
  *** block below is the certificate of this generation and
  *** the named path is no longer where it stands.

  the block:
    2318 |  (certification) final state, as written
    2319 |    active equations 6 of 30 in the inventory
    2320 |    tolerances (convergence study, contract sections 9 and 10; the hydrodynamic values were
    2321 |    re-anchored on the self-consistent JFNK residual after B5a):
    2322 |      hydrodynamic mass  3.00E-12, momentum  1.00E-08, energy  1.00E-06
    2323 |      the mass value is a FLOOR: where the cell's own rounding of the flux difference stands above it, 10.0 times that floor is the tolerance
    2324 |      of that cell (cert_tol_mass_at), and the verdict on the mass row is taken cell by cell
    2325 |      closure  1.00E-06, level  1.00E-06; and IN THE WIND, at r >= 1.20E+00, carrier  1.00E-05
    2326 |      and elemental transport  1.00E-05; the species rows of the cells below that radius are
    2327 |      REPORTED AND DO NOT GATE: the layer r < 1.10E+00, whose element fluxes
    2328 |      are not conserved (spread 9.8), and the band above it, whose discretization (4.6e-4 at the
    2329 |      binding radius) stands above any tolerance one could set there
    2330 |      the run's own "Resid tol" is  4.00E-06, which is the SOLVER's target and not a certification tolerance
    2331 |    hydrodynamic mass row                                evaluated       max= 6.925E-01  vol= 5.745E-02  cell=9  tol= 1.2E-07  ABOVE   
    2332 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
    2333 |         verdict at cell 11 (r= 1.00213E+00):  5.880E-01 against  1.2E-07, distance  5.075E+06 -- the rounding anchor of that cell
    2334 |         the largest measure sits at cell 9, whose own tolerance is  1.5E-07 there
    2335 |    hydrodynamic momentum row                            evaluated       max= 7.774E-07  vol= 2.295E-08  cell=495  tol= 1.0E-08  ABOVE   
    2336 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
    2337 |    hydrodynamic energy row                              evaluated       max= 1.949E+00  vol= 3.179E-01  cell=7  tol= 1.0E-06  ABOVE   
    2338 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
    2339 |    elemental transport He/H partition                   evaluated       max= 3.538E-02  vol= 1.140E-02  cell=389  tol= 1.0E-05  ABOVE   
    2340 |         scale: sum of the row's own terms [g cm^-3 s^-1]; floor 1e-20 rho X_base/(R0/v0)
    2341 |         gated at r >= 1.200E+00 against  1.00E-05:  3.538E-02 at cell 389; the cells below that radius are reported only
    2342 |    level balance He 2^3S                                evaluated       max= 9.338E-21  vol= 2.289E-21  cell=491  tol= 1.0E-06  within  
    2343 |         scale: the row's own turnover rate; dimensionless, scale 1
    2344 |    eliminated-species closure System_HeH_TR_metals      evaluated       max= 4.512E-16  vol= 1.662E-16  cell=480  tol= 1.0E-06  within  
    2345 |         scale: the row's own turnover rate; dimensionless, scale 1
    2346 |    cells without a chemical root: 0 (acceptance class 4 or 6)
    2347 |    mass closure of the composition, max |sum_i f_i A_i - 1| =  6.692E-16 at cell 289, r =   1.70732
    2348 |      reported only: the sum is one by the definition of the mass fractions, so this is the arithmetic of the state and no equation of it
    2349 |    validity states (B1a section 4):
    2350 |      active unvalidated physics: not produced
    2351 |      out-of-domain closure activations: 0
    2352 |      rejected trials with no adopted contribution: not produced
    2353 |      unbudgeted accepted corrections: 0
    2354 |      specified external reservoirs: not produced
    2355 |    attempts (not validity states): 0 energy-update and 0 conduction cell(s) reached the lower bracket end
    2356 |    NOT CERTIFIED: 4 entry/entries of the inventory refuse it
    2357 |      hydrodynamic mass row: row measure  5.880E-01 above  1.2E-07 at cell 11
    2358 |      hydrodynamic momentum row: row measure  7.774E-07 above  1.0E-08 at cell 495
    2359 |      hydrodynamic energy row: row measure  1.949E+00 above  1.0E-06 at cell 7
    2360 |      elemental transport He/H partition: gated row measure  3.538E-02 above  1.0E-05 at cell 389 (a wind cell)
    2361 |  
    2362 |  (write_output/eq) max |sum(O I levels)/n(O I) - 1| =  2.93E-16
    2363 |  (EXHALE_main) flux gate: accepted 6.5390E-03   as written 6.5390E-03
    2364 |  (EXHALE_main) flux spread by window: r>=1.03 1.4156E-01   r>=1.10 1.8037E-02  (reporting only; the gate is the r>=r_flux window above)
    2365 |  (EXHALE_main) Restart intent: stationary -- the stationary solve returned info = 1; output written.
    2366 |  (certification) the state written is a relaxation snapshot: the run made no stationary claim, so
    2367 | Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
    2368 |    it is written certified=F with the reason "no stationary claim" and the run exits 0.

run.log of this case
  path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/run.log
  lines             1740
  md5 MEASURED      4a82d4fcc667ced98feb59fe0d74553c
  written           2026-09-20T08:01:52
  runs in the file  1 (a run begins at line 2)
  is the log the manifest names     yes
  holds this generation's block     NO
  boundary model of its restart     characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3   [line 32]

  published generations whose certificate stands in this log,
  located by the certificate text (a pass certification of a
  state the run then moved away from is not one):
    NONE. No pass recorded in this log ended in a state that
    was published as a generation of this case.

  the binary that wrote this log:
    `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_75d55d9d.x`, md5 `75d55d9d4fd0e748cd01d6e35713e34c`
    [/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/not_solved.md, written by models/run_case.sh for the run that produced this log]

  the last iterate this log records:
    the last outer pass   [line 1683]
      (EXHALE_main) outer pass 4: hydro info=2, worst gated species row  3.29E-02 of  1.0E-05 at cell 387 (elemental transport He/H partition), mass  1.16E+00, momentum  1.69E-07, energy  1.10E+00, omega 0.250, trust 1.0E-02,  391.37 s
    the last inner solve verdict   [line 1655]
      (JFNK) done info=2 ||R||=  1.163E+00  flux spread= 9.356E-03  non-monotone accepts=24  solve=4
    the last residual gate   [line 1682]
      (JFNK) gate NOT met: flux spread 9.356E-03 >= 2.000E-05

  the last lines of the file:
    1738 |  (JFNK) cells outside the tolerance of their row, of 500: mass 500 (worst cell 437), momentum 254 (500), energy 500 (1)
    1739 |  (JFNK) forcing term 1.00E-01: the cycle reached  9.548E-07 in 1 product(s) of 40; the requested linear tolerance was reached
    1740 |  (JFNK) it  12  ||R||=  1.171E+00  ||Fs||2=  5.63E-04  lam= 1.00E+00  dtau= 1.18E-03  gm=  1  worst r=  1.000  worst row: mass of cell 1  solve=5

the identity of the run that wrote the case `run.log`
  [/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/not_solved.md, written 2026-09-20T08:02:36]
  A log never names its own executable, so this is where the
  identity of a run that published nothing is kept. It belongs to
  the LOG-ONLY table of rev9 section 4.1 and never to a row that
  also carries a stored state.
    case               `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7`
    ending class       `stopped_by_wall_ceiling`
    solver verdict     info = none printed
    outer passes       4 of the 40 this run allowed
    pseudo-time start  1.0
    binary             `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_75d55d9d.x`, md5 `75d55d9d4fd0e748cd01d6e35713e34c`
    host, threads      `lart4`, OMP_NUM_THREADS=8
    seed               /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7/states/g0005_20260919T212143Z_12918db7  tier3:models/atomic_photochem_gj1132x0.10_kzzprofile  HeH=9.71107  target=9.71107  dlog10=0.0000  candidates=562  (a certified case at another XUV normalization)  compat=same system (profile, read from the state)
    earlier copies of this record, one per stopped run:
      /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/not_solved_20260919192839.md
      /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/not_solved_20260919192903.md
      /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/not_solved_20260920080236.md

the publication test, MEASURED here
  Hydro_ioniz.txt
    output/ md5       (none)   (written (none))
    generation md5    c19b3df722fb19c91662c633bc5f1bc1
  Ion_species.txt
    output/ md5       (none)   (written (none))
    generation md5    9f04c20d6c365f87067c081ec2054b53
  verdict             NOT DECIDABLE: `output/` does not hold both halves of a state, so it cannot be compared with the stored generation

  *** The last iterate of the case `run.log` (line 1683) stands AFTER every
  *** line of that log, none of which is the certificate of a published
  *** generation of this case, so no state was written from it. `output/`
  *** holds no state pair, so the publication test could not be made. It
  *** exists ONLY in that log: it cannot be reloaded, it has no generation, no
  *** manifest and no certificate of its own, and it may not share a row of
  *** any table with a stored state (rev9 section 4.1).

7. FIELDS THIS RECORD COULD NOT FILL
------------------------------------
 1. the tolerance-normalized verdict of the row 'hydrodynamic momentum row'
    the certificate states the row measure and its tolerance but prints
    neither a `verdict at cell` line nor a `gated at r >=` line for it
 2. the tolerance-normalized verdict of the row 'hydrodynamic energy row'
    the certificate states the row measure and its tolerance but prints
    neither a `verdict at cell` line nor a `gated at r >=` line for it
 3. a recorded md5 for the log /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/run.log
    the manifest records no `source_md5`, so the file on disk cannot be
    shown to be the file the certificate was read from
 4. the publication test for g0004_20260919T094145Z_259fe9c3
    `output/` of this case does not hold both Hydro_ioniz.txt and
    Ion_species.txt
```

### 3.3 `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13`

```
==============================================================================
STAGE A0 IDENTITY RECORD
==============================================================================
case            atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13
case directory  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13
generation      g0002_20260919T004834Z_d59ec9ff
selected by     the reference `latest_complete` of state_index.json
latest_complete g0002_20260919T004834Z_d59ec9ff   (state_index.json)
latest_certified (none)  (state_index.json)
                the index names NO certified generation: no state
                of this case carries a certificate that covers its
                own bytes.
assembled       2026-09-21T08:14:57 by models/identity_record.py
classification  HISTORICAL AUDIT (PLAN_20260920_rev9 section 3). Every
                stored value below is READ; every md5 this tool
                recomputed is MEASURED. No binary was run and no
                state was evaluated, so nothing here is Mode R.

1. THE STATE PAIR
-----------------
directory       /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/states/g0002_20260919T004834Z_d59ec9ff

  Hydro_ioniz.txt
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/states/g0002_20260919T004834Z_d59ec9ff/Hydro_ioniz.txt
    bytes             93748
    md5 MEASURED      741dc1415efcebfb4f2e13e818ac6b5a
    md5 READ          741dc1415efcebfb4f2e13e818ac6b5a   [states/g0002_20260919T004834Z_d59ec9ff/manifest.json `components`]
    verdict           matches the recorded md5

  Ion_species.txt
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/states/g0002_20260919T004834Z_d59ec9ff/Ion_species.txt
    bytes             447309
    md5 MEASURED      60ab4afeff02fcb97e5ab37c39edb754
    md5 READ          60ab4afeff02fcb97e5ab37c39edb754   [states/g0002_20260919T004834Z_d59ec9ff/manifest.json `components`]
    verdict           matches the recorded md5

  the header of Hydro_ioniz.txt, which the md5 covers:
    # EXHALE schema 2
    # columns r[Rp] rho[mH/cm3] v[cm/s] p[cgs] T[K] heat[erg/cm3/s] cool[erg/cm3/s]
    # rows 504: 2 ghost cells at each end; physical cells are rows 3 to 502
    # coupling: sec_ion=T sec_ion_step=0 recon=WENO3 certified=F cert_reason=failing_entries mode=init
    # provenance: git=3c73905ca8a7 tree=dirty run=2026-09-19T09:48:34
    # provenance: ck_input=174102422691385346 ck_base=-1 ck_metals=-1
    # provenance: recon=WENO3 base_bc=characteristic carrier=none restart_schema=3 resid_def=145 N=500
    # boundary_model characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
    # boundary_reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp] 1.0000000001000000E+000 1.0000000000000000E+000 3.3088953495373719E-001 1.0000000000000000E+000
    # restart_schema 1
    # reservoir He/H 2.1299999999999999E+00
    # species_columns 34 r HI HII HeI HeII HeIII HeITR CI CII CIII OI OII OIII NI NII NIII MgI MgII MgIII SiI SiII SiIII CaI CaII CaIII NaI NaII KI KII SI SII FeI FeII FeIII
    # grid N 500 R0[cm] 1.1273716464000001E+09 r_min[Rp] 1.0001933187778382E+00 r_max[Rp] 2.9031224061195065E+01 mode Mixed
    # constants set IAU2015+CODATA2018 RJ[cm] 7.1492000000000000E+09 kB[erg/K] 1.3806490000000000E-16
    # options He23S=T metals=F eos_metals=T mol=F molbase=F oxychem=F carrier=F carrier_newton=F iontrans=F he_diff=T he_metal_diff=F sec_ion=T caloric_mono=F excH=F base_ir=F mol_ir=F mol_heat=T visc=F cond=F jlya=0 wellbal=T
    # t_phys[s] 0.0000000000000000E+00
    # source git=3c73905ca8a7 tree=dirty restart_input=provenance_unknown

2. ORIGIN
---------
generation id       g0002_20260919T004834Z_d59ec9ff   [states/<gen>/manifest.json]
iteration_phase     evaluate   [manifest.json]
published_at        2026-09-19T09:48:35   [manifest.json]
state_written       2026-09-19T09:48:34   [manifest.json]
published_by        models/publish_state.py on lart4   [manifest.json]
source_directory    runs/r20260919T004834Z_15338/eval/output   [manifest.json]
run_id              (none)   [manifest.json]
certification       NOT CERTIFIED   [manifest.json `certification.status`]
state claim         certified=False cert_reason=failing_entries mode=init   [manifest]

parent              g0001_20260917T172702Z_f98147ad
  source            (none)
  compatibility     stated by the run that took it

seed                [provenance/g0002_20260919T004834Z_d59ec9ff/provenance_recovered*.json (RECOVERED, not recorded by the run)]
  established       True
  kind              a published generation
  case              atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13
  generation        g0001_20260917T172702Z_f98147ad
  directory         (none)
  chosen by         the manifest of this generation names a parent generation of this same case, which is the state the pass was started from (MODELS.md section 9.2)
  transformation    (none)
  confidence        established
  evidence          atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/states/g0002_20260919T004834Z_d59ec9ff/manifest.json
  note              generation g0001_20260917T172702Z_f98147ad of this same case

ending class        evaluate   [manifest.json `ending.class`]
ending reason       the certified generation g0001_20260917T172702Z_f98147ad handed back and measured as it stands under the D9 step 3 catalog binary; an evaluation has no solve ending
index ending_class  evaluate   [state_index.json]
solver verdict      info = (none), ||R|| = (none)   [manifest.json `ending.solver`]

the case `ENDING` file (/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/ENDING, written 2026-09-20T07:02:34):
    stopped_by_wall_ceiling: stopped by the wall ceiling of 30m that models/run_campaign.sh applies (SIGTERM 1800 s after the campaign started the case), in the wind pass after outer pass 5 had completed; that pass had written no complete state of its own, so nothing new is published
    NOTE: this file is the ending of the LAST run of the case
    directory. It is not in general the ending of the run that
    published this generation, which is the `ending` block above.

3. THE OPERATOR
---------------
binary path         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_7670f310.x   [manifest.json `source_identity`]
  on disk now       yes
  md5 READ          7670f31031fb4db91d27b44cb0da6f70   [manifest.json]
  md5 MEASURED      7670f31031fb4db91d27b44cb0da6f70
  source manifest   BINARY_MANIFEST_7670f31031fb.txt
  its md5 READ      0376f5510808d64d32af864a9cf55742
  its md5 MEASURED  0376f5510808d64d32af864a9cf55742

model options       He23S=T metals=F eos_metals=T mol=F molbase=F oxychem=F carrier=F carrier_newton=F iontrans=F he_diff=T he_metal_diff=F sec_ion=T caloric_mono=F excH=F base_ir=F mol_ir=F mol_heat=T visc=F cond=F jlya=0 wellbal=T   [manifest.json `model_identity.options`]
reservoir           He/H 2.1299999999999999E+00

boundary model of the STORED state:
  state header      characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
                    [# boundary_model of Hydro_ioniz.txt, part of the hashed bytes]
  manifest          characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
  boundary reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp] 1.0000000001000000E+000 1.0000000000000000E+000 3.3088953495373719E-001 1.0000000000000000E+000

boundary model a log of this case names for its own restart:
  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/run.log:20
    (load_IC) boundary model of the restart: characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3 (this run's).

  *** THE TWO ARE DIFFERENT EXPERIMENTS.
  *** the stored state was written under
  ***   characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
  *** a log of this case ran under
  ***   characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3
  *** A run that rebuilds the boundary from the physical column
  *** and its own reservoir is the CURRENT operator on an OLD
  *** state (rev9 section 3); its residual is not the residual
  *** the stored state was written with, and the two may not
  *** share a row of any table.

4. THE CONFIGURATION
--------------------
  input.inp                      matches the recorded md5
    md5 MEASURED      c355d9827b63e92b3d2f711bf73e26a7
    md5 READ          c355d9827b63e92b3d2f711bf73e26a7   [manifest.json `configuration_identity`]
  base.inp                       absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  metals.inp                     absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  opacity.inp                    absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  lower_atmosphere_profile.dat   absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]

  spectrum
    named in input    ../../../sed/lhs1140_sed_gj1132_at_b_xuv0p01.txt   [input.inp `Spectrum file:`]
    resolved to       /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/sed/lhs1140_sed_gj1132_at_b_xuv0p01.txt
    resolution rule   relative to the case directory, which is the working directory of a run: the binary reads ./input.inp and opens the string as written
    md5 MEASURED      7da5fa41a0706feacac74b7d245fbec3
    md5 READ          7da5fa41a0706feacac74b7d245fbec3   [manifest.json]
    verdict           matches the recorded md5

  EXHALE_resolved.out
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/EXHALE_resolved.out
    md5 MEASURED      ef3fe44862164944a0f83ad19e6a96fb
    written           2026-09-20T06:32:41
    keys              38
    route keys it DOES carry:
      well_balanced              T
      carrier_transport          F
      carrier_in_newton          F
      carrier_newton_on_stall    F
      ionization_transport       F
      oxygen_chemistry           F
    route keys it does NOT carry, taken from elsewhere:
      - the boundary model identity (the state file header and the manifest carry it; the resolved writer never emits it, write_setup_report.f90:906-1120)
      - the numerical flux and the reconstruction method (stated in prose by EXHALE_setup.out)
      - the residual assembly selector EXHALE_RESID_QUAD, which is an environment variable read at hydrodynamic_rows.f90:177 and is in no resolved record

  EXHALE_setup.out (the resolved route, in prose)
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/EXHALE_setup.out
    md5 MEASURED      595740450f26e21ec9f52f3c1cdd9efa
    written           2026-09-20T06:32:41
    numerical flux    ROE   [line 46]
    reconstruction method PLM   [line 47]
    well balanced     on   [line 48]
    base boundary     characteristic condition at the face r_edg(0)   [line 53]

    NOTE: EXHALE_resolved.out and EXHALE_setup.out of this case
    were written AFTER this generation was published, so they are
    the record of a LATER run in the same directory. The input
    md5 comparison above is what ties the configuration to this
    generation; the route lines are read from the later report
    and are the route of that run.

5. THE ROW REGISTRY
-------------------
certificate         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/states/g0002_20260919T004834Z_d59ec9ff/certification.txt
  md5 MEASURED      a05585139d5090f197f993e34b7a23f1
  active equations 6 of 30 in the inventory   [line 2]

  rows evaluated (6):
    hydrodynamic mass row                         ABOVE  max= 9.589E-01 at cell 2 (r = 1.00039), tol= 3.7E-08   [line 14]
    hydrodynamic momentum row                     ABOVE  max= 5.539E-06 at cell 1 (r = 1.00019), tol= 1.0E-08   [line 17]
    hydrodynamic energy row                       ABOVE  max= 1.112E+00 at cell 1 (r = 1.00019), tol= 1.0E-06   [line 19]
    elemental transport He/H partition            within max= 8.535E-01 at cell 2 (r = 1.00039), tol= 1.0E-05   [line 21]
    level balance He 2^3S                         within max= 1.221E-20 at cell 464 (r = 15.96711), tol= 1.0E-06   [line 24]
    eliminated-species closure System_HeH_TR      within max= 1.475E-16 at cell 464 (r = 15.96711), tol= 1.0E-06   [line 26]

  verdict           NOT CERTIFIED: 3 entry/entries of the inventory refuse it   [line 38]

  refusing rows (3):

    hydrodynamic mass row: row measure  9.589E-01 above  3.7E-08 at cell 2
      row max         9.589E-01 at cell 2 (r = 1.00039, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       3.7E-08
      volume measure  1.014E-01
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      verdict cell    2 (r = 1.00039E+00 READ, 1.00039 MEASURED): 9.589E-01 against 3.7E-08, distance 2.610E+07
                      [line 16]

    hydrodynamic momentum row: row measure  5.539E-06 above  1.0E-08 at cell 1
      row max         5.539E-06 at cell 1 (r = 1.00019, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-08
      volume measure  4.259E-07
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      the certificate prints no tolerance-normalized verdict
      cell for this row (see part 7)

    hydrodynamic energy row: row measure  1.112E+00 above  1.0E-06 at cell 1
      row max         1.112E+00 at cell 1 (r = 1.00019, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-06
      volume measure  2.144E-01
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      the certificate prints no tolerance-normalized verdict
      cell for this row (see part 7)

6. THE LOG BLOCK
----------------
the log the manifest names for this generation
  path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/runs/r20260919T004834Z_15338/pp.log
  named in          states/g0002_20260919T004834Z_d59ec9ff/manifest.json `certification.source`
  on disk now       yes
  md5 MEASURED      404960269c7bc6e7da817cb629e28ec3
  md5 READ          (none)
  written           2026-09-19T09:48:35

where the certificate of this generation actually stands, located by
its own text and not by a file name:
  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/runs/r20260919T004834Z_15338/pp.log
  lines 169 to 213 (45 lines, the whole of certification.txt)

  the block:
     169 |  (certification) FINAL STATE AS WRITTEN -- this run is NOT a certified stationary solution
     170 |    active equations 6 of 30 in the inventory
     171 |    tolerances (convergence study, contract sections 9 and 10; the hydrodynamic values were
     172 |    re-anchored on the self-consistent JFNK residual after B5a):
     173 |      hydrodynamic mass  3.00E-12, momentum  1.00E-08, energy  1.00E-06
     174 |      the mass value is a FLOOR: where the cell's own rounding of the flux difference stands above it, 10.0 times that floor is the tolerance
     175 |      of that cell (cert_tol_mass_at), and the verdict on the mass row is taken cell by cell
     176 |      closure  1.00E-06, level  1.00E-06; and IN THE WIND, at r >= 1.20E+00, carrier  1.00E-05
     177 |      and elemental transport  1.00E-05; the species rows of the cells below that radius are
     178 |      REPORTED AND DO NOT GATE: the layer r < 1.10E+00, whose element fluxes
     179 |      are not conserved (spread 9.8), and the band above it, whose discretization (4.6e-4 at the
     180 |      binding radius) stands above any tolerance one could set there
     181 |      the run's own "Resid tol" is  5.00E-05, which is the SOLVER's target and not a certification tolerance
     182 |    hydrodynamic mass row                                evaluated       max= 9.589E-01  vol= 1.014E-01  cell=2  tol= 3.7E-08  ABOVE   
     183 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     184 |         verdict at cell 2 (r= 1.00039E+00):  9.589E-01 against  3.7E-08, distance  2.610E+07 -- the rounding anchor of that cell
     185 |    hydrodynamic momentum row                            evaluated       max= 5.539E-06  vol= 4.259E-07  cell=1  tol= 1.0E-08  ABOVE   
     186 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     187 |    hydrodynamic energy row                              evaluated       max= 1.112E+00  vol= 2.144E-01  cell=1  tol= 1.0E-06  ABOVE   
     188 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     189 |    elemental transport He/H partition                   evaluated       max= 8.535E-01  vol= 2.781E-02  cell=2  tol= 1.0E-05  within  
     190 |         scale: sum of the row's own terms [g cm^-3 s^-1]; floor 1e-20 rho X_base/(R0/v0)
     191 |         gated at r >= 1.200E+00 against  1.00E-05:  6.745E-06 at cell 393; the cells below that radius are reported only
     192 |    level balance He 2^3S                                evaluated       max= 1.221E-20  vol= 1.601E-21  cell=464  tol= 1.0E-06  within  
     193 |         scale: the row's own turnover rate; dimensionless, scale 1
     194 |    eliminated-species closure System_HeH_TR             evaluated       max= 1.475E-16  vol= 1.759E-17  cell=464  tol= 1.0E-06  within  
     195 |         scale: the row's own turnover rate; dimensionless, scale 1
     196 |    cells without a chemical root: 0 (acceptance class 4 or 6)
     197 |    mass closure of the composition, max |sum_i f_i A_i - 1| =  4.440E-16 at cell 207, r =   1.16838
     198 |      reported only: the sum is one by the definition of the mass fractions, so this is the arithmetic of the state and no equation of it
     199 |    validity states (B1a section 4):
     200 |      active unvalidated physics: not produced
     201 |      out-of-domain closure activations: 0
     202 |      rejected trials with no adopted contribution: not produced
     203 |      unbudgeted accepted corrections: 0
     204 |      specified external reservoirs: not produced
     205 |    attempts (not validity states): 0 energy-update and 0 conduction cell(s) reached the lower bracket end
     206 |    NOT CERTIFIED: 3 entry/entries of the inventory refuse it
     207 |      hydrodynamic mass row: row measure  9.589E-01 above  3.7E-08 at cell 2
     208 |      hydrodynamic momentum row: row measure  5.539E-06 above  1.0E-08 at cell 1
     209 |      hydrodynamic energy row: row measure  1.112E+00 above  1.0E-06 at cell 1
     210 |  
     211 |  (certification) the state was written in full; the run exits with status 2 because it is not certified.
     212 | Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
     213 | STOP 2

run.log of this case
  path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/run.log
  lines             2267
  md5 MEASURED      a466ea9baf9e393a2203b84d466fdc9a
  written           2026-09-20T07:02:29
  runs in the file  1 (a run begins at line 2)
  is the log the manifest names     no
  holds this generation's block     NO
  boundary model of its restart     characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3   [line 20]

  published generations whose certificate stands in this log,
  located by the certificate text (a pass certification of a
  state the run then moved away from is not one):
    NONE. No pass recorded in this log ended in a state that
    was published as a generation of this case.

  the binary that wrote this log:
    `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_75d55d9d.x`, md5 `75d55d9d4fd0e748cd01d6e35713e34c`
    [/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/not_solved.md, written by models/run_case.sh for the run that produced this log]

  the last iterate this log records:
    the last outer pass   [line 2152]
      (EXHALE_main) outer pass 5: hydro info=2, worst gated species row  4.38E-02 of  1.0E-05 at cell 331 (elemental transport He/H partition), mass  2.46E-01, momentum  3.29E-07, energy  1.22E-01, omega 0.500, trust 1.0E-02,  159.29 s
    the last inner solve verdict   [line 2124]
      (JFNK) done info=2 ||R||=  2.458E-01  flux spread= 1.956E-02  non-monotone accepts=12  solve=5
    the last residual gate   [line 2151]
      (JFNK) gate NOT met: flux spread 1.956E-02 >= 2.000E-05

  the last lines of the file:
    2265 |  (JFNK) cells outside the tolerance of their row, of 500: mass 500 (worst cell 402), momentum 5 (500), energy 497 (3)
    2266 |  (JFNK) forcing term 1.00E-01: the cycle reached  1.172E-08 in 1 product(s) of 40; the requested linear tolerance was reached
    2267 |  (JFNK) it  26  ||R||=  2.542E-01  ||Fs||2=  1.97E-05  lam= 1.00E+00  dtau= 6.47E-04  gm=  1  worst r=  1.261  worst row: mass of cell 232  solve=6

the identity of the run that wrote the case `run.log`
  [/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/not_solved.md, written 2026-09-20T07:02:34]
  A log never names its own executable, so this is where the
  identity of a run that published nothing is kept. It belongs to
  the LOG-ONLY table of rev9 section 4.1 and never to a row that
  also carries a stored state.
    case               `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13`
    ending class       `stopped_by_wall_ceiling`
    solver verdict     info = none printed
    outer passes       5 of the 40 this run allowed
    pseudo-time start  1.0
    binary             `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_75d55d9d.x`, md5 `75d55d9d4fd0e748cd01d6e35713e34c`
    host, threads      `lart4`, OMP_NUM_THREADS=8
    seed               /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13/states/g0005_20260919T212847Z_80bee6f1  tier3:models/atomic_scalar_gj1132x0.10_kzz1e9  HeH=2.13  target=2.13  dlog10=0.0000  candidates=101  (a certified case at another XUV normalization)  compat=same system (scalar, read from the state)

the publication test, MEASURED here
  Hydro_ioniz.txt
    output/ md5       741dc1415efcebfb4f2e13e818ac6b5a   (written 2026-09-19T09:48:34)
    generation md5    741dc1415efcebfb4f2e13e818ac6b5a
  Ion_species.txt
    output/ md5       60ab4afeff02fcb97e5ab37c39edb754   (written 2026-09-19T09:48:34)
    generation md5    60ab4afeff02fcb97e5ab37c39edb754
  verdict             the pair in `output/` is BITWISE the stored generation

  *** The last iterate of the case `run.log` (line 2152) stands AFTER every
  *** line of that log, none of which is the certificate of a published
  *** generation of this case, so no state was written from it. The pair in
  *** `output/` is bitwise the stored generation
  *** g0002_20260919T004834Z_d59ec9ff, which confirms that nothing was
  *** published after it. It exists ONLY in that log: it cannot be reloaded,
  *** it has no generation, no manifest and no certificate of its own, and it
  *** may not share a row of any table with a stored state (rev9 section 4.1).

7. FIELDS THIS RECORD COULD NOT FILL
------------------------------------
 1. the tolerance-normalized verdict of the row 'hydrodynamic momentum row'
    the certificate states the row measure and its tolerance but prints
    neither a `verdict at cell` line nor a `gated at r >=` line for it
 2. the tolerance-normalized verdict of the row 'hydrodynamic energy row'
    the certificate states the row measure and its tolerance but prints
    neither a `verdict at cell` line nor a `gated at r >=` line for it
 3. a recorded md5 for the log /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/runs/r20260919T004834Z_15338/pp.log
    the manifest records no `source_md5`, so the file on disk cannot be
    shown to be the file the certificate was read from
```

### 3.4 `molecular_scalar_gj1132_kzz1e9/HeH0.083`

```
==============================================================================
STAGE A0 IDENTITY RECORD
==============================================================================
case            molecular_scalar_gj1132_kzz1e9/HeH0.083
case directory  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083
generation      g0004_20260920T053734Z_79b42a03
selected by     the reference `latest_complete` of state_index.json
latest_complete g0004_20260920T053734Z_79b42a03   (state_index.json)
latest_certified (none)  (state_index.json)
                the index names NO certified generation: no state
                of this case carries a certificate that covers its
                own bytes.
assembled       2026-09-21T08:14:57 by models/identity_record.py
classification  HISTORICAL AUDIT (PLAN_20260920_rev9 section 3). Every
                stored value below is READ; every md5 this tool
                recomputed is MEASURED. No binary was run and no
                state was evaluated, so nothing here is Mode R.

1. THE STATE PAIR
-----------------
directory       /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/states/g0004_20260920T053734Z_79b42a03

  Hydro_ioniz.txt
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/states/g0004_20260920T053734Z_79b42a03/Hydro_ioniz.txt
    bytes             93830
    md5 MEASURED      fd882d8918052069f8ed229dd376f8c4
    md5 READ          fd882d8918052069f8ed229dd376f8c4   [states/g0004_20260920T053734Z_79b42a03/manifest.json `components`]
    verdict           matches the recorded md5

  Ion_species.txt
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/states/g0004_20260920T053734Z_79b42a03/Ion_species.txt
    bytes             499794
    md5 MEASURED      87288b03ad88a92b7b8b8b6c2b1bdd66
    md5 READ          87288b03ad88a92b7b8b8b6c2b1bdd66   [states/g0004_20260920T053734Z_79b42a03/manifest.json `components`]
    verdict           matches the recorded md5

  the header of Hydro_ioniz.txt, which the md5 covers:
    # EXHALE schema 2
    # columns r[Rp] rho[mH/cm3] v[cm/s] p[cgs] T[K] heat[erg/cm3/s] cool[erg/cm3/s]
    # rows 504: 2 ghost cells at each end; physical cells are rows 3 to 502
    # coupling: sec_ion=T sec_ion_step=0 recon=WENO3 certified=F cert_reason=no_stationary_claim mode=init
    # provenance: git=3c73905ca8a7 tree=dirty run=2026-09-20T14:37:34
    # provenance: ck_input=3499329486697119869 ck_base=1774340865516947195 ck_metals=-1
    # provenance: recon=WENO3 base_bc=characteristic carrier=transported restart_schema=3 resid_def=145 N=500
    # boundary_model characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3
    # boundary_reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp] 5.3832725770134127E-001 1.0000000000000000E+000 4.3847198334645321E-001 1.0000000000000000E+000
    # restart_schema 1
    # reservoir He/H 8.3000000000000004E-02
    # species_columns 38 r HI HII HeI HeII HeIII HeITR CI CII CIII OI OII OIII NI NII NIII MgI MgII MgIII SiI SiII SiIII CaI CaII CaIII NaI NaII KI KII SI SII FeI FeII FeIII H2 H2p H3p HeHp
    # grid N 500 R0[cm] 1.1273716464000001E+09 r_min[Rp] 1.0001933187778382E+00 r_max[Rp] 2.9031224061195065E+01 mode Mixed
    # constants set IAU2015+CODATA2018 RJ[cm] 7.1492000000000000E+09 kB[erg/K] 1.3806490000000000E-16
    # options He23S=T metals=F eos_metals=T mol=T molbase=T oxychem=F carrier=T carrier_newton=F iontrans=F he_diff=T he_metal_diff=F sec_ion=T caloric_mono=F excH=F base_ir=F mol_ir=F mol_heat=T visc=F cond=F jlya=0 wellbal=T
    # t_phys[s] 0.0000000000000000E+00
    # source git=3c73905ca8a7 tree=dirty restart_input=provenance_unknown

2. ORIGIN
---------
generation id       g0004_20260920T053734Z_79b42a03   [states/<gen>/manifest.json]
iteration_phase     evaluate   [manifest.json]
published_at        2026-09-20T14:37:36   [manifest.json]
state_written       2026-09-20T14:37:34   [manifest.json]
published_by        models/publish_state.py on lart4   [manifest.json]
source_directory    output   [manifest.json]
run_id              r20260920T053734Z_37865   [manifest.json]
certification       NOT CERTIFIED   [manifest.json `certification.status`]
state claim         certified=False cert_reason=no_stationary_claim mode=init   [manifest]

parent              g0002_20260916T220756Z_f9fca548
  source            (none)
  compatibility     stated by the run that took it

seed                [states/g0004_20260920T053734Z_79b42a03/manifest.json `seed`]
  established       True
  kind              internal continuation
  case              molecular_scalar_gj1132_kzz1e9/HeH0.083
  generation        g0002_20260916T220756Z_f9fca548
  directory         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/states/g0002_20260916T220756Z_f9fca548
  chosen by         the runner: this pass was started from the generation named here, a generation of this same case
  transformation    used as it stands: the pair was reloaded, with no mapping and no conversion
  md5 Hydro_ioniz.txt1ee5430a3c8560c40f017bc5eb618a8b
  md5 Ion_species.txt0e4aef56f5b58fbf878b7638d22dffed

ending class        evaluate   [manifest.json `ending.class`]
ending reason       the generation g0002_20260916T220756Z_f9fca548 handed back and measured as it stands; an evaluation has no solve ending
index ending_class  evaluate   [state_index.json]
solver verdict      info = (none), ||R|| = (none)   [manifest.json `ending.solver`]

the case `ENDING` file (/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/ENDING, written (none)):
    (the case records no ENDING file)

3. THE OPERATOR
---------------
binary path         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_75d55d9d.x   [manifest.json `source_identity`]
  on disk now       yes
  md5 READ          75d55d9d4fd0e748cd01d6e35713e34c   [manifest.json]
  md5 MEASURED      75d55d9d4fd0e748cd01d6e35713e34c
  source manifest   BINARY_MANIFEST_75d55d9d4fd0.txt
  its md5 READ      e8b5b83f3c8c82da8aae93000a25e75e
  its md5 MEASURED  e8b5b83f3c8c82da8aae93000a25e75e

model options       He23S=T metals=F eos_metals=T mol=T molbase=T oxychem=F carrier=T carrier_newton=F iontrans=F he_diff=T he_metal_diff=F sec_ion=T caloric_mono=F excH=F base_ir=F mol_ir=F mol_heat=T visc=F cond=F jlya=0 wellbal=T   [manifest.json `model_identity.options`]
reservoir           He/H 8.3000000000000004E-02

boundary model of the STORED state:
  state header      characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3
                    [# boundary_model of Hydro_ioniz.txt, part of the hashed bytes]
  manifest          characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3
  boundary reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp] 5.3832725770134127E-001 1.0000000000000000E+000 4.3847198334645321E-001 1.0000000000000000E+000

4. THE CONFIGURATION
--------------------
  input.inp                      matches the recorded md5
    md5 MEASURED      c3916719a26a2f5a7d1d7025a7cd71fa
    md5 READ          c3916719a26a2f5a7d1d7025a7cd71fa   [manifest.json `configuration_identity`]
  base.inp                       matches the recorded md5
    md5 MEASURED      1e33d197ead9269661dc6bb0d88c0e26
    md5 READ          1e33d197ead9269661dc6bb0d88c0e26   [manifest.json `configuration_identity`]
  metals.inp                     absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  opacity.inp                    absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  lower_atmosphere_profile.dat   absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]

  spectrum
    named in input    ../../../sed/lhs1140_sed_gj1132_at_b.txt   [input.inp `Spectrum file:`]
    resolved to       /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt
    resolution rule   relative to the case directory, which is the working directory of a run: the binary reads ./input.inp and opens the string as written
    md5 MEASURED      ac73a75c6396793435ed72e1cf34c427
    md5 READ          ac73a75c6396793435ed72e1cf34c427   [manifest.json]
    verdict           matches the recorded md5

  EXHALE_resolved.out
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/EXHALE_resolved.out
    md5 MEASURED      2027e7d8839e626335cac0dee4857af1
    written           2026-09-16T09:12:43
    keys              33
    route keys it DOES carry:
      well_balanced              T
      carrier_transport          T
      carrier_in_newton          F
      ionization_transport       F
      oxygen_chemistry           F
    route keys it does NOT carry, taken from elsewhere:
      - the boundary model identity (the state file header and the manifest carry it; the resolved writer never emits it, write_setup_report.f90:906-1120)
      - the numerical flux and the reconstruction method (stated in prose by EXHALE_setup.out)
      - the residual assembly selector EXHALE_RESID_QUAD, which is an environment variable read at hydrodynamic_rows.f90:177 and is in no resolved record

  EXHALE_setup.out (the resolved route, in prose)
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/EXHALE_setup.out
    md5 MEASURED      8471ed461e9d65764c8f19b5033c65a9
    written           2026-09-16T09:12:43
    numerical flux    HLLC   [line 52]
    reconstruction method PLM   [line 53]
    well balanced     on   [line 54]
    base boundary     characteristic condition at the face r_edg(0)   [line 58]

5. THE ROW REGISTRY
-------------------
certificate         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/states/g0004_20260920T053734Z_79b42a03/certification.txt
  md5 MEASURED      a153cc4915d8bfc89b6797977de5b0a1
  active equations 7 of 30 in the inventory   [line 2]

  rows evaluated (7):
    hydrodynamic mass row                         within max= 1.393E-08 at cell 1 (r = 1.00019), tol= 1.9E-08   [line 14]
    hydrodynamic momentum row                     within max= 4.854E-13 at cell 1 (r = 1.00019), tol= 1.0E-08   [line 17]
    hydrodynamic energy row                       within max= 6.461E-09 at cell 31 (r = 1.00599), tol= 1.0E-06   [line 19]
    carrier balance H2                            ABOVE  max= 1.744E-02 at cell 1 (r = 1.00019), tol= 1.0E-05   [line 21]
    elemental transport He/H partition            ABOVE  max= 4.656E-02 at cell 2 (r = 1.00039), tol= 1.0E-05   [line 24]
    level balance He 2^3S                         within max= 1.410E-19 at cell 459 (r = 14.71780), tol= 1.0E-06   [line 27]
    eliminated-species closure System_HeH_mol     within max= 3.331E-16 at cell 5 (r = 1.00097), tol= 1.0E-06   [line 29]

  verdict           NOT CERTIFIED: 2 entry/entries of the inventory refuse it   [line 52]

  refusing rows (2):

    carrier balance H2: gated row measure  8.150E-04 above  1.0E-05 at cell 499 (a wind cell)
      row max         1.744E-02 at cell 1 (r = 1.00019, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-05
      volume measure  2.983E-03
      scale           row_terms [cm^-3 s^-1], the sum of the row's own terms with the 1e-20 free-e
      gated at r >= 1.200E+00 against 1.00E-05: 8.150E-04 at cell 499 (r = 28.54893 MEASURED)
                      [line 23]

    elemental transport He/H partition: gated row measure  2.503E-04 above  1.0E-05 at cell 280 (a wind cell)
      row max         4.656E-02 at cell 2 (r = 1.00039, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-05
      volume measure  2.690E-04
      scale           sum of the row's own terms [g cm^-3 s^-1]; floor 1e-20 rho X_base/(R0/v0)
      gated at r >= 1.200E+00 against 1.00E-05: 2.503E-04 at cell 280 (r = 1.60443 MEASURED)
                      [line 26]

6. THE LOG BLOCK
----------------
the log the manifest names for this generation
  path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/runs/r20260920T053734Z_37865/pp.log
  named in          states/g0004_20260920T053734Z_79b42a03/manifest.json `certification.source`
  on disk now       yes
  md5 MEASURED      678b3f1a74c26f0c8c817324a422ae88
  md5 READ          678b3f1a74c26f0c8c817324a422ae88
  verdict           matches the recorded md5
  written           2026-09-20T14:37:35

where the certificate of this generation actually stands, located by
its own text and not by a file name:
  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/runs/r20260920T053734Z_37865/pp.log
  lines 101 to 196 (96 lines, the whole of certification.txt)

  the block:
     101 |  (certification) work state of the loaded restart (Restart intent: stationary)
     102 |    active equations 7 of 30 in the inventory
     103 |    tolerances (convergence study, contract sections 9 and 10; the hydrodynamic values were
     104 |    re-anchored on the self-consistent JFNK residual after B5a):
     105 |      hydrodynamic mass  3.00E-12, momentum  1.00E-08, energy  1.00E-06
     106 |      the mass value is a FLOOR: where the cell's own rounding of the flux difference stands above it, 10.0 times that floor is the tolerance
     107 |      of that cell (cert_tol_mass_at), and the verdict on the mass row is taken cell by cell
     108 |      closure  1.00E-06, level  1.00E-06; and IN THE WIND, at r >= 1.20E+00, carrier  1.00E-05
     109 |      and elemental transport  1.00E-05; the species rows of the cells below that radius are
     110 |      REPORTED AND DO NOT GATE: the layer r < 1.10E+00, whose element fluxes
     111 |      are not conserved (spread 9.8), and the band above it, whose discretization (4.6e-4 at the
     112 |      binding radius) stands above any tolerance one could set there
     113 |      the run's own "Resid tol" is  5.00E-05, which is the SOLVER's target and not a certification tolerance
     114 |    hydrodynamic mass row                                evaluated       max= 1.393E-08  vol= 2.497E-10  cell=1  tol= 1.9E-08  within  
     115 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     116 |         verdict at cell 1 (r= 1.00019E+00):  1.393E-08 against  1.9E-08, distance  7.175E-01 -- the rounding anchor of that cell
     117 |    hydrodynamic momentum row                            evaluated       max= 4.854E-13  vol= 2.468E-14  cell=1  tol= 1.0E-08  within  
     118 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     119 |    hydrodynamic energy row                              evaluated       max= 6.461E-09  vol= 9.207E-10  cell=31  tol= 1.0E-06  within  
     120 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     121 |    carrier balance H2                                   evaluated       max= 1.744E-02  vol= 2.983E-03  cell=1  tol= 1.0E-05  ABOVE   
     122 |         scale: row_terms [cm^-3 s^-1], the sum of the row's own terms with the 1e-20 free-e
     123 |         gated at r >= 1.200E+00 against  1.00E-05:  8.150E-04 at cell 499; the cells below that radius are reported only
     124 |    elemental transport He/H partition                   evaluated       max= 4.656E-02  vol= 2.690E-04  cell=2  tol= 1.0E-05  ABOVE   
     125 |         scale: sum of the row's own terms [g cm^-3 s^-1]; floor 1e-20 rho X_base/(R0/v0)
     126 |         gated at r >= 1.200E+00 against  1.00E-05:  2.503E-04 at cell 280; the cells below that radius are reported only
     127 |    level balance He 2^3S                                evaluated       max= 1.410E-19  vol= 2.941E-20  cell=459  tol= 1.0E-06  within  
     128 |         scale: the row's own turnover rate; dimensionless, scale 1
     129 |    eliminated-species closure System_HeH_mol            evaluated       max= 3.331E-16  vol= 4.400E-17  cell=5  tol= 1.0E-06  within  
     130 |         scale: the row's own turnover rate; dimensionless, scale 1
     131 |    face mass flux r_f^2 (rho v)_f of this state (reported beside the rows above; it gates nothing):
     132 |      operator: WENO3, numerical flux: HLLC, well balanced: T
     133 |      in units of the wind-window mean of rho v r^2, which is  1.11873E-07
     134 |        minimum over the faces  9.998913223006088E-01 at face 1
     135 |        maximum over the faces  9.998913362282768E-01 at face 0
     136 |        base face  9.998913362282768E-01, offset from the window mean -1.08664E-04
     137 |        a face value and the mean of a cell-centered product are not the same quantity, so an offset of the size of the
     138 |        wind's own discretization is not a mass leak
     139 |      base contact direction: read from the wind window, which agrees with the base face flux (w_rev = 0.0)
     140 |    cells without a chemical root: 0 (acceptance class 4 or 6)
     141 |    mass closure of the composition, max |sum_i f_i A_i - 1| =  6.014E-16 at cell 489, r =  24.14200
     142 |      reported only: the sum is one by the definition of the mass fractions, so this is the arithmetic of the state and no equation of it
     143 |    validity states (B1a section 4):
     144 |      active unvalidated physics: not produced
     145 |      out-of-domain closure activations: 334
     146 |        H3+ cooling domain records (informational; they do not invalidate):
     147 |          collider below the table: 272, temperature below the fits: 0, above them: 0, outside the non-LTE range: 62
     148 |      rejected trials with no adopted contribution: not produced
     149 |      unbudgeted accepted corrections: 0
     150 |      specified external reservoirs: not produced
     151 |    attempts (not validity states): 0 energy-update and 0 conduction cell(s) reached the lower bracket end
     152 |    NOT CERTIFIED: 2 entry/entries of the inventory refuse it
     153 |      carrier balance H2: gated row measure  8.150E-04 above  1.0E-05 at cell 499 (a wind cell)
     154 |      elemental transport He/H partition: gated row measure  2.503E-04 above  1.0E-05 at cell 280 (a wind cell)
     155 |  
     156 |  
     157 |  (EXHALE_main) the two answers of this evaluation:
     158 |    original claim of the file: certified=F cert_reason=no_stationary_claim
     159 |    original claim: NONE MADE -- the file states no stationary claim, so nothing is reproduced or refused
     160 | Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
     161 |    work state verdict: NOT CERTIFIED -- the certification report above names the entries that refuse it
     162 |    the pair written into the state: certified=F cert_reason=no_stationary_claim
     163 |      steps: 0 accepted of 0 attempted
     164 |      initialization / continuation: no physical elapsed time is claimed
     165 | 
     166 |    ATTEMPTED-STEP CONTROLLER
     167 |      outer attempts: 0, of which 0 were refused at the adoption boundary
     168 |      steps that needed at least one retry: 0
     169 |      hydro stage attempts (the nested positivity retry, a DIFFERENT count): 0 in 0 halvings
     170 |      ioniz-eq [marching states]: 2 equilibrium sweep(s)
     171 |      ioniz-eq molecular hybrd1 info: 1 1014, 2 0, 3 0, 4 0, 5 0, 0 0
     172 |      ioniz-eq acceptance: 1008 converged root(s) (max res  6.04E-17), 0 root(s) without solver convergence (max res  0.00E+00), 0 projected/handback root(s) (max res  0.00E+00)
     173 |      ioniz-eq residual decades (<=1e-16 .. >=1e-1) converged:   1014 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     174 |      flux correction: 0 positivity repair call(s), 0 face flux(es) dropped to first order, 0 of them in accepted steps
     175 |      coupled source step:   0.00 passes per coupled step on average over 0 of them, worst 0 (step 0)
     176 |      coupled source step: every step reached its fixed point within the pass cap
     177 |      coupled source step: 0 exit(s) on the estimated error, 0 on the last increment
     178 |      chemistry mass closure: worst |sum n_s m_s - rho|/rho = 0.000E+00 at step 0
     179 |      formation reservoir change, worst |du_form| / |dt (heat-cool)| = 0.000E+00
     180 |  (write_cool_breakdown_eq) max |sum(channels)/cool - 1| =  3.07E-16
     181 |  (EXHALE_main.f90) Starting the post processing routine..
     182 |  (post_process_adv) advection correction: 41 of 503 cells kept at the equilibrium ionization.
     183 |  (post_process_adv) energy correction: 0 of 502 cells kept the run temperature (thermal Damkohler > 1); largest Damkohler  1.33E-01 at r =  1.1205 Rp.
     184 |  (post_process_adv) conditional correction: mass row |R_1|/s_1 of the input state up to  1.01E+00 at r =  1.0000 Rp; above  1.00E-02, the fraction a corrected row is accurate to, in 3 of 503 cells, of which 0 the ionization conditions do not already refuse.
     185 |  (post_process_adv) stationarity: the mass row was assembled from the primitive state handed in, which is the run's state re-formed and not its conserved variables.
     186 |  (post_process_adv) energy correction: enthalpy flux of the mass-flux divergence / other terms up to  3.64E+03 at r =  1.0010 Rp; dominant in 11 of 503 cells (reported, not a refusal; level  1.00E+00).
     187 |  (post_process_adv) T solve: 1 of 502 cells fell back to eq T (non-positive or out-of-band root).
     188 |  (post_process_adv) energy equation: 7 cell solves stopped improving at a residual already at the cancellation floor of their own terms; the iterate is the root and was kept.
     189 |  (EXHALE_main.f90) Post processing routine done.
     190 |   
     191 |  ----- Results -----
     192 |   
     193 |  ---> 2D approximate method: Mdot
     194 |  ---> Log10 of steady-state Mdot =  7.39 g/s
     195 |  (EXHALE_main) Restart intent: stationary evaluate -- the loaded state was measured, the work state and
     196 |    the profiles derived from it were written, and no step and no solve were taken.

run.log of this case
  path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/run.log
  lines             8821
  md5 MEASURED      337102ed5fa45a65efcd4c9158c64af1
  written           2026-09-16T11:33:29
  runs in the file  1 (a run begins at line 2)
  is the log the manifest names     no
  holds this generation's block     NO

  published generations whose certificate stands in this log,
  located by the certificate text (a pass certification of a
  state the run then moved away from is not one):
    g0001_20260916T023329Z_32b34e57   lines 8769 to 8821

  the binary that wrote this log:
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7f.x
    md5 c2e9c9990b9f14f1be8cd77abca68945
    [states/g0001_20260916T023329Z_32b34e57/manifest.json `source_identity`, the generation this log's last solve published]

  the last iterate this log records:
    the last outer pass   [line 8682]
      (EXHALE_main) outer pass 40: hydro info=0, worst gated species row  2.97E-03 of  1.0E-05 at cell 500 (carrier balance H2), mass  1.51E-08, momentum  4.74E-13, energy  7.16E-09, omega 0.125, trust 1.3E-03,   98.53 s
    the last inner solve verdict   [line 8655]
      (JFNK) done info=0 ||R||=  1.369E-08  flux spread= 4.879E-12  non-monotone accepts=2

  the last lines of the file:
    8819 |  (certification) the state written is a relaxation snapshot: the run made no stationary claim, so
    8820 |    it is written certified=F with the reason "no stationary claim" and the run exits 0.
    8821 | Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL

the publication test, MEASURED here
  Hydro_ioniz.txt
    output/ md5       fd882d8918052069f8ed229dd376f8c4   (written 2026-09-20T14:37:34)
    generation md5    fd882d8918052069f8ed229dd376f8c4
  Ion_species.txt
    output/ md5       87288b03ad88a92b7b8b8b6c2b1bdd66   (written 2026-09-20T14:37:34)
    generation md5    87288b03ad88a92b7b8b8b6c2b1bdd66
  verdict             the pair in `output/` is BITWISE the stored generation

  *** The last solve of the case `run.log` ended in the published generation
  *** g0001_20260916T023329Z_32b34e57, whose certificate stands at lines 8769
  *** to 8821 of that log. Nothing in it is a log-only iterate.

7. FIELDS THIS RECORD COULD NOT FILL
------------------------------------
none: every field of parts 1 to 6 was found.
```

### 3.5 `molecular_photochem_gj1132_kzzprofile/HeH9`

```
==============================================================================
STAGE A0 IDENTITY RECORD
==============================================================================
case            molecular_photochem_gj1132_kzzprofile/HeH9
case directory  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9
generation      g0005_20260920T053700Z_f59de41d
selected by     the reference `latest_complete` of state_index.json
latest_complete g0005_20260920T053700Z_f59de41d   (state_index.json)
latest_certified (none)  (state_index.json)
                the index names NO certified generation: no state
                of this case carries a certificate that covers its
                own bytes.
assembled       2026-09-21T08:14:57 by models/identity_record.py
classification  HISTORICAL AUDIT (PLAN_20260920_rev9 section 3). Every
                stored value below is READ; every md5 this tool
                recomputed is MEASURED. No binary was run and no
                state was evaluated, so nothing here is Mode R.

1. THE STATE PAIR
-----------------
directory       /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9/states/g0005_20260920T053700Z_f59de41d

  Hydro_ioniz.txt
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9/states/g0005_20260920T053700Z_f59de41d/Hydro_ioniz.txt
    bytes             93894
    md5 MEASURED      348a6ab88a84b1b967d7582f14f9a723
    md5 READ          348a6ab88a84b1b967d7582f14f9a723   [states/g0005_20260920T053700Z_f59de41d/manifest.json `components`]
    verdict           matches the recorded md5

  Ion_species.txt
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9/states/g0005_20260920T053700Z_f59de41d/Ion_species.txt
    bytes             499875
    md5 MEASURED      e39c8734708d77894f578fe1df25283c
    md5 READ          e39c8734708d77894f578fe1df25283c   [states/g0005_20260920T053700Z_f59de41d/manifest.json `components`]
    verdict           matches the recorded md5

  the header of Hydro_ioniz.txt, which the md5 covers:
    # EXHALE schema 2
    # columns r[Rp] rho[mH/cm3] v[cm/s] p[cgs] T[K] heat[erg/cm3/s] cool[erg/cm3/s]
    # rows 504: 2 ghost cells at each end; physical cells are rows 3 to 502
    # coupling: sec_ion=T sec_ion_step=0 recon=WENO3 certified=F cert_reason=no_stationary_claim mode=init
    # provenance: git=3c73905ca8a7 tree=dirty run=2026-09-20T14:37:00
    # provenance: ck_input=4570321104784328650 ck_base=-1 ck_metals=-1
    # provenance: recon=WENO3 base_bc=characteristic carrier=transported restart_schema=3 resid_def=145 N=500
    # boundary_model characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3
    # boundary_reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp] 9.5014933560114534E-001 1.0000000000000000E+000 2.5853531227461574E-001 1.0000000000000000E+000
    # restart_schema 1
    # reservoir He/H 9.0102683952787981E+00 C/H 2.7780253208852597E-04 O/H 8.8626356925852883E-07 N/H 8.1907420604337587E-05
    # species_columns 38 r HI HII HeI HeII HeIII HeITR CI CII CIII OI OII OIII NI NII NIII MgI MgII MgIII SiI SiII SiIII CaI CaII CaIII NaI NaII KI KII SI SII FeI FeII FeIII H2 H2p H3p HeHp
    # grid N 500 R0[cm] 1.1593574793765156E+09 r_min[Rp] 1.0001933187778382E+00 r_max[Rp] 2.9031224061195065E+01 mode Mixed
    # constants set IAU2015+CODATA2018 RJ[cm] 7.1492000000000000E+09 kB[erg/K] 1.3806490000000000E-16
    # options He23S=T metals=T eos_metals=T mol=T molbase=T oxychem=F carrier=T carrier_newton=F iontrans=F he_diff=T he_metal_diff=F sec_ion=T caloric_mono=F excH=F base_ir=F mol_ir=F mol_heat=T visc=F cond=F jlya=0 wellbal=T
    # t_phys[s] 0.0000000000000000E+00
    # source git=3c73905ca8a7 tree=dirty restart_input=provenance_unknown

2. ORIGIN
---------
generation id       g0005_20260920T053700Z_f59de41d   [states/<gen>/manifest.json]
iteration_phase     evaluate   [manifest.json]
published_at        2026-09-20T14:37:02   [manifest.json]
state_written       2026-09-20T14:37:00   [manifest.json]
published_by        models/publish_state.py on lart4   [manifest.json]
source_directory    output   [manifest.json]
run_id              r20260920T053659Z_37074   [manifest.json]
certification       NOT CERTIFIED   [manifest.json `certification.status`]
state claim         certified=False cert_reason=no_stationary_claim mode=init   [manifest]

parent              g0003_20260917T222518Z_c55563e3
  source            (none)
  compatibility     stated by the run that took it

seed                [states/g0005_20260920T053700Z_f59de41d/manifest.json `seed`]
  established       True
  kind              internal continuation
  case              molecular_photochem_gj1132_kzzprofile/HeH9
  generation        g0003_20260917T222518Z_c55563e3
  directory         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9/states/g0003_20260917T222518Z_c55563e3
  chosen by         the runner: this pass was started from the generation named here, a generation of this same case
  transformation    used as it stands: the pair was reloaded, with no mapping and no conversion
  md5 Hydro_ioniz.txtd5816d608561f81814e40fbc1abc47f3
  md5 Ion_species.txt953b49d34d6da140015ac8bce41248c8

ending class        evaluate   [manifest.json `ending.class`]
ending reason       the generation g0003_20260917T222518Z_c55563e3 handed back and measured as it stands; an evaluation has no solve ending
index ending_class  evaluate   [state_index.json]
solver verdict      info = (none), ||R|| = (none)   [manifest.json `ending.solver`]

the case `ENDING` file (/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9/ENDING, written (none)):
    (the case records no ENDING file)

3. THE OPERATOR
---------------
binary path         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_75d55d9d.x   [manifest.json `source_identity`]
  on disk now       yes
  md5 READ          75d55d9d4fd0e748cd01d6e35713e34c   [manifest.json]
  md5 MEASURED      75d55d9d4fd0e748cd01d6e35713e34c
  source manifest   BINARY_MANIFEST_75d55d9d4fd0.txt
  its md5 READ      e8b5b83f3c8c82da8aae93000a25e75e
  its md5 MEASURED  e8b5b83f3c8c82da8aae93000a25e75e

model options       He23S=T metals=T eos_metals=T mol=T molbase=T oxychem=F carrier=T carrier_newton=F iontrans=F he_diff=T he_metal_diff=F sec_ion=T caloric_mono=F excH=F base_ir=F mol_ir=F mol_heat=T visc=F cond=F jlya=0 wellbal=T   [manifest.json `model_identity.options`]
reservoir           He/H 9.0102683952787981E+00 C/H 2.7780253208852597E-04 O/H 8.8626356925852883E-07 N/H 8.1907420604337587E-05

boundary model of the STORED state:
  state header      characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3
                    [# boundary_model of Hydro_ioniz.txt, part of the hashed bytes]
  manifest          characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3
  boundary reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp] 9.5014933560114534E-001 1.0000000000000000E+000 2.5853531227461574E-001 1.0000000000000000E+000

4. THE CONFIGURATION
--------------------
  input.inp                      matches the recorded md5
    md5 MEASURED      c0fd159ac7a8fbf52c7ba271b3209d9e
    md5 READ          c0fd159ac7a8fbf52c7ba271b3209d9e   [manifest.json `configuration_identity`]
  base.inp                       absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  metals.inp                     absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  opacity.inp                    absent in this case, and the manifest records none: the case does not use this file
    md5 MEASURED      (none)
    md5 READ          (none)   [manifest.json `configuration_identity`]
  lower_atmosphere_profile.dat   matches the recorded md5
    md5 MEASURED      8b3c9aea5834d51e836571bf082c5050
    md5 READ          8b3c9aea5834d51e836571bf082c5050   [manifest.json `configuration_identity`]

  spectrum
    named in input    ../../../sed/lhs1140_sed_gj1132_at_b.txt   [input.inp `Spectrum file:`]
    resolved to       /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt
    resolution rule   relative to the case directory, which is the working directory of a run: the binary reads ./input.inp and opens the string as written
    md5 MEASURED      ac73a75c6396793435ed72e1cf34c427
    md5 READ          ac73a75c6396793435ed72e1cf34c427   [manifest.json]
    verdict           matches the recorded md5

  EXHALE_resolved.out
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9/EXHALE_resolved.out
    md5 MEASURED      0c8626a1202a1cd2091e82f6eb3501b1
    written           2026-09-18T07:19:44
    keys              47
    route keys it DOES carry:
      well_balanced              T
      carrier_transport          T
      carrier_in_newton          F
      carrier_newton_on_stall    F
      ionization_transport       F
      oxygen_chemistry           F
    route keys it does NOT carry, taken from elsewhere:
      - the boundary model identity (the state file header and the manifest carry it; the resolved writer never emits it, write_setup_report.f90:906-1120)
      - the numerical flux and the reconstruction method (stated in prose by EXHALE_setup.out)
      - the residual assembly selector EXHALE_RESID_QUAD, which is an environment variable read at hydrodynamic_rows.f90:177 and is in no resolved record

  EXHALE_setup.out (the resolved route, in prose)
    path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9/EXHALE_setup.out
    md5 MEASURED      c68717396550b1441e2f8152d5c00ec5
    written           2026-09-18T07:19:44
    numerical flux    HLLC   [line 54]
    reconstruction method PLM   [line 55]
    well balanced     on   [line 56]
    base boundary     characteristic condition at the face r_edg(0)   [line 61]

5. THE ROW REGISTRY
-------------------
certificate         /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9/states/g0005_20260920T053700Z_f59de41d/certification.txt
  md5 MEASURED      7224ee1346769e78efdb3642f14a87bc
  active equations 7 of 30 in the inventory   [line 2]

  rows evaluated (7):
    hydrodynamic mass row                         ABOVE  max= 1.041E+00 at cell 2 (r = 1.00039), tol= 3.4E-10   [line 14]
    hydrodynamic momentum row                     ABOVE  max= 1.450E-03 at cell 1 (r = 1.00019), tol= 1.0E-08   [line 17]
    hydrodynamic energy row                       ABOVE  max= 1.021E+00 at cell 1 (r = 1.00019), tol= 1.0E-06   [line 19]
    carrier balance H2                            ABOVE  max= 1.000E+00 at cell 432 (r = 9.56781), tol= 1.0E-05   [line 21]
    elemental transport He/H partition            ABOVE  max= 9.903E-01 at cell 2 (r = 1.00039), tol= 1.0E-05   [line 24]
    level balance He 2^3S                         within max= 4.551E-19 at cell 354 (r = 3.19896), tol= 1.0E-06   [line 27]
    eliminated-species closure System_HeH_mol_metal within max= 3.258E-14 at cell 354 (r = 3.19896), tol= 1.0E-06   [line 29]

  verdict           NOT CERTIFIED: 5 entry/entries of the inventory refuse it   [line 53]

  refusing rows (5):

    hydrodynamic mass row: row measure  1.041E+00 above  3.4E-10 at cell 2
      row max         1.041E+00 at cell 2 (r = 1.00039, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       3.4E-10
      volume measure  6.449E-02
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      verdict cell    2 (r = 1.00039E+00 READ, 1.00039 MEASURED): 1.041E+00 against 3.4E-10, distance 3.102E+09
                      [line 16]

    hydrodynamic momentum row: row measure  1.450E-03 above  1.0E-08 at cell 1
      row max         1.450E-03 at cell 1 (r = 1.00019, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-08
      volume measure  6.009E-05
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      the certificate prints no tolerance-normalized verdict
      cell for this row (see part 7)

    hydrodynamic energy row: row measure  1.021E+00 above  1.0E-06 at cell 1
      row max         1.021E+00 at cell 1 (r = 1.00019, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-06
      volume measure  3.075E-01
      scale           residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
      the certificate prints no tolerance-normalized verdict
      cell for this row (see part 7)

    carrier balance H2: gated row measure  1.000E+00 above  1.0E-05 at cell 432 (a wind cell)
      row max         1.000E+00 at cell 432 (r = 9.56781, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-05
      volume measure  1.836E-01
      scale           row_terms [cm^-3 s^-1], the sum of the row's own terms with the 1e-20 free-e
      gated at r >= 1.200E+00 against 1.00E-05: 1.000E+00 at cell 432 (r = 9.56781 MEASURED)
                      [line 23]

    elemental transport He/H partition: gated row measure  1.494E-04 above  1.0E-05 at cell 217 (a wind cell)
      row max         9.903E-01 at cell 2 (r = 1.00039, MEASURED by indexing Hydro_ioniz.txt)
      tolerance       1.0E-05
      volume measure  2.307E-02
      scale           sum of the row's own terms [g cm^-3 s^-1]; floor 1e-20 rho X_base/(R0/v0)
      gated at r >= 1.200E+00 against 1.00E-05: 1.494E-04 at cell 217 (r = 1.20069 MEASURED)
                      [line 26]

6. THE LOG BLOCK
----------------
the log the manifest names for this generation
  path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9/runs/r20260920T053659Z_37074/pp.log
  named in          states/g0005_20260920T053700Z_f59de41d/manifest.json `certification.source`
  on disk now       yes
  md5 MEASURED      0dde68001cc3be7ae094668776d50579
  md5 READ          0dde68001cc3be7ae094668776d50579
  verdict           matches the recorded md5
  written           2026-09-20T14:37:02

where the certificate of this generation actually stands, located by
its own text and not by a file name:
  /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9/runs/r20260920T053659Z_37074/pp.log
  lines 105 to 204 (100 lines, the whole of certification.txt)

  the block:
     105 |  (certification) work state of the loaded restart (Restart intent: stationary)
     106 |    active equations 7 of 30 in the inventory
     107 |    tolerances (convergence study, contract sections 9 and 10; the hydrodynamic values were
     108 |    re-anchored on the self-consistent JFNK residual after B5a):
     109 |      hydrodynamic mass  3.00E-12, momentum  1.00E-08, energy  1.00E-06
     110 |      the mass value is a FLOOR: where the cell's own rounding of the flux difference stands above it, 10.0 times that floor is the tolerance
     111 |      of that cell (cert_tol_mass_at), and the verdict on the mass row is taken cell by cell
     112 |      closure  1.00E-06, level  1.00E-06; and IN THE WIND, at r >= 1.20E+00, carrier  1.00E-05
     113 |      and elemental transport  1.00E-05; the species rows of the cells below that radius are
     114 |      REPORTED AND DO NOT GATE: the layer r < 1.10E+00, whose element fluxes
     115 |      are not conserved (spread 9.8), and the band above it, whose discretization (4.6e-4 at the
     116 |      binding radius) stands above any tolerance one could set there
     117 |      the run's own "Resid tol" is  4.00E-06, which is the SOLVER's target and not a certification tolerance
     118 |    hydrodynamic mass row                                evaluated       max= 1.041E+00  vol= 6.449E-02  cell=2  tol= 3.4E-10  ABOVE   
     119 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     120 |         verdict at cell 2 (r= 1.00039E+00):  1.041E+00 against  3.4E-10, distance  3.102E+09 -- the rounding anchor of that cell
     121 |    hydrodynamic momentum row                            evaluated       max= 1.450E-03  vol= 6.009E-05  cell=1  tol= 1.0E-08  ABOVE   
     122 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     123 |    hydrodynamic energy row                              evaluated       max= 1.021E+00  vol= 3.075E-01  cell=1  tol= 1.0E-06  ABOVE   
     124 |         scale: residual_row_scale, the row's largest term; floor 1e-300 (numerical only)
     125 |    carrier balance H2                                   evaluated       max= 1.000E+00  vol= 1.836E-01  cell=432  tol= 1.0E-05  ABOVE   
     126 |         scale: row_terms [cm^-3 s^-1], the sum of the row's own terms with the 1e-20 free-e
     127 |         gated at r >= 1.200E+00 against  1.00E-05:  1.000E+00 at cell 432; the cells below that radius are reported only
     128 |    elemental transport He/H partition                   evaluated       max= 9.903E-01  vol= 2.307E-02  cell=2  tol= 1.0E-05  ABOVE   
     129 |         scale: sum of the row's own terms [g cm^-3 s^-1]; floor 1e-20 rho X_base/(R0/v0)
     130 |         gated at r >= 1.200E+00 against  1.00E-05:  1.494E-04 at cell 217; the cells below that radius are reported only
     131 |    level balance He 2^3S                                evaluated       max= 4.551E-19  vol= 4.133E-20  cell=354  tol= 1.0E-06  within  
     132 |         scale: the row's own turnover rate; dimensionless, scale 1
     133 |    eliminated-species closure System_HeH_mol_metal      evaluated       max= 3.258E-14  vol= 2.584E-16  cell=354  tol= 1.0E-06  within  
     134 |         scale: the row's own turnover rate; dimensionless, scale 1
     135 |    face mass flux r_f^2 (rho v)_f of this state (reported beside the rows above; it gates nothing):
     136 |      operator: WENO3, numerical flux: HLLC, well balanced: T
     137 |      in units of the wind-window mean of rho v r^2, which is  5.88765E-07
     138 |        minimum over the faces -1.550777432922951E+01 at face 1
     139 |        maximum over the faces  1.764584401504327E+00 at face 36
     140 |        base face -7.364251565210126E+00, offset from the window mean -8.36425E+00
     141 |        a face value and the mean of a cell-centered product are not the same quantity, so an offset of the size of the
     142 |        wind's own discretization is not a mass leak
     143 |      base contact direction: read from the matched face velocity (the window has no standing; w_rev = 1.0)
     144 |    cells without a chemical root: 0 (acceptance class 4 or 6)
     145 |    mass closure of the composition, max |sum_i f_i A_i - 1| =  1.037E-15 at cell 440, r =  10.85005
     146 |      reported only: the sum is one by the definition of the mass fractions, so this is the arithmetic of the state and no equation of it
     147 |    validity states (B1a section 4):
     148 |      active unvalidated physics: not produced
     149 |      out-of-domain closure activations: 1851
     150 |        H3+ cooling domain records (informational; they do not invalidate):
     151 |          collider below the table: 1083, temperature below the fits: 0, above them: 384, outside the non-LTE range: 384
     152 |      rejected trials with no adopted contribution: not produced
     153 |      unbudgeted accepted corrections: 0
     154 |      specified external reservoirs: not produced
     155 |    attempts (not validity states): 0 energy-update and 0 conduction cell(s) reached the lower bracket end
     156 | Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
     157 |    NOT CERTIFIED: 5 entry/entries of the inventory refuse it
     158 |      hydrodynamic mass row: row measure  1.041E+00 above  3.4E-10 at cell 2
     159 |      hydrodynamic momentum row: row measure  1.450E-03 above  1.0E-08 at cell 1
     160 |      hydrodynamic energy row: row measure  1.021E+00 above  1.0E-06 at cell 1
     161 |      carrier balance H2: gated row measure  1.000E+00 above  1.0E-05 at cell 432 (a wind cell)
     162 |      elemental transport He/H partition: gated row measure  1.494E-04 above  1.0E-05 at cell 217 (a wind cell)
     163 |  
     164 |  
     165 |  (EXHALE_main) the two answers of this evaluation:
     166 |    original claim of the file: certified=F cert_reason=no_stationary_claim
     167 |    original claim: NONE MADE -- the file states no stationary claim, so nothing is reproduced or refused
     168 |    work state verdict: NOT CERTIFIED -- the certification report above names the entries that refuse it
     169 |    the pair written into the state: certified=F cert_reason=no_stationary_claim
     170 |      steps: 0 accepted of 0 attempted
     171 |      initialization / continuation: no physical elapsed time is claimed
     172 | 
     173 |    ATTEMPTED-STEP CONTROLLER
     174 |      outer attempts: 0, of which 0 were refused at the adoption boundary
     175 |      steps that needed at least one retry: 0
     176 |      hydro stage attempts (the nested positivity retry, a DIFFERENT count): 0 in 0 halvings
     177 |      ioniz-eq [marching states]: 3 equilibrium sweep(s)
     178 |      ioniz-eq molecular hybrd1 info: 1 1520, 2 0, 3 0, 4 0, 5 0, 0 0
     179 |      ioniz-eq acceptance: 1512 converged root(s) (max res  2.24E-12), 0 root(s) without solver convergence (max res  0.00E+00), 0 projected/handback root(s) (max res  0.00E+00)
     180 |      ioniz-eq residual decades (<=1e-16 .. >=1e-1) converged:   1480 28 9 1 2 0 0 0 0 0 0 0 0 0 0 0
     181 |      flux correction: 0 positivity repair call(s), 0 face flux(es) dropped to first order, 0 of them in accepted steps
     182 |      coupled source step:   0.00 passes per coupled step on average over 0 of them, worst 0 (step 0)
     183 |      coupled source step: every step reached its fixed point within the pass cap
     184 |      coupled source step: 0 exit(s) on the estimated error, 0 on the last increment
     185 |      chemistry mass closure: worst |sum n_s m_s - rho|/rho = 0.000E+00 at step 0
     186 |      formation reservoir change, worst |du_form| / |dt (heat-cool)| = 0.000E+00
     187 |  (write_output/eq) max |sum(O I levels)/n(O I) - 1| =  2.86E-16
     188 |  (write_cool_breakdown_eq) max |sum(channels)/cool - 1| =  4.64E-16
     189 |  (EXHALE_main.f90) Starting the post processing routine..
     190 |  (post_process_adv) advection correction: 20 of 503 cells kept at the equilibrium ionization.
     191 |  (post_process_adv) energy correction: 0 of 502 cells kept the run temperature (thermal Damkohler > 1); largest Damkohler  1.41E-01 at r = 23.7421 Rp.
     192 |  (post_process_adv) conditional correction: mass row |R_1|/s_1 of the input state up to  1.04E+00 at r =  1.0004 Rp; above  1.00E-02, the fraction a corrected row is accurate to, in 20 of 503 cells, of which 14 the ionization conditions do not already refuse.
     193 |  (post_process_adv) stationarity: the mass row was assembled from the primitive state handed in, which is the run's state re-formed and not its conserved variables.
     194 |  (post_process_adv) energy correction: enthalpy flux of the mass-flux divergence / other terms up to  7.12E+02 at r =  1.0004 Rp; dominant in 17 of 503 cells (reported, not a refusal; level  1.00E+00).
     195 |  (post_process_adv) T solve: 1 of 502 cells fell back to eq T (non-positive or out-of-band root).
     196 |  (write_output/ad) max |sum(O I levels)/n(O I) - 1| =  2.92E-16
     197 |  (EXHALE_main.f90) Post processing routine done.
     198 |   
     199 |  ----- Results -----
     200 |   
     201 |  ---> 2D approximate method: Mdot
     202 |  ---> Log10 of steady-state Mdot =  7.93 g/s
     203 |  (EXHALE_main) Restart intent: stationary evaluate -- the loaded state was measured, the work state and
     204 |    the profiles derived from it were written, and no step and no solve were taken.

run.log of this case
  path              /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH9/run.log
  lines             495
  md5 MEASURED      57cd8c64f866a5fd7a1ce59b36849d9f
  written           2026-09-18T07:25:18
  runs in the file  1 (a run begins at line 2)
  is the log the manifest names     no
  holds this generation's block     NO

  published generations whose certificate stands in this log,
  located by the certificate text (a pass certification of a
  state the run then moved away from is not one):
    g0003_20260917T222518Z_c55563e3   lines 439 to 495

  the binary that wrote this log:
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/EXHALE_3146d11b.x
    md5 3146d11b4090306dcea75bb9718edd22
    [states/g0003_20260917T222518Z_c55563e3/manifest.json `source_identity`, the generation this log's last solve published]

  the last iterate this log records:
    the last outer pass   [line 412]
      (EXHALE_main) outer pass 1: REFUSED -- the element composition relaxation found no admissible advance and its entry composition was restored, so no further pass could differ from this one.
    the last inner solve verdict   [line 373]
      (JFNK) done info=2 ||R||=  3.167E-01  flux spread= 1.406E-02  non-monotone accepts=21  solve=1
    the last residual gate   [line 401]
      (JFNK) gate NOT met: flux spread 1.406E-02 >= 2.000E-05

  the last lines of the file:
     493 |  (EXHALE_main) Restart intent: stationary -- the stationary solve returned info = 1; output written.
     494 |  (certification) the state written is a relaxation snapshot: the run made no stationary claim, so
     495 |    it is written certified=F with the reason "no stationary claim" and the run exits 0.

the publication test, MEASURED here
  Hydro_ioniz.txt
    output/ md5       348a6ab88a84b1b967d7582f14f9a723   (written 2026-09-20T14:37:00)
    generation md5    348a6ab88a84b1b967d7582f14f9a723
  Ion_species.txt
    output/ md5       e39c8734708d77894f578fe1df25283c   (written 2026-09-20T14:37:00)
    generation md5    e39c8734708d77894f578fe1df25283c
  verdict             the pair in `output/` is BITWISE the stored generation

  *** The last solve of the case `run.log` ended in the published generation
  *** g0003_20260917T222518Z_c55563e3, whose certificate stands at lines 439
  *** to 495 of that log. Nothing in it is a log-only iterate.

7. FIELDS THIS RECORD COULD NOT FILL
------------------------------------
 1. the tolerance-normalized verdict of the row 'hydrodynamic momentum row'
    the certificate states the row measure and its tolerance but prints
    neither a `verdict at cell` line nor a `gated at r >=` line for it
 2. the tolerance-normalized verdict of the row 'hydrodynamic energy row'
    the certificate states the row measure and its tolerance but prints
    neither a `verdict at cell` line nor a `gated at r >=` line for it
```

## 4. Every field the tool could not fill, and why

### `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7`

- **the tolerance-normalized verdict of the row 'hydrodynamic momentum row'.** The certificate states the row measure and its tolerance but prints neither a `verdict at cell` line nor a `gated at r >=` line for it
- **the tolerance-normalized verdict of the row 'hydrodynamic energy row'.** The certificate states the row measure and its tolerance but prints neither a `verdict at cell` line nor a `gated at r >=` line for it
- **a recorded md5 for the log /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7/runs/r20260919T004843Z_15406/pp.log.** The manifest records no `source_md5`, so the file on disk cannot be shown to be the file the certificate was read from

### `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7`

- **the tolerance-normalized verdict of the row 'hydrodynamic momentum row'.** The certificate states the row measure and its tolerance but prints neither a `verdict at cell` line nor a `gated at r >=` line for it
- **the tolerance-normalized verdict of the row 'hydrodynamic energy row'.** The certificate states the row measure and its tolerance but prints neither a `verdict at cell` line nor a `gated at r >=` line for it
- **a recorded md5 for the log /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/run.log.** The manifest records no `source_md5`, so the file on disk cannot be shown to be the file the certificate was read from
- **the publication test for g0004_20260919T094145Z_259fe9c3.** `output/` of this case does not hold both Hydro_ioniz.txt and Ion_species.txt
- **The log the manifest names is no longer where the certificate
  stands.** The manifest names `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/run.log`; the block is in `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/attempt_1/run.log` at lines 2318 to 2368. The field is filled, and by a different file from the one the record points at.

### `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13`

- **the tolerance-normalized verdict of the row 'hydrodynamic momentum row'.** The certificate states the row measure and its tolerance but prints neither a `verdict at cell` line nor a `gated at r >=` line for it
- **the tolerance-normalized verdict of the row 'hydrodynamic energy row'.** The certificate states the row measure and its tolerance but prints neither a `verdict at cell` line nor a `gated at r >=` line for it
- **a recorded md5 for the log /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13/runs/r20260919T004834Z_15338/pp.log.** The manifest records no `source_md5`, so the file on disk cannot be shown to be the file the certificate was read from

### `molecular_scalar_gj1132_kzz1e9/HeH0.083`

None: every field of parts 1 to 6 was found.

### `molecular_photochem_gj1132_kzzprofile/HeH9`

- **the tolerance-normalized verdict of the row 'hydrodynamic momentum row'.** The certificate states the row measure and its tolerance but prints neither a `verdict at cell` line nor a `gated at r >=` line for it
- **the tolerance-normalized verdict of the row 'hydrodynamic energy row'.** The certificate states the row measure and its tolerance but prints neither a `verdict at cell` line nor a `gated at r >=` line for it

### What the five have in common

Three kinds of gap account for every entry above.

1. **The certificate prints a tolerance-normalized verdict cell only for
   the mass row.** The momentum and energy rows carry a row measure, a
   cell and a tolerance, and no `verdict at cell` line, because the
   tolerance anchor `cert_tol_mass_at`, whose verdict is taken cell by
   cell, exists for the mass row alone. The gated species rows carry a
   `gated at r >=` line instead,
   which the records report in that place. So for the momentum and energy
   rows the distance in units of the tolerance is not in the stored
   record and would have to be computed from the two numbers that are.
2. **A manifest written before the log md5 was carried records
   `source_md5: null`.** Where that is so, the log on disk cannot be shown
   to be the log the certificate was read from. The two molecular
   manifests do record it and both MEASURED equal. The three atomic ones
   do not.
3. **A case whose `output/` holds no state pair cannot take the
   publication test.** That is
   `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7`, whose `output/`
   holds only the `_IC` halves and an element flux profile.

## 5. What these records show that rev9 does not say

### 5.1 `EXHALE_resolved.out` does carry the derived carrier keys

The brief for this work, following rev9 section 3, states that the
machine-readable resolved record does not carry the derived route keys
`carrier_transport` and `Coupled carrier solve`. It does.
`write_setup_report.f90:993-995` writes `carrier_transport`,
`carrier_in_newton` and `carrier_newton_on_stall` into
`EXHALE_resolved.out` (INSPECTED), and all three stand in the file of
every one of the five cases (MEASURED by reading them). The derivation at
`input_read.f90:2193`,
`if (.not. carrier_transport_stated) carrier_transport = thereis_oxychem`,
is therefore VISIBLE in the resolved record: a case that transports
carriers although its input never states it says so there.

What the resolved record genuinely does not carry, and where each is taken
from instead:

- the **boundary model identity**: the state file header carries it on its
  own `# boundary_model` line, which the md5 of the state covers, and the
  manifest repeats it in `model_identity.boundary_model`;
- the **numerical flux** and the **reconstruction method**: stated in prose
  by `EXHALE_setup.out`;
- the **residual assembly selector** `EXHALE_RESID_QUAD`: an environment
  variable read at `hydrodynamic_rows.f90:177` and written into no
  configuration record at all, so a run that used it leaves no trace of
  that choice in any file beside the case.

### 5.2 The named certification source of the photochem checkpoint

`atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/states/g0004_20260919T094145Z_259fe9c3/manifest.json` names the case-level
`run.log` as `certification.source` and records no `source_md5`. That file
has since been rewritten by the stopped run of 2026-09-20 08:01:52, later
than the generation's `published_at` of 2026-09-19T18:41:50. The 51 lines
of the stored certificate stand in `attempt_1/run.log` at lines 2318 to
2368 (MEASURED by matching the text), not in the named file. A record that
followed the manifest path would read a different run's log.

### 5.3 A `certification.txt` that identifies nothing

Two generations of the molecular cases, `g0003_20260920T052722Z_a562871f`
of `molecular_scalar_gj1132_kzz1e9/HeH0.083` and
`g0004_20260920T052648Z_bae406ce` of
`molecular_photochem_gj1132_kzzprofile/HeH9`, carry a `certification.txt`
of exactly two lines: the sentence a relaxation snapshot is written with.
It states no inventory and no row, and that same sentence ends every such
run in the tree, so it cannot identify one block of one log. The tool
refuses to search for it rather than reporting a match.

### 5.4 A log never names its own executable

The identity of the run behind a `run.log` is nowhere in the log. Where
that run published a generation, the identity is in that generation's
manifest. Where it published nothing, the only record is `not_solved.md`,
which `models/run_case.sh` writes beside the case with the binary, the
host, the thread count and the seed. That is how the binary `75d55d9d` of
the three log-only iterates of table 2.2 is established, and it is the
only place it is written down. Note that `not_solved.md` is also used for
hand-written memos in the two molecular cases, where it carries prose and
no identity table.

