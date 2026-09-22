# What Cloudy c25.00 does with its H2 data, and what of it EXHALE should take

Item L7g, read-only study. Nothing was
built, nothing was run, no EXHALE source file was edited, and nothing in
`~/CLOUDY/c25.00/` was touched. Cloudy is third-party code under its own license
(`/home/kiseon/CLOUDY/c25.00/license.txt`); it is quoted below only in short
fragments with file and line, for reference, and no line of it is copied into
EXHALE.

Every claim about what Cloudy does carries the file and line of the routine that
does it. Every number is labeled MEASURED (computed here) or READ (taken from a
source file, a data file header or a document). A statement about a published
paper that Cloudy cites is labeled as Cloudy's citation unless the paper itself
sits in `../references/` and was read.

---

## 0. The one-paragraph answer

Cloudy has no formation heat fraction and no quenching fraction. It solves the
ro-vibrational level populations of the X state in full statistical equilibrium,
puts newly formed molecules into named (v, J) levels through normalized
distributions, and then obtains the gas heating as the NET of downward minus
upward collisional rates summed over every level pair and every collider, from
the same rate coefficients the level solver uses (`mole_h2.cpp:2211-2265`).
What EXHALE approximates by `D0 * [f_trans + (1 - f_trans) * f_quench]` is in
Cloudy the exact consequence of the level solution. That is the architectural
difference, and it is the reason Cloudy needs a nascent distribution for every
formation channel while EXHALE does not.

The price Cloudy pays is that it has no nascent distribution for most channels
either. Where it lacks one it spreads the source over the ladder in the LOCAL
LTE proportions (`mole_h2.cpp:175-190`), which by construction returns zero net
collisional heat, i.e. the same answer EXHALE's `f_quench = 1` gives in this
layer. Three-body association, H3+ + e and H2+ + e all fall in that class. So
Cloudy does not answer the question the inventory left open; it answers a
different, adjacent set of questions well, and four of its data sets are
directly useful to EXHALE.

---

## 1. Which data files the big H2 model reads, and how

### 1.1 The table

All files live in `~/CLOUDY/c25.00/data/h2/`. "Reader" is the routine that opens
the file; "wired in" is the line that names it.

| data file | contents | reader | file name wired in at | source as the file header states it (READ) |
|---|---|---|---|---|
| `energy_X.dat` | (v, J, E[cm^-1]) of the X ladder; this call also DEFINES the level list | `diatomics::H2_ReadEnergies`, `mole_h2_io.cpp:702-790`, driven from `mole_h2_io.cpp:633-700`, called at `mole_h2_create.cpp:138` | `mole_h2_io.cpp:706` | "Taken from Table 1 of Komasa et al. (2011)" |
| `energy_B.dat`, `energy_C_plus.dat`, `energy_C_minus.dat`, `energy_B_primed.dat`, `energy_D_plus.dat`, `energy_D_minus.dat` | the six excited singlet states | same routine, `nelec = 1..6` | `mole_h2_io.cpp:707-712` | per file |
| `energy_dissoc.dat` | ONE dissociation energy per electronic state, cm^-1 | `H2_ReadDissocEnergies`, `mole_h2_io.cpp:792-850`, called at `mole_h2_create.cpp:135` | `mole_h2_io.cpp:794` | "Sharp, T. E., 1971, Atomic Data, 2, 119"; X value 36118.11 cm^-1 |
| `transprob_X.dat` | quadrupole + magnetic dipole A values inside X | `H2_ReadTransprob`, `mole_h2_io.cpp:403-630`, called at `mole_h2_create.cpp:561` | `mole_h2_io.cpp:407` | "Wolniewicz, L., Simbotin, I., and Dalgarno, A., 1998, ApJS, 115, 293-313" |
| `transprob_B.dat` and the five siblings | electronic band A values | same | same block | per file |
| `dissprob_B.dat`, `dissprob_C_plus.dat`, `dissprob_C_minus.dat`, `dissprob_B_primed.dat`, `dissprob_D_plus.dat`, `dissprob_D_minus.dat` | for each (v, J) of each excited state: the spontaneous DISSOCIATION probability AND the kinetic energy released, in eV | `H2_ReadDissprob`, `mole_h2_io.cpp:853-939`, called at `mole_h2_create.cpp:404` | `mole_h2_io.cpp:856-862` | "Abgrall, H., Roueff, E., & Drira, I. 2000, A&AS, 141, 297-300" |
| `hminus_deposit.dat` | state-resolved rate coefficients for H- + H -> H2(v, J) + e, as -log10, at 7 temperatures 10 to 10000 K | `H2_Read_hminus_distribution`, `mole_h2_io.cpp:941-1035`, called at `mole_h2_create.cpp:378` | `mole_h2_io.cpp:951` | "Launay, J.M., Le Dourneuf, M., & Zeippen, 1991, A&A, 252, 842-852" |
| `coll_rates_H_15.dat` (DEFAULT for H) | de-excitation rate coefficients, 50 temperatures | `H2_CollidRateRead`, `mole_h2_coll.cpp:163-199`, called at `mole_h2_create.cpp:453` | `h2.cpp:169` | "Lique, F., 2015, MNRAS, 453, 810" |
| `coll_rates_H_07.dat`, `coll_rates_H_99.dat` | alternatives for H | same | `parse_atom_h2.cpp:202, 195` | Wrathmall, Gusdorf & Flower 2007; Le Bourlot, Pineau des Forets & Flower 1999 |
| `coll_rates_He_ORNL.dat` (DEFAULT for He) | de-excitation rates, 41 temperatures 1 to 10000 K, covering the WHOLE ladder to v = 14 | same | `h2.cpp:172` | "Lee, T. G. et al. 2007, ApJ, in preparation" |
| `coll_rates_He_LeBourlot.dat` | alternative for He | same | `parse_atom_h2.cpp:242` | Le Bourlot et al. 1999 |
| `coll_rates_H2ortho_ORNL.dat`, `coll_rates_H2para_ORNL.dat` (DEFAULT for the H2 perturber) | PURE ROTATIONAL rates inside v = 0 only | same | `h2.cpp:175, 178` | "Wan et al. 2018, ApJ, 862, 132" |
| `coll_rates_H2ortho_LeBourlot.dat`, `coll_rates_H2para_LeBourlot.dat` | alternatives, ro-vibrational to v = 3 | same | `parse_atom_h2.cpp:276, 303` | Le Bourlot et al. 1999 |
| `coll_rates_Hp.dat` (DEFAULT for H+) | rotational rates inside v = 0, 6 temperatures | same | `h2.cpp:181` | "Gerlich, D., 1990, J. Chem. Phys., 92, 2377-2388" |
| `cont_diss.dat` | continuum photodissociation cross sections | `Read_Mol_Diss_cross_sections`, called at `mole_h2_create.cpp:675`, file named at `mole_dissociate.cpp:18`; only when `mole_global.lgStancil` | | |
| `lte_cooling.dat` | H2 line cooling per molecule in LTE against temperature, for the SMALL model | `H2_Read_LTE_cooling_per_H2`, `mole_h2_io.cpp:1930-1978` | `mole_h2_io.cpp:1932` | generated by Cloudy itself, see section 3.4 |

The X ladder Cloudy builds from `energy_X.dat` has **303 bound levels**
(MEASURED: 303 data records, levels per v of 32, 31, 30, 28, 26, 24, 23, 21,
19, 17, 15, 13, 11, 8, 5 for v = 0 to 14). EXHALE's ladder has 302 (READ,
inventory section 3). The topmost Cloudy level, (14,4) at 36118.0431 cm^-1,
carries the in-file comment "#Energy level might be uncertain" and sits
0.07 cm^-1 below the X dissociation energy of `energy_dissoc.dat`
(36118.11 cm^-1). **No quasi-bound level is carried**: every level of
`energy_X.dat` lies below the dissociation energy, and the string "quasi" does
not occur anywhere in `source/*.cpp` outside an unrelated helium comment
(MEASURED by grep).

### 1.2 How a collision file is read

`H2_CollidRateRead` (`mole_h2_coll.cpp:163-199`) checks a magic number on line 1
(`110416` for every H2 collision file), then hands the open file to the generic
`ReadCollisionRateTable` (`atmdat.cpp:67-180`) together with a callback
`diatomics::GetIndices` (`mole_h2_coll.cpp:201-232`). `ReadCollisionRateTable`
skips comment lines beginning `#` or `!`, reads the first non-comment line as
the temperature row, allocates `collrates[ipHi][ipLo][nTemp]` over ALL level
pairs and zeroes it (`atmdat.cpp:100-112`), then reads one record per
transition. `GetIndices` converts `v' J' v J` into energy-sorted level indices
and returns `(-1, -1)` for any transition whose levels are outside the adopted
model, with the comment (`mole_h2_coll.cpp:208-210`):

> "check that we actually included the levels in the model representation /
> depending on the potential surface, the collision data set / may not agree
> with our adopted model - skip those"

Such records are silently dropped (`atmdat.cpp:125-127`). A transition the file
does not carry simply stays zero.

### 1.3 Ortho and para

Ortho/para is a property of the level, not of the file. It is assigned from the
parity of J plus an electronic-state offset at `mole_h2_create.cpp:187-207`,
and the statistical weight is `g = (2J+1) * (3 if ortho else 1)`
(`mole_h2_create.cpp:210-215`), cited in the code to "Shull, J.M., & Beckwith,
S., 1982, ARAA, 20, 163-188" (Cloudy's citation).

The H2 perturber is split into TWO colliders, ortho and para, with separate data
files and separate densities taken from the solved populations rather than from
an assumed 3:1 ratio (`mole_h2.cpp:216-224`):

```
collider_density[2] = h2.ortho_density_f;
collider_density[3] = h2.para_density_f;
collider_density[4] = dense.xIonDense[ipHYDROGEN][1];
collider_density[4] += (realnum)findspecieslocal("H3+")->den;
```

Note the last line: **H3+ is lumped into the proton collider** and quenches with
the Gerlich (1990) proton rates. The ortho/para ratio is itself a convergence
criterion of the level solution, default tolerance 1e-2 (`mole_h2.cpp:941`,
`mole_h2.cpp:1541`). Collisions that change ortho into para are ON by default
(`h2.cpp:45`, `lgH2_ortho_para_coll_on = true`) and can be killed
(`mole_h2_coll.cpp:124-127`).

### 1.4 Temperature interpolation of the collision rates

`InterpCollRate` (`atmdat.cpp:188-215`). It is **linear in T and linear in the
rate coefficient, not logarithmic**, through `linint`, and it is **clamped at
both ends**: below the first tabulated temperature it returns the first column,
above the last it returns the last column (`atmdat.cpp:194-201`). There is no
extrapolation anywhere. A commented-out block at `mole_h2_coll.cpp:70-87`
records an earlier exponential extension below 100 K, now disabled behind
`fixit("test that this does not matter with new rate table interpolation")`.

EXHALE's `h2_vibrational_relaxation.f90` holds its three tables at their end
values as well (READ, the module's VALIDITY block), so on the clamping the two
codes agree; on the interpolation EXHALE should state which it uses, because
linear-in-T interpolation of a rate that varies by four orders across a decade
of temperature is a real error at the low end of a coarse grid such as the
9-point Le Bourlot one.

### 1.5 How Cloudy fills a transition the file lacks: the g-bar rule

`H2_CollidRateEvalOne` (`mole_h2_coll.cpp:96-132`). If the interpolated rate is
exactly zero, and `lgColl_gbar` is on (DEFAULT true, `h2.cpp:31`), and the
transition does not change ortho into para, the rate is replaced by
`GbarRateCoeff` (`mole_h2_coll.cpp:133-158`):

```
// these are fits to the existing collision data
// used to create g-bar rates
...
// the fit is log(K)=y_0+a*((x)^b), where K is the rate coefficient,
// and x is the energy in wavenumbers
```

with one triple of coefficients per collider and the energy difference floored
at 100 cm^-1. **It has no temperature dependence at all.** The stated source is
Cloudy's own fit to its own collision data, not a publication; no `>>refer` line
accompanies it.

MEASURED, arithmetic on the coefficients read at `mole_h2_coll.cpp:137-144`,
evaluated at the (1,0) to (0,0) spacing 4161.17 cm^-1 (which is EXHALE's v = 1
origin of 0.5159 eV, READ from the inventory section 3):

| collider | g-bar rate [cm^3 s^-1], any T |
|---|---|
| H | 2.445e-15 |
| He | 6.073e-17 |
| ortho-H2 | 3.167e-17 |
| para-H2 | 1.333e-17 |
| H+ | 1.418e-10 |

**This matters for EXHALE's open H2-H2 question.** Cloudy's DEFAULT H2 perturber
files (Wan et al. 2018) are purely rotational inside v = 0, so in a default
Cloudy run **every vibrational H2-H2 quenching rate is the temperature-independent
g-bar value above, 1.3e-17 to 3.2e-17 cm^3 s^-1**. Against the Le Bourlot
thermal v = 1 to v' = 0 rate the inventory measured, 9.47e-16 at 808 K and
3.19e-14 at 1527 K (READ, inventory section 2.4), the g-bar is 70 to 2400 times
SMALLER. So Cloudy's default is not a usable check on EXHALE's H2-H2
coefficient, and the Le Bourlot files, which Cloudy ships and EXHALE already
uses, remain the best available.

### 1.6 How the collider data set is chosen, and why

Defaults are set in `diatoms_init` (`h2.cpp:156-199`) and overridden by
`DATABASE H2` commands parsed in `parse_atom_h2.cpp:170-310`. The defaults are
H = Lique 2015, He = ORNL, H2 ortho and para = ORNL, H+ = Gerlich 1990. The
source carries no comment giving a physical reason for any of these choices; the
only justification in the file is the label "(the default)" and, for the ORNL
sets, "Teck Lee et al. ApJ to be submitted" (`parse_atom_h2.cpp:227-228`,
`258-260`, `285-287`). This is worth saying plainly: **Cloudy's default He and
H2 sets are chosen over published alternatives without a stated reason.**

### 1.7 Levels with no data, including (14,4)

There is no repair rule for a level. A level with no collisional row simply has
zero collisional rates except through the g-bar rule of section 1.5, which fires
on every pair whose tabulated rate is zero and which therefore DOES connect
(14,4) and every other data-less level collisionally, at a temperature-independent
rate. Radiatively, `transprob_X.dat` carries no downward transition from (14,4)
at all (READ, inventory section 3, which measured 4661 X lines over 301 levels
with (0,0), (0,1) and (14,4) unconnected). In the level solver that level is
therefore held up only by its own source and sink terms and by g-bar collisions.

Cloudy's own data file flags the level as doubtful ("#Energy level might be
uncertain", `energy_X.dat` last record), and it lies 0.07 cm^-1 below
dissociation. **The judgment this supports for EXHALE: (14,4) is not a hole to
be filled, it is a level whose binding energy is smaller than the uncertainty on
it.** Dropping it, which is what Nesterenok et al. (2019) and Wolniewicz et al.
(1998) effectively do (READ, inventory section 3), is defensible and cheaper
than inventing a rate for it.

### 1.8 The He ORNL set and its publication status

The file header reads, verbatim and complete (READ,
`data/h2/coll_rates_He_ORNL.dat` line 3):

```
# >>refer   H2   collision    Lee, T. G. et al. 2007, ApJ, in preparation
```

The source comments repeat it as "Teck Lee et al. ApJ to be submitted"
(`parse_atom_h2.cpp:227-228`). MEASURED by ADS query, no paper matches that
description. Three published papers by the same group cover the same physics and
are plausible relatives, but none of them is stated anywhere in Cloudy to be the
origin of this file, and none of their abstracts says it distributes a
He + H2 rate set reaching v = 14:

- Lee, Balakrishnan, Forrey, Stancil, Shaw & Schultz (2008), ApJ 689, 1105,
  "Rotational Quenching Rate Coefficients for H2 in Collisions with H2 from 2 to
  10,000 K" (H2-H2, rotational).
- Mack, Clark, Forrey, Balakrishnan, Lee & Stancil (2006), PhRvA 74, 052718,
  "Cold He+H2 collisions near dissociation".
- Ohlinger, Forrey, Lee & Stancil (2007), PhRvA 76, 042712, "H2 dissociation due
  to collisions with He".

The inventory already carries Paolini, Ohlinger & Forrey (2011) from the same
group. **No number from `coll_rates_He_ORNL.dat` is taken here either**, and the
request list of section 5.4 keeps the inventory's request item open with these
three bibcodes attached as leads.

---

## 2. H2 formation heating and the nascent distribution

### 2.1 There is no formation heating term

MEASURED by grep over `source/*.cpp`: Cloudy has no variable, term or heating
slot for the energy released by H2 formation into the gas. The only place the
4.48 eV appears as an energy deposit is **into the grain**, not into the gas
(`grains.cpp:4229-4231`):

```
/* energy deposited into grain by formation of a single H2 molecule, in eV,
 * >>refer	grain	physics	Takahashi J., Uehara H., 2001, ApJ, 561, 843 */
const double H2_FORMATION_GRAIN_HEATING[H2_TOP] = { 0.20, 0.4, 1.72 };
```

used at `grains.cpp:4483-4489` for the three grain materials. The rest of the
binding energy reaches the gas ONLY as internal excitation of the molecule,
through the nascent (v, J) distribution, and is converted to gas heat (or
radiated) by the level solution. This is the architectural point of section 0.

### 2.2 The three formation channels that DO carry a nascent distribution

All three are assembled in `mole_H2_form` (`mole_h2_form.cpp:14-187`), which
fills `H2_X_formation[v][J]` in cm^-3 s^-1.

**(a) Grain formation**, `mole_h2_form.cpp:43-79`, using the normalized
distribution `H2_X_grain_formation_distribution[grain type][v][J]` built once in
`mole_h2_create.cpp:675-870`. Four options, selected by `hmi.chGrainFormPump`,
DEFAULT `'T'` (`hmi.cpp:102`):

| key | distribution | source as the code cites it |
|---|---|---|
| `'T'` (default) | Takahashi 2001: a vibrational distribution `H2_vib_dist` plus, for each v, a rotational distribution that is the mean of a Gaussian (their eq. 6) and a thermal form (their eq. 7), weighted by `(2J+1)` and the nuclear spin factor | "Takahashi, Junko, 2001, ApJ, 561, 254-263" and "Takahashi, J., & Uehara, H., 2001, ApJ, 561, 843-857" (`mole_h2_create.cpp:735-736`, `mole_h2_create.cpp:20-22`) |
| `'D'` | `g * (1+v) * exp(-E/T_form)` with `T_H2_FORM = 50000 K` | "Draine, B.T., & Bertoldi, F., 1996, ApJ, 468, 269-289", their eq. 19 (`mole_h2_create.cpp:691-694`) |
| `'t'` | a thermal distribution at `T_H2_FORM = 17329 K`, i.e. 1.5 eV | "Le Bourlot, J, 1991, A&A, 242, 235" (`mole_h2_create.cpp:830-834`) |
| `' '` | no formation pumping: everything into (0,0) and (0,1) | (`mole_h2_create.cpp:866-880`) |

The Takahashi branch carries the full energy partition explicitly
(`mole_h2_create.cpp:24-25` and `723-724`):

```
static double XVIB[H2_TOP] = { 0.70 , 0.60 , 0.20 };
static double Xdust[H2_TOP] = { 0.04 , 0.10 , 0.40 };
...
double Xrot[H2_TOP] = { 0.14 , 0.15 , 0.15 };
double Xtrans[H2_TOP] = { 0.12 , 0.15 , 0.25 };
```

three grain materials, four shares. `EH2_eval` (`mole_h2_create.cpp:31-52`)
subtracts the grain share and returns the energy left to the molecule; `Erot` is
then the rotational share of what remains after the vibrational level is fixed
(`mole_h2_create.cpp:751`):

```
Erot = (EH2 - Ev) * Xrot[ipH2] / (Xrot[ipH2] + Xtrans[ipH2]);
```

**The translational share Xtrans is used ONLY to set that ratio.** MEASURED by
grep: `Xtrans` occurs at `mole_h2_create.cpp:724` and `751` and nowhere else in
the whole source. **Cloudy computes the translational share of the formation
energy and then throws it away; it is never deposited as gas heat.** For grain
formation at Xtrans = 0.12 to 0.25 that is a real omission, and EXHALE's
`f_trans` treatment, which deposits the prompt share, is the better of the two
on this point. (EXHALE has no grains, so this is a statement about Cloudy, not a
channel EXHALE must reproduce.)

**(b) The H- route**, `mole_h2_form.cpp:85-160`. The distribution is the
Launay, Le Dourneuf & Zeippen (1991) state-resolved set of
`hminus_deposit.dat`, read as -log10 rate coefficients at seven temperatures
(`mole_h2_io.cpp:1008-1012`) and interpolated LINEARLY IN log10(T) between
bracketing points, clamped below 10 K and above 10000 K
(`mole_h2_form.cpp:90-113`, `131-133`). The interpolation is of the
DISTRIBUTION, normalized to unity, and the assertion at
`mole_h2_form.cpp:160` enforces that normalization to 1e-4. The same
distribution is used for the BACK reaction H2(v, J) + e -> H + H-
(`mole_h2_form.cpp:140`), which is detailed balance applied state by state.

**(c) H2+ + H -> H2 + H+**, `mole_h2_form.cpp:162-169`, verbatim:

```
/* >>chng 03 feb 10, add this population process */
/* H2+ + H => H2 + H+,
 * >>refer	H2	population	Krstic, P.S., preprint
 * all goes into v=4 but no J information, assume into J = 0 */
```

Everything into (4,0). The cited source is a preprint (Cloudy's citation, not
verified here).

### 2.3 Three-body association: Cloudy has the reaction and NO nascent distribution

The reactions exist (`mole_reactions.cpp:1974-1975`):

```
newreact("H,H,H=>H2,H","bh2g_dis_h",1.,0.,0.); /* back rate, three body recombination, 2H + S => H_2 + S */
newreact("H,H,H2=>H2,H2","bh2g_dis_h2",1.,0.,0.);
```

plus the H2* versions at `mole_reactions.cpp:2041`. The rate is obtained by
**detailed balance from the collisional dissociation rate**
(`mole_reactions.cpp:1343-1345`, `1359-1361`): `bh2g_dis_h2 = rh2g_dis_h2 *
hmi.rel_pop_LTE_H2g`. And the collisional dissociation rate itself, when the big
model supplies it, is built on a made-up coefficient
(`mole_h2_coll.cpp:34-47`, verbatim):

```
/* this is a guess of the collisional dissociation rate coefficient -
 * will be multiplied by the sum of densities of all colliders
 * except H2*/
double energy = H2_DissocEnergies[0] - states[ipHi].energy().WN();
ASSERT( energy > 0. );
/* we made this up - Boltzmann factor times rough coefficient */
H2_coll_dissoc_rate_coef[iVibHi][iRotHi] =
	1e-14f * (realnum)sexp(energy/phycon.te_wn) * lgColl_dissoc_coll;
/* collisions with H2 - pre coefficient changed from 1e-8
 * (from umist) to 1e-11 as per extensive discussion with Phillip Stancil */
H2_coll_dissoc_rate_coef_H2[iVibHi][iRotHi] =
	1e-11f * (realnum)sexp(energy/phycon.te_wn) * lgColl_dissoc_coll;
```

It IS state-resolved in the sense that the Boltzmann factor uses the level's own
binding energy, so the inverse association rate is automatically largest into the
highest levels, which is the qualitative shape Orel (1987), Esposito & Capitelli
(2009) and Paolini et al. (2011) give (READ, inventory section 1.6). But the
prefactor is stated in the code to be invented, so **no number from Cloudy's
three-body channel can be used, and its nascent distribution is not an
independent source.** For the level solver, three-body association is not in
`H2_X_formation` at all (MEASURED: `mole_h2_form.cpp` fills that array from
grains, H- and H2+ + H only), so it reaches the ladder through the residual
route of section 3.2.

### 2.4 The verdict for EXHALE's question

Cloudy does not supply what the inventory asked for. It has no nascent
distribution for gas-phase three-body association, it invents the rate, and it
discards the translational share it computes for the one channel where it does
partition the energy. **The inventory's recommendation of a bounded
two-parameter form with `f_trans = 0` stands, and Cloudy is evidence for it
rather than against it**: the most detailed H2 model in general use also places
newly formed molecules in internal levels and lets the collisional network decide
what becomes heat, which is exactly the physics `f_quench` approximates.

---

## 3. The level-resolved cascade itself

### 3.1 The solver

Full statistical equilibrium over the whole X ladder plus the six excited
electronic states, split into two pieces and iterated to convergence. The driver
is `H2_LevelPops` (`mole_h2.cpp:896-1930`), whose loop is
`mole_h2.cpp:1165-1930`:

1. `SolveExcitedElectronicLevels` (`mole_h2.cpp:1934-2043`), the B, C, B', D
   states.
2. `SolveSomeGroundElectronicLevels` (`mole_h2.cpp:2045-2180`): the levels of X
   ABOVE the matrix block are solved by downward back-substitution, starting at
   the highest level and working down, with the "in" rates taken from the
   previous iteration's populations (`mole_h2.cpp:2087-2092`). Gauss-Seidel,
   in other words.
3. `H2_Level_low_matrix` (`mole_h2.cpp:472-893`): the lowest `nXLevelsMatrix`
   levels are solved directly. **The default is 70** (`h2.cpp:126`, with the
   change history "04 oct 05, make default 30 levels / 04 dec 23, make default
   70 levels"); `DATABASE H2 MATRIX` sets it (`parse_atom_h2.cpp:330`), and a
   negative value means all of X (`mole_h2_create.cpp:247-250`). The matrix is
   assembled with radiative, pumped and collisional rates INCLUDING the exchange
   with levels above the block (`mole_h2.cpp:624-658`, `687-740`) and handed to
   the generic `Atom_LevelN` (`mole_h2.cpp:825-846`), which solves it by LU
   through `solve_system`, falling back to a rearranged Gaussian elimination
   `gthsolve` when that fails (`atom_leveln.cpp:569-596`).
4. Convergence is tested on the populations in both absolute and relative terms,
   on the Solomon dissociation rate, on the heating, and on the ortho/para ratio
   (`mole_h2.cpp:930-941` sets the tolerances; the tests are at
   `mole_h2.cpp:1440-1560`).

So: **a full 303-level statistical equilibrium, not a reduced set, with a
70-level direct block and the rest by iterated back-substitution.** Negative
populations out of the matrix are caught and renormalized with a printed warning
(`mole_h2.cpp:849-864`).

### 3.2 What happens to a formation channel with no state assignment

`H2_X_sink_and_source` (`mole_h2.cpp:52-190`). The explicitly state-resolved
sources of section 2.2 are added first (`mole_h2.cpp:82`), the total source from
the chemistry network is then compared against them, and the REMAINDER is spread
over the ladder in LTE proportions (`mole_h2.cpp:175-190`):

```
H2_X_source[ipHi] += source_left * H2_populations_LTE[iElec][iVib][iRot]*rpop_lte;
```

Sinks are treated the same way, spread uniformly in s^-1 over the levels
(`mole_h2.cpp:126-135`). **This is Cloudy's answer for three-body association,
H3+ + e -> H2 + H, and every other channel with no nascent distribution: deposit
it thermally.** A thermal deposit produces zero net collisional heating by
construction, so in Cloudy those channels contribute nothing to the gas energy
balance through the ladder. That is the same answer as EXHALE's `f_quench = 1`
limit, reached by a different route, and it is a defensible default rather than
a correct one.

### 3.3 How the ladder enters the gas energy balance

`H2_Cooling` (`mole_h2.cpp:2182-2297`) produces two numbers:

- `HeatDiss` (`mole_h2.cpp:2197-2205`): the Solomon-process dissociation heating,
  the sum over ALL states of `Pop * H2_dissprob[e][v][J] * H2_disske[e][v][J]`,
  i.e. the population of each excited-state level times its spontaneous
  dissociation probability times **the kinetic energy that level's dissociation
  releases**, both READ per (v, J) from the `dissprob_*.dat` files of Abgrall,
  Roueff & Drira (2000). `cool_eval.cpp:152-156` describes the result as
  "typically 0.25 - 0.5 eV per event".
- `HeatDexc` (`mole_h2.cpp:2211-2265`): the NET collisional heat of the X ladder.
  For every ordered pair of levels and every one of the five colliders it forms
  the downward rate from `CollRateCoeff` times the collider density, gets the
  upward rate by detailed balance from the same coefficient
  (`mole_h2.cpp:2238-2242`), multiplies each by the transition energy, and
  accumulates `heatone - coolone`. It also accumulates the temperature
  derivative for the thermal solver (`mole_h2.cpp:2249-2251`).

Both are handed to the energy balance in `cool_eval.cpp:148-166`, entered as
heating slot 17 (dissociation) and slot 8 (collisions within X), with the
negative part of the latter entered as the coolant `H2cX`
(`cool_eval.cpp:222-232`). During the initial temperature search `HeatDexc` is
forced to zero, with the reason given at `mole_h2.cpp:2287-2292`: far from
solution it is the small difference of two very large numbers and the noise
prevented convergence in dense cosmic-ray-heated clouds.

**This is the single most transferable idea in the whole model, and section 5
item 1 takes it up.** It replaces both of EXHALE's fractions with one expression
in the quantities EXHALE would already have if it solved the ladder, and it is
automatically correct in both the high-density and low-density limits and
everywhere between, with no `n_cr` and no `f_quench`.

Grain de-excitation of the ladder exists as a separate channel, sending every
(v, J) to (0,0) or (0,1) with spin conserved (`mole_h2.cpp:2110-2125`), DEFAULT
OFF (`h2.cpp:40`, `lgH2_grain_deexcitation = false`).

### 3.4 The small model, and `lte_cooling.dat`

The small model replaces the ladder with a fitted cooling function. The
prescription is Glover & Abel (2008), MNRAS 388, 1627, whose paper IS in
`../references/Glover_2008_MNRAS_388_1627.pdf`. `README-h2-lte.md` states the
scheme (READ, verbatim):

> "The small H2 model computes the H2 line cooling according to the prescription
> of Glover & Abel 2008, MNRAS, 388, 1627.  The total cooling at intermediate
> densities is computed according to their eqn (39), which combines the cooling
> at low densities and at LTE."

In the source: the low-density limit is a sum of fitted cooling functions, one
for each collider and separately for ortho and para H2, over the colliders H,
H2 (all four ortho/para combinations), He, H+ and e-, each a polynomial in log10(T/1000 K)
from their eq. 27 with their Tables 1 to 7 (`cool_eval.cpp:726-1250`). The LTE
limit is `interpolate_LTE_Cooling` (`mole_h2_etc.cpp:497-515`), a linear
interpolation in T over `lte_cooling.dat`. The two are combined at
`cool_eval.cpp:1334` as the harmonic mean `cooling_low * cooling_high /
(cooling_low + cooling_high)`, their eq. 39.

`lte_cooling.dat` is **generated by Cloudy itself**: `LTE_Cooling_per_H2`
(`mole_h2_etc.cpp:438-491`) runs the BIG model's line list with LTE level
populations, sums `Aul * E * n_LTE` over the X band only, and the README gives
the command that writes the table. So the small model's LTE limit is the big
model's own optically thin line cooling, tabulated. That is worth noting because
**it is exactly what EXHALE's `molecular_infrared_cooling.f90` already does**
(READ, inventory section 3: Boltzmann populations from the caloric EOS partition
function, optically thin line sum, tabulated against temperature). EXHALE's IR
cooling is therefore the Cloudy small-model LTE branch, computed in-house from a
better line list, and what EXHALE lacks against Cloudy's small model is the
low-density branch and the harmonic combination.

The small model's dissociation and vibrational heating come from one of four
fitted prescriptions selected by `hmi.chH2_small_model_type`
(`cool_eval.cpp:174-205`): Tielens & Hollenbach 1985, Burton, Hollenbach &
Tielens 1990, Bertoldi & Draine 1996, or Elwert, the default (Cloudy's
citations, not verified here). EXHALE's present `h2_vibrational_heat_fraction`
sits in this family: it is the Hollenbach & McKee / Burton et al. form.

### 3.5 Lyman-Werner pumping and the dissociated fraction

In the big model there is **no self-shielding function**. Every Lyman and Werner
line is transferred individually, with its own optical depth accumulated across
the zone and an escape probability (`mole_h2.cpp:388-470`, with the comment at
`mole_h2.cpp:411-413` that "lines become self-shielding surprisingly quickly").
The dissociated fraction is not a constant either: it is the per-level
spontaneous dissociation probability of the excited state, `dissprob_*.dat`
(Abgrall, Roueff & Drira 2000), applied to the solved excited-state populations,
and the heat it returns is the kinetic energy of that level from the same files
(section 3.3).

EXHALE, by contrast, holds the fragment kinetic energy at a constant 0.4 eV
(READ, `src/modules/lower_atmosphere/lyman_werner.f90` lines 349-351 and 405)
and uses a self-shielding table. The constant is inside Cloudy's stated 0.25 to
0.5 eV range, so it is not wrong; what it lacks is the dependence on which
levels are actually pumped. See section 5 item 2.

---

## 4. Dissociative recombination and the other formation channels

### 4.1 No product states are assigned

MEASURED by reading `mole_reactions.cpp`: the dissociative recombination
reactions are declared with bare products and no state information:

```
mole_reactions.cpp:2001:  newreact("e-,H3+=>H2,H","hmrate",2.5e-8,-0.3,0.);
mole_reactions.cpp:2002:  newreact("e-,H3+=>H,H,H","hmrate",7.5e-8,-0.3,0.);
mole_reactions.cpp:2007:  newreact("H2+,e-=>H,H","hmrate",1.6e-8,-0.43,0.);
```

The H2 produced by `e- + H3+` is the chemical species `H2`, and it reaches the
ladder through the LTE residual route of section 3.2. **Cloudy does not book the
vibrationally hot H2 fragment of H3+ dissociative recombination.** The branching
ratio it uses, 2.5e-8 against 7.5e-8, is 0.25 two-body to 0.75 three-body
(MEASURED, arithmetic on the two declarations), against the 0.70 +/- 0.07
three-body branching the inventory READ from Kokoouline, Greene & Esry (2001)
and against EXHALE's own 0.70/0.30 (READ, inventory section 5). So EXHALE's
branching is closer to the measured value than Cloudy's.

There is also no HeH+ dissociative recombination in the H2 part of the network
(MEASURED by grep for `HeH+` with `e-` in `mole_reactions.cpp`: no match), so
Cloudy says nothing about EXHALE's R16.

### 4.2 And the reaction heats are not booked at all

`t_mole_local::chem_heat` (`mole_reactions.cpp:4403-4500`) computes the network's
chemical heating from the formation enthalpies of reactants minus products,
skipping any reaction with a PHOTON or CRPHOT (on the ground that the photon
accounts for the energy difference, `mole_reactions.cpp:4430-4440`) and any
grain-catalyzed reaction (handled in the grain physics, `mole_reactions.cpp:4442-4446`).
Its own comment states the status (`mole_reactions.cpp:4475-4481`, verbatim):

> "This is only a subset and, thusfar, not actually used in getting the total
> heating.  Tests with pdr_leiden_hack_f1.in show that this heating rate can be
> up to 10% of the total heating"

and the call site confirms it (`cool_eval.cpp:307-313`):

```
fixit("test and enable chemical heating by default");
#if 0
	double chemical_heating = mole.chem_heat();
	thermal.setHeating(0,29, MAX2(0.,chemical_heating) );
	...
#endif
```

**Cloudy's chemical network heating is disabled.** EXHALE's
`molecular_reaction_heat.f90` therefore does something Cloudy does not do at all,
and there is no Cloudy result to compare the R5, R6, R7 and R16 corrections
against. The inventory's settlement of those channels from the primary
literature (Giusti-Suzor et al. 1983, Takagi 2002, Kokoouline et al. 2001,
Strasser et al. 2001, Guberman 1994, READ, inventory section 5) remains the only
basis, and it is a better one than Cloudy offers.

---

## 5. What is worth taking over

Judged list, most valuable first. Each entry gives what Cloudy does, whether it
is physically better than what EXHALE does or plans, the published source and
whether that source is in `../references/`, and the cost.

### 5.1 Item 1: the net collisional heat of the ladder, in place of both fractions

**(a)** `mole_h2.cpp:2211-2265`: gas heating from the vibrational ladder is
`sum over level pairs and colliders of (downward rate - upward rate) * dE`, with
the upward rate from the downward one by detailed balance, both from the same
tabulated coefficients that drive the level solver.

**(b) Physically better than what EXHALE does, unambiguously, IF EXHALE ever
solves level populations.** It has no `n_cr`, no `f_quench`, no v = 1 stand-in
for the whole ladder, and it is exact in both limits. It is NOT better than what
EXHALE does today, because EXHALE has no level populations to put in it: without
a solved ladder the expression cannot be evaluated. So this is the form to
adopt at the moment a cascade becomes buildable, and the inventory has already
established that it is not buildable now (three of the four ingredients missing,
READ, inventory section 0). **Recommendation: record it as the target form in
`h2_vibrational_relaxation.f90`, and do not build it.**

**(c)** No publication; it is the standard statistical-equilibrium expression.
Nothing to request.

**(d)** Cost: the cascade itself, which the inventory prices as blocked, not
expensive. The heat expression on top of a solved ladder is half a day.

### 5.2 Item 2: state-resolved Lyman-Werner dissociation heat

**(a)** `dissprob_*.dat` (Abgrall, Roueff & Drira 2000) carries, for every (v, J)
of each of the six excited singlet states, both the spontaneous dissociation
probability and the kinetic energy released in eV;
`mole_h2.cpp:2197-2205` sums `Pop * dissprob * disske` over them.

**(b) Better than EXHALE's constant 0.4 eV** (READ,
`src/modules/lower_atmosphere/lyman_werner.f90:349-351, 405`), but by how much is
unmeasured. Cloudy's own summary of the range is "typically 0.25 - 0.5 eV per
event" (`cool_eval.cpp:154-156`), so EXHALE's constant sits near the middle and
the error is bounded at roughly a factor 1.25 either way on a term that is itself
a small part of the molecular heat budget. **Not worth taking as a table**; worth
taking as a bound to be written at the EXHALE code site, replacing the present
bare assertion with "0.25 to 0.5 eV, Abgrall et al. 2000 per (v, J), 0.4 eV
adopted".

**(c)** Abgrall, Roueff & Drira (2000), A&AS 141, 297. **IN
`../references/Abgrall_2000AAS_141_297.pdf`.** The kinetic-energy column should
be checked against the paper before the bound is written into the code, which is
one PDF read and is not done here.

**(d)** One hour for the comment; a day if the level-resolved table is ever wanted.

### 5.3 Item 3: the Glover & Abel (2008) low-density branch of the H2 line cooling

**(a)** `cool_eval.cpp:726-1250` and `1260-1335`: cooling per H2 in the
low-density limit as fitted polynomials per collider and per ortho/para
combination, combined with the LTE limit by the harmonic mean of their eq. 39.

**(b) Better than EXHALE's `molecular_infrared_cooling.f90` in one specific
respect and one only.** EXHALE's module is the LTE branch: Boltzmann populations,
optically thin line sum, tabulated (READ, inventory section 3). That is exactly
Cloudy's `cooling_high`. What EXHALE has no equivalent of is `cooling_low`, the
regime where the levels are collisionally starved. **In the LHS 1140 b molecular
layer that regime is not reached**: the inventory measured the collider sum at
cell 1 to give `q = 8.0e+02 s^-1` against `A_tot,max = 5.594e-06 s^-1`, i.e.
`1 - f_quench = 7.0e-09` (READ, inventory section 2.5), so the ladder is in LTE
by eight orders and `cooling_low` would dominate the harmonic mean nowhere.
**So: not needed for this layer, and the reason it is not needed is already
measured.** It becomes relevant only for a thinner layer or a lower base
pressure, and the harmonic-mean combination is the cheap way to be right in both
limits at once when that day comes.

**(c)** Glover & Abel (2008), MNRAS 388, 1627. **IN
`../references/Glover_2008_MNRAS_388_1627.pdf`.** The fit coefficients are in
their Tables 1 to 7 and would be transcribed from the paper, not from Cloudy.

**(d)** Two days including the tests, if ever wanted.

### 5.4 Item 4: the He collision set that covers the whole ladder, still blocked

**(a)** `coll_rates_He_ORNL.dat`, the default He collider (`h2.cpp:172`), covers
300 of EXHALE's 302 levels with downward rows including all 36 levels at
v >= 11 (READ, inventory section 2.3).

**(b) It is the only set in hand that could make a helium cascade buildable, and
it cannot be used**, because its stated source is "Lee, T. G. et al. 2007, ApJ,
in preparation" (READ, the file header) and no such paper exists (MEASURED, ADS).
Section 1.8 adds three published leads from the same group. Until one of them is
shown to be the origin, nothing is taken from it. **Unchanged from the
inventory's verdict; what this memo adds is the three bibcodes.**

**(c)** Request list below.

**(d)** Zero until the request is answered.

### 5.5 Item 5: the LTE spread of unassigned formation channels

**(a)** `mole_h2.cpp:175-190`: any formation channel with no nascent
distribution is deposited over the ladder in local LTE proportions.

**(b) It is a clean default and worth recording as such, but it is NOT better
than what EXHALE does**, because the two are the same statement seen from
opposite ends. A thermal deposit yields zero net collisional heat, which is
`f_quench = 1` with `f_trans = 0`, which is what EXHALE's term evaluates to here.
The value of knowing this is argumentative: **the most detailed H2 model in
general use makes the same assumption EXHALE makes, and makes it silently.**
Worth one sentence at the EXHALE code site.

**(c)** No publication.

**(d)** Ten minutes.

### 5.6 What is worth NOT taking

- **The g-bar rule** (`mole_h2_coll.cpp:133-158`). A temperature-independent
  rate from an unpublished fit, which in a default Cloudy run supplies EVERY
  vibrational H2-H2 quenching rate at 70 to 2400 times below the Le Bourlot
  value (MEASURED, section 1.5). EXHALE's three measured tables are better than
  this everywhere they reach, and where they do not reach, holding at the end
  value is more honest than a constant.
- **The invented collisional dissociation coefficient** and therefore the
  three-body association rate that detailed balance builds on it
  (`mole_h2_coll.cpp:41-47`, "we made this up"). EXHALE's R15 rate must not be
  compared against it.
- **Cloudy's H3+ dissociative recombination branching**, 0.25/0.75
  (`mole_reactions.cpp:2001-2002`), which disagrees with the measured
  0.70 +/- 0.07 that EXHALE already uses.
- **Discarding the translational share** of the formation energy
  (section 2.2, `Xtrans` computed and never deposited). EXHALE's prompt
  deposit is the better treatment and the inventory has measured that it
  thermalizes with at least two orders of margin (READ, inventory section 4).
- **The disabled chemical heating** (`cool_eval.cpp:307-313`). EXHALE's ledger
  is ahead of Cloudy here, not behind it.

### 5.7 Request list

Papers that would have to be obtained for anything above to move. Each entry
says what must be checked in it and what is blocked without it.

| # | reference | bibcode | what to check | blocked without it |
|---|---|---|---|---|
| 1 | Lee, T. G., et al., the origin of `coll_rates_He_ORNL.dat`. Three published leads, none stated by Cloudy to be the origin: Mack, Clark, Forrey, Balakrishnan, Lee & Stancil (2006), PhRvA 74, 052718, "Cold He+H2 collisions near dissociation"; Ohlinger, Forrey, Lee & Stancil (2007), PhRvA 76, 042712, "H2 dissociation due to collisions with He"; Lee, Balakrishnan, Forrey, Stancil, Shaw & Schultz (2008), ApJ 689, 1105 | `2006PhRvA..74e2718M`, `2007PhRvA..76d2712O`, `2008ApJ...689.1105L` | whether any of them publishes the He + H2(v, J) de-excitation rate set that reaches v = 14, and over what temperature range | the only all-ladder collisional set stays unusable; a helium cascade stays unbuildable even in principle. This is the inventory's request item 2, kept open with leads attached |
| 2 | Takahashi, Junko (2001), ApJ 561, 254, "The Ortho/Para Ratio of H2 Newly Formed on Dust Grains"; and Takahashi, J., & Uehara, H. (2001), ApJ 561, 843, "H2 Emission Spectra with New Formation Pumping Models" | `2001ApJ...561..254T`, `2001ApJ...561..843T` | the four-way energy partition (their Xvib, Xrot, Xtrans, Xdust) and whether any of it transfers to a gas-phase channel | nothing in EXHALE. **These are grain-formation papers and EXHALE has no grains.** Requested only if the four-way partition is wanted as a template for the three-body one; LOW priority, and the inventory's gas-phase sources are the right ones |
| 3 | Launay, J. M., Le Dourneuf, M., & Zeippen, C. J. (1991), A&A 252, 842 | `1991A&A...252..842L` | the state-resolved H- + H -> H2(v, j) + e distribution and its temperature dependence | nothing in EXHALE today. EXHALE has no H- route in the molecular network; requested only if one is ever added. LOW priority |
| 4 | Komasa, J., et al. (2011), JCTC 7, 3105, "Quantum Electrodynamics Effects in Rovibrational Spectra of Molecular Hydrogen" | `2011JCTC....7.3105K` | the X-state level energies, Cloudy's source for `energy_X.dat` | nothing. EXHALE's ladder agrees with Cloudy's to the level (READ, inventory sections 2.3 and 3); requested only if the 302 against 303 difference ever needs adjudicating. LOW priority |

Items 2 to 4 are recorded for completeness and are NOT recommended for request
now. **Only item 1 blocks anything.**

---

## What EXHALE should and should not take

EXHALE should take one thing from Cloudy now, and it is a sentence rather than
a routine: the record, at the site of `h2_vibrational_heat_fraction`, that the
form this fraction approximates is the net of downward minus upward collisional
rates summed over the ladder, that Cloudy c25.00 evaluates that net exactly from
a full 303-level statistical equilibrium (`mole_h2.cpp:2211-2265`), and that the
two agree in this layer because the ladder here is in LTE by eight orders. That
converts the fraction from an approximation with an asserted validity into an
approximation with a named exact form behind it and a measured distance from it.
Alongside it belong three smaller corrections of the same kind: the
Lyman-Werner fragment kinetic energy stated as the 0.25 to 0.5 eV range of
Abgrall et al. (2000) with 0.4 eV adopted rather than as a bare constant; the
note that Cloudy deposits unassigned formation channels thermally, which is
EXHALE's `f_trans = 0`, `f_quench = 1` by another name; and the note that
EXHALE's infrared cooling is the LTE branch of the Glover & Abel (2008)
combination whose low-density branch this layer does not need. None of this
changes a number.

EXHALE should not take Cloudy's numbers. Its default H2-H2 vibrational rates are
a temperature-independent unpublished fit two to three orders below the measured
values; its three-body association rate rests on a coefficient the source
describes as made up; its H3+ dissociative recombination branching is 0.25/0.75
against the measured 0.70/0.30 EXHALE already carries; its chemical network
heating is switched off; and it computes the translational share of the grain
formation energy and then discards it. On each of these five EXHALE is ahead,
and on the one place where Cloudy holds something EXHALE wants, the He collision
set that covers the whole ladder, the set cannot be cited. The cascade the
inventory declared unbuildable remains unbuildable, and Cloudy is now a second,
independent reason to believe the bounded two-parameter form is the right thing
to keep.
