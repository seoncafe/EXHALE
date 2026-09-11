# A well-balanced flux difference for the near-hydrostatic layer

PLAN_20260909_rev1 item N37, decision 23 (a) of
`docs/To_be_determined_by_user_20260906.md` (decided by the user 2026-09-10).
Written before the code, as the item requires. The option it specifies is
`Well balanced:` and it is default off; nothing here makes it default.

## 1. What the published method is, and which of its two variants we take

Two papers define the method. Neither publisher PDF could be downloaded from
this machine (A&A serves `aanda.org` behind a bot filter that answers 403 to
`curl` and to the fetch tool alike; the JCP paper is closed at Elsevier), so
what was read in full is the authors' own institutional version of each, saved
as `../references/Kappeli2016_AA587_A94_SAMreport2015-40.pdf` (ETH SAM Research
Report 2015-40, November 2015, the A&A submission text) and
`../references/Kappeli2014_JCP259_199_SAMreport2013-05.pdf` (ETH SAM Research
Report 2013-05). Section and equation numbers below are those of the reports;
the published articles are still to be obtained and the numbers checked
against them.

**Käppeli and Mishra 2014 (JCP 259, 199), the isentropic local equilibrium.**
Within cell i the equilibrium is defined by holding the cell's own specific
entropy s_i and integrating the mechanical balance in the form
h + phi = const (their eq. 2.7), so the subcell equilibrium enthalpy is
h_0,i(x) = h_i + phi_i - phi(x); the equation of state then returns p_0,i(x)
and rho_0,i(x) from (s_i, h_0,i(x)) (their eqs. 2.9 to 2.12). Both the density
and the pressure are extrapolated to the faces from this equilibrium (their
eq. 2.13, and eq. 2.25 with the second-order perturbation added), and the
momentum source is the face-pressure difference of the cell's own equilibrium,
S_rho_v,i = [p_0,i(x_i+1/2) - p_0,i(x_i-1/2)]/dx (their eqs. 2.14 and 2.26),
which their theorem 1 shows is a second-order consistent discretization of
-rho dphi/dx and which, combined with a numerical flux that is consistent, makes
the flux difference equal the source exactly at the discrete equilibrium.

**Käppeli and Mishra 2016 (A&A 587, A94), arbitrary entropy stratification.**
The equilibrium is the same mechanical balance, but the subcell integral of
eq. (13) is evaluated with a CONSTANT density inside the cell (their eq. 14) and
a piecewise linear potential (their eq. 15), giving the pressure extrapolation
    p_0,i(x_i-1/2) = p_i + rho_i (phi_i - phi_i-1)/dx dx/2,
    p_0,i(x_i+1/2) = p_i - rho_i (phi_i+1 - phi_i)/dx dx/2       (their eq. 16)
and the non-uniform-mesh form in their appendix A, eq. (A.4). No thermodynamic
variable other than the cell's own pressure and density enters, so no thermal
equilibrium is assumed and "our reconstruction preserves discrete hydrostatic
equilibria with arbitrary temperature or entropy stratification" (their section
2.1.1). The discrete equilibrium it preserves is their eq. (18),
(p_i+1 - p_i)/dx = -(rho_i + rho_i+1)/2 (phi_i+1 - phi_i)/dx, which is exactly
the statement that the two neighboring extrapolations agree at the shared face.
Density and velocity keep the standard reconstruction, so at equilibrium the
face is a stationary contact discontinuity and the only requirement on the flux
is that it resolves one exactly, their eq. (10),
F([rho_L,0,p],[rho_R,0,p]) = [0,p,0] (HLLC and Roe do). The source stays the
centered difference of their eq. (11).

**Second order.** In both papers the reconstruction is applied to the
DEPARTURE: the equilibrium of cell i is extrapolated to the neighboring cell
CENTERS, the perturbation data p_1,i(x_i±1) = p_i±1 - p_0,i(x_i±1) is formed
(2016 eqs. 25 and 26; 2014 eq. 2.23), p_1,i(x_i) = 0 by construction, and the
standard limited slope is applied to that data (2016 eq. 27). At a discrete
equilibrium the perturbation data vanishes identically, so the reconstruction
returns the equilibrium and the scheme is exact; away from it the perturbation
is reconstructed to second order and "the perturbation part is not assumed to be
small" (2016 section 2.1.3).

**What the papers measure.** Their exactness test is the hydrostatic atmosphere
of section 3.1: an isentropic and an isothermal column in a constant field,
initialized by solving the DISCRETE equilibrium (their eq. 42) with a Newton
iteration, evolved for two sound crossing times at N = 32 to 2048, reporting
the L1 norm of the difference between the initial and final pressure (their
fig. 2). The well-balanced scheme sits at 1e-14 to 1e-16 at every N while the
standard scheme falls off as a power of the mesh width. Wave propagation on top
of the equilibrium (their section 3.1.2) is the second measurement, and the
2014 paper adds the same pair in three dimensions (its sections 3.2.1, 3.2.2).

**Our choice: the 2016 constant-density equilibrium**, not the 2014 isentrope,
for three reasons, all properties of our layer and not of theirs. (i) Our base
layer has a composition front and a temperature gradient in it, so an isentrope
through cell j does not pass through cell j+1; the equilibrium the 2014
reconstruction preserves would be one our column is not on, and the two
neighboring extrapolations would disagree at O(dx^2) instead of at rounding.
The 2016 equilibrium preserves an arbitrary stratification, which is what our
column is. (ii) The isentropic form needs an EOS inversion in every cell at
every residual evaluation (`continue_hydrostatic_isentrope` runs a Newton
iteration with the caloric EOS and the molecular mixture); the constant-density
form is two multiplications. (iii) The isentropic form would have to decide
whose composition the extrapolation carries across a front; the constant-density
form asks nothing about the composition at all. The ghost construction of
`base_boundary.f90` stays on the isentrope and is untouched: it continues the
solution BELOW the domain, where the entropy of cell 1 is the physical
statement, and it is not part of the interior pair.

## 2. The operator we have

`RK_rhs` (`src/modules/time_step/RK_rhs.f90`) forms, for cell j with faces at
r_edg(j-1) and r_edg(j), A+ = r_+^2, A- = r_-^2, dV = (r_+^3 - r_-^3)/3:

    dF(1,j) = (A+ F+_1 - A- F-_1)/dV
    dF(2,j) = (A+ F+_2 - A- F-_2)/dV        [ + (p_+ - p_-)/dr under WENO3 ]
    dF(3,j) = (A+ F+_3 - A- F-_3 + dF3p)/dV

with `source` (`src/modules/states/Source.f90`) returning
S(2,j) = -0.5 (rho_L + rho_R)(phi_i(j) - phi_i(j-1))/dr, plus, under PLM,
the geometric term (A+ - A-) p_j/dV, where under PLM the pressure sits INSIDE
F_2 (`Phys_flux`) and under WENO3 it does not. The marching stage is
u <- u - dt (dF - S) and the stationary residual is R = dF - S - (heat - cool).

At the base the terms of the momentum row are individually O(p) A/dV and cancel
to O(rho g); N33 measured the consequence, a rounding floor
epsilon x (face state) x r^2/dV, 1e6 to 2e9 times epsilon ||F|| on the mass and
energy rows of cells 1 to 58, with r^2/dV = 5.1e3 at a base cell.

## 3. The design, term by term

Write phi_up(j) = phi_i(j) - phi_c(j) and phi_dn(j) = phi_c(j) - phi_i(j-1),
both positive, from the production arrays `Gphi_i`, `Gphi_c`.

**(a) The local equilibrium of cell j.** At rest, constant density rho_j,
through the cell's own (rho_j, p_j):

    p_eq,j(r) = p_j - rho_j (phi(r) - phi_c(j)),

so the two face values and the two neighbor-center values are

    P_up(j) = p_j - rho_j phi_up(j)          at r_edg(j)
    P_dn(j) = p_j + rho_j phi_dn(j)          at r_edg(j-1)

and, continued THROUGH the face to a neighboring center (see (b)),

    p_eq,j(r_j+1) = P_up(j) - rho_j+1 phi_dn(j+1)
    p_eq,j(r_j-1) = P_dn(j) + rho_j-1 phi_up(j-1).

This is the 2016 eq. (16) with the potential evaluated at the true face radius
rather than by the linear staggered interpolation of their eq. (15): our
potential is an analytic function of r, known at faces and at centers alike
(`set_gravity_grid`), so the interpolation their eq. (15) needs (they allow for
a potential known only at cell centers) is not needed and is not used. On a
non-uniform mesh this is their appendix A with the exact potential difference
in place of the weighted one.

**(b) The reconstruction carries the departure.** The perturbation data of
cell j is the neighboring cell's pressure measured against the equilibrium
continued to it THROUGH THE FACE, with cell j's density up to the face and
the neighbor's beyond it,

    d_plus(j)  = (p_j+1 - p_j) + rho_j phi_up(j) + rho_j+1 phi_dn(j+1)
    d_minus(j) = (p_j-1 - p_j) - rho_j phi_dn(j) - rho_j-1 phi_up(j-1)

with the value at the cell's own center identically zero. Continuing cell j's
own constant density over the whole gap to the neighboring center instead
would leave data of size (rho_j - rho_j+1) x (potential difference across
half a cell), which does NOT vanish on the discrete equilibrium: the scheme
would then be second order and not exact (MEASURED: with that form the WENO3
scheme leaves 2.5e-5 of the momentum weight on the discrete column at N = 250
instead of 3.9e-14). The form above is the spherical, non-uniform-mesh
counterpart of the average density in the 2016 paper's eq. (26), and
d_plus(j) is exactly the face mismatch dp_eq(j) of (c) below. The limited slope
(PLM) or the WENO3 stencil is applied to (d_minus(j), 0, d_plus(j)) exactly as
it is applied to (p_j-1, p_j, p_j+1) today, including the same limiter, the
same smoothness indicators and the same volume-share coefficients, and returns
the two face DEPARTURES q_L(j) (upper face) and q_R(j) (lower face). The face
states handed to the Riemann solver are

    WL(3,j)   = P_up(j) + q_L(j),      WR(3,j-1) = P_dn(j) + q_R(j),

and density and velocity are reconstructed exactly as today (2016 section
2.1.3: no equilibrium reconstruction is applied to the density). At a discrete
equilibrium d_plus and d_minus vanish, the slope is zero, and the face
pressures are the equilibrium ones.

**(c) The jumps, and why they are formed from small numbers.** The two
extrapolations at face f (between cells f and f+1) differ by

    dp_eq(f) = P_dn(f+1) - P_up(f)
             = (p_f+1 - p_f) + rho_f+1 phi_dn(f+1) + rho_f phi_up(f),

a difference of the two cell pressures (exact in binary floating point whenever
the two lie within a factor of two of each other, Sterbenz's lemma, which holds
throughout a resolved layer) plus two terms of the size of the hydrostatic
pressure drop across a cell. It is the discrete equilibrium condition of the
scheme, our form of their eq. (18), and it vanishes on it. The pressure jump
the Riemann solver uses is then

    dp(f) = dp_eq(f) + q_R(f) - q_L(f),

every term of which is of the size of the departure and not of the state. This
is the point of the whole change on the mass and energy rows: today the same
number is formed as pR - pL from two O(p) face states, so it carries their
rounding. The density jump and the velocity jump keep the standard form, as in
the paper: at equilibrium the density jump is the stationary contact the flux
must resolve.

**(d) The pressure force, cancelled analytically.** Let q_up(f) = p_out(f) -
P_up(f) and q_dn(f) = p_out(f) - P_dn(f+1) be the face pressure the flux
routine returns, measured from each side's own equilibrium. Both are formed
from dp_eq, q_L and q_R and never as a difference of two O(p) numbers: for the
Roe and LLF fluxes p_out = (p_L + p_R)/2 and q_up = (q_L + q_R + dp_eq)/2, for
HLLC p_out is one side's pressure and q_up is that side's departure. Then, for
the PLM form of the momentum equation,

    (A+ (F+_2) - A- (F-_2))/dV - (A+ - A-) p_j/dV + rho_j [A+ phi_up + A- phi_dn]/dV
      = (A+ (F+_2 - P_up) - A- (F-_2 - P_dn))/dV
      = (A+ (Fk+_2 + q_up(j)) - A- (Fk-_2 + q_dn(j-1)))/dV,

because A+ P_up - A- P_dn - (A+ - A-) p_j = -rho_j (A+ phi_up + A- phi_dn)
IDENTICALLY, the cell's own p_j cancelling term by term. Fk_2 is the momentum
flux WITHOUT the pressure, which is what `Phys_flux` already builds under
WENO3. So with the key on the row is assembled in the last form, the gravity
source is the bracket on the left, and `source` returns S(2,j) = 0: the
equilibrium's flux difference and its source have cancelled in the algebra, and
what remains is a sum of departures. The same for the WENO3 form,

    (A+ Fk+_2 - A- Fk-_2)/dV + [q_up(j) - q_dn(j-1) - rho_j (phi_i(j) - phi_i(j-1))]/dr,

where the last term is P_up(j) - P_dn(j) written out and is again exactly the
gravity source it replaces. The mass and energy rows keep their formulas: at
rest their equilibrium flux is zero exactly, because the reconstructed velocity
of a cell at rest is zero exactly (a limited slope of a constant is zero), the
Roe dissipation of the mass row is dp/(2 a) with the dp of (c), and the
energy row's gravitational term dF3p rides on the mass flux and vanishes with
it.

**What is NOT changed.** The reconstruction of density and velocity and the
limiter itself; the Riemann solvers and their wave-speed estimates and entropy
fix (the O(1) face states they are evaluated at are unchanged in meaning); the
boundary construction of `base_boundary.f90`, which is already an isentropic
continuation; the low-Mach damping; the positivity limiter, which acts on the
assembled face state; the energy source; the species rows, which ride on the
mass flux.

**(e) The key.** `Well balanced: True|False`, default False, K46 of
`docs/input_schema.md` and the manual's key table. `write_setup_report` echoes
it, the resolved-input file carries it, and it is the 21st token of the restart
metadata's option field (`wellbal`), nameable in `Restart option change`: it
does not change how many unknowns a state has, so a state may be continued
across it once named.

**(f) Expected movement.** Every case moves with the key on: the pressure
reconstruction and the discrete gravity are different (both second-order
consistent, neither a subset of the other), and the pressure/gravity pair is in
every run. With the key off every path is bitwise the one the goldens hold.

**(g) For the code site.** The equilibrium is local to the cell and re-formed
at every evaluation from that cell's own (rho, p) and the grid's potential;
nothing is stored between evaluations and no reference state is assumed known.
The scheme is exact for the discrete equilibrium of (c) and second order for
departures from it.

## 4. Where this design departs from the papers, and why

1. Spherical geometry. Their flux difference is (F+ - F-)/dx; ours carries the
   areas A± and the cell volume, and the PLM form of the momentum equation
   carries a geometric pressure term. The analytic cancellation of section 3(d)
   is therefore ours and not theirs: their eq. (17) is the Cartesian case
   A+ = A- = 1 of it.
2. The source is the face-pressure difference of the cell's own equilibrium
   (the 2014 paper's eq. 2.26) and not the centered difference of the 2016
   paper's eq. (11). The two agree to second order, but only the first cancels
   the equilibrium flux difference in the ALGEBRA rather than in floating
   point, which is what this item is for.
3. The departure is carried as a separate small number through the jump and the
   pressure force (section 3(c) and 3(d)) instead of being re-formed by
   subtracting two assembled face states. Their papers state the property up to
   machine precision and do not distinguish these; for us the distinction is
   the whole point.
4. WENO3. Their second-order scheme uses a limited linear slope; our WENO3 form
   applies the same equilibrium/departure split to the third-order stencil.
   Nothing in the argument uses the linearity of the reconstruction, only that
   it returns zero for zero data.
5. The caloric EOS with a molecular mixture. Nothing in the constant-density
   equilibrium refers to the EOS, which is why this variant was chosen.
6. `Numerical flux: LLF` does not resolve a stationary contact (their eq. 10),
   so with LLF the scheme is NOT well balanced; the key still applies the
   equilibrium reconstruction and source, and the residual keeps the LLF
   dissipation of the density jump. The same holds at a face the positivity
   repair replaces with a first-order LLF flux.

## 5. What this cannot reach, stated before it is measured

The cell pressure p_j is itself a rounded O(1) double, so it moves in steps of
one unit in its last place, about 1.9e-16 x p, as the state moves continuously.
Every jump above is formed from p_j and p_j+1, so its step size, and therefore
the step size of the interface flux and of the row, cannot fall below
epsilon x p x r^2/dV: that is exactly N33's arithmetic bound, and it is the
bound and not the measured floor that this design attacks. What the design
removes is the FACTOR between the measured floor and that bound (N33 measured 6
to 12) and, far more importantly, the O(dx) spurious flux that the standard
scheme puts into a hydrostatic column, which is 1e12 times larger than either.
On the exactness measurement (a column at the discrete equilibrium) the
prediction is the rounding level; on N33's second-difference floor the
prediction is a factor of order ten, that is about one decade, NOT the three
decades the item's acceptance names. If that is what the measurement says, the
finding is that the floor is set by the representation of the pressure in the
cell average and not by the flux assembly's grouping, and no rearrangement of
the operator reaches it: only a formulation whose UNKNOWN is the departure
would, which is a different item.
