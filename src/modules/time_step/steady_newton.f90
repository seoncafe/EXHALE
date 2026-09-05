      module steady_newton
      ! Steady-state residual in vector form, F(Y) = 0, for a direct
      ! (Newton / pseudo-transient-continuation) solve. Increment (ii)-2:
      ! pack/unpack + newton_residual. The Jacobian and the PTC driver are
      ! added in later increments.
      !
      ! Unknowns: the hydro conserved variables on the PHYSICAL cells
      ! j = 1..N (neq = 3*N). The ghost cells (1-Ng..0 lower, N+1..N+Ng
      ! upper) are NOT unknowns: Apply_BC fills them from the interior each
      ! residual evaluation (lower = fixed rho_bc / p = 1+dp_bc with the
      ! v-valve max(v_1,0); upper = zero-gradient / WENO3 extrapolation).
      !
      ! Ionization/temperature are eliminated LOCALLY inside the residual:
      ! given Y, ioniz_eq re-solves the equilibrium fractions and returns the
      ! consistent heat/cool, so F(Y) is a pure-hydro residual of size 3*N
      ! (the "method A" choice from ATES_sundials_solver_plan.md).

      use global_parameters
      use species_table,  only: n_mion, isp_H2
      use Conversion,      only: U_to_W, U_to_W_interior
      use composition,     only: get_species_densities, comp_T_from_p
      use BC_Apply,        only: Apply_BC
      use ionization_equilibrium, only: ioniz_eq, ioniz_eq_ledger,      &
                 set_ioniz_eq_sweep_state_kind, ieq_state_marching,      &
                 ieq_state_steady_iterate, ieq_state_steady_candidate,   &
                 finite_real, keep_background_of_adopted_state,          &
                 keep_background_of_best_iterate,                        &
                 adopt_background_of_best_iterate,                       &
                 install_background_of_adopted_state
      use excited_hydrogen,       only: excited_H_update
      use steady_residual_mod,    only: assemble_residual,             &
                                        residual_row_scale,            &
                                        relnorm_over_cells,            &
                                        state_scales_of_cell,          &
                                        flux_spread_of_state,          &
                                        flux_spread_above_radius,      &
                                        steady_gates_met,              &
                                        write_residual_breakdown
      use diffusive_photochemistry, only: carrier_steady_residual,     &
                                        carrier_column_scale,          &
                                        carrier_row_scale_H2,          &
                                        carrier_headroom,              &
                                        n_carrier_max, ic_H2

      implicit none
      private
      public :: neq_newton, pack_U, unpack_U, newton_residual,         &
                eval_residual, frozen_residual, build_banded_jac,      &
                band_matvec, band_matvec_transpose, kl_jac, ku_jac,    &
                ncolor_jac, solve_steady_ptc,                          &
                jv_product, solve_steady_jfnk, set_base_fix,           &
                gate_rnorm_accepted, gate_fspread_accepted,            &
                set_carrier_unknown, carrier_unknown_on

      ! Banded-Jacobian geometry. WENO3 reconstruction reaches +/-2 cells,
      ! and within a cell all 3 hydro variables couple, so in the flat
      ! ordering Y(3*(j-1)+k) a column couples to rows within
      ! |row-col| <= 3*2 + 2 = 8. Half-bandwidths kl = ku = 8; a graph
      ! coloring with stride ncolor = kl+ku+1 = 17 makes same-color columns'
      ! row supports disjoint, so the FROZEN-radiation residual (strictly
      ! local) is captured exactly by 17 finite-difference probes.
      ! The two gate values of the state the LAST steady solve returned, kept
      ! so that the caller can re-measure them on the state it finally writes
      ! and prove the two are the same (EXHALE_main, after write_output).
      ! Nothing between the solve's return and the write touches u -- the
      ! refresh at EXHALE_main lines 1599-1607 recomputes rho, v, p FROM u and
      ! only moves the composition -- so the check is an assertion, not a
      ! correction; it exists because "the state the user reads" and "the state
      ! the gate accepted" must be the same object and nothing was saying so
      ! (section 133.6).
      real*8 :: gate_rnorm_accepted   = -1.0d0
      real*8 :: gate_fspread_accepted = -1.0d0

      ! WITH THE CARRIER UNKNOWN THE GEOMETRY CHANGES, so these are set once
      ! per solve by set_carrier_unknown rather than fixed at compile time.
      !
      !   nvar = 3 (hydro alone):   |row-col| <= 3*2 + 2 = 8,  17 colors.
      !   nvar = 4 (+ n(H2)):       |row-col| <= 4*2 + 3 = 11, 23 colors.
      !
      ! The +/-2 is the WENO3 reach; the carrier row itself reaches only
      ! +/-1, because its second-order correction is a frozen source inside
      ! the model (carrier_advection_correction), so the hydro rows still set
      ! the bandwidth and the extra variable only widens the stride.
      integer :: nvar_jac = 3
      integer :: kl_jac = 8, ku_jac = 8
      integer :: ncolor_jac = 17
      ! J^T J of a band of half-width kl_jac has half-width 2*kl_jac, and it
      ! is what the damped Gauss-Newton (Levenberg-Marquardt) step below
      ! factors when the Newton direction ascends the merit.
      integer :: kl_jac_normal = 16
      integer :: ku_jac_normal = 16
      ! LAPACK general-band leading dimension for that matrix.
      integer :: ldab_normal   = 49

      ! HOW A PROBE STATE THAT CANNOT BE DESCRIBED IS ANSWERED
      ! (docs/Update_EXHALE.md section 121).
      !
      ! The Jacobian columns and the Krylov products are directional
      ! derivatives OF THE RESIDUAL AT Y. The point Y + h*v at which they are
      ! sampled is not physics: h is ours to choose, and v is a Krylov
      ! direction with no physical meaning at all. So when the sample point
      ! leaves the set of states the code can describe -- a balance row of the
      ! chemical network goes to NaN or overflows, and the equilibrium sweep
      ! reports it -- the answer is to sample CLOSER TO Y, not to accept the
      ! non-finite residual and not to declare the derivative undefined. The
      ! step is halved and the probe repeated, at most this many times; over
      ! six halvings the finite-difference step falls by 64 and stays far
      ! above the round-off floor of the residual.
      !
      ! If even the smallest of those steps is inadmissible, the boundary of
      ! the describable states runs through Y itself along v, and there is no
      ! derivative to sample. Then the derivative is NOT invented: the Krylov
      ! cycle stops at the directions it did build (a shorter GMRES cycle is
      ! still a valid approximate solve, and the line search judges the step
      ! it produces), and a Jacobian color that cannot be probed leaves its
      ! columns at zero, which degrades the PRECONDITIONER of those unknowns
      ! to plain relaxation and costs iterations, not correctness. Both
      ! events are counted and printed.
      integer, parameter :: n_probe_step_halvings = 6

      ! A SECOND WAY A PROBE STATE CANNOT BE DESCRIBED, and the same answer.
      ! The residual's energy row contains heat - cool, and both are evaluated
      ! on the composition the equilibrium sweep accepted. When that
      ! composition is not a certified root of the network -- acceptance class
      ! 4, kept only by the relaxation amnesty, and class 5, the constrained
      ! continuation -- the residual is a number computed on a state the
      ! chemistry does not describe, and its value carries no information
      ! about how far the hydro is from steady.
      !
      ! The test is therefore a COMPARISON, not a threshold: a trial is
      ! refused when it leaves MORE cells uncertified than the iterate it is
      ! being measured against. That is the statement the line search actually
      ! needs -- it compares a trial's merit with the iterate's, and the
      ! comparison means something only while the trial's chemistry is no
      ! worse -- and at a fully certified iterate it is the same thing as
      ! refusing any uncertified trial at all.
      !
      ! WHY NOT THE ABSOLUTE FORM. Refusing every trial with an uncertified
      ! cell was written first and measured: on the hot Uranus hand-off the
      ! iterate ITSELF carries 11 such cells, so all 17 Jacobian colors and
      ! every Krylov direction were refused, the solve had no Newton model at
      ! all (`||grad merit|| = 0`) and aborted at iteration 1. A rule that
      ! rejects the neighbourhood of the point it is standing on cannot be
      ! used to leave that point.
      !
      ! It is applied only where `admissible` is asked for, so the marching
      ! path and the current iterate are untouched.
      !
      ! Grounds (docs/Update_EXHALE.md section 138; measurement P49). On the
      ! hot Uranus hand-off the solve spent 66 iterations driving the layer at
      ! 1.19-1.23 R_p from 3.5e-14 to 7.7e-16 g cm^-3 and from 900 K to
      ! 65,000 K while ||R|| improved from 8.2e-2 to 2.8e-2. What held that
      ! state up was cell 259, accepted class 4 with a reaction residual of
      ! 2.4e-3, whose H3+ infrared cooling -- 93% of the blocking energy row --
      ! came from an H3+ population that cannot exist at 65,000 K. The
      ! residual was measuring a composition, not a wind.

      ! Frozen-base option: the residual of the first nfix_base physical
      ! cells is replaced by the anchor row F_j = Y_j - Yfix_j, pinning them
      ! to the warm-start state, so that Newton solves for the WIND on top of
      ! a given base state. Set via set_base_fix before a solve.
      ! NB: it was previously assumed that the cell-1 momentum row cannot be
      ! satisfied at all, the lower BC (hard-pinned rho_bc/p_bc ghosts + the
      ! v-valve max(v,0)) leaving it non-zero at any steady interior solution.
      ! Direct measurement on the HD 209458 b hand-off state does not support
      ! that: R(2,1) is a four-order cancellation of gravity against the
      ! pressure gradient whose remainder is 6.5e-5 of the gravity term, and a
      ! 2.9 ppm change of the ghost pressure drives it to zero
      ! (docs/newton_scaling_and_base_wall.md). The option is kept as an
      ! experiment, not as a remedy for a base that cannot be solved.
      ! Trial states this solve refused because their chemistry was not
      ! certified (see eval_residual). Counted here rather than in an ioniz_eq
      ! ledger, so that what a probe state did stays out of the run's own
      ! acceptance statistics. Reset at the start of each solve.
      integer :: n_trial_uncertified   = 0
      integer :: uncertified_cell_max  = 0
      real*8  :: uncertified_res_max   = 0.0d0
      ! Cells the LAST equilibrium sweep accepted without certifying them as a
      ! root, and the same count for the state the solve currently holds. The
      ! test is the comparison of the two (eval_residual).
      integer :: n_uncertified_last    = 0
      integer :: n_uncertified_iterate = 0

      integer :: nfix_base = 0
      real*8, allocatable :: Yfix_base(:)

      ! THE CARRIER ROW OF THE LAST FULL RESIDUAL EVALUATION (section 139).
      ! With n(H2) among the unknowns the state has a fourth row per cell,
      ! and it is not a hydrodynamic one: it is assembled by
      ! diffusive_photochemistry, in cm^-3 s^-1, on the frozen background the
      ! same evaluation's equilibrium sweep left. These hold its
      ! volume-weighted measure and its worst cell for the state
      ! eval_residual last saw, and again for the BEST iterate -- kept
      ! separately because the last evaluation of a solve is a REJECTED
      ! trial, not the state the solve returns.
      real*8  :: carrier_relnorm_last = 0.0d0
      real*8  :: carrier_cellmax_last = 0.0d0
      integer :: carrier_cell_worst   = 0
      real*8  :: carrier_relnorm_best = 0.0d0
      real*8  :: carrier_cellmax_best = 0.0d0
      integer :: carrier_cell_best    = 0
      ! Trial states refused because a cell's n(H2) left the element budget
      ! or went negative. The budget is a CONSTRAINT, not a row: clamping it
      ! would put a kink in the residual (the defect of section 138), and a
      ! state outside it is not one the write-back can describe.
      integer :: n_trial_headroom  = 0
      integer :: n_headroom_last    = 0
      integer :: n_headroom_iterate = 0
      ! Full residual evaluations of a solve, the unit its cost is counted
      ! in: each one is an equilibrium sweep over the column.
      integer :: n_resid_eval      = 0

      ! Explicit interfaces for the external LAPACK banded-LU routines used by
      ! the direct/preconditioned Newton solves below (double-precision,
      ! general band form). The array dummies are assumed-size so that the
      ! call sites, which pass the right-hand side as a rank-1 vector
      ! (nrhs = 1), match without any argument change; this only gives the
      ! compiler kind/rank/intent information and does not alter the call.
      interface
         subroutine dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
            integer,          intent(in)    :: m, n, kl, ku, ldab
            real*8,           intent(inout) :: ab(ldab,*)
            integer,          intent(out)   :: ipiv(*)
            integer,          intent(out)   :: info
         end subroutine dgbtrf

         subroutine dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
            character(1),     intent(in)    :: trans
            integer,          intent(in)    :: n, kl, ku, nrhs, ldab, ldb
            real*8,           intent(in)    :: ab(ldab,*)
            integer,          intent(in)    :: ipiv(*)
            real*8,           intent(inout) :: b(*)
            integer,          intent(out)   :: info
         end subroutine dgbtrs
      end interface

      contains

      ! ------------------------------------------------------!

      subroutine set_carrier_unknown(on)
      ! Put n(H2) among the Newton unknowns, or take it out, and set the
      ! band geometry that goes with it. Called once before a solve; the
      ! 3-unknown path is untouched when `on` is false, which is what keeps
      ! every run without the carrier transport bit-for-bit unchanged.
      logical, intent(in) :: on
      if (on) then
         nvar_jac = 4;  kl_jac = 11;  ku_jac = 11
      else
         nvar_jac = 3;  kl_jac = 8;   ku_jac = 8
      endif
      ncolor_jac    = kl_jac + ku_jac + 1
      kl_jac_normal = 2*kl_jac
      ku_jac_normal = 2*kl_jac
      ldab_normal   = 2*kl_jac_normal + ku_jac_normal + 1
      end subroutine set_carrier_unknown

      ! ------------------------------------------------------!

      logical function carrier_unknown_on()
      carrier_unknown_on = (nvar_jac .eq. 4)
      end function carrier_unknown_on

      ! ------------------------------------------------------!

      subroutine set_base_fix(Y, nfix)
      ! Anchor the first nfix physical cells at their current values.
      real*8, dimension(nvar_jac*N), intent(in) :: Y
      integer, intent(in) :: nfix
      nfix_base = max(0, min(nfix, N-10))
      if (.not. allocated(Yfix_base)) allocate(Yfix_base(nvar_jac*N))
      Yfix_base = Y
      end subroutine set_base_fix

      ! ------------------------------------------------------!

      subroutine cell_state_scales(Y, D)
      ! Diagonal scaling of the Newton system and of the line-search merit:
      ! the characteristic scale of each conserved unknown in its own cell.
      ! The three numbers come from state_scales_of_cell (steady_residual.f90).
      ! Since section 143 they are NOT what the convergence measure divides by:
      ! each residual row is divided by the largest term that row itself
      ! contains (residual_row_scale), so the merit and the acceptance test are
      ! no longer one expression apart. That is deliberate and it was measured
      ! -- scaling this system by those quantities as well, the mass row's
      ! varying by a factor 9 across the first two cells, left the molecular
      ! hot-Uranus solve with no descent direction after 179 iterations against
      ! 9 for the measure-only build (docs/p54_base_layer_mass_flux.md section
      ! 10.4). The unknowns keep the scales below; the rows do not.
      !
      ! It sets the finite-difference step sizes, the scaled Newton system
      ! D^-1 J D, and the line-search merit ||D^-1 F||_2.
      real*8, dimension(nvar_jac*N), intent(in)  :: Y
      real*8, dimension(nvar_jac*N), intent(out) :: D
      real*8  :: vv, cs
      integer :: j, i1, i2, i3
      do j = 1, N
         i1 = nvar_jac*(j-1) + 1;  i2 = i1 + 1;  i3 = i1 + 2
         call state_scales_of_cell(j, Y(i1), Y(i2), Y(i3),               &
                                   D(i1), D(i2), D(i3), vv, cs)
         ! THE CARRIER UNKNOWN'S COLUMN SCALE is its own magnitude with the
         ! floor diffusive_photochemistry builds for it -- the same
         ! construction as the momentum unknown's |rho v| + rho c_s, and for
         ! the same reason. Why it is not the ROW scale, and what happens if
         ! it is, is recorded where the two are defined.
         if (nvar_jac .ge. 4)                                            &
            D(i3+1) = max(carrier_column_scale(j)/n0,                    &
                          1.0d-30*max(abs(Y(i3+1)), 1.0d0))
      enddo
      end subroutine cell_state_scales

      ! ------------------------------------------------------!

      subroutine cell_row_scales(Y, D, Dr)
      ! ROW scaling of the Newton system, beside the column scaling D.
      !
      ! For the three hydrodynamic rows it is IDENTICAL to D -- copied, not
      ! recomputed -- which is the single scaling the solve has always used.
      ! The carrier row is the one place the two must differ, and why is
      ! recorded at carrier_row_scale_H2.
      !
      ! IT IS A LOCAL OF THE SOLVE AND NOT A MODULE ARRAY, so that it has the
      ! same aliasing status as D and the solve owns its own copy. That alone
      ! does NOT make a three-unknown solve bit-identical to one that divides
      ! by D throughout -- measured, the two WASP-121b Newton cases moved by
      ! 1.2e-9 and 4.0e-8 -- so every site that divides by the row scale
      ! branches on nvar_jac and keeps the original expression, character for
      ! character, in the three-unknown branch. The run-time band geometry is
      ! NOT the cause of that difference and that was tested: a build with
      ! kl_jac and the rest back at compile time is bit-identical to this one.
      real*8, dimension(nvar_jac*N), intent(in)  :: Y, D
      real*8, dimension(nvar_jac*N), intent(out) :: Dr
      integer :: j
      Dr = D
      if (nvar_jac .lt. 4) return
      do j = 1, N
         Dr(nvar_jac*(j-1)+4) =                                           &
            max(carrier_row_scale_H2(j)/n0,                               &
                1.0d-30*max(abs(Y(nvar_jac*(j-1)+4)), 1.0d0))
      enddo
      end subroutine cell_row_scales

      ! ------------------------------------------------------!

      subroutine apply_base_fix(Y, Fvec)
      ! Replace the residual rows of the frozen base cells by anchor rows.
      real*8, dimension(nvar_jac*N), intent(in)    :: Y
      real*8, dimension(nvar_jac*N), intent(inout) :: Fvec
      integer :: j, k
      do j = 1, nfix_base
         do k = 1, nvar_jac
            Fvec(nvar_jac*(j-1)+k) = Y(nvar_jac*(j-1)+k) - Yfix_base(nvar_jac*(j-1)+k)
         enddo
      enddo
      end subroutine apply_base_fix

      ! ------------------------------------------------------!

      integer function neq_newton()
      ! Number of Newton unknowns = 3 hydro variables over physical cells.
      neq_newton = nvar_jac*N
      end function neq_newton

      ! ------------------------------------------------------!

      subroutine pack_U(u, Y)
      ! Physical-cell conserved variables -> flat unknown vector.
      ! Y(3*(j-1)+k) = u(j,k), j = 1..N, k = 1..3.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(nvar_jac*N),         intent(out) :: Y
      integer :: j, k
      do j = 1, N
         do k = 1, 3
            Y(nvar_jac*(j-1)+k) = u(k,j)
         enddo
      enddo
      end subroutine pack_U

      ! ------------------------------------------------------!

      subroutine unpack_U(Y, u)
      ! Flat unknown vector -> physical-cell conserved variables.
      ! Ghost cells are left untouched (Apply_BC sets them).
      real*8, dimension(nvar_jac*N),         intent(in)    :: Y
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u
      integer :: j, k
      do j = 1, N
         do k = 1, 3
            u(k,j) = Y(nvar_jac*(j-1)+k)
         enddo
      enddo
      end subroutine unpack_U

      ! ------------------------------------------------------!

      subroutine pack_carrier(u, f_sp, Y)
      ! Fill the carrier slot of the unknown vector from the state's own H2:
      ! Y(nvar*(j-1)+4) = n(H2)_j/n0, IN THE CODE'S OWN DENSITY UNITS and not
      ! in cm^-3. The Newton takes norms of the whole vector -- the Krylov
      ! step size is sqrt(eps)(1 + ||Y||)/||v|| -- so a component thirteen
      ! decades from the hydrodynamic ones would set that step by itself.
      real*8, dimension(3,1-Ng:N+Ng),         intent(in)    :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)    :: f_sp
      real*8, dimension(nvar_jac*N),          intent(inout) :: Y
      integer :: j
      if (nvar_jac .lt. 4) return
      do j = 1, N
         Y(nvar_jac*(j-1)+4) = f_sp(j,isp_H2)*u(1,j)
      enddo
      end subroutine pack_carrier

      ! ------------------------------------------------------!

      subroutine unpack_carrier(Y, u, f_sp)
      ! Inverse of pack_carrier, onto the physical cells only.
      real*8, dimension(nvar_jac*N),          intent(in)    :: Y
      real*8, dimension(3,1-Ng:N+Ng),         intent(in)    :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      integer :: j
      if (nvar_jac .lt. 4) return
      do j = 1, N
         f_sp(j,isp_H2) = Y(nvar_jac*(j-1)+4)/max(u(1,j), 1.0d-99)
      enddo
      end subroutine unpack_carrier

      ! ------------------------------------------------------!

      subroutine eval_residual(Y, f_sp_seed, f_sp, Fvec, heat, cool, n_part, &
                               admissible, may_be_adopted)
      ! Full steady residual F(Y) with local ionization-equilibrium
      ! elimination, AND the heat/cool it used (so the caller can FREEZE the
      ! radiation when building the banded Jacobian).
      !   unpack Y -> u(1:N); Apply_BC fills ghosts;
      !   (rho,v,p) -> densities, T; refresh excited-H; ioniz_eq -> heat,cool;
      !   Apply_BC again, so the ghosts carry THAT composition;
      !   R = assemble_residual(u, n_part, heat, cool); pack R(1:N) -> Fvec.
      ! State, then the composition of that state, then the ghosts, then the
      ! fluxes -- the order the marching loop has, so that the two are the same
      ! discrete operator (sections 141.6 and 144.1).
      ! THE SEED IS AN ARGUMENT AND THE RESULT IS ANOTHER ONE.  f_sp_seed is
      ! the composition the equilibrium sweep starts from -- the cell's stored
      ! composition, or whatever initial guess the caller means -- and f_sp is
      ! where the sweep's answer is written.  They are separate dummies, and
      ! the FIRST thing this routine does is overwrite f_sp from f_sp_seed, so
      ! whatever f_sp happened to hold on entry is never read.
      !
      ! WHY THAT IS THE INTERFACE.  The sweep is a Newton solve per cell and it
      ! starts from a guess; when the guess came from whatever the caller's
      ! workspace last held, F was not a function of Y.  MEASURED on the
      ! converged 1 microbar hot Uranus (section 146.3): evaluate F at a state,
      ! at two others, then at the first again through one shared workspace and
      ! 1406 of 1500 entries differ, by 3.0e-3 of the row scale against an
      ! acceptance tolerance of 1e-5; hold the seed and the three evaluations
      ! are bitwise identical.  Every call site used to defend against that by
      ! copying the adopted composition into a workspace first, five times, by
      ! convention.  The copy is now here, once, and a caller cannot get it
      ! wrong: it has to name the seed.
      !
      ! It is the same rule the chemistry of a single cell already keeps -- a
      ! cell's root may not depend on the cell solved before it, which is what
      ! seed_species_of_cell states in constrained_chemical_equilibrium.f90 --
      ! raised one level, from the cell to the sweep.
      !
      ! n_part = n_tot + n_e is returned too, so a caller that later builds a
      ! FROZEN-radiation residual can hand the same particle count back and
      ! keep T = p/n_part (hence the transport coefficients) consistent.
      real*8, dimension(nvar_jac*N),                  intent(in)    :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in)    :: f_sp_seed
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(out)   :: f_sp
      real*8, dimension(nvar_jac*N),                  intent(out)   :: Fvec
      real*8, dimension(1-Ng:N+Ng),            intent(out)   :: heat, cool
      real*8, dimension(1-Ng:N+Ng), optional,  intent(out)   :: n_part
      ! Whether the state Y can be described at all: .false. when the
      ! equilibrium sweep met a cell whose reaction residual is not a finite
      ! number, when any cell's composition was accepted WITHOUT being a
      ! certified root of the network (acceptance class 4, the relaxation
      ! amnesty, and class 5, the constrained continuation), or when a row of
      ! the assembled residual is not finite. A
      ! caller that is PROBING -- a Krylov direction, a Jacobian column, a
      ! line-search trial -- must not put such a residual into its Newton
      ! model; the marching loop never asks and is unaffected
      ! (docs/Update_EXHALE.md section 121).
      logical, optional,                       intent(out)   :: admissible
      ! Whether the caller could ADOPT this state. A line-search or damped
      ! Gauss-Newton candidate could; a Jacobian color or a Krylov
      ! directional-derivative sample could not -- it is a point where the
      ! residual is read, never a state anyone keeps. The uncertified-
      ! composition refusal below applies only to the first kind: a derivative
      ! sampled through a cell whose chemistry the sweep could not certify is
      ! a worse derivative, but refusing it leaves the solve with no Newton
      ! model at all (13 of 17 colors zeroed, measured on the hot Uranus
      ! hand-off), and a model built from an imperfect derivative is still a
      ! model. Defaults to .true., the strict reading.
      logical, optional,                       intent(in)    :: may_be_adopted

      real*8, dimension(3,1-Ng:N+Ng) :: u, W, R
      real*8, dimension(1-Ng:N+Ng)   :: rho, v, p, T
      real*8, dimension(1-Ng:N+Ng)   :: nhi, nhii, nhei, nheii, nheiii, nheiTR
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
      real*8, dimension(1-Ng:N+Ng)   :: ne, n_tot, eta
      real*8 :: rel_change
      type(ioniz_eq_ledger) :: sweep
      integer :: irow
      logical :: gate_uncertified, headroom_ok
      real*8, dimension(1:N,n_carrier_max) :: cres, cterms
      real*8  :: crc_max, crc_vol, tscale_code, nH2j
      integer :: jc4, jj4, crc_j, crc_ic

      n_resid_eval = n_resid_eval + 1
      ! The sweep starts from the seed the caller named, and from nothing else.
      f_sp = f_sp_seed
      u = 0.0d0
      call unpack_U(Y, u)
      ! THE CARRIER UNKNOWN ENTERS HERE, and it enters as the composition the
      ! equilibrium sweep is handed. With the carrier transport on, ioniz_eq
      ! already holds H2 fixed at the partition it is given (x_h2_fix) and
      ! solves everything else around it -- that is what made the Picard
      ! splitting possible in the first place. The only change the coupled
      ! solve makes is WHERE that number comes from: not from the last
      ! transport step, but from the Newton unknown. So the sweep, the
      ! particle count it returns, the temperature p/(n_tot+n_e) that follows
      ! from it, and every heating and cooling term are all functions of
      ! Y(4) as well as of the hydro rows, and the Jacobian sees it.
      !
      ! The ghosts are NOT set from Y: the base ghost's partition is incoming
      ! data (section 117 pins it from the handoff) and the upper ghosts are
      ! extrapolated. Only the physical cells are unknowns.
      headroom_ok = .true.
      n_headroom_last = 0
      if (nvar_jac .ge. 4) then
         do jj4 = 1, N
            jc4  = nvar_jac*(jj4-1) + 4
            nH2j = Y(jc4)
            if (nH2j .lt. 0.0d0) headroom_ok = .false.
            if (nH2j*n0 .gt. carrier_headroom(jj4))                      &
               n_headroom_last = n_headroom_last + 1
            f_sp(jj4,isp_H2) = nH2j/max(u(1,jj4), 1.0d-99)
         enddo
         ! A COMPARISON, NOT A THRESHOLD, for the reason section 138 gives
         ! for the chemistry gate: the headroom is frozen at the iterate and
         ! the iterate itself can sit a cell or two outside it, and a rule
         ! that rejects the neighbourhood of the point it stands on cannot
         ! be used to leave that point. A negative density has no such
         ! excuse and is refused outright.
         headroom_ok = headroom_ok .and.                                 &
                       (n_headroom_last .le. n_headroom_iterate)
      endif
      ! Composition of the INTERIOR of Y, evaluated before the ghosts are
      ! filled: Apply_BC reads n_part_cell1 for the continuous-temperature base
      ! ghost (T(1) = p(1)/n_part_cell1) and that global is written by
      ! get_species_densities. Without this call it still holds the value left
      ! by the PREVIOUS residual evaluation, so F would depend on the previous Y
      ! as well as on Y -- and a finite-difference Jacobian column would then
      ! mix two states. Cost is one array pass; ioniz_eq below dominates.
      ! u(1,:) is already the density, so no U_to_W here: the ghosts are still
      ! zero at this point and U_to_W would divide by them.
      rho = u(1,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,         &
                                 nheiii,nheiTR,nm,ne,n_tot)
      call Apply_BC(u)           ! fill ghosts from the interior

      call U_to_W(u, W)
      rho = W(1,:);  v = W(2,:);  p = W(3,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,         &
                                 nheiii,nheiTR,nm,ne,n_tot)
      call comp_T_from_p(p, n_tot, ne, T)
      if (use_excited_H) call excited_H_update(T,rho,f_sp,v,rel_change)
      call ioniz_eq(T,rho,f_sp,heat,cool,eta,sweep)
      ! Refresh the particle count from the equilibrium fractions, exactly as
      ! the marching loop does before its transport stage, so the residual's
      ! T = p/(n_tot+n_e) is the same temperature the marching step diffuses.
      ! Only the transport terms read it; with them off nothing changes.
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,         &
                                 nheiii,nheiTR,nm,ne,n_tot)
      ! THE GHOSTS ARE WRITTEN AGAIN, WITH THE COMPOSITION THE FLUXES WILL
      ! USE. The Apply_BC above ran before the sweep, because U_to_W needs a
      ! ghost density to divide by and ioniz_eq needs a ghost pressure; the
      ! conservative ghost state it left therefore carried the composition of
      ! the PREVIOUS evaluation, while Reconstruct and RK_rhs below read it
      ! alongside interior cells carrying this one. Under a constant adiabatic
      ! index the primitive-to-conservative round trip returns the pressure the
      ! boundary condition asked for and the two cannot disagree; under the
      ! caloric EOS of section 141 they can, by as much as the ghost
      ! composition moved in one sweep -- largest during relaxation, which is
      ! when the solver is trying to make progress (section 141.6, item (AB)).
      !
      ! This is now the order the MARCHING LOOP has always had, and that is the
      ! point: there Apply_BC follows the composition refresh, so the ghosts a
      ! step's fluxes read were written with that step's composition. The
      ! residual and the marching right-hand side have to be the same discrete
      ! operator or they cannot have the same fixed point (section 144.1).
      !
      ! It costs nothing but the call: Apply_BC writes the ghosts and returns
      ! the interior unchanged, bit for bit, which is what made a second call
      ! possible at all.
      !
      ! WHAT IS STILL ONE SWEEP BEHIND, and deliberately. The ghost cells' own
      ! composition was solved by the sweep above at the ghost pressure the
      ! FIRST Apply_BC wrote, so the EOS used here is that composition and not
      ! one re-solved at the pressure just written. Iterating the pair to a
      ! ghost fixed point would make this residual something the marching loop
      ! is not, and the requirement above is the stronger one.
      call Apply_BC(u)

      call assemble_residual(u, n_tot + ne, heat, cool, R)
      call pack_R(R, Fvec)
      ! The carrier row of the same state: the steady continuity equation of
      ! n(H2), assembled by the module that owns it, on the background this
      ! evaluation's own sweep just left. Its natural units are cm^-3 s^-1
      ! while the hydrodynamic rows are per code time, so it is converted by
      ! the code's own time scale -- otherwise the merit would be comparing
      ! two different clocks.
      if (nvar_jac .ge. 4) then
         call carrier_steady_residual(rho, v, T, f_sp, crc_max, crc_j,   &
                                      crc_ic, rvol=crc_vol,              &
                                      res_out=cres, terms_out=cterms)
         tscale_code = R0/(v0*n0)
         do jj4 = 1, N
            Fvec(nvar_jac*(jj4-1)+4) = cres(jj4,ic_H2)*tscale_code
         enddo
         carrier_relnorm_last = crc_vol
         carrier_cellmax_last = crc_max
         carrier_cell_worst   = crc_j
      endif
      if (nfix_base .gt. 0) call apply_base_fix(Y, Fvec)
      if (present(n_part)) n_part = n_tot + ne
      ! Recorded on EVERY evaluation, so that a caller which does not ask for
      ! admissibility -- the one that evaluates the current iterate -- can
      ! still read off what that iterate's own chemistry cost.
      n_uncertified_last = sweep%acc_n(4) + sweep%acc_n(5)

      if (present(admissible)) then
         ! Two independent statements, both required. The sweep's own count
         ! catches a state that broke a balance row of the chemical network,
         ! which is where a probe state fails first and where the diagnosis
         ! lies; the row scan catches anything non-finite that reached the
         ! residual by another route (an overflowed flux, a cooling term).
         gate_uncertified = .true.
         if (present(may_be_adopted)) gate_uncertified = may_be_adopted
         ! The element budget is a CONSTRAINT on the state, so it is tested
         ! here and not imposed by a clamp inside the row. A trial whose
         ! n(H2) exceeds the hydrogen its cell has, or goes negative, is not
         ! a state the write-back could describe; refusing it is the same
         ! treatment the positivity of rho and p already gets in the line
         ! search.
         if (.not. headroom_ok) n_trial_headroom = n_trial_headroom + 1
         admissible = (sweep%n_nonfinite .eq. 0) .and. headroom_ok        &
                      .and. ((.not. gate_uncertified)                     &
                             .or. (n_uncertified_last                     &
                                   .le. n_uncertified_iterate))
         if (gate_uncertified .and.                                       &
             n_uncertified_last .gt. n_uncertified_iterate) then
            n_trial_uncertified  = n_trial_uncertified + 1
            uncertified_cell_max = max(uncertified_cell_max,              &
                                       n_uncertified_last)
            uncertified_res_max  = max(uncertified_res_max,               &
                 max(sweep%acc_resmax(4), sweep%acc_resmax(5)))
         endif
         if (admissible) then
            do irow = 1, nvar_jac*N
               if (.not. finite_real(Fvec(irow))) then
                  admissible = .false.;  exit
               endif
            enddo
         endif
      endif

      end subroutine eval_residual

      ! ------------------------------------------------------!

      subroutine newton_residual(Y, f_sp_seed, f_sp, Fvec)
      ! Thin wrapper: full residual, discarding the heat/cool diagnostics.
      ! Seed in, composition out, for the reason eval_residual gives.
      real*8, dimension(nvar_jac*N),                  intent(in)    :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in)    :: f_sp_seed
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(out)   :: f_sp
      real*8, dimension(nvar_jac*N),                  intent(out)   :: Fvec
      real*8, dimension(1-Ng:N+Ng) :: heat, cool
      call eval_residual(Y, f_sp_seed, f_sp, Fvec, heat, cool)
      end subroutine newton_residual

      ! ------------------------------------------------------!

      subroutine frozen_residual(Y, n_part, heat, cool, Fvec)
      ! FROZEN-radiation residual: the hydro flux/gravity residual at Y with
      ! heat/cool held FIXED (no ioniz_eq, no column-density recompute). This
      ! is strictly local (WENO3 stencil) and hence exactly banded, so the
      ! colored finite-difference Jacobian of THIS residual is exact. It is
      ! the approximate Jacobian / preconditioner for the inexact Newton: the
      ! weakly non-local radiation response is left to the outer iteration.
      real*8, dimension(nvar_jac*N),       intent(in)  :: Y
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: n_part, heat, cool
      real*8, dimension(nvar_jac*N),       intent(out) :: Fvec
      real*8, dimension(3,1-Ng:N+Ng) :: u, R
      u = 0.0d0
      call unpack_U(Y, u)
      call Apply_BC(u)
      call assemble_residual(u, n_part, heat, cool, R)
      call pack_R(R, Fvec)
      if (nfix_base .gt. 0) call apply_base_fix(Y, Fvec)
      end subroutine frozen_residual

      ! ------------------------------------------------------!

      subroutine pack_R(R, Fvec)
      ! Pack residual array (physical cells) into the flat vector.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: R
      real*8, dimension(nvar_jac*N),         intent(out) :: Fvec
      integer :: j, k
      do j = 1, N
         do k = 1, 3
            Fvec(nvar_jac*(j-1)+k) = R(k,j)
         enddo
      enddo
      end subroutine pack_R

      ! ------------------------------------------------------!

      subroutine build_banded_jac(Y, n_part, heat, cool, ab)
      ! Colored finite-difference banded Jacobian of the FROZEN residual,
      ! stored in LAPACK general-band form for dgbtrf/dgbtrs:
      !   ab(kl+ku+1 + i - j, j) = J(i,j),  i in [j-ku, j+kl]
      ! ldab = 2*kl+ku+1. Columns sharing a color (stride ncolor = kl+ku+1)
      ! have disjoint row supports, so one probe per color recovers their
      ! band entries with no cross-contamination.
      real*8, dimension(nvar_jac*N),               intent(in)  :: Y
      real*8, dimension(1-Ng:N+Ng),         intent(in)  :: n_part, heat, cool
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N), intent(out) :: ab
      real*8, dimension(nvar_jac*N) :: F0, Fp, Yp, dYc
      integer :: neq, color, jcol, irow, ilo, ihi
      real*8  :: sqeps

      neq   = nvar_jac*N
      sqeps = sqrt(epsilon(1.0d0))
      ab    = 0.0d0

      call frozen_residual(Y, n_part, heat, cool, F0)

      do color = 1, ncolor_jac
         Yp  = Y
         dYc = 0.0d0
         do jcol = color, neq, ncolor_jac
            dYc(jcol) = sqeps*max(abs(Y(jcol)), 1.0d0)
            Yp(jcol)  = Y(jcol) + dYc(jcol)
         enddo
         call frozen_residual(Yp, n_part, heat, cool, Fp)
         do jcol = color, neq, ncolor_jac
            ilo = max(1,   jcol - ku_jac)
            ihi = min(neq, jcol + kl_jac)
            do irow = ilo, ihi
               ab(kl_jac+ku_jac+1 + irow - jcol, jcol) =                &
                    (Fp(irow) - F0(irow))/dYc(jcol)
            enddo
         enddo
      enddo

      end subroutine build_banded_jac

      ! ------------------------------------------------------!

      subroutine build_banded_jac_full(Y, f_sp_base, ab, n_color_unresolved)
      ! Colored-FD banded Jacobian of the FULL residual (eval_residual, i.e.
      ! including the local ionization-equilibrium + heat/cool response). The
      ! frozen version omits the dominant local d(heat-cool)/dE coupling, so
      ! Newton has no descent direction even near the solution; this version
      ! captures it. The weakly non-local column-density coupling slightly
      ! contaminates same-color columns (treated as preconditioner error).
      ! Each probe restores f_sp from the base copy (ioniz_eq mutates it).
      real*8, dimension(nvar_jac*N),                  intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N), intent(out) :: ab
      ! Colors whose probe stayed inadmissible down to the smallest step, so
      ! their columns are left at zero (see n_probe_step_halvings above).
      integer, intent(out) :: n_color_unresolved
      real*8, dimension(nvar_jac*N) :: F0, Fp, Yp, dYc, Dsc
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng) :: heat, cool
      integer :: neq, color, jcol, irow, ilo, ihi, ish
      real*8  :: sqeps, hstep
      logical :: probe_ok

      neq   = nvar_jac*N
      sqeps = sqrt(epsilon(1.0d0))
      ab    = 0.0d0
      n_color_unresolved = 0

      ! FD column steps RELATIVE to the scale for each unknown (a flat
      ! max(|Y|,1) floor gives wind cells, ~1e-7 of the base in code
      ! units, order-unity relative kicks and garbage columns).
      call cell_state_scales(Y, Dsc)

      call eval_residual(Y, f_sp_base, fwork, F0, heat, cool,             &
                         admissible=probe_ok,                             &
                         may_be_adopted=.false.)
      if (.not. probe_ok) then
         ! The base point of the differences is itself indescribable. No
         ! column of this Jacobian means anything, so none is built and the
         ! caller sees every color unresolved.
         n_color_unresolved = ncolor_jac
         return
      endif

      do color = 1, ncolor_jac
         hstep    = sqeps
         probe_ok = .false.
         do ish = 0, n_probe_step_halvings
            Yp  = Y;  dYc = 0.0d0
            do jcol = color, neq, ncolor_jac
               dYc(jcol) = hstep*Dsc(jcol)
               Yp(jcol)  = Y(jcol) + dYc(jcol)
            enddo
            call eval_residual(Yp, f_sp_base, fwork, Fp, heat, cool,      &
                               admissible=probe_ok,                       &
                               may_be_adopted=.false.)
            if (probe_ok) exit
            hstep = 0.5d0*hstep
         enddo
         if (.not. probe_ok) then
            n_color_unresolved = n_color_unresolved + 1
            cycle
         endif
         do jcol = color, neq, ncolor_jac
            ilo = max(1,   jcol - ku_jac)
            ihi = min(neq, jcol + kl_jac)
            do irow = ilo, ihi
               ab(kl_jac+ku_jac+1 + irow - jcol, jcol) =                &
                    (Fp(irow) - F0(irow))/dYc(jcol)
            enddo
         enddo
      enddo

      end subroutine build_banded_jac_full

      ! ------------------------------------------------------!

      subroutine band_matvec(ab, x, y)
      ! y = J*x for the band-stored Jacobian ab (same layout as above).
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N), intent(in)  :: ab
      real*8, dimension(nvar_jac*N),                    intent(in)  :: x
      real*8, dimension(nvar_jac*N),                    intent(out) :: y
      integer :: neq, irow, jcol, jlo, jhi
      neq = nvar_jac*N
      do irow = 1, neq
         y(irow) = 0.0d0
         jlo = max(1,   irow - kl_jac)
         jhi = min(neq, irow + ku_jac)
         do jcol = jlo, jhi
            y(irow) = y(irow) + ab(kl_jac+ku_jac+1 + irow - jcol, jcol)*x(jcol)
         enddo
      enddo
      end subroutine band_matvec

      ! ------------------------------------------------------!

      subroutine band_matvec_transpose(ab, x, y)
      ! y = J^T*x for the band-stored Jacobian ab (same layout as above).
      ! The transpose is needed because the gradient of the least-squares
      ! merit 1/2*||F||^2 is J^T F, and the Newton/PTC step does not use it.
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N), intent(in)  :: ab
      real*8, dimension(nvar_jac*N),                    intent(in)  :: x
      real*8, dimension(nvar_jac*N),                    intent(out) :: y
      integer :: neq, irow, jcol, jlo, jhi
      neq = nvar_jac*N
      y = 0.0d0
      do irow = 1, neq
         jlo = max(1,   irow - kl_jac)
         jhi = min(neq, irow + ku_jac)
         do jcol = jlo, jhi
            y(jcol) = y(jcol) + ab(kl_jac+ku_jac+1 + irow - jcol, jcol)*x(irow)
         enddo
      enddo
      end subroutine band_matvec_transpose

      ! ------------------------------------------------------!

      subroutine normal_equations_band(ab, bn)
      ! bn = J^T J in band form, from the band-stored J. A band of
      ! half-width kl_jac squared has half-width 2*kl_jac, and the product is
      ! symmetric, so bn is stored in the same LAPACK general-band layout
      ! with kl = ku = kl_jac_normal.
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N),           intent(in)  :: ab
      real*8, dimension(ldab_normal, nvar_jac*N),           intent(out) :: bn
      integer :: neq, i, j, k, klo, khi
      real*8  :: acc
      neq = nvar_jac*N
      bn  = 0.0d0
      do j = 1, neq
         do i = max(1, j-kl_jac_normal), min(neq, j+ku_jac_normal)
            acc = 0.0d0
            klo = max(1,   max(i-ku_jac, j-ku_jac))
            khi = min(neq, min(i+kl_jac, j+kl_jac))
            do k = klo, khi
               acc = acc + ab(kl_jac+ku_jac+1+k-i, i)*ab(kl_jac+ku_jac+1+k-j, j)
            enddo
            bn(kl_jac_normal+ku_jac_normal+1+i-j, j) = acc
         enddo
      enddo
      end subroutine normal_equations_band

      ! ------------------------------------------------------!

      subroutine levenberg_marquardt_step(bn, grad_merit, mu_damp, dZ, ok)
      ! One damped Gauss-Newton (Levenberg-Marquardt) step of the scaled
      ! model:  (J^T J + mu_damp*diag(J^T J)) dZ = -grad_merit,
      ! grad_merit = J^T s.  The matrix is symmetric positive definite for
      ! every mu_damp > 0 whenever diag(J^T J) > 0, so
      ! grad^T dZ = -grad_merit^T (.)^-1 grad_merit < 0: the step is a
      ! descent direction of the merit for ANY damping. mu_damp is the trust
      ! parameter -- mu_damp -> infinity recovers the scaled steepest-descent
      ! direction, mu_damp -> 0 the Gauss-Newton step. The Marquardt scaling
      ! (mu_damp times the diagonal, not mu_damp times the identity) makes
      ! the damping independent of the units of the individual unknowns.
      real*8, dimension(ldab_normal, nvar_jac*N), intent(in)  :: bn
      real*8, dimension(nvar_jac*N),                    intent(in)  :: grad_merit
      real*8,                                    intent(in)  :: mu_damp
      real*8, dimension(nvar_jac*N),                    intent(out) :: dZ
      logical,                                   intent(out) :: ok
      real*8,  allocatable :: bf(:,:)
      integer, allocatable :: ip(:)
      integer :: neq, ld, j, lpinfo
      neq = nvar_jac*N;  ld = ldab_normal
      allocate(bf(ld,neq), ip(neq))
      bf = bn
      do j = 1, neq
         bf(kl_jac_normal+ku_jac_normal+1, j) =                          &
              bf(kl_jac_normal+ku_jac_normal+1, j)                        &
            + mu_damp*abs(bn(kl_jac_normal+ku_jac_normal+1, j))
      enddo
      call dgbtrf(neq, neq, kl_jac_normal, ku_jac_normal, bf, ld, ip, lpinfo)
      if (lpinfo .ne. 0) then
         dZ = 0.0d0;  ok = .false.
         deallocate(bf, ip);  return
      endif
      dZ = -grad_merit
      call dgbtrs('N', neq, kl_jac_normal, ku_jac_normal, 1, bf, ld, ip,  &
                  dZ, neq, lpinfo)
      ok = (lpinfo .eq. 0)
      if (.not. ok) dZ = 0.0d0
      deallocate(bf, ip)
      end subroutine levenberg_marquardt_step

      ! ------------------------------------------------------!

      subroutine levenberg_marquardt_descent(Y, F, f_sp_base, D, Drow,    &
                                             ab, f2,                      &
                                             dY, f2_new, mu_used, gnorm,  &
                                             found)
      ! Find a step that LOWERS the merit when the Newton/PTC direction
      ! raises it.
      !
      ! The merit the solve minimizes is m(x) = 1/2*||s||^2 with s = F/D and
      ! x the scaled unknowns (Y = D*x); its gradient is A^T s with
      ! A = D^-1 J D, which the band-stored ab already holds. The Newton/PTC
      ! direction solves (I/dtau + A) dZ = -s -- an implicit Euler step, not a
      ! minimizer of m -- and on the He/H = 0.3 hand-off state it points
      ! UPHILL at every damping, including the dtau -> 0 limit dZ = -s, which
      ! is uphill exactly when s^T A s < 0 (docs/Update_EXHALE.md section
      ! 126). The damped Gauss-Newton family above descends for every mu, so
      ! what is left is choosing mu, and that is a one-dimensional search
      ! along the Levenberg-Marquardt curve: walk mu down by decades from a
      ! near-steepest-descent value while the TRUE merit keeps improving, and
      ! return the best step seen.
      !
      ! The trials are measured with weno_mode = 0, the same true merit the
      ! backtracking line search uses, so the two are directly comparable.
      real*8, dimension(nvar_jac*N),                    intent(in)  :: Y, F, D
      real*8, dimension(nvar_jac*N),                    intent(in)  :: Drow
      real*8, dimension(1-Ng:N+Ng,n_species),    intent(in)  :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1, nvar_jac*N), intent(in)  :: ab
      real*8,                                    intent(in)  :: f2
      real*8, dimension(nvar_jac*N),                    intent(out) :: dY
      real*8,                                    intent(out) :: f2_new
      real*8,                                    intent(out) :: mu_used
      ! Norm of the merit gradient in the scaled space, ||A^T s||. Reported so
      ! that a search that finds nothing can be told apart from a stationary
      ! point, where there is nothing to find.
      real*8,                                    intent(out) :: gnorm
      logical,                                   intent(out) :: found
      real*8, allocatable :: bn(:,:)
      real*8, dimension(nvar_jac*N) :: s, grad_merit, dZ, dYtry, Ytry, Ftry
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng)   :: heat, cool
      real*8, dimension(3,1-Ng:N+Ng) :: utry, Wtry
      real*8  :: mu_damp, f2t, f2prev
      integer :: im, saved_mode
      logical :: okstep, okres
      ! Damping sweep: mu_start is far enough on the steepest-descent side
      ! that the first step is short and safe, and n_mu decades reach the
      ! Gauss-Newton end. Measured on the He/H = 0.3 hand-off state, the merit
      ! falls monotonically from mu_damp = 1e2 down to 1e-4 (0.514 -> 0.260)
      ! and rises again at 1e-6.
      real*8,  parameter :: mu_start = 1.0d2
      integer, parameter :: n_mu     = 10
      dY = 0.0d0;  f2_new = f2;  mu_used = 0.0d0;  found = .false.
      saved_mode = weno_mode
      allocate(bn(ldab_normal, nvar_jac*N))
      if (nvar_jac .ge. 4) then
         s = F/Drow
      else
         s = F/D
      endif
      call band_matvec_transpose(ab, s, grad_merit)
      gnorm = sqrt(sum(grad_merit*grad_merit))
      call normal_equations_band(ab, bn)
      f2prev = huge(1.0d0)
      mu_damp = mu_start
      do im = 1, n_mu
         call levenberg_marquardt_step(bn, grad_merit, mu_damp, dZ, okstep)
         if (okstep) then
            dYtry = D*dZ
            Ytry  = Y + dYtry
            call unpack_U(Ytry, utry)
            call U_to_W_interior(utry, Wtry)
            if (minval(Wtry(1,1:N)) .gt. 0.0d0 .and.                     &
                minval(Wtry(3,1:N)) .gt. 0.0d0) then
               weno_mode = 0
               call eval_residual(Ytry, f_sp_base, fwork, Ftry, heat,     &
                                  cool,                                   &
                                  admissible=okres)
               if (nvar_jac .ge. 4) then
                  f2t = sqrt(sum((Ftry/Drow)**2))
               else
                  f2t = sqrt(sum((Ftry/D)**2))
               endif
               if (okres .and. f2t .lt. f2_new) then
                  dY = dYtry;  f2_new = f2t;  mu_used = mu_damp;  found = .true.
               endif
               ! The merit along the LM curve falls, reaches a minimum and
               ! rises again; once it has turned, smaller mu only overshoots
               ! further, so the sweep stops at the turn.
               if (okres .and. f2t .gt. f2prev) exit
               if (okres) f2prev = f2t
            endif
         endif
         mu_damp = 0.1d0*mu_damp
      enddo
      weno_mode = saved_mode
      deallocate(bn)
      end subroutine levenberg_marquardt_descent

      ! ------------------------------------------------------!


      ! ------------------------------------------------------!

      subroutine resid_relnorm(F, u, rc, rnorm)
      ! THE CONVERGENCE MEASURE of the steady solve: the relative residual of
      ! each row over the WHOLE physical column [1:N], every cell divided by
      ! residual_row_scale (steady_residual.f90) -- |u| for the mass and
      ! energy rows, |rho v| for momentum in the wind and the gravitational
      ! force density for momentum below the escape radius, where the layer
      ! is quasi-hydrostatic and |rho v| carries no information.
      !
      ! The wind [j_min:N] and the layer [1:j_min-1] are normed SEPARATELY
      ! and combined by the LARGER of the two. Until section 127 this measure
      ! looked at [j_min:N] alone, i.e. at r > r_esc, and the marching du it
      ! takes over from still does; the column between the base and r_esc was
      ! then controlled by nothing, and two solves of one hand-off state both
      ! returned info = 0 on winds differing by 0.22 dex in the rate, the
      ! WORSE of them scoring the better ||R||. Combining by the maximum,
      ! rather than by a volume- or cell-count-weighted sum, is what keeps
      ! that from recurring: the outer region carries almost all of
      ! sum r^2 dr, so any sum lets it average the inner column away.
      real*8, dimension(nvar_jac*N),         intent(in)  :: F
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3),           intent(out) :: rc
      real*8,                         intent(out) :: rnorm
      real*8, dimension(3,1-Ng:N+Ng) :: Rres
      real*8, dimension(3)           :: rc_wind, rc_layer
      integer :: j, k
      Rres = 0.0d0
      do j = 1, N
         do k = 1, 3
            Rres(k,j) = F(nvar_jac*(j-1)+k)
         enddo
      enddo
      call relnorm_over_cells(Rres, u, j_min, N,       rc_wind)
      call relnorm_over_cells(Rres, u, 1,     j_min-1, rc_layer)
      rc    = max(rc_wind, rc_layer)
      rnorm = maxval(rc)
      ! THE CARRIER ROW IS NOT FOLDED IN HERE, and that is the same
      ! judgment section 133 made about the flux gate: one number cannot say
      ! two things. `Resid tol` is calibrated against a row scale that BOUNDS
      ! each hydrodynamic row's largest term, and converged states sit at
      ! 1e-6 on it; the carrier row is measured against the terms
      ! themselves, where a converged row sits near 1e-3. Taking the maximum
      ! would hold the carrier row to a tolerance meant for a different
      ! scale. It is a third gate instead (steady_gates_met).
      end subroutine resid_relnorm

      ! ------------------------------------------------------!

      subroutine resid_relnorm_below_escape(F, u, rc, j_worst, r_worst, &
                                            k_worst)
      ! The relative residual of the layer [1:j_min-1] BELOW the escape
      ! radius on its own, and the cell in it that carries the largest
      ! cell-wise scaled residual.
      !
      ! Since section 127 this layer is one of the two halves of the
      ! convergence measure (resid_relnorm returns the larger of it and the
      ! wind), so what this adds to that number is WHERE in the layer the
      ! residual sits -- which the volume-weighted sum averages away -- and
      ! how the layer stands against the wind at the same iterate.
      !
      ! The row scales are residual_row_scale's: |u| for mass and energy,
      ! and for momentum the gravitational force density
      !    s_grav(j) = |u(1,j)| * |Gphi_i(j) - Gphi_i(j-1)| / dr_j(j)
      ! (the same discrete potential difference the source term uses;
      ! Source.f90), so rc(2) is the fractional violation of hydrostatic
      ! balance in the layer.
      real*8, dimension(nvar_jac*N),         intent(in)  :: F
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3),           intent(out) :: rc
      real*8,                         intent(out) :: r_worst
      integer,                        intent(out) :: j_worst, k_worst
      real*8, dimension(3,1-Ng:N+Ng) :: Rres
      integer :: j, k
      real*8  :: cellrel, amx, uscale

      rc = 0.0d0;  j_worst = 0;  k_worst = 0;  r_worst = 0.0d0
      if (j_min .le. 1) return          ! no cells below the escape radius

      Rres = 0.0d0
      do j = 1, j_min-1
         do k = 1, 3
            Rres(k,j) = F(nvar_jac*(j-1)+k)
         enddo
      enddo
      call relnorm_over_cells(Rres, u, 1, j_min-1, rc)

      amx = -1.0d0
      do j = 1, j_min-1
         do k = 1, 3
            uscale  = residual_row_scale(k, j, u)
            cellrel = abs(F(nvar_jac*(j-1)+k))/max(uscale, 1.0d-30)
            if (cellrel .gt. amx) then
               amx = cellrel;  j_worst = j;  k_worst = k
            endif
         enddo
      enddo
      if (j_worst .gt. 0) r_worst = r(j_worst)

      end subroutine resid_relnorm_below_escape

      ! ------------------------------------------------------!

      subroutine write_resid_below_escape(tag, F, u)
      ! One line for the layer below the escape radius, which is one of the
      ! two halves of the convergence measure (see resid_relnorm_below_escape
      ! and resid_relnorm). Silent when the layer is empty. (F, u) must be a
      ! consistent pair: the caller passes the residual of the state
      ! currently held in u. The three numbers are the volume-weighted row
      ! values: mass and energy relative to |u| (rates), momentum relative to
      ! the gravitational force density (fractional violation of hydrostatic
      ! balance in the layer). The trailing cell is where the largest
      ! cell-wise scaled residual of the layer sits.
      character(*),                   intent(in) :: tag
      real*8, dimension(nvar_jac*N),         intent(in) :: F
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      real*8  :: rc_below(3), r_worst
      integer :: j_worst, k_worst
      if (j_min .le. 1) return
      call resid_relnorm_below_escape(F, u, rc_below, j_worst, r_worst,   &
                                      k_worst)
      write(*,'(A,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,I0,A,F8.4,A,I0)')     &
           tag, ' below r_esc [1:',j_min-1,                               &
           '] |R|: mass=',rc_below(1), ' mom/grav=',rc_below(2),          &
           ' energy=',rc_below(3),                                        &
           '  max cell j=',j_worst,' r=',r_worst,' k=',k_worst
      end subroutine write_resid_below_escape

      ! ------------------------------------------------------!

      subroutine reset_uncertified_trial_count
      n_trial_uncertified  = 0
      uncertified_cell_max = 0
      uncertified_res_max  = 0.0d0
      n_uncertified_iterate = 0
      end subroutine reset_uncertified_trial_count

      ! ------------------------------------------------------!

      subroutine write_uncertified_trial_count(tag)
      ! What the solve refused, and how far outside a root those states were.
      ! Silent when every probe and trial landed on a certified composition.
      character(*), intent(in) :: tag
      if (n_trial_uncertified .le. 0) return
      write(*,'(A,A,I0,A,I0,A,ES10.3)') tag,                             &
           ' trial states refused for an uncertified composition: ',      &
           n_trial_uncertified, '; worst cell count ', uncertified_cell_max,&
           ', worst reaction residual ', uncertified_res_max
      end subroutine write_uncertified_trial_count

      ! ------------------------------------------------------!

      subroutine solve_steady_ptc(u, f_sp, resid_tol, maxit, dtau0, info)
      ! NOT AVAILABLE WITH THE CARRIER UNKNOWN. This route builds its
      ! Jacobian from frozen_residual, which by construction does no
      ! equilibrium sweep and therefore has no carrier row; giving it one
      ! would mean assembling that row from a background nothing in this
      ! routine refreshes. The coupled solve is the JFNK route, and the
      ! caller is told so rather than being given a silently 3-row answer.
      ! Pseudo-transient-continuation inexact Newton solve of the steady
      ! residual F(Y)=0. Per iteration:
      !   F = eval_residual(Y)                       (full, nonlocal radiation)
      !   J = build_banded_jac (frozen radiation)    (banded approx / precond)
      !   (I/dtau + J) dY = -F   via dgbtrf/dgbtrs    (LAPACK banded LU)
      !   Y <- Y + lam*dY        (backtracking line search + positivity)
      !   dtau <- dtau * rnorm_old/rnorm_new (SER ramp; cut on failure)
      ! As dtau->inf this is Newton; small dtau behaves like explicit
      ! relaxation, giving the robust startup PTC is designed for.
      real*8, dimension(3,1-Ng:N+Ng),         intent(inout) :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8,  intent(in)  :: resid_tol, dtau0
      integer, intent(in)  :: maxit
      integer, intent(out) :: info

      integer :: neq, ldab, iter, ls, lpinfo, jc
      real*8, allocatable :: Y(:), F(:), Ftry(:), dY(:), Ytry(:)
      real*8, allocatable :: ab(:,:), abf(:,:)
      integer, allocatable :: ipiv(:)
      real*8, dimension(1-Ng:N+Ng)            :: heat0, cool0
      real*8, dimension(3,1-Ng:N+Ng)          :: utry, Wtry
      real*8, dimension(1-Ng:N+Ng,n_species)  :: f_sp_j
      real*8  :: rnorm, rnorm_try, rc(3), dtau, lam, f2, f2_try
      real*8  :: fspread
      integer :: n_jac_unresolved
      logical :: ok, try_ok

      if (nvar_jac .ge. 4) then
         write(*,'(A)') ' (PTC) the carrier unknown is a JFNK-only'//&
              ' route; this solve is not attempted'
         info = 2;  return
      endif
      neq  = nvar_jac*N
      ldab = 2*kl_jac + ku_jac + 1
      allocate(Y(neq), F(neq), Ftry(neq), dY(neq), Ytry(neq))
      allocate(ab(ldab,neq), abf(ldab,neq), ipiv(neq))
      dtau = dtau0
      info = 1
      call reset_uncertified_trial_count

      call pack_U(u, Y)
      ! The sweeps below evaluate states the SOLVER invents, not states the
      ! run holds; tag them so that the acceptance ledgers, the non-root
      ! streak and its stop keep the two apart (section 121).
      call set_ioniz_eq_sweep_state_kind(ieq_state_steady_iterate)
      ! Seed from the adopted composition, write the swept one into the work
      ! array, then ADOPT it. The copy is an adoption and not a seeding
      ! convention: the seed is named in the call.
      call eval_residual(Y, f_sp, f_sp_j, F, heat0, cool0)
      f_sp = f_sp_j
      n_uncertified_iterate = n_uncertified_last
      call keep_background_of_adopted_state
      call resid_relnorm(F, u, rc, rnorm)
      f2 = sqrt(sum(F*F))                ! smooth line-search merit
      write(*,'(A,ES11.3,A,3ES10.2,A,ES10.2)') ' (PTC) start ||R||=',rnorm, &
           '  R(m,p,E)=',rc,'  ||F||2=',f2
      call write_resid_below_escape(' (PTC)', F, u)

      do iter = 1, maxit
         if (steady_gates_met(rnorm, u, resid_tol, fspread,        &
                          nvar_jac .ge. 4, carrier_relnorm_last)) then
            info = 0;  exit
         endif

         ! Full-residual banded Jacobian (includes the local source
         ! derivative d(heat-cool)/dE), then the PTC system (I/dtau + J).
         call set_ioniz_eq_sweep_state_kind(ieq_state_steady_candidate)
         call build_banded_jac_full(Y, f_sp, ab, n_jac_unresolved)
         if (n_jac_unresolved .gt. 0)                                    &
            write(*,'(A,I0,A,I0,A)') ' (PTC) Jacobian: ',                &
                 n_jac_unresolved, ' of ', ncolor_jac,                   &
                 ' color(s) had no admissible probe; their columns are'//&
                 ' zero'
         abf = ab
         do jc = 1, neq
            abf(kl_jac+ku_jac+1, jc) = abf(kl_jac+ku_jac+1, jc) + 1.0d0/dtau
         enddo
         dY = -F
         call dgbtrf(neq, neq, kl_jac, ku_jac, abf, ldab, ipiv, lpinfo)
         if (lpinfo .ne. 0) then
            write(*,'(A,I0,A)') ' (PTC) dgbtrf info=',lpinfo,            &
                 ' (singular); cutting dtau'
            dtau = max(dtau*0.25d0, dtau0*1.0d-3);  cycle
         endif
         call dgbtrs('N', neq, kl_jac, ku_jac, 1, abf, ldab, ipiv,      &
                     dY, neq, lpinfo)

         ! Backtracking line search on the SMOOTH merit ||F||_2 (the
         ! Newton direction reduces this; the component-wise max-relative
         ! rnorm is non-smooth and rejected valid steps). Armijo condition
         ! with positivity (rho>0, p>0).
         lam = 1.0d0;  ok = .false.
         do ls = 1, 20
            Ytry = Y + lam*dY
            call unpack_U(Ytry, utry)
            call U_to_W_interior(utry, Wtry)   ! ghosts of utry are not set
            if (minval(Wtry(1,1:N)) .gt. 0.0d0 .and.                    &
                minval(Wtry(3,1:N)) .gt. 0.0d0) then
               call eval_residual(Ytry, f_sp, f_sp_j, Ftry, heat0, cool0, &
                                  admissible=try_ok)
               f2_try = sqrt(sum(Ftry*Ftry))
               ! A trial the code cannot describe is not a descent
               ! candidate whatever its merit evaluates to; it is rejected
               ! and the step is halved, exactly as a negative density is.
               if (try_ok .and. f2_try .lt. (1.0d0 - 1.0d-4*lam)*f2) then
                  ok = .true.;  exit
               endif
            endif
            lam = 0.5d0*lam
         enddo

         if (ok) then
            Y = Ytry;  F = Ftry;  f_sp = f_sp_j
            n_uncertified_iterate = n_uncertified_last
            ! The last residual evaluation was this accepted trial, so the
            ! frozen background now describes the state just adopted.
            call keep_background_of_adopted_state
            call unpack_U(Y, u);  call Apply_BC(u)
            call resid_relnorm(F, u, rc, rnorm)
            ! SER ramp on the merit ratio, scaled by the accepted step.
            dtau = min(dtau*max(lam,0.1d0)*(f2/max(f2_try,1.0d-30)),     &
                       1.0d14*dtau0)
            f2 = f2_try
         else
            ! No descent found: shrink the pseudo-time step (more
            ! relaxation-like) but keep it >= the explicit-stable dtau0.
            dtau = max(dtau*0.25d0, dtau0)
         endif

         write(*,'(A,I4,A,ES11.3,A,ES10.2,A,ES10.2,A,ES9.2)') ' (PTC) it', &
              iter,'  ||R||=',rnorm,'  ||F||2=',f2,'  dtau=',dtau,'  lam=',lam
      enddo

      call unpack_U(Y, u);  call Apply_BC(u)
      call resid_relnorm(F, u, rc, rnorm)
      call set_ioniz_eq_sweep_state_kind(ieq_state_marching)
      call install_background_of_adopted_state
      write(*,'(A,I0,A,ES11.3)') ' (PTC) done info=',info,' ||R||=',rnorm
      call write_uncertified_trial_count(' (PTC)')
      call write_resid_below_escape(' (PTC)', F, u)

      gate_rnorm_accepted = rnorm;  gate_fspread_accepted = fspread
      deallocate(Y,F,Ftry,dY,Ytry,ab,abf,ipiv)
      end subroutine solve_steady_ptc

      ! ------------------------------------------------------!

      subroutine jv_product(Y, F0, f_sp_base, v, Jv, ok)
      ! Matrix-free Jacobian-vector product J*v by a forward directional
      ! finite difference of the FULL residual (captures the non-local
      ! radiation coupling the banded Jacobian omits):
      !   J*v ~= ( F(Y + eps*v) - F0 ) / eps,   F0 = F(Y).
      ! eps is the standard scaled step. f_sp is restored from the base
      ! copy (eval_residual mutates it).
      real*8, dimension(nvar_jac*N),                  intent(in)  :: Y, F0, v
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in)  :: f_sp_base
      real*8, dimension(nvar_jac*N),                  intent(out) :: Jv
      ! .false. when no step along v, down to the smallest of
      ! n_probe_step_halvings, lands on a state the code can describe: then
      ! there is no directional derivative to return and Jv is meaningless
      ! (it is set to zero so that nothing non-finite can escape) -- the
      ! caller must stop using v, not use the value.
      logical,                                 intent(out) :: ok
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(nvar_jac*N) :: Fp
      real*8, dimension(1-Ng:N+Ng) :: heat, cool
      real*8  :: vn, eps
      integer :: ish
      Jv = 0.0d0
      ok = .true.
      vn = sqrt(sum(v*v))
      if (vn .le. 0.0d0) return
      eps = sqrt(epsilon(1.0d0))*(1.0d0 + sqrt(sum(Y*Y)))/vn
      do ish = 0, n_probe_step_halvings
         call eval_residual(Y + eps*v, f_sp_base, fwork, Fp, heat, cool,   &
                            admissible=ok,                                 &
                         may_be_adopted=.false.)
         if (ok) exit
         eps = 0.5d0*eps
      enddo
      if (.not. ok) return
      Jv = (Fp - F0)/eps
      end subroutine jv_product

      ! ------------------------------------------------------!

      subroutine pgmres(Y, F0, f_sp_base, D, Drow, abf, ipiv, idtau, b, x, m, &
                        rtol, gm_iters, gm_truncated, Ax)
      ! Right-preconditioned GMRES(m), single cycle, for the JFNK step, in
      ! the DIAGONALLY SCALED space:
      !   solve  A_z x = b,   A_z v = v*idtau + D^-1 * J*(D v),
      ! (J*(Dv) matrix-free; b and x are scaled quantities; the caller maps
      ! dY = D*x back). Preconditioner M = factored banded SCALED system
      ! (I/dtau + D^-1 J_banded D) supplied in abf.
      real*8, dimension(nvar_jac*N),  intent(in)  :: Y, F0, b, D, Drow
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in) :: abf
      integer, dimension(nvar_jac*N), intent(in)  :: ipiv
      real*8,  intent(in)  :: idtau, rtol
      integer, intent(in)  :: m
      real*8, dimension(nvar_jac*N), intent(out) :: x
      ! The SCALED operator applied to the step it returns, A_z x, taken from
      ! the Arnoldi relation A_z M^-1 V_k = V_{k+1} Hbar_k rather than from a
      ! fresh product, so it costs no residual evaluation. The trust region
      ! needs it for the predicted reduction, and needs it to be the SAME
      ! model the Krylov step was built from. Hbar is the Hessenberg BEFORE
      ! the Givens rotations, which is why it is kept separately below.
      ! OPTIONAL: the three-unknown route never asks for it.
      real*8, dimension(nvar_jac*N), intent(out), optional :: Ax
      integer, intent(out) :: gm_iters
      ! .true. when the cycle stopped early because the next Krylov
      ! direction could not be sampled at all (jv_product returned ok =
      ! .false.). The step is then built from the subspace that WAS
      ! constructed; with not one direction available it is the zero step,
      ! and the line search reports no descent, which is the honest outcome.
      logical, intent(out) :: gm_truncated

      integer :: neq, i, j, kk, lpinfo, kv_set
      real*8, allocatable :: V(:,:), Hs(:,:), Hb(:,:), gg(:), cs(:), sn(:), yy(:)
      real*8, allocatable :: z(:), w(:), u(:)
      real*8 :: beta, hij, nrm, denom, tmp
      logical :: jv_ok

      neq = nvar_jac*N
      allocate(V(neq,m+1), Hs(m+1,m), gg(m+1), cs(m), sn(m), yy(m))
      allocate(z(neq), w(neq), u(neq))
      ! Allocated last and on its own, so that every array this cycle has
      ! always had keeps the address it has always had (section 139.9).
      allocate(Hb(m+1,m))
      Hs = 0.0d0; gg = 0.0d0; cs = 0.0d0; sn = 0.0d0; x = 0.0d0
      yy = 0.0d0
      ! V is deliberately NOT zeroed: touching an array the cycle has always
      ! had is what moved two atomic cases by 1e-9 in section 139.9. The Ax
      ! sum below skips any column the Arnoldi never wrote, which is only
      ! ever a column whose coefficient is exactly zero anyway.
      Hb = 0.0d0;  kv_set = 1
      if (present(Ax)) Ax = 0.0d0
      gm_truncated = .false.

      ! r0 = b - A*0 = b
      beta = sqrt(sum(b*b))
      gm_iters = 0
      if (beta .le. 0.0d0) then
         deallocate(V,Hs,Hb,gg,cs,sn,yy,z,w,u);  return
      endif
      V(:,1) = b/beta
      gg(1)  = beta

      do j = 1, m
         gm_iters = j
         ! z = M^{-1} V(:,j)   (right preconditioning)
         z = V(:,j)
         call dgbtrs('N', neq, kl_jac, ku_jac, 1, abf, 2*kl_jac+ku_jac+1, &
                     ipiv, z, neq, lpinfo)
         ! w = A_z z = z*idtau + D^-1 * J*(D z)
         call jv_product(Y, F0, f_sp_base, D*z, w, jv_ok)
         if (.not. jv_ok) then
            ! No admissible sample of the residual along this direction:
            ! stop the cycle here rather than feed the Arnoldi recursion a
            ! product that does not exist.
            gm_iters     = j - 1
            gm_truncated = .true.
            exit
         endif
         ! THE THREE-UNKNOWN BRANCH IS THE ORIGINAL EXPRESSION, CHARACTER
         ! FOR CHARACTER. Drow holds exactly D's values when there is no
         ! carrier row, so the two branches are the same arithmetic -- but
         ! not the same instruction sequence, and a reduction summed in a
         ! different order differs in its last bits. A change with no
         ! intended effect has to be byte-identical (section 139.9).
         if (nvar_jac .ge. 4) then
            w = w/Drow + idtau*(D/Drow)*z
         else
            w = w/D + idtau*z
         endif
         ! Modified Gram-Schmidt
         do i = 1, j
            Hs(i,j) = sum(w*V(:,i))
            w = w - Hs(i,j)*V(:,i)
         enddo
         Hs(j+1,j) = sqrt(sum(w*w))
         if (Hs(j+1,j) .gt. 0.0d0) then
            V(:,j+1) = w/Hs(j+1,j)
            kv_set   = j + 1
         endif
         ! The Arnoldi column as built, kept before the rotations overwrite it.
         Hb(1:j+1,j) = Hs(1:j+1,j)
         ! Apply previous Givens rotations to column j
         do i = 1, j-1
            tmp       =  cs(i)*Hs(i,j) + sn(i)*Hs(i+1,j)
            Hs(i+1,j) = -sn(i)*Hs(i,j) + cs(i)*Hs(i+1,j)
            Hs(i,j)   =  tmp
         enddo
         ! New Givens rotation to zero Hs(j+1,j)
         denom  = sqrt(Hs(j,j)**2 + Hs(j+1,j)**2)
         if (denom .le. 0.0d0) denom = 1.0d0
         cs(j)  = Hs(j,j)/denom
         sn(j)  = Hs(j+1,j)/denom
         Hs(j,j)   = cs(j)*Hs(j,j) + sn(j)*Hs(j+1,j)
         Hs(j+1,j) = 0.0d0
         gg(j+1) = -sn(j)*gg(j)
         gg(j)   =  cs(j)*gg(j)
         if (abs(gg(j+1)) .le. rtol*beta) exit
      enddo

      ! Back-substitute H(1:kk,1:kk) yy = gg(1:kk)
      kk = gm_iters
      do i = kk, 1, -1
         tmp = gg(i)
         do j = i+1, kk
            tmp = tmp - Hs(i,j)*yy(j)
         enddo
         yy(i) = tmp/Hs(i,i)
      enddo
      ! u = V(:,1:kk) yy ;  x = M^{-1} u  (undo right preconditioning)
      u = 0.0d0
      do i = 1, kk
         u = u + yy(i)*V(:,i)
      enddo
      call dgbtrs('N', neq, kl_jac, ku_jac, 1, abf, 2*kl_jac+ku_jac+1,  &
                  ipiv, u, neq, lpinfo)
      x = u

      ! A_z x = A_z M^-1 (V_kk yy) = V_{kk+1} Hbar(1:kk+1,1:kk) yy.
      if (present(Ax)) then
         do j = 1, kk
            do i = 1, min(j+1, kv_set)
               if (Hb(i,j) .ne. 0.0d0) Ax = Ax + (Hb(i,j)*yy(j))*V(:,i)
            enddo
         enddo
      endif

      deallocate(V,Hs,Hb,gg,cs,sn,yy,z,w,u)
      end subroutine pgmres

      ! ------------------------------------------------------!

      ! DOGLEG GEOMETRY.
      !
      ! The classical Powell dogleg: the path runs from the origin to the
      ! Cauchy point sU (the minimizer of the linear model along the steepest
      ! descent direction) and from there to the Newton point sN, and the step
      ! is where that path first meets the ball of radius delta. Returned as
      ! the two coefficients of s = aU*sU + aN*sN, because the trust region
      ! needs A*s and gets it from A*sU and A*sN by linearity rather than from
      ! another matrix-vector product.
      !
      ! Three cases, in the order they are tested:
      !   ||sN|| <= delta            -> the Newton point itself
      !   ||sU|| >= delta            -> the Cauchy direction cut to the ball
      !   otherwise                  -> the crossing of the second leg
      subroutine dogleg_step(sU, sN, delta, aU, aN, on_boundary)
      real*8, dimension(nvar_jac*N), intent(in)  :: sU, sN
      real*8,                        intent(in)  :: delta
      real*8,                        intent(out) :: aU, aN
      logical,                       intent(out) :: on_boundary
      real*8, dimension(nvar_jac*N) :: d
      real*8 :: nU, nN, a, b, c, disc, tau
      nN = sqrt(sum(sN*sN))
      nU = sqrt(sum(sU*sU))
      on_boundary = .true.
      if (nN .le. delta) then
         aU = 0.0d0;  aN = 1.0d0;  on_boundary = (nN .ge. delta)
         return
      endif
      if (nU .ge. delta) then
         ! the Cauchy leg already leaves the ball
         if (nU .gt. 0.0d0) then
            aU = delta/nU
         else
            aU = 0.0d0
         endif
         aN = 0.0d0
         return
      endif
      ! second leg: sU + tau (sN - sU), tau in [0,1], ||.|| = delta
      d = sN - sU
      a = sum(d*d)
      b = 2.0d0*sum(sU*d)
      c = sum(sU*sU) - delta*delta
      if (a .le. 0.0d0) then
         aU = 1.0d0;  aN = 0.0d0;  return
      endif
      disc = b*b - 4.0d0*a*c
      if (disc .lt. 0.0d0) disc = 0.0d0
      tau = (-b + sqrt(disc))/(2.0d0*a)
      tau = min(1.0d0, max(0.0d0, tau))
      aU = 1.0d0 - tau
      aN = tau
      end subroutine dogleg_step

      ! ------------------------------------------------------!

      ! A SCALED TRUST-REGION STEP FOR THE FOUR-UNKNOWN COUPLED SOLVE.
      !
      ! WHY. The pseudo-transient line search controls the step by a
      ! pseudo-time and then asks the merit whether the result is acceptable.
      ! Measured (docs/Update_EXHALE.md section 142), on the coupled hot
      ! Uranus that control fails in a specific way: the step it builds is an
      ! ASCENT direction of the merit in the model's own arithmetic from outer
      ! iteration 5 onward, because the I/dtau term's contribution to the
      ! merit slope carries no fixed sign under the two-sided scaling; and
      ! raising dtau until the slope turns negative does not help either,
      ! because the arc over which the raised direction descends is far
      ! shorter than the step. Both failures are failures of STEP LENGTH
      ! control, which is what a trust region is.
      !
      ! WHAT. The model is the same scaled linear model the Krylov cycle is
      ! built from,
      !
      !     m(s) = 1/2 || r0 + A s ||^2 ,   r0 = F/Drow ,   A = Drow^-1 J D ,
      !
      ! minimized over ||s|| <= delta by a dogleg between the Cauchy point and
      ! the Krylov point. A*sN comes from the Arnoldi relation for free; A*sU
      ! costs one matrix-free product; every other A*s on the dogleg path is a
      ! linear combination of those two, so the predicted reduction is exact
      ! for the model and needs no further evaluation.
      !
      ! The gradient A^T r0 is taken from the BANDED A, because a transpose
      ! is not available matrix-free. That makes the Cauchy direction an
      ! approximation while the Newton direction and every predicted reduction
      ! are not; the ray test below is what stops a bad Cauchy direction from
      ! being accepted on a model it does not describe.
      !
      ! ADMISSIBILITY IS PART OF THE STEP, not a test after it: the step is
      ! scaled back by theta until the trial state has positive rho and p, a
      ! non-negative carrier density, a carrier density inside the hydrogen
      ! its cell has, and a chemistry the equilibrium sweep can certify --
      ! which is exactly what eval_residual's "admissible" reports. The
      ! predicted reduction is recomputed for the scaled-back step, so the
      ! ratio compares like with like.
      !
      ! THE MODEL IS REJECTED, not merely disbelieved, when its directional
      ! derivative disagrees with a central finite difference of the true
      ! merit along the same ray. That is the check the ray scans of section
      ! 142 showed is needed: at the stagnating iterates the merit turns over
      ! inside a step, and a model that does not see the turn cannot be used
      ! to size one.
      subroutine trust_region_step(Y, F, f_sp, D, Drow, ab, abf, ipiv,     &
                                   idtau, gm_m, gm_rtol, f2, delta,        &
                                   dY, f2_new, pred, actual, ratio, snorm, &
                                   n_theta_cuts, gmit, gm_truncated,       &
                                   model_ok, accepted, Yacc, Facc, f_sp_acc, &
                                   f2_ref_out)
      real*8, dimension(nvar_jac*N),          intent(in)    :: Y, F, D, Drow
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)    :: f_sp
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(in)    :: ab
      real*8, dimension(2*kl_jac+ku_jac+1,nvar_jac*N), intent(inout) :: abf
      integer, dimension(nvar_jac*N),         intent(inout) :: ipiv
      real*8,  intent(in)    :: idtau, gm_rtol, f2
      integer, intent(in)    :: gm_m
      real*8,  intent(inout) :: delta
      real*8, dimension(nvar_jac*N), intent(out) :: dY
      ! The accepted trial, handed back whole: the caller adopts Y, F and the
      ! species array together, and the LAST residual evaluation this routine
      ! makes is of exactly this state, so the frozen background the caller
      ! keeps describes the state it adopts.
      real*8, dimension(nvar_jac*N),          intent(out) :: Yacc, Facc
      real*8, dimension(1-Ng:N+Ng,n_species), intent(out) :: f_sp_acc
      real*8,  intent(out) :: f2_new, pred, actual, ratio, snorm
      ! The merit of the iterate measured in the TRIAL's mode. Reported beside
      ! the caller's own f2 so that the gap between the two evaluation modes
      ! is visible rather than hidden inside the ratio.
      real*8,  intent(out) :: f2_ref_out
      integer, intent(out) :: n_theta_cuts, gmit
      logical, intent(out) :: gm_truncated, model_ok, accepted

      real*8, dimension(nvar_jac*N) :: r0, g, Ag, sN, AsN, sU, AsU, s, As
      real*8, dimension(nvar_jac*N) :: Ytry, Ftry, Jv
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng)   :: heat, cool
      real*8, dimension(3,1-Ng:N+Ng) :: utry, Wtry
      real*8  :: aU, aN, nAg, tau_c, theta, f2t, mslope, fslope
      real*8  :: f2p, f2m, hray, gnorm, f2_ref
      integer :: it, saved_mode
      logical :: okres, on_boundary, ok_state
      ! Acceptance threshold on the reduction ratio, and the radius bounds.
      ! eta = 0.1 is the standard "accept anything that reduces the merit by a
      ! tenth of what the model promised"; the radius floor exists so that a
      ! solve that cannot move says so instead of grinding the radius to zero.
      real*8,  parameter :: eta_accept  = 1.0d-1
      real*8,  parameter :: delta_floor = 1.0d-10
      integer, parameter :: n_theta_max = 12
      ! Ray test: the fraction of the step the finite difference is taken
      ! over, and how far the two slopes may differ before the model is
      ! refused. 1e-3 of the step is inside the arc the section 142 ray scans
      ! measured as model-valid; a factor 4 is loose enough that ordinary
      ! finite-difference noise does not reject a good model.
      real*8,  parameter :: ray_h      = 1.0d-3
      real*8,  parameter :: ray_factor = 4.0d0
      character(len=8) :: tr_trace_env
      logical :: tr_trace
      real*8  :: idtau_leg

      call get_environment_variable('EXHALE_TR_TRACE', tr_trace_env)
      tr_trace = (trim(tr_trace_env) .eq. '1')
      ! The dogleg's Newton leg is solved WITHOUT the pseudo-transient shift
      ! (A sN = -r0, not (idtau I + A) sN = -r0), because the dogleg's
      ! monotone-model premise requires the leg to be the minimizer of the
      ! SAME model the predicted reduction is measured in.  With the shifted
      ! leg the model residual at full length is exactly the shift,
      ! idtau ||(D/Drow) sN||, which sat above ||r0|| and made the model
      ! predict a rise on half of all iterations, so the radius walked
      ! halve-halve-double downward and could never carry a long step
      ! (verified by trace and by a Krylov-tolerance ladder that changed
      ! nothing, section 159.3).  The banded preconditioner keeps its shift:
      ! a preconditioner need not be exact.  EXHALE_TR_SHIFTED_LEG=1 restores
      ! the shifted leg for measurement.
      call get_environment_variable('EXHALE_TR_SHIFTED_LEG', tr_trace_env)
      idtau_leg = 0.0d0
      if (trim(tr_trace_env) .eq. '1') idtau_leg = idtau
      dY = 0.0d0;  f2_new = f2;  accepted = .false.;  model_ok = .true.
      f2_ref = f2
      Yacc = Y;  Facc = F;  f_sp_acc = f_sp
      pred = 0.0d0;  actual = 0.0d0;  ratio = 0.0d0;  snorm = 0.0d0
      n_theta_cuts = 0;  gmit = 0;  gm_truncated = .false.
      saved_mode = weno_mode

      ! THE REFERENCE MERIT IS MEASURED IN THE SAME MODE AS THE TRIALS.
      ! The caller's f2 is the merit of this iterate with the WENO3 weights
      ! and the limited slope FROZEN at it (mode 1); every trial is measured
      ! with them recomputed (mode 0), because that is the residual the solve
      ! is driving to zero (section 126). Those are different functions, and a
      ! reduction ratio taken across them is not a ratio of anything: measured
      ! on the coupled hot Uranus the two differ by enough that a step which
      ! lowers the mode-0 merit is followed by a mode-1 re-evaluation that
      ! raises it, so the region grows on a reduction the next iteration does
      ! not see. One extra evaluation at the iterate fixes it, and both the
      ! predicted and the actual reduction are then taken about the same
      ! baseline. The MODEL is still the Krylov one; only its constant term is
      ! the residual actually being reduced.
      weno_mode = 0
      call eval_residual(Y, f_sp, fwork, Ftry, heat, cool)
      r0 = Ftry/Drow
      f2_ref = sqrt(sum(r0*r0))

      ! --- the Krylov (Newton) point, and A applied to it, for free ---
      call pgmres(Y, F, f_sp, D, Drow, abf, ipiv, idtau_leg, -r0, sN,    &
                  gm_m, gm_rtol, gmit, gm_truncated, AsN)
      if (gmit .le. 0) then
         model_ok = .false.
         f2_ref_out = f2_ref
         weno_mode = saved_mode
         return
      endif
      ! A_z = idtau*(D/Drow) + A, so remove the pseudo-transient part to get
      ! the Jacobian action the model is written in.
      AsN = AsN - idtau_leg*(D/Drow)*sN
      ! The model residual at the FULL Newton leg is, for an exact leg,
      ! r0 + A sN = -idtau (D/Drow) sN: the pseudo-transient shift, not zero.
      ! Printed against ||r0|| because their ratio decides the sign of the
      ! full-length predicted reduction.
      if (tr_trace) write(*,'(A,ES11.3,A,ES11.3,A,ES11.3)')               &
           ' (JFNK) [TR-trace] ||r0||=', sqrt(sum(r0*r0)),                &
           '  idtau*||(D/Drow)sN||=', idtau*sqrt(sum(((D/Drow)*sN)**2)),  &
           '  ||r0+AsN||=', sqrt(sum((r0+AsN)*(r0+AsN)))

      ! --- the Cauchy point ---
      call band_matvec_transpose(ab, r0, g)
      gnorm = sqrt(sum(g*g))
      if (gnorm .le. 0.0d0) then
         sU = 0.0d0;  AsU = 0.0d0
      else
         call jv_product(Y, F, f_sp, D*g, Jv, ok_state)
         if (.not. ok_state) then
            model_ok = .false.
            f2_ref_out = f2_ref
            weno_mode = saved_mode
            return
         endif
         Ag  = Jv/Drow
         nAg = sum(Ag*Ag)
         if (nAg .le. 0.0d0) then
            sU = 0.0d0;  AsU = 0.0d0
         else
            tau_c = (gnorm*gnorm)/nAg
            sU  = -tau_c*g
            AsU = -tau_c*Ag
         endif
      endif

      ! --- the dogleg, then admissibility, then the ratio ---
      do it = 1, n_theta_max
         call dogleg_step(sU, sN, delta, aU, aN, on_boundary)
         if (tr_trace) write(*,'(A,I3,A,ES11.3,A,L1)')                    &
              ' (JFNK) [TR-trace] trial', it, '  delta_in=', delta,        &
              '  on_boundary=', on_boundary
         s  = aU*sU  + aN*sN
         As = aU*AsU + aN*AsN
         snorm = sqrt(sum(s*s))
         if (snorm .le. 0.0d0) then
            model_ok = .false.;  exit
         endif
         ! Scale back until the trial is a state the code can describe. The
         ! model is linear, so the predicted reduction of the scaled-back step
         ! is recomputed rather than scaled.
         theta = 1.0d0
         ok_state = .false.
         do while (theta .ge. 1.0d0/2.0d0**n_theta_max)
            Ytry = Y + D*(theta*s)
            call unpack_U(Ytry, utry)
            call U_to_W_interior(utry, Wtry)
            if (minval(Wtry(1,1:N)) .gt. 0.0d0 .and.                     &
                minval(Wtry(3,1:N)) .gt. 0.0d0) then
               weno_mode = 0
               call eval_residual(Ytry, f_sp, fwork, Ftry, heat, cool,    &
                                  admissible=okres)
               if (okres) then
                  f2t = sqrt(sum((Ftry/Drow)**2))
                  ok_state = .true.
                  exit
               endif
            endif
            theta = 0.5d0*theta
            n_theta_cuts = n_theta_cuts + 1
         enddo
         if (.not. ok_state) then
            ! nothing on this dogleg is describable; shrink the region
            delta = 0.25d0*delta
            if (delta .lt. delta_floor) then
               model_ok = .false.;  exit
            endif
            cycle
         endif
         s     = theta*s
         As    = theta*As
         snorm = theta*snorm
         pred   = 0.5d0*(sum(r0*r0) - sum((r0 + As)*(r0 + As)))
         actual = 0.5d0*(f2_ref*f2_ref - f2t*f2t)
         if (pred .le. 0.0d0) then
            ! a model that promises no reduction is not a model to step on
            model_ok = .false.
            if (tr_trace) write(*,'(A,I3,A,ES11.3,A,ES11.3)')             &
                 ' (JFNK) [TR-trace] trial', it, '  pred<=0  pred=', pred, &
                 '  snorm=', snorm
            delta = 0.25d0*snorm
            if (delta .lt. delta_floor) exit
            cycle
         endif
         ratio = actual/pred

         ! --- the ray test: does the model's slope describe the function? ---
         mslope = sum(r0*As)
         hray   = ray_h
         weno_mode = 0
         Ytry = Y + D*(hray*s)
         call eval_residual(Ytry, f_sp, fwork, Ftry, heat, cool,          &
                            admissible=okres)
         f2p = 0.5d0*sum((Ftry/Drow)**2)
         Ytry = Y - D*(hray*s)
         call eval_residual(Ytry, f_sp, fwork, Ftry, heat, cool,          &
                            admissible=okres)
         f2m = 0.5d0*sum((Ftry/Drow)**2)
         fslope = (f2p - f2m)/(2.0d0*hray)
         model_ok = (mslope*fslope .gt. 0.0d0) .and.                      &
                    (abs(fslope) .le. ray_factor*abs(mslope)) .and.       &
                    (abs(mslope) .le. ray_factor*abs(fslope))
         if (.not. model_ok) then
            delta = 0.25d0*snorm
            exit
         endif

         ! --- radius update and acceptance ---
         if (ratio .lt. 0.25d0) then
            delta = 0.25d0*snorm
         else if (ratio .gt. 0.75d0 .and. on_boundary) then
            delta = 2.0d0*delta
         endif
         if (ratio .gt. eta_accept) then
            ! Re-evaluate the accepted trial so that it is the LAST state this
            ! routine evaluated: the caller's keep_background_of_adopted_state,
            ! carrier_relnorm_last and the acceptance ledgers all read what the
            ! previous evaluation left behind, and the ray test above has been
            ! evaluating points that are not the step.
            weno_mode = 0
            Yacc = Y + D*s
            call eval_residual(Yacc, f_sp, f_sp_acc, Facc, heat, cool,    &
                               admissible=okres)
            if (okres) then
               f2t = sqrt(sum((Facc/Drow)**2))
               dY = D*s;  f2_new = f2t;  accepted = .true.
            else
               ! it was admissible a moment ago and is not now; refuse it
               model_ok = .false.
               delta = 0.25d0*snorm
            endif
         endif
         exit
      enddo
      f2_ref_out = f2_ref
      weno_mode = saved_mode
      end subroutine trust_region_step

      ! ------------------------------------------------------!

      subroutine solve_steady_jfnk(u, f_sp, resid_tol, maxit, dtau0,    &
                                   gm_m, info)
      ! Pseudo-transient-continuation Jacobian-free Newton-Krylov solve.
      ! Per outer iteration:
      !   F  = eval_residual(Y)                          (full residual)
      !   M  = factored banded (I/dtau + J_banded)       (preconditioner)
      !   solve (I/dtau + J) dY = -F  by right-precond. GMRES (J*v matrix-
      !                                free; captures non-local radiation)
      !   Y <- Y + lam*dY   (||F||_2 line search + positivity)
      !   dtau <- SER ramp
      real*8, dimension(3,1-Ng:N+Ng),         intent(inout) :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8,  intent(in)  :: resid_tol, dtau0
      integer, intent(in)  :: maxit, gm_m
      integer, intent(out) :: info

      integer :: neq, ldab, iter, ls, lpinfo, jc, gmit
      integer :: jj, kk, jworst, kworst, irow, ilo, ihi
      real*8  :: amx
      real*8, allocatable :: Y(:), F(:), Ftry(:), dY(:), Ytry(:)
      real*8, allocatable :: dZ(:), D(:), dY_lm(:), Drow(:)
      real*8, allocatable :: ab(:,:), abf(:,:)
      integer, allocatable :: ipiv(:)
      real*8, dimension(1-Ng:N+Ng)            :: heat0, cool0
      real*8, dimension(3,1-Ng:N+Ng)          :: utry, Wtry
      real*8, dimension(1-Ng:N+Ng,n_species)  :: f_sp_j, f_sp_best
      ! Receives each hand-back sweep's answer: the seed and the result of
      ! eval_residual may not be the same array (section 147), so the fixed
      ! point of section 155 is iterated through this one explicitly.
      real*8, dimension(1-Ng:N+Ng,n_species)  :: f_sc_out
      real*8  :: rnorm, rc(3), dtau, lam, f2, f2_try, idtau, amx0
      real*8  :: f2hist(5), f2ref, rnorm_best
      ! Flux gate: the spread of rho v r^2 over r >= r_flux, of the current
      ! iterate and of the best one kept, and whether each meets BOTH gates.
      real*8  :: fspread, fspread_best
      real*8  :: fspread_103, fspread_110, dum_fmean
      logical :: gates_now, gates_best
      real*8  :: rn_prev
      integer :: it_sc
      character(len=32) :: env_sc
      integer :: n_no_descent, n_nonmonotone_accepts
      integer :: n_since_best
      logical :: monotone_search
      ! Stagnation limit: consecutive outer iterations in which the line
      ! search found NO acceptable step at all. Measured on the HD 209458 b
      ! hand-off state (docs/newton_scaling_and_base_wall.md): runs that go on
      ! to converge never string more than 4 such iterations together, runs
      ! that are truly stuck string 34 or more.
      integer, parameter :: n_no_descent_max = 12
      ! P51. STAGNATION OF THE BEST ITERATE, and what to do about it.
      !
      ! The non-monotone (Grippo) line search compares a trial against the
      ! WORST merit of the last five iterates, so it accepts sideways steps.
      ! That is what lets the solve leave a bad neighbourhood, and it is also
      ! what lets it circle: measured on the coupled hot Uranus hand-off, the
      ! best iterate is reached at outer iteration 10 and then 490 further
      ! iterations produce no improvement of ||R|| at all while the search
      ! accepts 61 sideways steps.
      !
      ! So when the best iterate has not improved for n_stall_best
      ! iterations the solve GOES BACK TO IT and switches the acceptance to
      ! monotone Armijo -- descent against the current merit rather than
      ! against the worst of five. If that still buys nothing in another
      ! n_stall_best iterations, the solve stops and says it STAGNATED,
      ! which is a different statement from "no descent direction exists".
      !
      ! WHERE 20 COMES FROM. Over every solve of the P50 campaign that ended
      ! info = 0 on a configuration that survived measurement, the largest
      ! gap between two improvements of the best iterate is 6 outer
      ! iterations (both WASP-121b Newton cases); on the coupled hot Uranus
      ! it is 3. Two accepted solves show gaps of 68 and 85, and both are
      ! configurations that were measured and rejected. So 20 stands a
      ! factor 3.3 above anything a surviving route has needed, and 40 --
      ! the point at which the solve gives up -- stands a factor 6.7 above
      ! it and below every stall observed (tails of 278 to 729).
      integer, parameter :: n_stall_best = 20
      logical :: ok, try_ok, gm_truncated, carrier_ok
      ! --- D2: the scaled trust region on the four-unknown branch ---
      ! Enabled by EXHALE_TRUST_REGION=1. When off, not one statement of the
      ! path below runs and the three-unknown route never reaches it at all.
      logical :: use_tr
      real*8  :: tr_delta, tr_pred, tr_actual, tr_ratio, tr_snorm, tr_dmax
      ! Relative tolerance of the Krylov solve that builds the dogleg's
      ! Newton leg.  The dogleg falls monotonically along its path only
      ! when that leg is the model's own minimizer, so a loose leg makes
      ! the model predict an increase at full step length (section 159.3).
      ! Overridable so the dependence can be measured; the default is the
      ! value the section 153 draft used.
      real*8  :: tr_gm_rtol
      real*8  :: tr_f2ref, tr_f2_end, tr_f2_oldscale
      integer :: tr_cuts, n_tr_accept, n_tr_reject, n_tr_model_bad
      logical :: tr_model_ok, tr_accept
      character(len=24) :: tr_env
      ! Damped Gauss-Newton escape from an ascending Newton direction.
      logical :: lm_tried, lm_found
      real*8  :: f2_lm, mu_lm, gnorm_lm
      ! Probe states of THIS solve that could not be described: Jacobian
      ! colors left at zero, and GMRES cycles cut short because the next
      ! Krylov direction had no admissible sample (section 121).
      integer :: n_jac_unresolved, n_jac_unresolved_tot, n_gm_truncated
      real*8, allocatable :: Ybest(:)

      neq  = nvar_jac*N
      ldab = 2*kl_jac + ku_jac + 1
      allocate(Y(neq), F(neq), Ftry(neq), dY(neq), Ytry(neq))
      allocate(dZ(neq), D(neq), Ybest(neq), dY_lm(neq))
      allocate(ab(ldab,neq), abf(ldab,neq), ipiv(neq))
      ! Allocated last and on its own, so that every array the three-unknown
      ! solve has always had keeps the address it has always had.
      allocate(Drow(neq))
      dtau = dtau0;  info = 1

      n_jac_unresolved_tot = 0;  n_gm_truncated = 0
      use_tr = .false.
      call get_environment_variable('EXHALE_TRUST_REGION', tr_env)
      if (trim(tr_env) .eq. '1' .and. nvar_jac .ge. 4) use_tr = .true.
      n_tr_accept = 0;  n_tr_reject = 0;  n_tr_model_bad = 0
      tr_f2_end = -1.0d0
      tr_delta = -1.0d0
      tr_dmax  = -1.0d0
      tr_gm_rtol = 1.0d-1
      call get_environment_variable('EXHALE_TR_GMRTOL', tr_env)
      if (len_trim(tr_env) .gt. 0) read(tr_env,*) tr_gm_rtol
      call get_environment_variable('EXHALE_TR_DELTA0', tr_env)
      if (len_trim(tr_env) .gt. 0) read(tr_env,*) tr_delta
      call reset_uncertified_trial_count

      call pack_U(u, Y)
      call pack_carrier(u, f_sp, Y)
      n_trial_headroom = 0;  n_resid_eval = 0
      ! From here to the end of the solve the equilibrium sweeps are the
      ! steady solver's, not the run's. The current iterate is a state the
      ! run holds; everything the solver evaluates to BUILD its Newton model
      ! -- Jacobian columns, Krylov products, line-search trials -- is not,
      ! and is tagged separately so that the acceptance ledgers, the non-root
      ! streak and its stop are not written by states no one adopted
      ! (docs/Update_EXHALE.md section 121).
      call set_ioniz_eq_sweep_state_kind(ieq_state_steady_iterate)
      ! Seed from the adopted composition, write the swept one into the work
      ! array, then ADOPT it. The copy is an adoption and not a seeding
      ! convention: the seed is named in the call.
      call eval_residual(Y, f_sp, f_sp_j, F, heat0, cool0)
      f_sp = f_sp_j
      n_uncertified_iterate = n_uncertified_last
      call keep_background_of_adopted_state
      call resid_relnorm(F, u, rc, rnorm)
      call cell_state_scales(Y, D)
      if (nvar_jac .ge. 4) call cell_row_scales(Y, D, Drow)
      ! Line-search merit: the scaled 2-norm over the WHOLE domain 1..N. It has
      ! covered every physical cell since it was written, and it stays that
      ! way: restricting it to the wind -- a volume-weighted RMS of F/D over
      ! [j_min:N] -- was tried on 2026-08-11 and again on 2026-08-13, and both
      ! times it broke solves that converge with the whole-domain norm:
      ! photo_deep_secion_cont goes info = 0 at ||R|| = 7.28e-4 in 4 outer
      ! iterations against info = 2 at 3.25e-3 in 15, and examples/15
      ! (HD 209458 b, molecular) goes info = 0 at 9.54e-6 in 272 against
      ! info = 2 at 1.17e-4 in 30. The Newton step is computed from the full
      ! residual over ALL rows, so a merit that ignores most of those rows
      ! rejects the steps that step actually takes.
      !
      ! The merit and rnorm now cover the same cells, but they remain
      ! DIFFERENT FUNCTIONALS: a 2-norm of F/D, D being the cell's own signal
      ! scale rho(|v|+c_s) in the momentum row (cell_state_scales), against a
      ! maximum over rows and over the two regions of volume-weighted ratios
      ! taken on the scale each region's physics sets (residual_row_scale).
      ! So the iterate with the smallest ||R|| still need not be the last one
      ! -- handled at the end of this routine, not by changing the merit
      ! (docs/newton_scaling_and_base_wall.md section 4,
      ! docs/hd209_metal_stagnation.md).
      if (nvar_jac .ge. 4) then
         f2 = sqrt(sum((F/Drow)**2))  ! merit in the SCALED space
      else
         f2 = sqrt(sum((F/D)**2))     ! merit in the SCALED space
      endif
      f2hist = f2                  ! non-monotone line-search memory
      rnorm_best = rnorm;  Ybest = Y;  f_sp_best = f_sp;  n_no_descent = 0
      n_since_best = 0;  monotone_search = .false.
      carrier_relnorm_best = carrier_relnorm_last
      carrier_cellmax_best = carrier_cellmax_last
      carrier_cell_best    = carrier_cell_worst
      gates_best = steady_gates_met(rnorm, u, resid_tol, fspread_best,   &
                                    nvar_jac .ge. 4, carrier_relnorm_last)
      gates_now  = gates_best;  fspread = fspread_best
      call keep_background_of_best_iterate
      n_nonmonotone_accepts = 0
      write(*,'(A,ES11.3,A,ES10.2)') ' (JFNK) start ||R||=',rnorm,        &
           '  ||Fs||2=',f2
      call write_resid_below_escape(' (JFNK)', F, u)

      do iter = 1, maxit
         if (steady_gates_met(rnorm, u, resid_tol, fspread,        &
                          nvar_jac .ge. 4, carrier_relnorm_last)) then
            info = 0;  exit
         endif
         idtau = 1.0d0/dtau

         ! Diagonal scaling for this outer iteration.
         call cell_state_scales(Y, D)
         if (nvar_jac .ge. 4) call cell_row_scales(Y, D, Drow)

         ! FREEZE the WENO3 weights at the current iterate: one mode-1
         ! residual evaluation stores the smoothness factors and yields the
         ! frozen-weights F; the inner evaluations that build the Newton model
         ! (Jacobian probes, J*v) then reuse them (mode 2), so the inner
         ! problem excludes the strongly nonlinear weight response (the
         ! standard lagged-weights remedy for FV steady solves). The line
         ! search does NOT use them -- see there.
         weno_mode = 1
         call set_ioniz_eq_sweep_state_kind(ieq_state_steady_iterate)
         call eval_residual(Y, f_sp, f_sp_j, F, heat0, cool0)
         ! The reference the trial admissibility test compares against: how
         ! many cells THIS iterate's own chemistry left uncertified. Refreshed
         ! here, at the top of every outer iteration, from the iterate itself.
         n_uncertified_iterate = n_uncertified_last
         call keep_background_of_adopted_state
         call set_ioniz_eq_sweep_state_kind(ieq_state_steady_candidate)
         weno_mode = 2
         if (nvar_jac .ge. 4) then
            f2 = sqrt(sum((F/Drow)**2)) ! TRUE merit of the current iterate
         else
            f2 = sqrt(sum((F/D)**2))    ! TRUE merit of the current iterate
         endif

         ! Grippo non-monotone reference: the worst true merit of the last 5
         ! outer iterates, the current one included. Recording it here, and not
         ! only when a step is accepted, is what keeps the reference from
         ! lagging behind the state the line search actually starts from.
         f2hist = (/ f2hist(2:5), f2 /)

         ! Banded preconditioner of the SCALED system:
         ! M = I/dtau + D^-1 J_banded D, factored.
         call build_banded_jac_full(Y, f_sp, ab, n_jac_unresolved)
         n_jac_unresolved_tot = n_jac_unresolved_tot + n_jac_unresolved
         ! The pseudo-transient term is the identity of the UNSCALED
         ! system, so under a two-sided scaling it is diag(D/Drow), not the
         ! identity -- the same wherever the two scalings are, which is
         ! every row but the carrier's. The three-unknown branch is kept as
         ! the original text so that it stays byte-identical.
         if (nvar_jac .ge. 4) then
            do jc = 1, neq
               ilo = max(1,   jc - ku_jac)
               ihi = min(neq, jc + kl_jac)
               do irow = ilo, ihi
                  ab(kl_jac+ku_jac+1 + irow - jc, jc) =                  &
                       ab(kl_jac+ku_jac+1 + irow - jc, jc)*D(jc)/Drow(irow)
               enddo
            enddo
            abf = ab
            do jc = 1, neq
               abf(kl_jac+ku_jac+1, jc) = abf(kl_jac+ku_jac+1, jc)       &
                                        + idtau*D(jc)/Drow(jc)
            enddo
         else
            do jc = 1, neq
               ilo = max(1,   jc - ku_jac)
               ihi = min(neq, jc + kl_jac)
               do irow = ilo, ihi
                  ab(kl_jac+ku_jac+1 + irow - jc, jc) =                  &
                       ab(kl_jac+ku_jac+1 + irow - jc, jc)*D(jc)/D(irow)
               enddo
            enddo
            abf = ab
            do jc = 1, neq
               abf(kl_jac+ku_jac+1, jc) = abf(kl_jac+ku_jac+1, jc) + idtau
            enddo
         endif
         call dgbtrf(neq, neq, kl_jac, ku_jac, abf, ldab, ipiv, lpinfo)
         if (lpinfo .ne. 0) then
            dtau = max(dtau*0.25d0, dtau0);  cycle
         endif

         ! ---- D2: the scaled trust region, four-unknown branch only ----
         if (use_tr) then
            ! The initial radius is the MEASURED model-valid arc, not a raw
            ! norm: the ray scans of section 142 put the whole of the merit's
            ! descent inside t = 1e-2 of the Newton step and the rise dominant
            ! by t = 0.1, so the region starts at a hundredth of the first
            ! Newton step's length and is then controlled by the ratio.
            if (tr_delta .lt. 0.0d0) then
               call pgmres(Y, F, f_sp, D, Drow, abf, ipiv, idtau, -F/Drow, &
                           dZ, gm_m, 1.0d-1, gmit, gm_truncated)
               tr_delta = 1.0d-2*sqrt(sum(dZ*dZ))
               tr_dmax  = 1.0d2*tr_delta
               write(*,'(A,ES11.3,A,ES11.3)') ' (JFNK) [TR] initial radius'//&
                    ' from the model-valid arc:', tr_delta,                 &
                    '  (1e-2 of the first Newton step,', sqrt(sum(dZ*dZ))
               write(*,'(A,ES9.2)') ' (JFNK) [TR] Krylov tolerance of the'//&
                    ' dogleg Newton leg:', tr_gm_rtol
            endif
            ! IS THE MERIT THE SAME FUNCTION IT WAS LAST ITERATION? Drow is
            ! rebuilt from the state at the top of every outer iteration, so
            ! ||F/Drow|| is a different functional each time. A trust region
            ! compares a predicted with an actual reduction of ONE function;
            ! if the function itself moves between iterations the comparison
            ! is only valid inside an iteration. This prints the size of that
            ! move: the merit of the state just adopted, measured with the
            ! scaling it was adopted under, against the same state measured
            ! with the scaling of the iteration that follows.
            if (tr_f2_end .gt. 0.0d0)                                      &
               write(*,'(A,ES11.3,A,ES11.3,A,F9.3)')                       &
                    ' (JFNK) [TR]      merit across the iteration boundary:'//&
                    ' old scaling', tr_f2_end, '   new scaling', f2,        &
                    '   factor', f2/max(tr_f2_end,1.0d-99)
            call trust_region_step(Y, F, f_sp, D, Drow, ab, abf, ipiv,     &
                                   idtau, gm_m, tr_gm_rtol, f2, tr_delta,  &
                                   dY, f2_try, tr_pred, tr_actual,         &
                                   tr_ratio, tr_snorm, tr_cuts, gmit,      &
                                   gm_truncated, tr_model_ok, tr_accept,   &
                                   Ytry, Ftry, f_sp_j, tr_f2ref)
            if (tr_dmax .gt. 0.0d0) tr_delta = min(tr_delta, tr_dmax)
            if (gm_truncated) n_gm_truncated = n_gm_truncated + 1
            ok  = tr_accept
            lam = 1.0d0
            tr_f2_end = -1.0d0
            if (tr_accept) then
               tr_f2_end   = f2_try
               n_tr_accept = n_tr_accept + 1
            else
               n_tr_reject = n_tr_reject + 1
            endif
            if (.not. tr_model_ok) n_tr_model_bad = n_tr_model_bad + 1
            ! One line per iteration, accepted or not: radius, predicted and
            ! actual reduction, their ratio, and the step actually taken.
            write(*,'(A,I4,A,ES10.3,A,ES11.3,A,ES11.3,A,F9.3,A,ES10.3,'//  &
                    'A,I2,A,L1,A,L1)')                                     &
                 ' (JFNK) [TR] it', iter, '  delta=', tr_delta,            &
                 '  pred=', tr_pred, '  actual=', tr_actual,               &
                 '  ratio=', tr_ratio, '  ||s||=', tr_snorm,               &
                 '  cuts=', tr_cuts, '  model_ok=', tr_model_ok,           &
                 '  accepted=', tr_accept
            write(*,'(A,ES11.3,A,ES11.3,A,F8.3)')                          &
                 ' (JFNK) [TR]      merit of the iterate: frozen-mode',    &
                 f2, '   trial-mode', tr_f2ref, '   ratio', f2/max(tr_f2ref,1.0d-99)
         else
         ! GMRES solve of the scaled PTC system; unscale the step.
         if (nvar_jac .ge. 4) then
            call pgmres(Y, F, f_sp, D, Drow, abf, ipiv, idtau, -F/Drow,  &
                        dZ, gm_m, 1.0d-1, gmit, gm_truncated)
         else
         call pgmres(Y, F, f_sp, D, Drow, abf, ipiv, idtau, -F/D, dZ,     &
                     gm_m, 1.0d-1, gmit, gm_truncated)
         endif
         if (gm_truncated) n_gm_truncated = n_gm_truncated + 1
         dY = D*dZ
         endif

         ! Scaled ||F/D||_2 backtracking line search with positivity, tested
         ! against the non-monotone reference above.
         !
         ! The trial states are evaluated with weno_mode = 0, i.e. with the
         ! smoothness weights recomputed at the trial state, so what decides
         ! acceptance is the residual the solve is driving to zero. Evaluating
         ! the trials with the weights frozen at Y instead -- the quantity the
         ! inner Newton model minimizes -- is a different measure: on the
         ! HD 189733 b hand-off state the two differ by a median factor 2.2 at
         ! lam = 1, and 7 of the 14 steps the frozen test accepted RAISED the
         ! true residual, by up to 9x (docs/newton_scaling_and_base_wall.md
         ! §10). Mode 0 leaves the stored weights alone, so the frozen model of
         ! this outer iteration survives the search.
         ! Armijo after a stagnation restart: descent against the merit of
         ! the state the search starts from, not against the worst of the
         ! last five. Until the restart happens this is the original line.
         if (monotone_search) then
            f2ref = f2
         else
            f2ref = maxval(f2hist)
         endif
         if (.not. use_tr) then
         lam = 1.0d0;  ok = .false.;  f2_try = huge(1.0d0)
         weno_mode = 0
         ! A cycle that could sample NO Krylov direction returns the zero
         ! step. The zero step must not go through the line search: the
         ! non-monotone test compares the trial against the WORST merit of
         ! the last five iterates, so Ytry = Y passes it whenever the current
         ! iterate is not that worst one, and the iteration would be recorded
         ! as a descent step that moves nothing -- resetting the no-descent
         ! counter and letting the solve spin to maxit instead of aborting.
         ! There is no Newton model this iteration; that is a no-descent
         ! iteration, and it is counted as one.
         if (gmit .le. 0) then
            lam = 0.0d0
            write(*,'(A,I4,A)') ' (JFNK) it',iter,'  no Krylov direction'// &
                 ' could be sampled; no Newton model this iteration'
         endif
         do ls = 1, 20
            if (gmit .le. 0) exit
            Ytry = Y + lam*dY
            call unpack_U(Ytry, utry)
            call U_to_W_interior(utry, Wtry)   ! ghosts of utry are not set
            ! The carrier unknown gets the same treatment rho and p get: a
            ! trial that puts a negative H2 density into a cell is rejected
            ! before the residual is asked for it.
            carrier_ok = .true.
            if (nvar_jac .ge. 4) then
               do jj = 1, N
                  if (Ytry(nvar_jac*(jj-1)+4) .lt. 0.0d0)                &
                     carrier_ok = .false.
               enddo
            endif
            if (carrier_ok .and.                                        &
                minval(Wtry(1,1:N)) .gt. 0.0d0 .and.                    &
                minval(Wtry(3,1:N)) .gt. 0.0d0) then
               call eval_residual(Ytry, f_sp, f_sp_j, Ftry, heat0, cool0, &
                                  admissible=try_ok)
               if (nvar_jac .ge. 4) then
                  f2_try = sqrt(sum((Ftry/Drow)**2))
               else
                  f2_try = sqrt(sum((Ftry/D)**2))
               endif
               ! A trial state the code cannot describe is rejected here,
               ! not left to the merit comparison: a non-finite merit fails
               ! the test by accident of IEEE ordering, and a trial whose
               ! chemistry broke in only a few cells can still produce a
               ! finite, smaller merit and be adopted (section 121).
               if (try_ok .and.                                          &
                   f2_try .lt. (1.0d0 - 1.0d-4*lam)*f2ref) then
                  ok = .true.
                  if (f2_try .ge. (1.0d0 - 1.0d-4*lam)*f2)              &
                       n_nonmonotone_accepts =                    &
                            n_nonmonotone_accepts + 1
                  exit
               endif
            endif
            lam = 0.5d0*lam
         enddo

         ! WHEN BACKTRACKING FINDS NOTHING AND SHORTENING THE STEP CANNOT
         ! HELP. The search above only shortens the Newton/PTC direction. If
         ! that direction points uphill, every backtrack points uphill too --
         ! the excess merit falls linearly with lam and never changes sign --
         ! and with dtau already at its floor the next outer iteration would
         ! rebuild the same preconditioner, the same Krylov direction and the
         ! same search, and repeat this iteration bit for bit until the
         ! stagnation counter fires (11 identical repetitions, measured on the
         ! He/H = 0.3 hand-off state). What is wrong there is the DIRECTION,
         ! so a direction that is guaranteed downhill is built instead
         ! (docs/Update_EXHALE.md section 126). While dtau is still above its
         ! floor the reduction below does change the next iteration, and the
         ! PTC ramp is left to do its work.
         endif   ! .not. use_tr -- the trust region replaces the search
         ! The damped Gauss-Newton escape is what the trust region REPLACES:
         ! it looks for a descent direction after the fact, where the region
         ! controls the step by model agreement before the fact. Running both
         ! would make the reduction ratio meaningless, so it is skipped.
         lm_tried = .false.;  lm_found = .false.
         if ((.not. use_tr) .and. dtau .le. dtau0*(1.0d0 + 1.0d-12)) then
            ! Two ways for the Newton/PTC direction to leave the iterate where
            ! it is, both of them repeatable bit for bit at this dtau: the
            ! backtracking search rejects every step, or it accepts one that
            ! does not lower the true merit -- the non-monotone test compares
            ! against the worst of the last five iterates, so a sideways step
            ! passes it and the solve can circle for hundreds of iterations
            ! (measured: 487 such acceptances in 500 iterations on the
            ! He/H = 0.3 hand-off state).
            lm_tried = (.not. ok)
            if (ok) lm_tried = (f2_try .ge. f2)
         endif
         if (lm_tried) then
            call levenberg_marquardt_descent(Y, F, f_sp, D, Drow, ab, f2,&
                            dY_lm, f2_lm, mu_lm, gnorm_lm, lm_found)
            if (lm_found .and. f2_lm .ge. f2_try) lm_found = .false.
            if (lm_found) then
               ! BACKTRACK ALONG THE DAMPED GAUSS-NEWTON DIRECTION, in the
               ! coupled route only. The escape used to evaluate the full
               ! step and nothing else, which is enough while the only way
               ! it can fail is to ascend -- the direction is built to
               ! descend, so a shorter one descends too. With the carrier
               ! unknown a second failure mode appears and it is not about
               ! descent: the full step lowers the merit (2.28 -> 1.50,
               ! measured) but puts 73 cells' molecular equilibrium outside
               ! the physical simplex, so the trial is refused and the solve
               ! aborts with a descent direction in hand. Shortening the
               ! step is what that asks for. The 3-unknown route keeps the
               ! single evaluation it has always had, so no atomic solve
               ! moves.
               if (nvar_jac .lt. 4) then
               Ytry = Y + dY_lm
               weno_mode = 0
               call eval_residual(Ytry, f_sp, f_sp_j, Ftry, heat0, cool0, &
                                  admissible=try_ok)
               f2_try = sqrt(sum((Ftry/D)**2))
               ok  = try_ok
               lam = 1.0d0
               else
                  lam = 1.0d0
                  do ls = 1, 12
                     Ytry = Y + lam*dY_lm
                     carrier_ok = .true.
                     do jj = 1, N
                        if (Ytry(nvar_jac*(jj-1)+4) .lt. 0.0d0)          &
                           carrier_ok = .false.
                     enddo
                     if (carrier_ok) then
                        weno_mode = 0
                        call eval_residual(Ytry, f_sp, f_sp_j, Ftry,      &
                                           heat0,                        &
                                           cool0, admissible=try_ok)
                        f2_try = sqrt(sum((Ftry/Drow)**2))
                     else
                        try_ok = .false.;  f2_try = huge(1.0d0)
                     endif
                     if (try_ok .and. f2_try .lt. f2) exit
                     lam = 0.5d0*lam
                  enddo
                  ok  = try_ok .and. f2_try .lt. f2
               endif
               write(*,'(A,I4,A,ES9.2,A,ES10.3,A,ES10.3)')               &
                    ' (JFNK) it',iter,'  damped Gauss-Newton escape,'//  &
                    ' mu=',mu_lm,'  ||Fs||2',f2,' ->',f2_try
            else if (.not. ok) then
               ! Nothing downhill along either direction. With the damped
               ! Gauss-Newton family exhausted this is a stationary point of
               ! the merit to within what the banded model can see, so the
               ! gradient norm is what says whether that is true.
               write(*,'(A,I4,A,ES10.3)') ' (JFNK) it',iter,'  no descent'// &
                    ' along the Newton or the damped Gauss-Newton'//     &
                    ' direction; ||grad merit||=', gnorm_lm
            endif
         endif

         if (ok) then
            Y = Ytry;  F = Ftry;  f_sp = f_sp_j
            ! The last residual evaluation was this accepted trial, so the
            ! frozen background now describes the state just adopted.
            call keep_background_of_adopted_state
            call unpack_U(Y, u);  call Apply_BC(u)
            call resid_relnorm(F, u, rc, rnorm)
            dtau = min(dtau*max(lam,0.1d0)*(f2/max(f2_try,1.0d-30)),     &
                       1.0d14*dtau0)
            f2 = f2_try
         else
            dtau = max(dtau*0.25d0, dtau0)
         endif

         ! Keep the best iterate seen (by the convergence measure) so that a
         ! failed solve returns it rather than wherever it stopped.
         ! Keep the best iterate by ||R||, but never prefer one that fails
         ! the flux gate over one that passes it: the returned state is judged
         ! on BOTH gates, so the iterate kept has to be the best ACCEPTABLE
         ! one wherever an acceptable one has been seen.
         gates_now = steady_gates_met(rnorm, u, resid_tol, fspread,        &
                          nvar_jac .ge. 4, carrier_relnorm_last)
         if ((gates_now .and. .not. gates_best) .or.                      &
             ((gates_now .eqv. gates_best) .and. rnorm .lt. rnorm_best)) then
            rnorm_best = rnorm;  Ybest = Y;  f_sp_best = f_sp
            fspread_best = fspread;  gates_best = gates_now
            call keep_background_of_best_iterate
            n_since_best = 0
            carrier_relnorm_best = carrier_relnorm_last
            carrier_cellmax_best = carrier_cellmax_last
            carrier_cell_best    = carrier_cell_worst
         else
            n_since_best = n_since_best + 1
         endif

         ! P51: the best iterate has stopped moving. Go back to it and stop
         ! accepting sideways steps; if that buys nothing either, stop.
         if (n_since_best .ge. n_stall_best) then
            if (.not. monotone_search) then
               Y = Ybest;  f_sp = f_sp_best;  rnorm = rnorm_best
               call adopt_background_of_best_iterate
               call unpack_U(Y, u)
               call unpack_carrier(Y, u, f_sp)
               call Apply_BC(u)
               monotone_search = .true.
               n_since_best    = 0
               write(*,'(A,I0,A,ES11.3,A)') ' (JFNK) best iterate has not'//&
                    ' improved in ', n_stall_best, ' iterations;'//        &
                    ' restarting from it (||R||=', rnorm_best,            &
                    ') with a monotone line search'
            else
               info = 2
               write(*,'(A,I0,A)') ' (JFNK) STAGNATED: no improvement of'//&
                    ' the best iterate in ', 2*n_stall_best,              &
                    ' iterations, monotone search included -- stopping'
               exit
            endif
         endif

         ! Stagnation: the failure mode of the non-monotone search is that it
         ! stops finding ANY acceptable step and the iterate no longer moves.
         ! Count consecutive such iterations. (A watchdog on rnorm instead
         ! cannot work here: the solver minimizes the merit, not rnorm, and
         ! the opening pseudo-transient legitimately raises rnorm for ~30
         ! iterations while it repairs the sub-sonic region -- see
         ! docs/newton_scaling_and_base_wall.md.)
         if (ok) then
            n_no_descent = 0
         else
            n_no_descent = n_no_descent + 1
         endif
         if (lm_tried .and. .not. ok) then
            ! The escape above was reached, which means repeating this
            ! iteration would reproduce it exactly. There is nothing left to
            ! try, so the solve stops here instead of counting to
            ! n_no_descent_max identical iterations.
            info = 2
            write(*,'(A)') ' (JFNK) no descent direction exists for the'// &
                 ' banded model at this state -- aborting'
            exit
         endif
         if (n_no_descent .ge. n_no_descent_max) then
            info = 2
            write(*,'(A,I0,A)') ' (JFNK) line search found no descent '// &
                 'step in ', n_no_descent_max, ' consecutive iterations'// &
                 ' -- aborting'
            exit
         endif

         ! Locate the cell carrying the largest scaled residual |F/D|, i.e.
         ! the term the merit is actually dominated by, over the whole domain
         ! 1..N. (Normalizing by max_j|u(k,j)| instead reports the base cell
         ! almost always, because the base holds the global maximum of |rho v|
         ! while carrying no wind.)
         amx = -1.0d0;  jworst = 1;  kworst = 1
         do jj = 1, N
            do kk = 1, nvar_jac
               if (nvar_jac .ge. 4) then
                  amx0 = abs(F(nvar_jac*(jj-1)+kk))                      &
                       / Drow(nvar_jac*(jj-1)+kk)
               else
                  amx0 = abs(F(nvar_jac*(jj-1)+kk))                      &
                       / D(nvar_jac*(jj-1)+kk)
               endif
               if (amx0 .gt. amx) then
                  amx = amx0
                  jworst = jj;  kworst = kk
               endif
            enddo
         enddo
         write(*,'(A,I4,A,ES11.3,A,ES10.2,A,ES9.2,A,I3,A,I4,A,F7.3,A,I1)') &
              ' (JFNK) it',iter,'  ||R||=',rnorm,'  ||Fs||2=',f2,        &
              '  lam=',lam,'  gm=',gmit,'  worst j=',jworst,            &
              ' r=',r(jworst),' k=',kworst
      enddo

      weno_mode = 0                 ! restore default reconstruction
      call set_ioniz_eq_sweep_state_kind(ieq_state_marching)

      ! What the solve could not describe. Both are properties of the
      ! iterates this solve visited, not of the returned state, and neither
      ! is an error by itself: a zeroed Jacobian color only degrades the
      ! preconditioner, and a truncated cycle only shortens the Krylov
      ! subspace. They are printed because a solve that ends info = 2 with
      ! either of them non-zero failed for a reason the iteration counts do
      ! not show.
      if (use_tr) write(*,'(A,I0,A,I0,A,I0,A,ES10.3)')                    &
           ' (JFNK) [TR] steps accepted ', n_tr_accept, ', rejected ',    &
           n_tr_reject, ', model refused by the ray test ',               &
           n_tr_model_bad, '; final radius', tr_delta
      if (n_jac_unresolved_tot .gt. 0 .or. n_gm_truncated .gt. 0)        &
         write(*,'(A,I0,A,I0,A)') ' (JFNK) probe states without an '//   &
              'admissible residual: ', n_jac_unresolved_tot,             &
              ' Jacobian color(s) zeroed, ', n_gm_truncated,             &
              ' GMRES cycle(s) truncated'

      ! Residual of the layer below the escape radius, resolved into its three
      ! rows and its worst cell -- the half of ||R|| that the single number
      ! reports only when it is the larger one. Printed HERE, before the
      ! best-iterate restore below, because this is the last point at which F
      ! and u are a consistent pair: the restore replaces Y (hence u) by the
      ! best iterate without re-evaluating F, and the line is not worth a
      ! further residual evaluation. So it refers to the LAST iterate visited,
      ! which is also the returned state whenever no restore happens.
      call write_resid_below_escape(' (JFNK)', F, u)

      ! Return the best iterate SEEN, not the last one visited, and judge
      ! convergence on it. The line search minimizes a merit; the solve is
      ! accepted on ||R||. These are different functionals, so the iterate with
      ! the smallest ||R|| is not in general the last one, and a run that
      ! reached ||R|| < resid_tol at some iterate and then wandered off has
      ! nonetheless produced a state that satisfies the code's own convergence
      ! criterion -- reporting it as a failure throws that state away and sends
      ! the caller back to time-marching (measured on examples/16: a 3.958e-5
      ! iterate discarded at iteration 2, docs/hd209_metal_stagnation.md).
      call unpack_U(Y, u);  call Apply_BC(u)
      gates_now = steady_gates_met(rnorm, u, resid_tol, fspread,        &
                          nvar_jac .ge. 4, carrier_relnorm_last)
      if ((gates_best .and. .not. gates_now) .or.                         &
          ((gates_best .eqv. gates_now) .and. rnorm_best .lt. rnorm)) then
         Y = Ybest;  f_sp = f_sp_best;  rnorm = rnorm_best
         fspread = fspread_best
         call adopt_background_of_best_iterate
         write(*,'(A,ES11.3)') ' (JFNK) returning best iterate, '//      &
              '||R||=', rnorm
      endif
      call unpack_U(Y, u)
      call unpack_carrier(Y, u, f_sp)
      call Apply_BC(u)
      ! ONE FINAL MEASUREMENT OF THE STATE ACTUALLY RETURNED, coupled route
      ! only. ||R|| is carried from the iterate that set rnorm_best, and its
      ! carrier half was evaluated with the limited slopes frozen at a
      ! different outer iterate; after a P51 restart the two can be several
      ! iterations apart. Re-evaluating here costs one equilibrium sweep and
      ! makes the number reported, the number the gate is applied to, and
      ! the state written to disk the same object. The three-unknown route
      ! keeps the residual it has always reported.
      if (nvar_jac .ge. 4) then
         weno_mode = 0
         f_sp_best = f_sp
         call eval_residual(Y, f_sp_best, f_sp_j, F, heat0, cool0)
         f_sp_best = f_sp_j
         f_sp = f_sp_best
         call resid_relnorm(F, u, rc, rnorm)
         carrier_relnorm_best = carrier_relnorm_last
         carrier_cellmax_best = carrier_cellmax_last
         carrier_cell_best    = carrier_cell_worst
      endif
      ! THE RESIDUAL THE GATE READS IS THE RESIDUAL OF THE STATE THAT IS
      ! WRITTEN, at the composition that state carries. Until now the number
      ! tested was the one the line search produced: eval_residual(Ytry,
      ! f_sp_j) with f_sp_j copied from the PREVIOUS iterate, so it is the
      ! residual of the new hydro state under the old composition. The state
      ! handed back carries the composition the accepted trial left, and its
      ! residual there is 4 to 28 times larger (section 155). Re-evaluating at
      ! the state's own composition, and iterating the sweep at fixed Y until
      ! the composition stops moving, makes the number reported, the number
      ! the gate is applied to, and the state written to disk one object.
      !
      ! Iterating rather than evaluating once: the composition at fixed Y has
      ! its own fixed point, reached in about ten sweeps, and one evaluation
      ! is still short of it by up to a factor two (WASP-121b). The joint
      ! fixed point of hydro balance AND chemical equilibrium is what a steady
      ! state is.
      !
      ! DEFAULT ON (`resid_at_own_composition`, parameters.f90): a state that
      ! is accepted has to satisfy its own residual, so the burden of a reason
      ! is on turning this OFF, and the only reason is reproducing a
      ! pre-section-155 number. `Resid tol` values calibrated against the old
      ! measure are correspondingly loose under this one.
      !
      ! THE SEED AND THE RESULT ARE SEPARATE ARRAYS (section 147), so the fixed
      ! point is reached by feeding each sweep's answer back as the next
      ! sweep's seed explicitly. f_sp holds the composition of the state being
      ! handed back and is the seed of the first pass; f_sc_out receives each
      ! answer and is copied into f_sp once it stops moving, so the state, its
      ! composition and the number the gate reads are one object.
      call get_environment_variable('EXHALE_RESID_SELFCONSISTENT', env_sc)
      if (resid_at_own_composition .or. len_trim(env_sc) .gt. 0) then
         weno_mode = 0
         rn_prev = huge(1.0d0)
         do it_sc = 1, n_selfconsistent_max
            call eval_residual(Y, f_sp, f_sc_out, F, heat0, cool0)
            f_sp = f_sc_out
            call resid_relnorm(F, u, rc, rnorm)
            if (abs(rnorm - rn_prev) .le. 1.0d-3*max(rnorm,1.0d-300))    &
               exit
            rn_prev = rnorm
         enddo
         write(*,'(A,I0,A,ES11.3)') ' (JFNK) residual at the state''s own'//&
              ' composition, ', it_sc, ' sweeps: ||R||=', rnorm
      endif

      if (steady_gates_met(rnorm, u, resid_tol, fspread,        &
                          nvar_jac .ge. 4, carrier_relnorm_last)) info = 0
      ! The frozen background the carrier transport reads must describe the
      ! state handed back, not the last state this solve happened to evaluate
      ! (docs/Update_EXHALE.md section 121).
      call install_background_of_adopted_state
      gate_rnorm_accepted = rnorm;  gate_fspread_accepted = fspread
      write(*,'(A,I0,A,ES11.3,A,ES10.3,A,I0)') ' (JFNK) done info=',info, &
           ' ||R||=',rnorm,'  flux spread=',fspread,                      &
           '  non-monotone accepts=', n_nonmonotone_accepts
      ! Both norms on the state handed back, with the window contributions
      ! and the worst cell's terms (section 145).
      block
        real*8, dimension(3,1-Ng:N+Ng) :: Rp145
        integer :: jp145, kp145
        Rp145 = 0.0d0
        do jp145 = 1, N
           do kp145 = 1, 3
              Rp145(kp145,jp145) = F(nvar_jac*(jp145-1)+kp145)
           enddo
        enddo
        call write_residual_breakdown(Rp145, u, heat0, cool0,             &
             'JFNK state handed back')
      end block

      ! Where the non-flatness reaches, for reading only: the acceptance test
      ! is steady_gates_met on the r >= r_flux window alone and it is unchanged.
      call flux_spread_above_radius(u, 1.03d0, fspread_103, dum_fmean)
      call flux_spread_above_radius(u, 1.10d0, fspread_110, dum_fmean)
      write(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A)')                           &
           ' (JFNK) flux spread by window: r>=1.03', fspread_103,         &
           '   r>=1.10', fspread_110, '   r>=r_flux', fspread,            &
           '  (only the last is the gate)'
      if (nvar_jac .ge. 4) then
         write(*,'(A,ES10.2,A,ES10.2,A,I4,A,I0)') ' (JFNK) carrier row '//&
              'of the returned state: volume-weighted',                  &
              carrier_relnorm_best, ', worst cell',                      &
              carrier_cellmax_best, ' at j=', carrier_cell_best,         &
              '; trials refused for the element budget: ',               &
              n_trial_headroom
         write(*,'(A,I0,A,I0,A,I0,A)') ' (JFNK) cost: ', nvar_jac,       &
              ' unknowns per cell, ', ncolor_jac, ' colors, ',           &
              n_resid_eval, ' residual evaluations'
      endif
      call write_uncertified_trial_count(' (JFNK)')
      ! WHICH GATES ARE LEFT -- each one asked separately, because there are
      ! three of them and "the one left" was a guess about which.
      !
      ! The chain this replaces tested the residual and the flux gate TOGETHER
      ! first and then fell through to the carrier row, so a state that met the
      ! residual gate and missed the other two printed "the CARRIER gate is the
      ! one left" while the flux gate was also unmet, by a factor 4.9.
      ! Measured on the 1 microbar hot Uranus (docs/p23_transport_on_state.md
      ! section 6): ||R|| = 1.547e-6 against 1e-5, carrier row 7.535e-3 against
      ! 1e-3, flux spread 2.447e-2 against 5e-3 -- two gates left, one named.
      ! One line per unmet gate says what is true whichever combination occurs.
      if (info .ne. 0) then
         if (rnorm .ge. resid_tol)                                        &
            write(*,'(A,ES10.3,A,ES10.3)') ' (JFNK) gate NOT met: '//     &
                 'residual ||R|| =', rnorm, ' >=', resid_tol
         if (flux_spread_th .gt. 0.0d0 .and. fspread .ge. flux_spread_th) &
            write(*,'(A,ES10.3,A,ES10.3)') ' (JFNK) gate NOT met: '//     &
                 'flux spread', fspread, ' >=', flux_spread_th
         if (nvar_jac .ge. 4 .and.                                        &
             carrier_relnorm_best .ge. carrier_resid_th)                  &
            write(*,'(A,ES10.3,A,ES10.3)') ' (JFNK) gate NOT met: '//     &
                 'carrier row', carrier_relnorm_best,                     &
                 ' >=', carrier_resid_th
      endif
      deallocate(Y,F,Ftry,dY,Ytry,dZ,D,Ybest,dY_lm,Drow,ab,abf,ipiv)
      end subroutine solve_steady_jfnk

      ! End of module
      end module steady_newton
