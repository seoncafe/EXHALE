   module ionization_equilibrium
	! Evaluate the ionization structure and the heating and cooling functions for a given temperature

	use global_parameters
   use ion_cell_state, only: ieq_cell, ion_rates
   use species_table, only: n_mion, mion_fsp, n_melem, melem_i0,        &
                            melem_top, mion_stage,                       &
                            isp_H2, isp_H2p, isp_H3p, isp_HeHp,          &
                            isp_OH, isp_H2O, isp_CO, iel_O, iel_C
   use utils
   use utils_ion_eq
   use Cooling_Coefficients      ! eval_cool, recombination/ionization rates
   use System_HeH                ! Equilibrium equations
	use System_HeH_TR
	use System_HeH_mol            ! molecular network
	use System_HeH_mol_metals     ! merged molecular network + metals
	use lower_column, only: q_h2_equilibrium
	use composition, only: base_h2_nuclei_fraction,                       &
	                       base_h2_composition_imposed
	use lyman_werner_photodissociation, only: e_lw_fragment_erg
	! Oxygen chemistry (the A2 option, docs/a2_oxygen_option_design.md):
	! the FUV photolysis bands and the CO reservoir.
	use water_photolysis, only: n_fuv_band,                               &
	                            heat_per_water_dissociation,              &
	                            heat_per_hydroxyl_dissociation
	use oxygen_rates, only: co_equilibrium_density,                       &
	                        oxygen_chemical_equilibrium_fractions
   use System_HeH_metals
   use System_HeH_TR_metals      ! merged He-triplet + metals system
	! Constrained element-conserving chemical equilibrium in positive
	! species densities, tracked from the radiation-free molecular limit to
	! the cell's field (docs/supersonic_molecular_base.md section 11.5-B).
	use constrained_chemical_equilibrium, only:                           &
	                             equilibrium_from_molecular_limit
   use charge_exchange, only: cx_set_cell, cx_metal_base,  &   ! Huang Table 4 charge exchange
                              he_h_cx_rates,               &   ! He <-> H pair (group B)
                              cx_add_to_turnover               ! cx bound of a row's turnover
   use System_H
   use newton_solver, only: solve_ieq   ! Task 2: analytic-Jacobian Newton (+ hybrd1 fallback)
   use opacity_models            ! opacity_pT_factor for the 'P' model

   implicit none

	! molecular species densities (cols 1 H2, 2 H2+, 3 H3+, 4 HeH+;
	! zero unless thereis_mol).  Module state: written by the equilibrium
	! solve, read by write_output for the extra output columns.
	real*8, dimension(:,:), allocatable :: nmol_eq

	! Lyman-Werner photodissociation diagnostics, filled by the equilibrium
	! solve when a band flux is supplied and read by write_lyman_werner:
	! star-ward H2 column [cm^-2], the self-shielding factor the rate
	! carries (Richings, Schaye & Oppenheimer 2014; lyman_werner.f90
	! sec. 2), and the resulting dissociation rate [s^-1]. Untouched
	! (shielding 1, rate 0) when the run supplies no band flux.
	real*8, dimension(:), allocatable :: NH2_col_lw
	real*8, dimension(:), allocatable :: f_shield_lw
	real*8, dimension(:), allocatable :: k_lw_diss
	! Transmission of the Lyman-Werner LINES down to each cell's inner face,
	! 1 - A: the fraction of the 912-1110 A band that H2 has NOT already
	! taken, which the H2O and OH continua of the same interval multiply
	! their own attenuation by (water_photolysis.f90 sec. 3). Identically 1
	! without a Lyman-Werner band flux. a_lines_lw is the largest A the
	! column reached, i.e. the deepest cell's removed fraction; the run
	! reports once if the H2 column leaves the checked range of either of
	! the two self-shielding fits (lyman_werner.f90 sec. 2).
	real*8, dimension(:), allocatable :: tr_lines_lw
	! H2 photoionization rate of the last sweep [s^-1], kept for the H2
	! budget decomposition write_output prints: it is the one loss channel
	! of the H2 row that cannot be rebuilt from the output files.
	real*8, dimension(:), allocatable :: P_H2_eq
	real*8 :: a_lines_lw = 0.0d0
	logical :: warned_lw_fit_range = .false.
	! Largest H2 columns at which the two self-shielding fits of
	! lyman_werner.f90 were checked against an exact calculation.
	!
	! NH2_db96_max: Draine & Bertoldi (1996) sec. 5.2, the agreement is
	! excellent "even out to the largest column densities considered,
	! N2 = 3 x 10^21 cm^-2".  They state no upper limit on the fit itself,
	! so a deeper column is undemonstrated rather than outside a published
	! range.  This fit sets the band share A only.
	!
	! NH2_richings_max: Richings, Schaye & Oppenheimer (2014) sec. 3.2 give
	! an explicit range -- the fit "agrees with CLOUDY to within 30 per
	! cent at 100 K for NH2 < 10^21 cm^-2, and to within 60 per cent at
	! 5000 K for NH2 < 10^20 cm^-2".  This fit sets the photodissociation
	! RATE, and 10^21 cm^-2 is the looser of the two bounds, so it is the
	! one a molecular base crosses first.
	real*8, parameter :: NH2_db96_max     = 3.0d21
	real*8, parameter :: NH2_richings_max = 1.0d21

	! Oxygen-chemistry state (the A2 option; all zero unless
	! thereis_oxychem).  Written by the equilibrium solve, read by
	! write_output for output/Ion_species.txt, output/Oxygen_chemistry.txt
	! and output/FUV_bands.txt, and by the heat breakdown.
	!   nox_eq   cols 1 OH, 2 H2O, 3 CO [cm^-3]
	!   n_o1d_eq O(1D) density from its local steady state [cm^-3]
	!   NH2O_col, NOH_col   star-ward columns of the two FUV absorbers
	!   j_h2o_fuv, j_oh_fuv band-resolved photodissociation rates [s^-1]
	!   tau_fuv             band-resolved star-ward optical depth
	!   heat_fuv            FUV photolysis heating [erg cm^-3 s^-1]
	real*8, dimension(:,:), allocatable :: nox_eq
	real*8, dimension(:),   allocatable :: n_o1d_eq
	real*8, dimension(:),   allocatable :: NH2O_col, NOH_col
	real*8, dimension(:,:), allocatable :: j_h2o_fuv, j_oh_fuv, tau_fuv
	real*8, dimension(:),   allocatable :: heat_fuv

	! The coefficient state each cell was left in by the last equilibrium
	! sweep, kept so the carrier transport (diffusive_photochemistry) can
	! evaluate the SAME chemistry rows this module solves without rebuilding
	! the photoionization, recombination and collisional rates. It is a copy
	! of the ieq_cell each cell was solved with, which is why a field added
	! to that type reaches the transport operator with no second edit here.
	! bg_ready is false until the first sweep has filled it: a transport step
	! before that would run on zeros.
	type(ion_rates), dimension(:), allocatable :: bg_cell
	logical :: bg_ready = .false.

	! Run-wide totals of the atomic ionization root validation, reported once
	! at the end of the run (EXHALE_main) next to the Newton usage counters:
	! stored states rejected as a starting point, cell solves that needed a
	! second or third starting point, first roots that lay outside the
	! physical simplex, and cells where no starting point produced an
	! admissible root. All zero for a run that never leaves the simplex.
	integer, save :: ieq_n_reseed = 0, ieq_n_retry  = 0
	integer, save :: ieq_n_unphys = 0, ieq_n_noroot = 0

	! Run-wide total of molecular cells whose equilibrium roots all left the
	! physical simplex, so the closest one was clamped back onto the element
	! budget (ioniz_eq). Zero for a run whose molecular solve stays inside the
	! simplex everywhere, and zero for an atomic run.
	integer, save :: ieq_n_mol_clamped = 0

	! Run-wide histogram of the hybrd1 exit code of the molecular cell solves,
	! indexed by info (0 = improper input or iflag < 0, 1 = converged to tol,
	! 2 = iteration limit, 3 = xtol too small, 4/5 = no progress). Every
	! attempt of every molecular cell is counted, so the total exceeds the
	! cell count whenever a cell needs its second or third starting point.
	! Zero for an atomic run.
	integer, save :: ieq_n_mol_info(0:5) = 0

	! Acceptance tolerance on the normalized reaction residual of an
	! equilibrium state (docs/supersonic_molecular_base.md section 11.5-A;
	! docs/Update_EXHALE.md section 113). A state is accepted as a ROOT of
	! the ionization/chemical network only when it lies inside the element
	! bounds AND the largest reaction imbalance of any row, in units of the
	! row's turnover rate (normalized_reaction_residual), is at or below
	! this value; the solver exit code is neither sufficient nor necessary
	! (MINPACK info = 1 is an xtol statement about the step, and info = 4
	! routinely returns finished roots it cannot certify).
	!
	! The value is set from the measured populations of the full 8-case
	! regression matrix plus the He/H = 1 molecular arm (2026-08-31):
	! across the matrix every solver-converged accepted root sits at
	! res <= 7.8e-8 and every root accepted without solver convergence at
	! res <= 2.7e-7 -- the scale is hybrd1's xtol = sqrt(eps) ~ 1.5e-8
	! times a scaled-Jacobian norm of order 1-10 -- while the smallest
	! above-tolerance acceptance of the matrix (a cold-start base cell
	! recovering over the next sweeps) is 1.1e-3 and the poisoned events
	! of the He/H = 1 arm reach 2e-2..1.4e6. 1e-6 sits 3.7x above the
	! measured root tail, in a decade ([1e-6, 1e-5)) that is empty in
	! every matrix histogram, and 1100x below the matrix's smallest
	! non-root; the He/H = 1 arm additionally carries one borderline
	! accepted iterate at 2.1e-6, which this value classifies (and marks)
	! as a non-root rather than stretching the root band to cover it.
	real*8, parameter :: ieq_res_tol = 1.0d-6

	! Relaxation amnesty and its limit. The cold-start relaxation of a
	! healthy run genuinely passes through non-root acceptances and
	! recovers: measured on the regression matrix, the metals-on molecular
	! gates accept res ~ 28 states at step 0 in the shielded base and res
	! ~ 2.6e-3 states in one or two cells for up to 53 CONSECUTIVE sweeps
	! (mol_ir_bands, cell 260) before the solve lands on roots for the
	! rest of the run. An unconditional stop on the first non-root would
	! therefore kill working configurations, and no severity threshold
	! separates them (the healthy step-0 spike, 28, exceeds the poisoned
	! He/H = 1 step-0 event, 2e-2). What does separate a recovering
	! transient from a solution RESTING on a non-root is persistence: a
	! steady or pseudo-steady state re-evaluates the same cell every sweep,
	! so a wind built on a non-root fails the same cell indefinitely. A
	! cell may therefore carry a non-root acceptance -- loudly counted and
	! reported, never silent -- for at most this many consecutive sweeps
	! (~19x the largest healthy streak measured); one sweep more and the
	! run stops with full diagnostics (nonroot_equilibrium_stop).
	integer, parameter :: ieq_nonroot_streak_stop = 1000
	! Consecutive sweeps each cell has spent on a non-root acceptance
	! (reset to zero by any accepted root), and the run-wide peak.
	integer, allocatable, save :: ieq_nonroot_streak(:)
	integer, save :: ieq_nonroot_streak_peak = 0

	! Run-wide acceptance statistics of the He-branch equilibrium solves
	! (docs/Update_EXHALE.md section 113), reported once at the end of the
	! run: how many cell states were accepted as (1) solver-converged roots,
	! (2) roots without solver convergence, (3) projected/handback states
	! whose rechecked residual still marks a root, (4) NON-ROOT states
	! accepted under the relaxation amnesty, (5) roots of the constrained
	! element-conserving continuation solve
	! (constrained_chemical_equilibrium), with the largest normalized
	! reaction residual each class carried; and the residual-decade
	! histograms of the physical iterates (converged and not), which locate
	! the gap between roots and non-roots the tolerance sits in.
	integer, save :: ieq_acc_n(5) = 0
	real*8,  save :: ieq_acc_resmax(5) = 0.0d0
	integer, save :: ieq_hist_conv(0:15) = 0
	integer, save :: ieq_hist_uncv(0:15) = 0
	! Print budget of the acceptance report lines: a pathological
	! run must not flood the log; the counters above keep the population.
	integer, save :: ieq_acc_nprint = 0
	integer, parameter :: ieq_acc_nprint_max = 2000

	! Run-wide totals of the constrained element-conserving continuation
	! solve (constrained_chemical_equilibrium, section 11.5-B): molecular
	! cells promoted to it because the three fraction-layout starting points
	! produced no root, cells whose promoted candidate was accepted as a
	! root by the same judge, hybrd1 calls the continuation spent, and the
	! wall time it took. All zero for a run whose molecular cells all find a
	! root the ordinary way, and for an atomic run.
	integer, save :: ieq_n_cce_attempt = 0
	integer, save :: ieq_n_cce_root    = 0
	integer, save :: ieq_n_cce_solve   = 0
	real*8,  save :: ieq_cce_seconds   = 0.0d0

	contains

	subroutine ioniz_eq_allocate_arrays
	! Allocate the grid-sized module arrays once the number of cells N is
	! known; called from EXHALE_main right after input_read. The values are
	! the initializers the declarations used to carry.

	allocate(nmol_eq(1-Ng:N+Ng,4))
	allocate(bg_cell(1-Ng:N+Ng))
	allocate(ieq_nonroot_streak(1-Ng:N+Ng))
	ieq_nonroot_streak = 0
	allocate(NH2_col_lw(1-Ng:N+Ng), f_shield_lw(1-Ng:N+Ng),               &
	         k_lw_diss(1-Ng:N+Ng), tr_lines_lw(1-Ng:N+Ng),                &
	         P_H2_eq(1-Ng:N+Ng))
	allocate(nox_eq(1-Ng:N+Ng,3), n_o1d_eq(1-Ng:N+Ng),                    &
	         NH2O_col(1-Ng:N+Ng), NOH_col(1-Ng:N+Ng),                     &
	         j_h2o_fuv(1-Ng:N+Ng,n_fuv_band),                             &
	         j_oh_fuv(1-Ng:N+Ng,n_fuv_band),                              &
	         tau_fuv(1-Ng:N+Ng,n_fuv_band), heat_fuv(1-Ng:N+Ng))

	nmol_eq     = 0.0d0
	NH2_col_lw  = 0.0d0
	f_shield_lw = 1.0d0
	tr_lines_lw = 1.0d0
	k_lw_diss   = 0.0d0
	P_H2_eq     = 0.0d0
	nox_eq      = 0.0d0
	n_o1d_eq    = 0.0d0
	NH2O_col    = 0.0d0
	NOH_col     = 0.0d0
	j_h2o_fuv   = 0.0d0
	j_oh_fuv    = 0.0d0
	tau_fuv     = 0.0d0
	heat_fuv    = 0.0d0

	end subroutine ioniz_eq_allocate_arrays

	subroutine ioniz_eq(T_in,n_io,f_sp_io,heat_out,cool_out,q)
      	 		  
	integer :: j,im
	logical :: usednt                     ! Task 2: Newton-vs-fallback flag

	real*8, dimension(1-Ng:N+Ng),   intent(in) :: T_in
	! Density and species fractions are read on entry and overwritten with the
	! equilibrium result on exit (in place); callers must not alias them.
	real*8, dimension(1-Ng:N+Ng),   intent(inout) :: n_io
	real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp_io

	real*8, dimension(1-Ng:N+Ng) ::  T_K      ! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng) ::  nh,nhi,nhii,                   & ! Species densities
	                                 nhe,nhei,nheii,nheiii,nheiTR,  &
	                                 ne,n_in_dim,n_tot
	! Ionized fraction of the H+He nuclei, for the SvS85 secondary ionization.
	real*8, dimension(1-Ng:N+Ng) ::  xion
   ! Metal ion densities in canonical species_table order (col im maps
   ! to f_sp column mion_fsp(im)); used throughout in place of named
   ! scalars for each ion so the driver scales with the number of metals.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  nm
   ! Total density of each metal element (canonical element order),
   ! held constant across the cell sweep (sum of its three stages).
   real*8, dimension(1-Ng:N+Ng,n_melem) ::  nm_tot
   ! ELEMENT total of each metal element, i.e. every nucleus of it in the
   ! cell whatever carries it. It differs from nm_tot only for oxygen and
   ! carbon, and only with the oxygen chemistry on: nm_tot is what the metal
   ! ionization block is solved against (the FREE oxygen family, and the
   ! carbon not locked in CO), nm_el is the conserved element. Identical
   ! arrays without the option.
   real*8, dimension(1-Ng:N+Ng,n_melem) ::  nm_el
   ! Carbon monoxide density of each cell [cm^-3], from the CO <-> C + O
   ! chemical equilibrium, and the free oxygen family it leaves behind.
   real*8, dimension(1-Ng:N+Ng) ::  nCO_cell
   ! Row of the first oxygen-carrier unknown, and band loop index.
   integer :: iox, ib

   ! Metal recombination and collisional ionization rates for each ion from
   ! eval_cool (canonical order); bridged to the named rc*/a_ion_*
   ! scalars below for the params packing.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  rec_m,aion_m

   ! Photo ionization rates
   real*8, dimension(1-Ng:N+Ng) ::  P_HI,P_HeI,P_HeII,P_HeITR
   real*8, dimension(1-Ng:N+Ng) ::  P_H2      ! (molecular; zero unless mol)
   ! Metal photoionization rates for each ion (canonical order) from PH_heat.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  P_m
                       	
   ! Heating, cooling
   real*8, dimension(1-Ng:N+Ng) ::  heat,cool    
                 
   ! Recombination coefficients
   real*8, dimension(1-Ng:N+Ng) ::  rchiiB,rcheiiB,rcheiiiB,rcheiTR

   ! He recombination radiation -> H ionization coupling scratch
   ! (use_he_rec_coupling; zero-effect when off).
   real*8, dimension(1-Ng:N+Ng) ::  rcheiiB_hrc,dP_HI_hrc,dheat_hrc

	real*8, dimension(1-Ng:N+Ng) :: q13,q31a,q31b,Q31
	real*8 :: A31

   ! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng) ::  a_ion_HI,a_ion_HeI,a_ion_HeII,a_ion_HeITR

	! Equilibrium system setup
   real*8 :: tol,dpmpar
   real*8, dimension(60) :: params
   ! Each element's metal coefficients for ion_system_HeH_metals (canonical order),
   ! built per cell from the 2D rate arrays and handed to set_metal_coeffs.
   real*8, dimension(n_melem) :: meg_ntot,meg_g0,meg_g1,meg_b0,meg_b1, &
                                 meg_a1,meg_a2
   ! Highest stage per element handed to set_metal_coeffs (2 = three-stage,
   ! 1 = two-stage); see species_table::melem_top.
   integer, dimension(n_melem) :: meg_top
   integer :: i0,top
   ! Base index of the first metal element's X+ unknown in sys_x: 4 normally,
   ! 5 when the He triplet occupies x(4) (merged HeITR+metals system), and
   ! metal_row_base() = 8 or 9 in the molecular layout, where x(4..7) are the
   ! H nuclei bound into H2/H2+/H3+/HeH+ and x(8) the triplet.
   integer :: mbase
   ! Molecular solve (thereis_mol): count of cells whose roots all left the
   ! physical simplex, so the closest one was clamped back onto the element
   ! budget, and the local gas pressure [bar] the H2 dissociation
   ! equilibrium of the retry seed is evaluated at.
   integer :: n_mol_clamped
   ! Histogram of the molecular hybrd1 exit code over the sweep (see
   ! ieq_n_mol_info above); combined over threads by the reduction below.
   integer :: n_mol_info(0:5)
   real*8  :: pbar_loc

   ! Ionization solve: validation of the returned root.
   ! The equilibrium systems are polynomial and possess roots outside the
   ! physical simplex (negative stage fractions, or ionized stages of one
   ! element summing above its nucleus total); such a root makes the neutral
   ! density negative and is not a solution of the physical problem even when
   ! the algebraic residual vanishes. Each cell therefore tries up to three
   ! starting points and keeps the best admissible root; x_root_best holds it,
   ! ok_rank grades the attempt (0 none, 1 physical, 2 physical + converged),
   ! and the counters feed the one-line sweep summary.
   ! Largest N_eq over the layouts: molecular (7) + triplet (1) + the two
   ! oxygen carriers + two per metal element. The oxygen pair was added on
   ! 2026-08-30 with the oxygen chemistry, and leaving it out was a real
   ! out-of-bounds write, not a spare-capacity question: with every option
   ! on, N_eq = 30 against a bound of 28, and the copy of the root into
   ! x_root_best ran two elements past the end of a thread-private stack
   ! array. -O3 did not notice; run_fcheck.sh did.
   integer, parameter :: n_x_max = 10 + 2*n_melem
   integer :: iatt, info_ieq, ok_rank, best_rank
   integer :: n_ieq_reseed, n_ieq_retry, n_ieq_unphys, n_ieq_fail
   logical :: conv_ieq, phys_ieq
   real*8, dimension(n_x_max) :: x_root_best
   ! Molecular cell whose roots all left the physical simplex: how far the
   ! least-offending one lies outside it, so the closest to a state can be
   ! selected and clamped back onto the element budget.
   real*8  :: viol, viol_best
   ! Acceptance statistics of this sweep (section 113): counts and largest
   ! normalized reaction residual by acceptance class (1 = solver-converged
   ! root, 2 = root without solver convergence, 3 = projected/handback
   ! state that rechecks as a root, 4 = non-root under the relaxation
   ! amnesty, 5 = root of the constrained element-conserving continuation
   ! solve), and the residual-decade histograms of the physical iterates.
   integer :: n_acc(5), hist_conv(0:15), hist_uncv(0:15)
   real*8  :: acc_resmax(5)
   integer :: acc_class
   real*8  :: acc_res, res_att
   ! Root selection (section 113): residual of the kept rank-1/2 root, the
   ! first bounded iterate whose residual failed the tolerance (the state
   ! of last resort, matching the old rank-1 order), and whether any
   ! attempt stayed inside the simplex / supplied a clamp candidate.
   real*8  :: res_best, res_nonroot
   real*8, dimension(n_x_max) :: x_nonroot
   integer :: info_best, info_nonroot
   logical :: have_nonroot, have_clamp

   ! Promotion of a molecular cell with no root to the constrained
   ! element-conserving continuation solve (section 11.5-B): the candidate
   ! it returns, whether the continuation reached the cell's full radiation
   ! field, whether the same judge marked the candidate a root, its
   ! residual, and the cost -- hybrd1 calls and wall time. The clock counts
   ! are integer(8) because a nanosecond-resolution system_clock overflows
   ! the default integer in a few seconds.
   real*8, dimension(n_x_max) :: x_cce
   logical :: cce_full, cce_root
   integer :: n_cce_fs
   real*8  :: res_cce
   integer(8) :: clk_beg, clk_end, clk_rate
   ! Sweep totals of the same, combined over threads by the reduction below.
   integer :: n_cce_attempt, n_cce_root, n_cce_solve
   real*8  :: cce_seconds

   ! Output heating,cooling and absorbed energy
   real*8, dimension(1-Ng:N+Ng),intent(out) :: heat_out,cool_out,q

   !----------------------------------------------------------!      
   ! Global parameters
      
   ! Numerical tolerance for system solution
   tol = sqrt(dpmpar(1))

	!----------------------------------!
	
	! Preliminary profiles exctraction
	
	! Dimensional total number density profile and temperature
	n_in_dim = n_io*n0
	T_K      = T_in*T0
		
	! Extract species profiles
	nhi    = f_sp_io(:,1)*n_in_dim    ! HI
	nhii   = f_sp_io(:,2)*n_in_dim    ! HII
	if (thereis_He) then
		nhei   = f_sp_io(:,3)*n_in_dim    ! HeI
		nheii  = f_sp_io(:,4)*n_in_dim    ! HeII
		nheiii = f_sp_io(:,5)*n_in_dim    ! HeIII
		nheiTR = f_sp_io(:,6)*n_in_dim    ! HeITR
	else
		! Enforce condition of zero helium
		nhei   = 0.0
		nheii  = 0.0
		nheiii = 0.0
		nheiTR = 0.0
	endif

    ! Metal ion densities in canonical species_table order (col im maps to
    ! f_sp column mion_fsp(im)).
    do im = 1,n_mion
       nm(:,im) = f_sp_io(:,mion_fsp(im))*n_in_dim
    enddo

	! molecular species (zero when thereis_mol is off)
	if (thereis_mol) then
		nmol_eq(:,1) = f_sp_io(:,isp_H2)  *n_in_dim
		nmol_eq(:,2) = f_sp_io(:,isp_H2p) *n_in_dim
		nmol_eq(:,3) = f_sp_io(:,isp_H3p) *n_in_dim
		nmol_eq(:,4) = f_sp_io(:,isp_HeHp)*n_in_dim
	endif

	! oxygen carriers (zero when thereis_oxychem is off)
	if (thereis_oxychem) then
		nox_eq(:,1) = f_sp_io(:,isp_OH) *n_in_dim
		nox_eq(:,2) = f_sp_io(:,isp_H2O)*n_in_dim
		nox_eq(:,3) = f_sp_io(:,isp_CO) *n_in_dim
	endif

    ! Total density of each metal element (sum of its three stages, in the
    ! same neutral+singly+doubly order as the original nC=nci+ncii+nciii).
    do im = 1,n_melem
       i0 = melem_i0(im)
       if (melem_top(im) .ge. 2) then
          nm_tot(:,im) = nm(:,i0) + nm(:,i0+1) + nm(:,i0+2)
       else
          ! Two-stage element: neutral + singly ionized only.
          nm_tot(:,im) = nm(:,i0) + nm(:,i0+1)
       endif
    enddo
    nm_el    = nm_tot
    nCO_cell = 0.0d0

    ! ---- oxygen and carbon element totals, and the CO reservoir ----
    ! With the oxygen chemistry on, the O I / C I columns hold FREE ATOMIC
    ! oxygen and carbon: the rest of each element sits in OH, H2O and CO.
    ! The element total -- the quantity that is conserved and that the
    ! abundance melem_ab fixes -- is therefore the stage sum PLUS the
    ! molecular carriers, counted through bsp_nO / bsp_nC of the species
    ! table. Without this the closure
    !   n_O,tot = n(OI)+n(OII)+n(OIII)+n(OH)+n(H2O)+n(CO)
    ! is not closed and half the element leaks away one sweep at a time.
    !
    ! CO is then taken out of both elements at the CO <-> C + O chemical
    ! equilibrium of the local (n, T) and frozen there for the solve
    ! (decision D4: CO reacts with nothing in the audited set, but it holds
    ! 45-46% of the oxygen at every measured level, so a network that gave
    ! the whole oxygen abundance to the water family would over-supply the
    ! OH cycle by about a factor two). What is left is the FREE OXYGEN
    ! FAMILY shared by O I/II/III, OH and H2O, and the carbon left for the
    ! metal block's C I/II/III.
    if (thereis_oxychem) then
       nm_el(:,iel_O) = nm_tot(:,iel_O)                                  &
                      + nox_eq(:,1) + nox_eq(:,2) + nox_eq(:,3)
       nm_el(:,iel_C) = nm_tot(:,iel_C) + nox_eq(:,3)
       if (oxygen_transport) then
          ! CO is a TRANSPORTED reservoir: the carrier operator moved it,
          ! and this sweep takes the CO the cell has rather than the CO the
          ! local (n, T) would make. The difference is not cosmetic -- the
          ! equilibrium form re-forms CO in a cool outer wind, where a
          ! transported one would have been carried out of the molecular
          ! layer and quenched.
          nCO_cell = max(nox_eq(:,3), 0.0d0)
          nCO_cell = min(nCO_cell, nm_el(:,iel_C))
          nCO_cell = min(nCO_cell, nm_el(:,iel_O))
       else
          do j = 1-Ng,N+Ng
             nCO_cell(j) = co_equilibrium_density(nm_el(j,iel_C),        &
                                                  nm_el(j,iel_O), T_K(j))
          enddo
       endif
       nm_tot(:,iel_O) = max(nm_el(:,iel_O) - nCO_cell, 0.0d0)
       nm_tot(:,iel_C) = max(nm_el(:,iel_C) - nCO_cell, 0.0d0)
       nox_eq(:,3)     = nCO_cell
    endif
	
	! Total number densities (with molecules: H and He NUCLEI totals --
	! the molecular system conserves elements, and nmol_eq is zero otherwise)
	nh  = nhi  + nhii
	nhe = nhei + nheii + nheiii
	if (thereis_mol) then
		nh  = nh  + 2.0d0*(nmol_eq(:,1) + nmol_eq(:,2))                  &
		          + 3.0d0*nmol_eq(:,3) + nmol_eq(:,4)
		nhe = nhe + nmol_eq(:,4)
	endif
	! OH carries one H nucleus and H2O two (bsp_nH of the species table).
	! Written as its own statement so the sum above stays bit-for-bit the
	! one a run without the oxygen chemistry evaluates.
	if (thereis_oxychem) nh = nh + nox_eq(:,1) + 2.0d0*nox_eq(:,2)

	! Ionized fraction of the H+He nuclei, for the SvS85 secondary-ionization
	! partition (metals excluded; nheiTR is neutral and not in the numerator).
	xion = min(max((nhii + nheii + nheiii)/max(nh + nhe, 1.0d-99), 0.0d0), 1.0d0)

	! Free electron density (assuming overall neutrality; nm adds the
	! metal electrons under the eos_metals policy, nmol_eq the molecular-ion
	! electrons -- it is zero for an atomic run, so the sum is unchanged there)
	call calc_ne(nhii,nheii,nheiii,ne,nm,nmol_eq)

	! Total gas-particle density (electrons excluded), i.e. the density of
	! third bodies M for the three-body molecular reactions R12/R13/R15, and
	! -- with the electrons added back -- the gas pressure p = (n_tot+n_e)kT
	! of the chemical-equilibrium retry seed below.  It is NOT n_in_dim:
	! n_in_dim = rho/m_H is a MASS density in m_H units (calc_rho weights
	! each species by bsp_mass), so it over-counts the particles by the mean
	! particle mass -- 2.3x at an H2-rich base, up to 4x in a He-dominated
	! one.  Only the molecular path reads it, so an atomic run is unchanged.
	if (thereis_mol) then
		if (thereis_oxychem) then
			call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol_eq,   &
			               nox_eq)
		else
			call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol_eq)
		endif
	else
		n_tot = 0.0d0
	endif

	! Cell-by-cell pressure-broadening factor for the opacity ('P' model).
	! opacity_pT_factor returns 1.0 for all other models, so opa_pf=1
	! and the column densities are unchanged (bit-identical).
	do j = 1-Ng,N+Ng
		opa_pf(j) = opacity_pT_factor((nh(j)+nhe(j)+ne(j))*kb_erg*T_K(j))
	enddo

	!----------------------------------!

	! The FUV photon field of the molecular layer, for ALL of its absorbers
	! at once: H2 in the Lyman and Werner lines (Draine & Bertoldi 1996) and
	! H2O and OH in their continua (the A2 oxygen option).  Over 912-1110 A
	! the three share one beam, so they are solved together --
	! fuv_lw_photon_field in utils_ion_eq states why and how the shared beam
	! splits.  The columns use the same radial integration and opa_pf
	! weighting as every other absorber column, and the incoming (pre-solve)
	! densities, like the photoionization columns built inside PH_heat_HHe.
	call fuv_lw_photon_field(nmol_eq(:,1), nox_eq(:,2), nox_eq(:,1), T_K, &
	                         NH2_col_lw, NH2O_col, NOH_col,               &
	                         f_shield_lw, tr_lines_lw, tau_fuv,           &
	                         k_lw_diss, j_h2o_fuv, j_oh_fuv, a_lines_lw)
	if (thereis_mol .and. F_LW_star .gt. 0.0d0 .and.                      &
	    .not. warned_lw_fit_range) then
		if (maxval(NH2_col_lw) .gt. NH2_richings_max) then
			warned_lw_fit_range = .true.
			write(*,'(a,es9.2,a,es9.2,a,es9.2,a)') ' (ioniz_eq)'//        &
			  ' WARNING: the star-ward H2 column reaches ',               &
			  maxval(NH2_col_lw), ' cm^-2. The self-shielding fits are'// &
			  ' checked only to ', NH2_richings_max,                      &
			  ' cm^-2 (Richings, Schaye & Oppenheimer 2014, which sets'// &
			  ' the photodissociation rate) and ', NH2_db96_max,          &
			  ' cm^-2 (Draine & Bertoldi 1996, which sets the band'//     &
			  ' fraction the H2 lines remove), so the deepest cells are'//&
			  ' extrapolating at least the first of the two.'
		endif
	endif

	!----------------------------------!

    !---- Photoionization and photoheating ----!
      
	if (thereis_He) then
		if (thereis_mol) then
			call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm, xion,    &
			         P_HI,P_HeI,P_HeII,P_HeITR, P_m,             &
			         heat,q, nmol_eq(:,1),P_H2)
		else
      	call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm, xion,      &
      			     P_HI,P_HeI,P_HeII,P_HeITR, P_m,        &
      			     heat,q)
		endif
	else
		call PH_heat_H(nhi, xion, P_HI,heat,q)
		P_m = 0.0
	endif

	! excited-H H(n=2) feedback (zero unless use_excited_H). The Balmer
	! photoionization of H(n=2) adds an effective HI photoionization rate
	! [s^-1] (proton source), and the Balmer photoelectric (+ optional
	! collisional de-excitation) heating adds to the photoheating rate
	! [erg cm^-3 s^-1]. Both global arrays are filled from the previous
	! converged outer pass by excited_hydrogen::excited_H_update, so the
	! coupling is decoupled from the hydro sub-step. Default-off => no change.
	if (use_excited_H) then
		gph_ground_HI = P_HI         ! capture pure ground-state rate (pre-n=2)
		P_HI = P_HI + gph_balmer_HI
		heat = heat + heat_balmer
	endif


   !----------------------------------!
      
   ! Evaluate cooling rates and recombination/collisional 
   ! 	ionization rates
      
    ! nmol_eq gives the cooling the same electron density this routine
    ! balances the ionization against, and the H3+/H2 densities its infrared
    ! cooling channel needs (zero for an atomic run). Nothing is added to
    ! `cool` after this call: the marching temperature update rebuilds the
    ! cooling from eval_cool alone, so any term added here would be a source
    ! the marching relaxes without and the steady residual demands.
    call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,           &
			   	   	rchiiB,rcheiiB,rcheiiiB, rec_m,             &
			    	a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,         &
			    	cool, nheiTR=nheiTR, a_ion_HeITR=a_ion_HeITR,  &
			    	nmol=nmol_eq, nox=nox_eq)

	! Capture the ground-state H proton-budget coefficients on
	! every pass (the converged pass is the one read out by write_excited_H).
	if (use_excited_H) then
		cion_HI  = a_ion_HI          ! collisional ionization [cm^3 s^-1]
		arec_HII = rchiiB            ! case-B recombination  [cm^3 s^-1]
	endif

	! Charge-exchange rate coefficients are evaluated per cell below
	! (cx_set_cell) before each metal ionization solve.

	if (thereis_HeITR) then
		call HeITR_coeffs(T_K,rcheiTR,rcheiiB,A31,q13,q31a,q31b,Q31)
		! NOTE: rcheiiB is alpha1 from Oklopcic - being overwritten
	endif

	! He recombination radiation ionizing H I (Draine 2011 on-the-spot y/z;
	! default off => no change). Uses the pre-solve (lagged) densities, like
	! the other lagged rate terms; corrects the He II recombination coefficient
	! rcheiiB, adds an H I photoionization rate to P_HI, and adds photoelectron
	! heating to heat. The He II recombination *cooling* (rec_cool_HeII =
	! kT alpha_B) is left unchanged; the mismatch is <= y alpha_1 kT, negligible.
	if (use_he_rec_coupling .and. thereis_He) then
		call he_rec_coupling(T_K, nhi, nhei, nheii, nheiTR, ne,           &
		                     A31, q31a, q31b,                             &
		                     rcheiiB_hrc, dP_HI_hrc, dheat_hrc)
		rcheiiB = rcheiiB_hrc
		P_HI    = P_HI + dP_HI_hrc
		heat    = heat + dheat_hrc
	endif

	! Penning ionization heating: He(2^3S)+H0 -> He(1^1S)+H+ + e- releases the
	! electron kinetic energy e_th_HeI - e_th_HeTR - e_th_HI (= 6.2 eV) into the
	! gas. Q31 is the total ionization rate, so only its Penning branch
	! (f_penning_HeI23S) carries this exothermicity; the associative branch
	! ends in HeH+ and has a different one. Lagged (pre-solve) densities, like
	! every other channel above; nheiTR is the same array he_rec_coupling
	! already consumes.
	if (thereis_HeITR) heat = heat                                    &
	     + f_penning_HeI23S*nheiTR*nhi*Q31                             &
	       *(e_th_HeI - e_th_HeTR - e_th_HI)/erg2eV

	! Molecular Penning ionization heating: He(2^3S)+H2 -> He(1^1S)+H2+ + e-
	! releases the electron kinetic energy (e_th_HeI - e_th_HeTR) - e_th_H2
	! (= 24.6 - 4.80 - 15.4 = 4.4 eV) into the gas. Lagged (pre-solve)
	! densities; nmol_eq(:,1) is the neutral-H2 number density. ioniz_HeI23S_H2
	! (Cool_coeff.f90) is the total, scaled here to the Penning branch. Zero
	! unless a molecular run also tracks the triplet.
	if (thereis_mol .and. thereis_HeITR) heat = heat                  &
	     + f_penning_HeI23S*nheiTR*nmol_eq(:,1)*ioniz_HeI23S_H2(T_K)   &
	       *((e_th_HeI - e_th_HeTR) - e_th_H2)/erg2eV

	! Lyman-Werner photodissociation heating: H2 + hv -> H + H leaves the
	! fragment pair with about 0.4 eV of kinetic energy (Black & Dalgarno
	! 1977, ApJS 34, 405, p. 418). The 4.48 eV bond energy is paid by the
	! absorbed photon, not by the gas, so it is NOT a thermal sink of this
	! channel. Lagged (pre-solve) H2 density, like every channel above.
	if (thereis_mol .and. F_LW_star .gt. 0.0d0)                       &
	     heat = heat + k_lw_diss*nmol_eq(:,1)*e_lw_fragment_erg

	! FUV photolysis heating: each H2O or OH dissociation leaves the
	! fragments with the excess of the absorbed photon over the bond energy,
	! the same ledger the Lyman-Werner channel above uses (the bond energy
	! is paid by the photon, not by the gas, so it is not a thermal sink).
	! The H2 + O(1D) branch keeps its 1.96 eV of electronic excitation out
	! of this sum: that energy leaves as O(1D) and is released later, in the
	! O(1D) + H2 -> OH + H reaction. THAT exothermicity is NOT deposited by
	! this network -- an omission of the same kind as the missing thermal
	! dissociation sink of R12/R14, and of the same size (a few percent of
	! the photolysis heat), recorded here rather than hidden. Lagged
	! (pre-solve) densities, like every channel above.
	heat_fuv = 0.0d0
	if (thereis_oxychem) then
		do ib = 1,n_fuv_band
			heat_fuv = heat_fuv                                          &
			  + j_h2o_fuv(:,ib)*nox_eq(:,2)                              &
			    *heat_per_water_dissociation(ib)                         &
			  + j_oh_fuv(:,ib)*nox_eq(:,1)                               &
			    *heat_per_hydroxyl_dissociation(ib)
		enddo
		heat = heat + heat_fuv
	endif

   !----------------------------------!

   ! Ionization equilibrium system solution

	if (.not.thereis_He) then ! If no helium

		! Ionization solves in each cell are independent (the count>0 warm-start uses
		! this cell's own previous-step value), so the sweep is OpenMP-parallel
		! over cells. sys_x/wa/info are threadprivate (global_parameters); only
		! the subroutine-local scratch is private. count==0 runs serial (the if
		! clause) because its first-step warm-start reads the neighbour cell.
		!$omp parallel do default(shared) schedule(dynamic,8)                  &
		!$omp   private(params, usednt) if(count > 0)
		do j = N+Ng,1-Ng,-1

			! Lazily allocate this thread's threadprivate NL scratch.
			if (.not. allocated(sys_x)) allocate(sys_x(N_eq))
			if (.not. allocated(wa))    allocate(wa(lwa))

			! Ionization equilibrium system setup: named-field cell state
			! (Inc 4). System_H reads these; params is now only the MINPACK
			! transport argument (unread by the converted system).
			ieq_cell%P_HI     = P_HI(j)
			ieq_cell%rchiiB   = rchiiB(j)
			ieq_cell%nh       = nh(j)
			ieq_cell%a_ion_HI = a_ion_HI(j)

			 ! Initial guess
			if (count.le.0) then
				if(r(j).le.(1.5))then
					sys_x(1) = r(j)-0.5
				else
					sys_x(1) = 1.0
				endif     
			else

				sys_x(1) = nhii(j)/nh(j)

			endif

		 	! Analytic-Jacobian Newton (Task 2); hybrd1 fallback inside solve_ieq.
			call solve_ieq(ion_system_H,jac_system_H,N_eq,sys_x,    &
				      params,tol,wa,lwa,usednt)

			! Extract solution profiles
			nhi(j)    = nh(j)*(1.0 - sys_x(1))
			nhii(j)   = nh(j)*sys_x(1)

		enddo
		!$omp end parallel do

		nhei   = 0.0
		nheii  = 0.0
		nheiii = 0.0
		nheiTR = 0.0

	else  ! If there is helium

		! Metal unknowns start at sys_x(4) normally, but shift to sys_x(5)
		! when the He triplet occupies sys_x(4) (merged HeITR+metals system).
		! In the molecular layout the molecular unknowns own sys_x(4..7) and
		! the metals follow at metal_row_base() (single definition, shared
		! with the merged molecular+metals residual).
		if (thereis_mol) then
			mbase = metal_row_base()
			iox   = oxygen_row_base()
		else
			mbase = 4
			if (thereis_HeITR) mbase = 5
			iox   = 0
		endif

		! OpenMP-parallel cell sweep (see the no-He branch above). The cell-by-cell
		! metal coefficients (met_*, System_HeH_metals) and charge-exchange rates
		! (cx_kc, cx_metal_base) are threadprivate, so each thread keeps its own;
		! cx_metal_base is broadcast (copyin) and toggled 4<->5 per cell. All the
		! subroutine-local scratch is private. count==0 stays serial (neighbour
		! warm-start).
		! Failed-cell counters (combined via the reduction below); reset once
		! per equilibrium sweep. n_mol_clamped counts molecular cells whose roots
		! were all outside the physical simplex; the n_ieq_* counters are the
		! atomic root validation (see the declarations above).
		n_mol_clamped = 0
		n_mol_info(:) = 0
		n_ieq_reseed = 0
		n_ieq_retry  = 0
		n_ieq_unphys = 0
		n_ieq_fail   = 0
		n_acc(:)      = 0
		acc_resmax(:) = 0.0d0
		hist_conv(:)  = 0
		hist_uncv(:)  = 0
		n_cce_attempt = 0
		n_cce_root    = 0
		n_cce_solve   = 0
		cce_seconds   = 0.0d0
		!$omp parallel do default(shared) schedule(dynamic,8) copyin(cx_metal_base) &
		!$omp   private(params, usednt, i0, top, im, meg_ntot, meg_g0, meg_g1,      &
		!$omp           meg_b0, meg_b1, meg_a1, meg_a2, meg_top,                     &
		!$omp           pbar_loc, viol, viol_best,                                  &
		!$omp           iatt, info_ieq, ok_rank, best_rank, conv_ieq, phys_ieq,      &
		!$omp           x_root_best, acc_class, acc_res, res_att,                    &
		!$omp           res_best, res_nonroot, x_nonroot, info_best, info_nonroot,   &
		!$omp           have_nonroot, have_clamp,                                    &
		!$omp           x_cce, cce_full, cce_root, n_cce_fs, res_cce,                &
		!$omp           clk_beg, clk_end, clk_rate)                                  &
		!$omp   reduction(+:n_mol_clamped,n_mol_info,n_ieq_reseed,n_ieq_retry, &
		!$omp               n_ieq_unphys,n_ieq_fail,n_acc,hist_conv,hist_uncv, &
		!$omp               n_cce_attempt,n_cce_root,n_cce_solve,cce_seconds)  &
		!$omp   reduction(max:acc_resmax) if(count > 0)
		do j = N+Ng,1-Ng,-1

			! Lazily allocate this thread's threadprivate NL scratch.
			if (.not. allocated(sys_x))   allocate(sys_x(N_eq))
			if (.not. allocated(sys_sol)) allocate(sys_sol(N_eq))
			if (.not. allocated(wa))      allocate(wa(lwa))

			! System coefficients: named-field cell state (Inc 4). The He
			! equilibrium systems read these; params is now only the MINPACK
			! transport argument (unread by the converted systems).
			ieq_cell%P_HI       = P_HI(j)
			ieq_cell%P_HeI      = P_HeI(j)
			ieq_cell%P_HeII     = P_HeII(j)
			ieq_cell%rchiiB     = rchiiB(j)
			ieq_cell%rcheiiB    = rcheiiB(j)
			ieq_cell%rcheiiiB   = rcheiiiB(j)
			ieq_cell%nh         = nh(j)
			ieq_cell%nhe        = nhe(j)
			ieq_cell%a_ion_HI   = a_ion_HI(j)
			ieq_cell%a_ion_HeI  = a_ion_HeI(j)
			ieq_cell%a_ion_HeII = a_ion_HeII(j)

			! He <-> H charge-exchange rate coefficients (Huang Table 4 group
			! B): read by he_h_cx_fvec/he_h_cx_jac in every He system. Depends
			! only on T, so evaluate once per cell here (cheap). The residual
			! routines add nothing when he_h_charge_exchange is off, so this is
			! harmless (and bit-identical) in that case.
			call he_h_cx_rates(T_K(j), ieq_cell%kcx_He0_Hp,               &
			                           ieq_cell%kcx_Hep_H0)

			! Add more if HeITR is present
			if (thereis_HeITR) then
				ieq_cell%rcheiTR = rcheiTR(j)
				ieq_cell%A31     = A31
				ieq_cell%P_HeITR = P_HeITR(j)
				ieq_cell%q13     = q13(j)
				ieq_cell%q31a    = q31a(j)
				ieq_cell%q31b    = q31b(j)
				ieq_cell%Q31     = Q31(j)
				ieq_cell%a_ion_HeITR = a_ion_HeITR(j)
			endif

			! molecular cell state (System_HeH_mol layout)
			if (thereis_mol) then
				if (.not. thereis_HeITR) then
					ieq_cell%rcheiTR = 0.0d0     ! no triplet channels
					ieq_cell%A31     = 0.0d0
					ieq_cell%P_HeITR = 0.0d0
					ieq_cell%q13     = 0.0d0
					ieq_cell%q31a    = 0.0d0
					ieq_cell%q31b    = 0.0d0
					ieq_cell%Q31     = 0.0d0
					ieq_cell%a_ion_HeITR = 0.0d0
				endif
				ieq_cell%P_H2 = P_H2(j)
				P_H2_eq(j)    = P_H2(j)
				ieq_cell%k_LW = k_lw_diss(j)   ! 0 without a LW band flux
				ieq_cell%T_K  = T_K(j)
				ieq_cell%ntot = n_tot(j)        ! M for the 3-body rates
				! Compute the molecular rate coefficients that are invariant
				! across this cell's Newton solve (they depend only on T and
				! n_tot); the residual then reads them, like set_metal_coeffs.
				! n_tot, not n_in_dim: R13 and R15 fold the third body M into
				! their coefficient, so they need the same particle density
				! the R12 term of the residual reads (see the n_tot comment
				! above).
				call set_mol_coeffs(T_K(j), n_tot(j))
				! Oxygen-chemistry cell state: the free oxygen family the
				! carrier fractions are measured against, the frozen CO
				! reservoir, and the rate coefficients and band photolysis
				! rates of the carriers.
				if (thereis_oxychem) then
					ieq_cell%n_ofam = nm_tot(j,iel_O)
					ieq_cell%n_co   = nCO_cell(j)
					call set_oxygen_coeffs(T_K(j), j_h2o_fuv(j,:),        &
					                               j_oh_fuv(j,:))
				else
					ieq_cell%n_ofam = 0.0d0
					ieq_cell%n_co   = 0.0d0
				endif
				! Turnover scale of each balance row, from the coefficients
				! just built and this cell's densities; the residual divides
				! by it so the helium and molecular blocks reach hybrd1 with
				! the same weight.
				! photo_scale = 1: the cell's own radiation field. Only
				! the continuation of constrained_chemical_equilibrium
				! ever passes anything else, and it restores this on
				! every exit.
				call set_mol_turnover_rates(ne(j), 1.0d0)
				! The oxygen rows get theirs after that reset, and it also
				! adds the oxygen terms to the H2 row's bound.
				if (thereis_oxychem)                                  &
					call set_oxygen_turnover_rates(iox, 1.0d0)
				! Carrier partitions this cell does not own. Default: it
				! owns all of them and every row is a balance.
				ieq_cell%x_h2_fixed     = .false.
				ieq_cell%x_ox_fixed     = .false.
				ieq_cell%x_h2_fix       = 0.0d0
				ieq_cell%x_oh_fix       = 0.0d0
				ieq_cell%x_h2o_fix      = 0.0d0
				! With the carriers transported, their partition is not a
				! local root any more: the transport-chemistry solve owns
				! H2, OH and H2O and this sweep is handed the answer. The
				! remaining rows are then solved against it, which is what
				! keeps the ionization stages consistent with the
				! transported composition.
				! bg_ready OR a restart: on a cold start the very first
				! sweep is what INITIALIZES the carriers, from the local
				! equilibrium of the seeded state, and imposing a seed on
				! it would be imposing nothing. On a restart the loaded
				! state already carries them and re-solving it locally would
				! throw the transported partition away -- which is what the
				! round-trip gate G6 measures.
				if (thereis_oxychem .and. oxygen_transport                &
				    .and. (bg_ready .or. do_load_IC)) then
					ieq_cell%x_h2_fixed = .true.
					ieq_cell%x_ox_fixed = .true.
					ieq_cell%x_h2_fix  = 2.0d0*nmol_eq(j,1)/nh(j)
					ieq_cell%x_oh_fix  = nox_eq(j,1)                      &
					                   /max(nm_tot(j,iel_O),1.0d-30)
					ieq_cell%x_h2o_fix = nox_eq(j,2)                      &
					                   /max(nm_tot(j,iel_O),1.0d-30)
				endif
				! THE LOWER-BOUNDARY RESERVOIR. The ghost cells below the
				! base are not a piece of atmosphere this cell's radiation
				! field gets to determine: they are the gas flowing IN, and
				! its molecular partition was set in the lower atmosphere,
				! where the photodissociating field and the mixing that the
				! base cell cannot see both act. A shielded base cell left
				! to its own balance drives H2 to the fully molecular limit
				! -- measured at 0.9997 of the H nuclei against the 0.9251
				! the handoff states -- and the equation of state and the
				! species state then describe different gas at the same
				! ghost. So the handoff value is imposed here and the
				! ionization stages are solved against it.
				!
				! Both lower ghosts are pinned: they are the same reservoir,
				! and Rec_BC builds face states from each of them.
				!
				! Only when a handoff actually stated the value. Without one
				! the partition comes from the chemical-equilibrium fit,
				! which is this cell's own local estimate rather than
				! upstream information -- imposing it would be imposing the
				! solver's own answer on the solver.
				if (base_h2_composition_imposed() .and. j .le. 0) then
					ieq_cell%x_h2_fixed = .true.
					ieq_cell%x_h2_fix   = base_h2_nuclei_fraction()
				endif
				! Keep this cell's coefficient state for the carrier
				! transport operator, which evaluates the same rows on the
				! same background between sweeps. Each iteration writes its
				! own element, so the parallel sweep needs no guard.
				bg_cell(j) = ieq_cell
			endif

			! Each element's metal coefficients are handed to
			! ion_system_HeH_metals via set_metal_coeffs; the charge-
			! exchange rate coefficients are stored for this cell by
			! cx_set_cell (used inside the residual by cx_add_to_fvec).
			if (thereis_metals) then
				call cx_set_cell(T_K(j))

				! Build the metal coefficients for each element in canonical order
				! from the 2D rate arrays (col i0 = neutral, i0+1 = +,
				! i0+2 = ++) and store them for the residual.
				do im = 1,n_melem
					i0  = melem_i0(im)
					top = melem_top(im)
					meg_top(im)  = top
					meg_ntot(im) = nm_tot(j,im)
					meg_g0(im)   = P_m(j,i0)
					meg_b0(im)   = aion_m(j,i0)
					meg_a1(im)   = rec_m(j,i0+1)
					if (top .ge. 2) then
						meg_g1(im) = P_m(j,i0+1)
						meg_b1(im) = aion_m(j,i0+1)
						meg_a2(im) = rec_m(j,i0+2)
					else
						! No second ionization stage for this element.
						meg_g1(im) = 0.0d0
						meg_b1(im) = 0.0d0
						meg_a2(im) = 0.0d0
					endif
				enddo
				call set_metal_coeffs(n_melem, meg_ntot, meg_g0, meg_g1, &
				                    meg_b0, meg_b1, meg_a1, meg_a2,     &
				                    meg_top)

				! Metal rows of the molecular system get their turnover
				! scale too, after set_mol_turnover_rates has reset the
				! array and set_metal_coeffs has filled met_*.
				if (thereis_mol)                                     &
					call set_mol_metal_turnover_rates(ne(j), 1.0d0)
			endif

			! Initial guess
			if (count .eq. 0) then
				if (j .eq. N+Ng) then
					sys_x(1) = 1.0
					sys_x(2) = 1.0
					sys_x(3) = 1.0
					if (thereis_HeITR) sys_x(4) = 0.01
					if (thereis_mol) then
						! start fully ionized aloft; molecules negligible
						sys_x(4:7) = 1.0d-10
						if (thereis_HeITR) sys_x(8) = 0.01
						! Oxygen carriers likewise negligible aloft: the
						! oxygen there is atomic and ionized.
						if (thereis_oxychem) then
							sys_x(iox)   = 1.0d-12
							sys_x(iox+1) = 1.0d-12
						endif
					endif
					if (thereis_metals) then
						! Each metal starts fully singly ionized.
						do im = 1,n_melem
							sys_x(mbase + 2*(im-1))   = 1.0  ! X+  frac
							sys_x(mbase+1 + 2*(im-1)) = 0.0  ! X++ frac
						enddo
					endif
				else
					sys_x(1) = nhii(j+1)/nh(j+1)
					sys_x(2) = nheii(j+1)/nhe(j+1)
					sys_x(3) = nheiii(j+1)/nhe(j+1)
					if (thereis_HeITR) &
						sys_x(4) = nheiTR(j+1)/nhe(j+1)
					if (thereis_mol) then
						sys_x(4) = max(2.0d0*nmol_eq(j+1,1)/nh(j+1), 1.0d-10)
						sys_x(5) = max(2.0d0*nmol_eq(j+1,2)/nh(j+1), 1.0d-12)
						sys_x(6) = max(3.0d0*nmol_eq(j+1,3)/nh(j+1), 1.0d-12)
						sys_x(7) = max(nmol_eq(j+1,4)/nh(j+1), 1.0d-14)
						if (thereis_HeITR) sys_x(8) = nheiTR(j+1)/nhe(j+1)
						if (thereis_oxychem) then
							sys_x(iox)   = max(nox_eq(j+1,1)              &
							         /max(nm_tot(j+1,iel_O),1.0d-30),     &
							                   1.0d-14)
							sys_x(iox+1) = max(nox_eq(j+1,2)              &
							         /max(nm_tot(j+1,iel_O),1.0d-30),     &
							                   1.0d-14)
						endif
					endif
					if (thereis_metals) then
						do im = 1,n_melem
							i0 = melem_i0(im)
							sys_x(mbase+2*(im-1)) = nm(j+1,i0+1)        &
							                  /max(nm_tot(j+1,im),1.0d-30)
							if (melem_top(im) .ge. 2) then
								sys_x(mbase+1+2*(im-1)) = nm(j+1,i0+2)        &
								                  /max(nm_tot(j+1,im),1.0d-30)
							else
								sys_x(mbase+1+2*(im-1)) = 0.0d0
							endif
						enddo
					endif
				endif
			else
				sys_x(1) = nhii(j)/nh(j)
				sys_x(2) = nheii(j)/nhe(j)
				sys_x(3) = nheiii(j)/nhe(j)
				if (thereis_HeITR) sys_x(4) = nheiTR(j)/nhe(j)
				if (thereis_mol) then
					sys_x(4) = 2.0d0*nmol_eq(j,1)/nh(j)
					sys_x(5) = 2.0d0*nmol_eq(j,2)/nh(j)
					sys_x(6) = 3.0d0*nmol_eq(j,3)/nh(j)
					sys_x(7) = nmol_eq(j,4)/nh(j)
					if (thereis_HeITR) sys_x(8) = nheiTR(j)/nhe(j)
					if (thereis_oxychem) then
						sys_x(iox)   = nox_eq(j,1)                        &
						             /max(nm_tot(j,iel_O),1.0d-30)
						sys_x(iox+1) = nox_eq(j,2)                        &
						             /max(nm_tot(j,iel_O),1.0d-30)
					endif
				endif
				if (thereis_metals) then
					do im = 1,n_melem
						i0 = melem_i0(im)
						sys_x(mbase+2*(im-1)) = nm(j,i0+1)              &
						                  /max(nm_tot(j,im),1.0d-30)
						if (melem_top(im) .ge. 2) then
							sys_x(mbase+1+2*(im-1)) = nm(j,i0+2)              &
							                  /max(nm_tot(j,im),1.0d-30)
						else
							sys_x(mbase+1+2*(im-1)) = 0.0d0
						endif
					enddo
				endif

				! The guess just built is this cell's previous state. If that
				! state is not physical it must not seed the solve again: an
				! unphysical root fed back as the next step's guess reproduces
				! itself indefinitely (the self-sticking negative H II / Fe I
				! state at the base). Start from the local ionization balance
				! instead, which lies inside the simplex by construction. The
				! molecular network keeps its own chemical-equilibrium retry.
				if (.not.thereis_mol) then
					if (.not.ionization_fractions_physical(sys_x,N_eq,mbase)) then
						call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
						                                    mbase,ne(j))
						n_ieq_reseed = n_ieq_reseed + 1
					endif
				endif
			endif

		 	! Analytic-Jacobian Newton (Task 2); hybrd1 fallback inside solve_ieq.
			! The He metastable-triplet systems keep the MINPACK solve (no
			! analytic Jacobian written for the triplet kinetics).
			if (thereis_mol) then
				! Molecular network, with the trace metals solved in the same
				! system when they are present (System_HeH_mol_metals: the two
				! blocks share the free electron density, which the metals
				! dominate in the shielded molecular base).
				!
				! hybrd1 is bistable in its initial guess, because the cell
				! itself is: a molecular basin (the dense, optically-thick
				! base, where the true root is strongly molecular) and an
				! atomic basin (the wind above the H2 -> H front). From a
				! zero-molecular warm start the solve fails at the base; at
				! the front cells themselves neither warm start need land in
				! the right basin. One starting point in each basin is
				! therefore tried in turn:
				!   1  the guess built above (previous state / neighbour),
				!   2  molecular basin -- the H2 dissociation equilibrium of
				!      the local (p, T) with every element in its own
				!      ionization balance at the incoming n_e,
				!   3  atomic basin -- the same balance with no H2,
				! and an iterate is accepted as a ROOT only if it is physical
				! (every stage fraction >= 0, each element's tracked stages
				! summing to at most its nuclei) AND its normalized reaction
				! residual is within ieq_res_tol. The solver exit code is
				! neither sufficient nor necessary for that judgement (see
				! ieq_res_tol above).
				!
				! A cell that finds such a root takes it (a solver-converged
				! one immediately). A cell that does not keeps the first
				! bounded iterate as a MARKED non-root under the relaxation
				! amnesty (nonroot_streak_update), and when NONE of the
				! attempts even stayed inside the simplex it clamps the one
				! lying closest onto the element budget and RECHECKS the
				! reaction residual, which the projection changes. Every
				! outcome is computed at THIS cell's state; none of them is
				! the cell's previous composition, which would make the
				! equilibrium -- and the steady residual built on it -- a
				! function of the sequence of evaluations rather than of the
				! state being evaluated.
				best_rank = 0
				viol_best = 0.0d0
				res_best  = 0.0d0
				info_best = 0
				have_nonroot = .false.
				have_clamp   = .false.
				! Gas pressure p = (n_tot + n_e) kB T [bar], the same ideal-gas
				! law the EOS uses (comp_p_from_T).
				pbar_loc  = (n_tot(j) + ne(j))*kb_erg*T_K(j)/1.0d6
				do iatt = 1,3
					if (iatt .eq. 2) then
						call dissociation_ionization_balance_at_fixed_ne(  &
						                 sys_x,N_eq,mbase,ne(j),pbar_loc,T_K(j))
						n_ieq_retry = n_ieq_retry + 1
					else if (iatt .eq. 3) then
						call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
						                                    mbase,ne(j))
						n_ieq_retry = n_ieq_retry + 1
					endif

					if (thereis_metals) then
						! Metals appended above the molecular unknowns; point
						! charge exchange at their rows for this solve.
						cx_metal_base = mbase
						call hybrd1(ion_system_HeH_mol_metals,N_eq,sys_x,  &
						            sys_sol,tol,info,wa,lwa,params)
						cx_metal_base = 4
					else
						call hybrd1(ion_system_HeH_mol,N_eq,sys_x,sys_sol, &
						            tol,info,wa,lwa,params)
					endif
					conv_ieq = (info .eq. 1)
					n_mol_info(max(0,min(5,info))) =                     &
						n_mol_info(max(0,min(5,info))) + 1

					phys_ieq = ionization_fractions_physical(sys_x,N_eq,mbase)
					if (iatt .eq. 1 .and. .not.phys_ieq)                &
						n_ieq_unphys = n_ieq_unphys + 1
					! Reaction residual of a physical iterate (section 113):
					! the histograms measure where converged roots sit and
					! where the bounded-but-unconverged iterates do.
					res_att = 0.0d0
					if (phys_ieq) then
						res_att = normalized_reaction_residual(sys_x,     &
						                              N_eq,mbase,ne(j))
						if (conv_ieq) then
							hist_conv(residual_decade(res_att)) =         &
								hist_conv(residual_decade(res_att)) + 1
						else
							hist_uncv(residual_decade(res_att)) =         &
								hist_uncv(residual_decade(res_att)) + 1
						endif
					endif
					ok_rank = 0
					if (phys_ieq .and. res_att .le. ieq_res_tol) ok_rank = 1
					if (ok_rank .eq. 1 .and. conv_ieq) ok_rank = 2
					! Strict improvement only, so attempt 1 wins any tie.
					if (ok_rank .gt. best_rank) then
						best_rank = ok_rank
						x_root_best(1:N_eq) = sys_x(1:N_eq)
						res_best  = res_att
						info_best = info
					else if (best_rank .eq. 0) then
						if (phys_ieq) then
							! A bounded iterate that is not a root: the state
							! of last resort, first one kept (the old rank-1
							! order).
							if (.not. have_nonroot) then
								have_nonroot = .true.
								x_nonroot(1:N_eq) = sys_x(1:N_eq)
								res_nonroot  = res_att
								info_nonroot = info
							endif
						else
							! Outside the simplex: remember the iterate that
							! lies closest, to be clamped onto it below.
							viol = element_budget_violation(sys_x,N_eq,mbase)
							if (.not. have_clamp .or. viol .lt. viol_best) then
								have_clamp = .true.
								viol_best  = viol
								x_root_best(1:N_eq) = sys_x(1:N_eq)
							endif
						endif
					endif
					if (ok_rank .eq. 2) exit
				enddo

				! No starting point produced a root. Before the degraded
				! paths below -- the marked non-root of last resort and the
				! projection onto the element budget -- solve the CONSTRAINED
				! problem: positive species densities as unknowns, element
				! conservation as explicit rows (so that HeH+ and the oxygen
				! carriers satisfy BOTH element budgets they belong to at
				! once, which a projection cannot do), charge neutrality
				! identically, and a continuation in the radiation field from
				! the dense molecular limit instead of three unconnected
				! guesses (docs/supersonic_molecular_base.md section 11.5-B).
				!
				! Its result is a candidate, not an acceptance: it faces the
				! SAME judge every other candidate faces, bounds and the
				! normalized reaction residual against ieq_res_tol. Reaching
				! the full field is a statement about the continuation, not
				! about the cell's network.
				!
				! Only the molecular branch is promoted. The atomic matrix
				! cases carry no measured class-3/4 population (section 113),
				! and the constrained formulation is aimed at what the
				! molecular network has and the atomic one does not: species
				! that carry two elements at once, whose two budgets a
				! fraction layout can only satisfy one at a time.
				cce_root = .false.
				res_cce  = 0.0d0
				if (best_rank .eq. 0) then
					n_cce_attempt = n_cce_attempt + 1
					x_cce(1:N_eq) = 0.0d0
					call system_clock(count=clk_beg, count_rate=clk_rate)
					call equilibrium_from_molecular_limit(x_cce,N_eq,  &
					              mbase,iox,ne(j),pbar_loc,cce_full,   &
					              n_cce_fs)
					call system_clock(count=clk_end)
					if (clk_rate .gt. 0)                               &
						cce_seconds = cce_seconds                      &
						  + dble(clk_end - clk_beg)/dble(clk_rate)
					n_cce_solve = n_cce_solve + n_cce_fs
					if (cce_full) then
						phys_ieq =                                     &
						  ionization_fractions_physical(x_cce,N_eq,mbase)
						if (.not.phys_ieq) then
							! The conservation rows are residuals, so
							! they leave at most tolerance-level slack;
							! removing it is a round-off correction of
							! the solve, and the residual is measured
							! AFTER it, which is the section 11.5-A
							! recheck rule. A larger violation is not
							! round-off and is not clamped away.
							if (element_budget_violation(x_cce,N_eq,   &
							             mbase) .le. 1.0d-6) then
								call clamp_fractions_to_element_budget(&
								             x_cce,N_eq,mbase)
								phys_ieq =                             &
								  ionization_fractions_physical(x_cce, &
								                        N_eq,mbase)
							endif
						endif
						if (phys_ieq) then
							res_cce =                                  &
							  normalized_reaction_residual(x_cce,N_eq, &
							                        mbase,ne(j))
							if (res_cce .le. ieq_res_tol)              &
								cce_root = .true.
						endif
					endif
					if (cce_root) n_cce_root = n_cce_root + 1
				endif

				! Accept the state (section 113). Classes: 1 = solver-
				! converged root, 2 = root without solver convergence,
				! 3 = projected state whose RECHECKED residual still marks a
				! root, 4 = non-root under the relaxation amnesty
				! (nonroot_streak_update reports it and stops the run if it
				! persists), 5 = root of the constrained element-conserving
				! continuation solve just above.
				if (best_rank .eq. 2) then
					acc_class = 1
					acc_res   = res_att   ! sys_x is the attempt just measured
				else if (best_rank .eq. 1) then
					sys_x(1:N_eq) = x_root_best(1:N_eq)
					acc_class = 2
					acc_res   = res_best
					info      = info_best
				else if (cce_root) then
					sys_x(1:N_eq) = x_cce(1:N_eq)
					acc_class = 5
					acc_res   = res_cce
					! info stays the exit code of the last fraction-layout
					! attempt: it is what the promotion followed, and the
					! continuation's own status is not a solver code (its
					! cost is in the end-of-run constrained-solve line).
					call report_acceptance_event(                      &
						'constrained-continuation root accepted',      &
						j,count,r(j),T_K(j),info,viol_best,acc_res)
				else if (have_nonroot) then
					sys_x(1:N_eq) = x_nonroot(1:N_eq)
					acc_class = 4
					acc_res   = res_nonroot
					info      = info_nonroot
				else
					! Every attempt left the physical simplex: project the
					! closest one onto the element budget -- the violations
					! are the solver's own resolution of a root sitting on a
					! face (a fully dissociated or fully ionized species) --
					! and recheck the reaction residual, which the nonlinear
					! projection changes.
					sys_x(1:N_eq) = x_root_best(1:N_eq)
					call clamp_fractions_to_element_budget(sys_x,N_eq,mbase)
					n_mol_clamped = n_mol_clamped + 1
					acc_res   = normalized_reaction_residual(sys_x,N_eq,   &
					                                         mbase,ne(j))
					acc_class = 3
					if (acc_res .gt. ieq_res_tol) acc_class = 4
				endif
				n_acc(acc_class)      = n_acc(acc_class) + 1
				acc_resmax(acc_class) = max(acc_resmax(acc_class),acc_res)
				call nonroot_streak_update(acc_class,j,count,r(j),T_K(j),  &
				                           n_in_dim(j),ne(j),info,         &
				                           viol_best,acc_res,sys_x,N_eq)
			else
				! Atomic H/He (+ He 2^3S) (+ metals). Acceptance as a root
				! requires the same two conditions as the molecular block:
				! a physical state (every stage fraction non-negative, the
				! ionized stages of each element summing to at most its
				! nucleus total) AND a normalized reaction residual within
				! ieq_res_tol -- a vanishing solver status is neither
				! sufficient nor necessary. The cell is solved from up to
				! three starting points and the best root is kept:
				!   1  the guess built above (previous state / neighbour),
				!   2  the uncoupled ionization balance at the incoming n_e,
				!   3  the optically thick limit, all nuclei neutral.
				! Attempt 1 is accepted as it stands whenever it converges to
				! a physical root within tolerance, so a healthy cell takes
				! exactly the same solver path as before.
				best_rank = 0
				res_best  = 0.0d0
				info_best = 0
				have_nonroot = .false.
				do iatt = 1,3
					if (iatt .eq. 2) then
						call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
						                                    mbase,ne(j))
						n_ieq_retry = n_ieq_retry + 1
					else if (iatt .eq. 3) then
						sys_x(1:N_eq) = 0.0d0
					endif

					if (thereis_HeITR .and. thereis_metals) then
						! Merged He-triplet + metals: triplet at sys_x(4),
						! metals at sys_x(5..). The cell-by-cell metal
						! coefficients (set_metal_coeffs) and charge-exchange
						! rates (cx_set_cell) were already loaded above in the
						! thereis_metals block; here we only point charge
						! exchange at the shifted metal rows for the merged
						! solve.
						cx_metal_base = 5
						call hybrd1(ion_system_HeH_TR_metals,N_eq,sys_x,  &
						            sys_sol,tol,info,wa,lwa,params)
						cx_metal_base = 4
						conv_ieq = (info .eq. 1)
						info_ieq = info
					else if (thereis_HeITR) then
						call hybrd1(ion_system_HeH_TR,N_eq,sys_x,sys_sol, &
						            tol,info,wa,lwa,params)
						conv_ieq = (info .eq. 1)
						info_ieq = info
					else if (thereis_metals) then
						call solve_ieq(ion_system_HeH_metals,             &
						            jac_system_HeH_metals,                &
						            N_eq,sys_x,params,tol,wa,lwa,usednt,  &
						            info_ieq)
						conv_ieq = (info_ieq .eq. 1)
					else
						call solve_ieq(ion_system_HeH,jac_system_HeH,     &
						            N_eq,sys_x,params,tol,wa,lwa,usednt,  &
						            info_ieq)
						conv_ieq = (info_ieq .eq. 1)
					endif

					phys_ieq = ionization_fractions_physical(sys_x,N_eq,mbase)
					if (iatt .eq. 1 .and. .not.phys_ieq)                &
						n_ieq_unphys = n_ieq_unphys + 1
					! Reaction residual of a physical iterate (section 113);
					! same convention as the molecular block above.
					res_att = 0.0d0
					if (phys_ieq) then
						res_att = normalized_reaction_residual(sys_x,     &
						                              N_eq,mbase,ne(j))
						if (conv_ieq) then
							hist_conv(residual_decade(res_att)) =         &
								hist_conv(residual_decade(res_att)) + 1
						else
							hist_uncv(residual_decade(res_att)) =         &
								hist_uncv(residual_decade(res_att)) + 1
						endif
					endif
					ok_rank = 0
					if (phys_ieq .and. res_att .le. ieq_res_tol) ok_rank = 1
					if (ok_rank .eq. 1 .and. conv_ieq) ok_rank = 2
					! Strict improvement only, so attempt 1 wins any tie and
					! an already healthy cell is untouched.
					if (ok_rank .gt. best_rank) then
						best_rank = ok_rank
						x_root_best(1:N_eq) = sys_x(1:N_eq)
						res_best  = res_att
						info_best = info_ieq
					else if (best_rank .eq. 0 .and. phys_ieq              &
					         .and. .not. have_nonroot) then
						! A bounded iterate that is not a root: the state of
						! last resort, first one kept (the old rank-1 order).
						have_nonroot = .true.
						x_nonroot(1:N_eq) = sys_x(1:N_eq)
						res_nonroot  = res_att
						info_nonroot = info_ieq
					endif
					if (ok_rank .eq. 2) exit
				enddo

				! Accept the state; same classes as the molecular block. A
				! handback state's residual is measured under the FULL
				! coupled system -- the quantity that decides whether the
				! uncoupled balance is a root of it at all.
				if (best_rank .eq. 2) then
					acc_class = 1
					acc_res   = res_att   ! sys_x is the attempt just measured
				else if (best_rank .eq. 1) then
					sys_x(1:N_eq) = x_root_best(1:N_eq)
					acc_class = 2
					acc_res   = res_best
					info_ieq  = info_best
				else if (have_nonroot) then
					sys_x(1:N_eq) = x_nonroot(1:N_eq)
					acc_class = 4
					acc_res   = res_nonroot
					info_ieq  = info_nonroot
				else
					! Every attempt left the physical simplex. Rather than
					! propagate negative densities, hand back the uncoupled
					! ionization balance, which is admissible by
					! construction, and recheck its residual under the
					! coupled system.
					call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
					                                    mbase,ne(j))
					n_ieq_fail = n_ieq_fail + 1
					acc_res   = normalized_reaction_residual(sys_x,N_eq,   &
					                                         mbase,ne(j))
					acc_class = 3
					if (acc_res .gt. ieq_res_tol) acc_class = 4
				endif
				n_acc(acc_class)      = n_acc(acc_class) + 1
				acc_resmax(acc_class) = max(acc_resmax(acc_class),acc_res)
				call nonroot_streak_update(acc_class,j,count,r(j),T_K(j),  &
				                           n_in_dim(j),ne(j),info_ieq,     &
				                           0.0d0,acc_res,sys_x,N_eq)
			endif

			! Extract solution profiles
			if (thereis_mol) then
				! guard tiny negatives from the NL solve
				sys_x(1:N_eq) = max(sys_x(1:N_eq), 0.0d0)
				nhii(j)   = nh(j)*sys_x(1)
				nmol_eq(j,1) = 0.5d0*sys_x(4)*nh(j)
				nmol_eq(j,2) = 0.5d0*sys_x(5)*nh(j)
				nmol_eq(j,3) = sys_x(6)*nh(j)/3.0d0
				nmol_eq(j,4) = sys_x(7)*nh(j)
				if (thereis_oxychem) then
					nox_eq(j,1) = sys_x(iox)  *nm_tot(j,iel_O)
					nox_eq(j,2) = sys_x(iox+1)*nm_tot(j,iel_O)
					n_o1d_eq(j) = excited_oxygen_density(nmol_eq(j,1),    &
					                                     nox_eq(j,2))
				endif
				nhi(j)    = nh(j)*max(1.0d0 - sys_x(1) - sys_x(4)      &
				              - sys_x(5) - sys_x(6) - sys_x(7), 0.0d0)
				! The H nuclei bound in OH and H2O are not atomic H.
				if (thereis_oxychem) nhi(j) = max(nhi(j)                  &
				              - nox_eq(j,1) - 2.0d0*nox_eq(j,2), 0.0d0)
				nheii(j)  = nhe(j)*sys_x(2)
				nheiii(j) = nhe(j)*sys_x(3)
				nhei(j)   = max(nhe(j)*(1.0d0 - sys_x(2) - sys_x(3))   &
				              - nmol_eq(j,4), 0.0d0)
				if (thereis_HeITR) nheiTR(j) = nhe(j)*sys_x(8)
			else
			nhi(j)    = nh(j)*(1.0 - sys_x(1))
			nhii(j)   = nh(j)*sys_x(1)
			nhei(j)   = nhe(j)*(1.0 - sys_x(2) - sys_x(3))
			nheii(j)  = nhe(j)*sys_x(2)
			nheiii(j) = nhe(j)*sys_x(3)
			if (thereis_HeITR) nheiTR(j) = nhe(j)*sys_x(4)
			endif
			if (thereis_metals) then
				do im = 1,n_melem
					i0 = melem_i0(im)
					nm(j,i0+1) = nm_tot(j,im)*sys_x(mbase+2*(im-1))
					if (melem_top(im) .ge. 2) then
						nm(j,i0+2) = nm_tot(j,im)*sys_x(mbase+1+2*(im-1))
						nm(j,i0)   = nm_tot(j,im)                            &
						           *(1.0 - sys_x(mbase+2*(im-1)) - sys_x(mbase+1+2*(im-1)))
					else
						! Two-stage element: neutral = total - singly ionized.
						nm(j,i0)   = nm_tot(j,im)*(1.0 - sys_x(mbase+2*(im-1)))
					endif
				enddo
				! O I means FREE ATOMIC neutral oxygen with the oxygen
				! chemistry on: the closure above still holds the oxygen
				! bound in OH and H2O, which the residual removed from the
				! same quantity (System_HeH_mol_metals). Removing it here as
				! well is what keeps the written column and the solved one
				! the same object.
				if (thereis_oxychem)                                       &
					nm(j,melem_i0(iel_O)) = max(nm(j,melem_i0(iel_O))      &
					        - nox_eq(j,1) - nox_eq(j,2), 0.0d0)
			endif

		enddo
		!$omp end parallel do

		! The frozen background is now a complete sweep old at worst, so the
		! carrier transport operator may run.
		if (thereis_mol) bg_ready = .true.

		! One summary line per sweep when a molecular cell's roots were all
		! outside the physical simplex, so the closest one was clamped onto the
		! element budget.
		if (thereis_mol .and. n_mol_clamped .gt. 0) then
			write(*,'(A,I0,A)') ' (ioniz_eq) WARNING: every molecular '//     &
				'equilibrium root left the physical simplex at ',              &
				n_mol_clamped, ' cell(s), clamped onto the element budget'
		endif

		! One summary line per sweep whenever a state left the physical
		! simplex: how many cells had their stored state rejected as a
		! starting point, how many first roots were outside the simplex, and
		! how many ended on the ionization balance because no starting point
		! produced an admissible root. A cell that is merely retried because
		! the solver did not reach its tolerance is routine and stays silent
		! here; it is counted in the run-wide totals reported at the end.
		if (n_ieq_reseed + n_ieq_unphys + n_ieq_fail .gt. 0) then
			write(*,'(A,I0,A,I0,A,I0,A,I0,A)')                             &
				' (ioniz_eq) step ', count,                                &
				': ionization roots - ', n_ieq_reseed,                     &
				' stored state(s) rejected, ', n_ieq_unphys,               &
				' root(s) outside the simplex, no admissible root at ',    &
				n_ieq_fail, ' cell(s)'
		endif

		ieq_n_reseed = ieq_n_reseed + n_ieq_reseed
		ieq_n_retry  = ieq_n_retry  + n_ieq_retry
		ieq_n_unphys = ieq_n_unphys + n_ieq_unphys
		ieq_n_noroot = ieq_n_noroot + n_ieq_fail
		ieq_n_mol_clamped = ieq_n_mol_clamped + n_mol_clamped
		ieq_n_mol_info(:) = ieq_n_mol_info(:) + n_mol_info(:)
		ieq_acc_n(:)      = ieq_acc_n(:) + n_acc(:)
		ieq_acc_resmax(:) = max(ieq_acc_resmax(:), acc_resmax(:))
		ieq_hist_conv(:)  = ieq_hist_conv(:) + hist_conv(:)
		ieq_hist_uncv(:)  = ieq_hist_uncv(:) + hist_uncv(:)
		ieq_n_cce_attempt = ieq_n_cce_attempt + n_cce_attempt
		ieq_n_cce_root    = ieq_n_cce_root    + n_cce_root
		ieq_n_cce_solve   = ieq_n_cce_solve   + n_cce_solve
		ieq_cce_seconds   = ieq_cce_seconds   + cce_seconds
		ieq_nonroot_streak_peak = max(ieq_nonroot_streak_peak,             &
		                              maxval(ieq_nonroot_streak))

	endif

	
	! Density with atomic numbers (nm adds the metal mass under the
	! eos_metals policy)
   if (thereis_mol) then
      if (thereis_oxychem) then
         call calc_rho(nhi,nhii,nhei,nheii,nheiii,n_io,nm,nmol_eq,nox_eq)
      else
         call calc_rho(nhi,nhii,nhei,nheii,nheiii,n_io,nm,nmol_eq)
      endif
   else
      call calc_rho(nhi,nhii,nhei,nheii,nheiii,n_io,nm)
   endif

   ! Abundancies profiles
   f_sp_io(:,1) = nhi/n_io
   f_sp_io(:,2) = nhii/n_io
   f_sp_io(:,3) = nhei/n_io
   f_sp_io(:,4) = nheii/n_io
   f_sp_io(:,5) = nheiii/n_io
   f_sp_io(:,6) = nheiTR/n_io
   ! molecular abundances
   if (thereis_mol) then
      f_sp_io(:,isp_H2)   = nmol_eq(:,1)/n_io
      f_sp_io(:,isp_H2p)  = nmol_eq(:,2)/n_io
      f_sp_io(:,isp_H3p)  = nmol_eq(:,3)/n_io
      f_sp_io(:,isp_HeHp) = nmol_eq(:,4)/n_io
   else
      f_sp_io(:,isp_H2:isp_HeHp) = 0.0d0
   endif

   ! oxygen-carrier abundances
   if (thereis_oxychem) then
      f_sp_io(:,isp_OH)  = nox_eq(:,1)/n_io
      f_sp_io(:,isp_H2O) = nox_eq(:,2)/n_io
      f_sp_io(:,isp_CO)  = nox_eq(:,3)/n_io
   else
      f_sp_io(:,isp_OH:isp_CO) = 0.0d0
   endif

   ! Metal abundances in canonical order (col mion_fsp(im) of f_sp_io).
   do im = 1,n_mion
      f_sp_io(:,mion_fsp(im)) = nm(:,im)/n_io
   enddo

   ! Adimensional number density profile
	n_io = n_io/n0

   ! Adimensional heating and cooling rates
   heat_out = heat/q0
   cool_out = cool/q0
      
   ! Adjust value of pressure boundary condition (the base electron
   ! density in units of n0; with eos_metals the metal electrons are
   ! included, consistently with calc_ne). Molecular ions are deliberately
   ! omitted as trace electron donors: the base is nearly neutral, so the
   ! molecular-ion electrons are negligible in dp_bc.
   dp_bc = (nhii(1-Ng) + nheii(1-Ng) + 2.0*nheiii(1-Ng))/n0
   if (eos_include_metals .and. thereis_metals) then
      do im = 1,n_mion
         if (mion_stage(im) .gt. 0)                                     &
            dp_bc = dp_bc + dble(mion_stage(im))*nm(1-Ng,im)/n0
      enddo
   endif

   ! ...and the heavy-particle count of the same ghost, so that the pressure
   ! boundary condition and the species state are one description of one gas.
   !
   ! ntot_bc was resolved once at startup from the base H2 mixing ratio.
   ! That is a count of nuclei corrected for H2 binding; the species state is
   ! a count of actual particles, and with the molecular network solved they
   ! are not the same number -- the ghost also carries H2+, H3+ and HeH+, and
   ! any hydrogen the imposed partition leaves atomic. Measured before this
   ! was closed, on the molecular regression cases: the startup value was 6.4%
   ! above the species one, and since the ghost pressure is (ntot_bc+dp_bc)T0
   ! at the pinned density, the base the run actually marched on sat at
   ! 1213.3 K where T0 = 1140 K had been asked for. The isothermal base
   ! boundary condition was not isothermal.
   !
   ! Recomputed through calc_ntot rather than summed here, so the particle-
   ! count policy (what counts as one particle, and whether metal nuclei are
   ! in the budget) keeps its single definition.
   !
   ! Molecular runs only. With an atomic network every nucleus is its own
   ! particle and the two counts agree identically, so recomputing would only
   ! move the value by round-off; and a passive molecular base (EOS-only,
   ! species left atomic) states its H2 binding through ntot_bc alone, so
   ! taking the count from the atomic species there would silently undo it.
   if (thereis_mol) then
      if (thereis_oxychem) then
         call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol_eq,nox_eq)
      else
         call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol_eq)
      endif
      ntot_bc = n_tot(1-Ng)/n0
   endif

	! End of subroutine
	end subroutine ioniz_eq

	!----------------------------------!

	logical function ionization_fractions_physical(x,n,mbase) result(ok)
	! Is a root of the equilibrium system a physically admissible state?
	!
	! Every unknown of these systems is the fraction of one element's nuclei
	! found in one ionization stage, so a physical state has
	!   (i)  every fraction >= 0, and
	!   (ii) for each element, the tracked ionized stages summing to <= 1,
	!        the neutral stage being the remainder 1 - sum.
	! The residuals are polynomials in the fractions and do admit roots
	! outside this simplex; such a root gives a negative neutral density and
	! is not a solution of the physical problem, however small the residual.
	! Downstream it poisons the photoionization integrals (a negative
	! absorber column gives negative photoheating) and the element budget.
	!
	! The comparisons carry a small tolerance so that a root sitting on a
	! face of the simplex (a fully neutral or fully ionized element) is not
	! rejected for round-off; it is far below the violations this test is
	! meant to catch (fractions of ~1e-6 to ~1e-1 outside the simplex).
	!
	! Widening it to the solver's own xtol = sqrt(machine epsilon) was tried
	! and measured: it lets a cell stop on a root that sits just outside a
	! face instead of retrying from another starting point, and the retried
	! root is the better one. On examples/15 (HD 209458 b, molecular, cold)
	! the wider band turns info = 0 at ||R|| = 6.105e-6 in 89 outer
	! iterations into info = 2 at 1.112e-4 in 133, so the band stays at
	! round-off.
	!
	! Layout: x(1) = H II / H, x(2) = He II / He, x(3) = He III / He,
	! x(4) = He 2^3S / He when the triplet is tracked (x(8) in the molecular
	! layout, where x(4..7) are the H nuclei bound in H2, H2+, H3+, HeH+),
	! and the metal stages X+ / X++ from x(mbase) upwards, two per element.

	integer, intent(in) :: n, mbase
	real*8,  intent(in) :: x(n)
	real*8, parameter   :: ftol = 1.0d-10
	integer :: im, ix, iy
	real*8  :: s

	ok = .false.

	! Hydrogen nuclei: H II, plus the H bound in molecules where tracked,
	! plus the H bound in the oxygen carriers (OH one nucleus, H2O two).
	! Those two are fractions of the OXYGEN family, so they are converted
	! with this cell's n_O(free)/n_H, the same way HeH+ is converted into
	! the helium budget below. Leaving them out lets a root whose free
	! atomic H is negative pass the test, and the extraction then clips it
	! to zero and loses the nuclei.
	if (x(1) .lt. -ftol) return
	s = x(1)
	if (thereis_mol) then
		do ix = 4,7
			if (x(ix) .lt. -ftol) return
			s = s + x(ix)
		enddo
	endif
	if (thereis_oxychem .and. ieq_cell%nh .gt. 0.0d0) then
		iy = oxygen_row_base()
		s = s + (x(iy) + 2.0d0*x(iy+1))*ieq_cell%n_ofam/ieq_cell%nh
	endif
	if (s .gt. 1.0d0 + ftol) return

	! Helium nuclei: He II, He III, the 2^3S metastable (which the systems
	! carry as a separate level inside the neutral stage) and, in the
	! molecular layout, the He nucleus bound in HeH+. The molecular systems
	! close the neutral He as n_He(1 - x2 - x3) - n_HeH+, so a root leaving
	! that negative is inadmissible even when x2 + x3 <= 1; x(7) is HeH+ per
	! H nucleus, hence the n_H/n_He conversion from this cell's state.
	if (thereis_He) then
		if (x(2) .lt. -ftol) return
		if (x(3) .lt. -ftol) return
		s = x(2) + x(3)
		if (thereis_HeITR) then
			ix = 4
			if (thereis_mol) ix = 8
			if (x(ix) .lt. -ftol) return
			s = s + x(ix)
		endif
		if (thereis_mol .and. ieq_cell%nhe .gt. 0.0d0)                &
			s = s + x(7)*ieq_cell%nh/ieq_cell%nhe
		if (s .gt. 1.0d0 + ftol) return
	endif

	! Each metal element separately: X+ (+ X++ for the three-stage elements).
	! With the oxygen chemistry on the oxygen element carries two more
	! fractions, n_OH and n_H2O over the free oxygen family, and the four of
	! them share one budget: free atomic O is the remainder.
	if (thereis_metals) then
		do im = 1,n_melem
			ix = mbase + 2*(im-1)
			if (x(ix) .lt. -ftol) return
			s = x(ix)
			if (melem_top(im) .ge. 2) then
				if (x(ix+1) .lt. -ftol) return
				s = s + x(ix+1)
			endif
			if (thereis_oxychem .and. im .eq. iel_O) then
				iy = oxygen_row_base()
				if (x(iy)   .lt. -ftol) return
				if (x(iy+1) .lt. -ftol) return
				s = s + x(iy) + x(iy+1)
			endif
			if (s .gt. 1.0d0 + ftol) return
		enddo
	endif

	ok = .true.

	end function ionization_fractions_physical

	!----------------------------------!

	subroutine ionization_balance_at_fixed_ne(x,n,mbase,n_e)
	! Stage fractions of each element in its OWN ionization balance at a
	! given electron density: photoionization plus electron-impact ionization
	! against radiative recombination,
	!
	!    n_k (gamma_k + beta_k n_e) = alpha_{k+1} n_e n_{k+1},
	!
	! with the couplings between elements (charge exchange, and the
	! dependence of n_e on the unknowns themselves) dropped and n_e taken
	! from the incoming state. Writing u_k = gamma_k + beta_k n_e for the
	! rate out of stage k and d_k = alpha_k n_e for the rate back into stage
	! k-1, the three-stage solution is
	!
	!    (n_0, n_1, n_2) proportional to (d_1 d_2, u_0 d_2, u_0 u_1),
	!
	! non-negative and normalized to one, so the result always lies inside
	! the physical simplex whatever state it is asked to replace. It is used
	! as a starting point for the coupled solve, and as the state of last
	! resort for a cell where no starting point produced an admissible root.
	!
	! The He 2^3S fraction is the steady state of the metastable level at
	! those populations: recombination and collisional excitation feed it,
	! radiative decay, collisional de-excitation, photoionization and Penning
	! ionization on H0 drain it. It is capped by the neutral He fraction.
	!
	! Rates come from the same cell state the residuals read (ieq_cell and
	! the met_* coefficients of System_HeH_metals), so this is the incoming
	! cell's own physics, not a generic guess. Molecules are not part of this
	! balance: it is the molecule-free limit of the state, so in the molecular
	! layout the H2/H2+/H3+/HeH+ fractions x(4..7) are left at zero and the
	! metastable sits at x(8).

	integer, intent(in)  :: n, mbase
	real*8,  intent(in)  :: n_e
	real*8,  intent(out) :: x(n)
	integer :: itr
	real*8  :: u0,u1,d1,d2,w1,w2,s
	real*8  :: xneu, n_hi_loc

	x(1:n) = 0.0d0

	! Hydrogen
	u0 = ieq_cell%P_HI + ieq_cell%a_ion_HI*n_e
	d1 = ieq_cell%rchiiB*n_e
	s  = u0 + d1
	if (s .gt. 0.0d0) x(1) = u0/s

	if (thereis_He) then
		! Helium. Both He+ recombination channels return to neutral He: the
		! metastable capture rcheiTR is a branch of the recombination, so it
		! adds to the rate back into He I when the triplet is tracked.
		u0 = ieq_cell%P_HeI + ieq_cell%a_ion_HeI*n_e
		d1 = ieq_cell%rcheiiB*n_e
		if (thereis_HeITR) d1 = d1 + ieq_cell%rcheiTR*n_e
		u1 = ieq_cell%P_HeII + ieq_cell%a_ion_HeII*n_e
		d2 = ieq_cell%rcheiiiB*n_e
		w1 = u0*d2
		w2 = u0*u1
		s  = d1*d2 + w1 + w2
		if (s .gt. 0.0d0) then
			x(2) = w1/s
			x(3) = w2/s
		endif

		if (thereis_HeITR) then
			itr = 4
			if (thereis_mol) itr = 8
			xneu     = max(1.0d0 - x(2) - x(3), 0.0d0)
			n_hi_loc = max(1.0d0 - x(1), 0.0d0)*ieq_cell%nh
			s = ieq_cell%P_HeITR + ieq_cell%A31                          &
			  + n_hi_loc*ieq_cell%Q31                                    &
			  + (ieq_cell%q31a + ieq_cell%q31b                           &
			     + ieq_cell%a_ion_HeITR)*n_e
			if (s .gt. 0.0d0) x(itr) = min(                              &
			      n_e*(x(2)*ieq_cell%rcheiTR + xneu*ieq_cell%q13)/s, xneu)
		endif
	endif

	if (thereis_metals)                                                  &
		call metal_ionization_balance_at_fixed_ne(x,n,mbase,n_e)

	end subroutine ionization_balance_at_fixed_ne

	!----------------------------------!

	real*8 function element_budget_violation(x,n,mbase) result(viol)
	! How far a root of an equilibrium system lies outside the physically
	! allowed states: the largest of the negative stage fractions and of the
	! amounts by which one element's tracked stages exceed its nuclei, both
	! measured as fractions of the element. Zero for an admissible state.
	! Used to pick, among roots that are all inadmissible, the one closest to
	! a state (see clamp_fractions_to_element_budget). Same layout as
	! ionization_fractions_physical.

	integer, intent(in) :: n, mbase
	real*8,  intent(in) :: x(n)
	integer :: im, ix, iy
	real*8  :: s

	viol = 0.0d0
	do ix = 1,n
		viol = max(viol, -x(ix))
	enddo

	s = x(1)
	if (thereis_mol) s = s + x(4) + x(5) + x(6) + x(7)
	! The H nuclei bound in OH and H2O (per H nucleus), as in
	! ionization_fractions_physical.
	if (thereis_oxychem .and. ieq_cell%nh .gt. 0.0d0) then
		iy = oxygen_row_base()
		s = s + (x(iy) + 2.0d0*x(iy+1))*ieq_cell%n_ofam/ieq_cell%nh
	endif
	viol = max(viol, s - 1.0d0)

	if (thereis_He) then
		s = x(2) + x(3)
		if (thereis_HeITR) then
			ix = 4
			if (thereis_mol) ix = 8
			s = s + x(ix)
		endif
		! The He nucleus bound in HeH+ (x(7) is per H nucleus), as in
		! ionization_fractions_physical.
		if (thereis_mol .and. ieq_cell%nhe .gt. 0.0d0)                &
			s = s + x(7)*ieq_cell%nh/ieq_cell%nhe
		viol = max(viol, s - 1.0d0)
	endif

	if (thereis_metals) then
		do im = 1,n_melem
			ix = mbase + 2*(im-1)
			s  = x(ix)
			if (melem_top(im) .ge. 2) s = s + x(ix+1)
			! The oxygen budget also carries OH and H2O.
			if (thereis_oxychem .and. im .eq. iel_O) then
				iy = oxygen_row_base()
				s  = s + x(iy) + x(iy+1)
			endif
			viol = max(viol, s - 1.0d0)
		enddo
	endif

	end function element_budget_violation

	!----------------------------------!

	double precision function normalized_reaction_residual(x,n,mbase,       &
	                                                       n_e_ref) result(res)
	! Dimensionless reaction imbalance of a candidate equilibrium state: the
	! residual vector of the SAME system the cell solve used, evaluated at x,
	! with each balance row measured against its turnover scale -- the rate
	! at which the species that row balances can be produced or destroyed in
	! this cell (the section-93 convention: every species set to the whole of
	! its element). A root has res -> 0 to the accuracy of the solve; element
	! bounds and conservation alone say nothing about it, which is why the
	! acceptance test of ioniz_eq requires BOTH
	! (docs/supersonic_molecular_base.md section 11.5).
	!
	! The molecular systems already reach hybrd1 with their rows divided by
	! the turnover scale (set_mol_turnover_rates, set_mol_metal_turnover_-
	! rates), so their residual is used as returned, widened only by the
	! metal <-> H/He charge-exchange bound that the solver scale leaves out
	! of its rows: one contribution among several for the solver's path, but
	! in the wind, where charge exchange couples a trace metal's rows to the
	! whole H reservoir, the dominant term of the row's turnover -- judging a
	! candidate against a scale without it would reject roots. The atomic
	! systems are unscaled inside the solver (scaling would change the
	! MINPACK path and the byte-identical history for no physical need), so
	! the same style of scale is built here from this cell's own
	! coefficients (ieq_cell, met_*): photo- and collisional ionization,
	! recombination, charge exchange, and the He 2^3S kinetics where the
	! triplet is tracked. n_e_ref is the cell's incoming electron density,
	! as in set_mol_turnover_rates. A row whose scale vanishes carries no
	! reaction: its residual is exactly zero and the scale is left at 1, as
	! are the dimensionless identity rows (absent elements, the pinned X++
	! of a two-stage element, a transported carrier's x - x_fix row).
	!
	! Serves the He-branch acceptance only; the H-only system (System_H) has
	! a single analytic balance with a unique root in [0,1] and no
	! acceptance ladder.
	integer, intent(in) :: n, mbase
	real*8,  intent(in) :: x(n)
	real*8,  intent(in) :: n_e_ref
	real*8  :: fv(n), srow(n), cxb(n), par(60), el_tot(12)
	real*8  :: nH, nHe, ne, cx_heh, s
	integer :: iflag, i, e, ix

	par(:)    = 0.0d0                 ! transport argument only, unread
	iflag     = 1
	srow(1:n) = 1.0d0

	if (thereis_metals) then
		el_tot(:) = 0.0d0
		el_tot(1:met_nelem) = met_ntot(1:met_nelem)
		el_tot(11) = ieq_cell%nh      ! cx_H
		el_tot(12) = ieq_cell%nhe     ! cx_He
	endif

	if (thereis_mol) then
		if (thereis_metals) then
			cx_metal_base = mbase
			call ion_system_HeH_mol_metals(n,x,fv,iflag,par)
			! Charge-exchange bound of each row (absolute), converted to
			! the dimensionless factor (s_solver + s_cx)/s_solver each
			! already-scaled row is divided by. Rows charge exchange never
			! reaches keep exactly 1.
			cxb(1:n) = 0.0d0
			call cx_add_to_turnover(cxb, el_tot)
			cx_metal_base = 4
			do i = 1,n
				srow(i) = 1.0d0 + cxb(i)*mol_inv_turnover(i)
			enddo
		else
			call ion_system_HeH_mol(n,x,fv,iflag,par)
		endif
	else
		! Atomic layouts: raw balance rows [cm^-3 s^-1] of the system the
		! solve used, then the turnover scale of each row.
		if (thereis_HeITR .and. thereis_metals) then
			cx_metal_base = 5
			call ion_system_HeH_TR_metals(n,x,fv,iflag,par)
			cx_metal_base = 4
		else if (thereis_HeITR) then
			call ion_system_HeH_TR(n,x,fv,iflag,par)
		else if (thereis_metals) then
			call ion_system_HeH_metals(n,x,fv,iflag,par)
		else
			call ion_system_HeH(n,x,fv,iflag,par)
		endif

		nH  = ieq_cell%nh
		nHe = ieq_cell%nhe
		ne  = n_e_ref
		cx_heh = (ieq_cell%kcx_He0_Hp + ieq_cell%kcx_Hep_H0)*nH*nHe
		! (1) H+ : photo- and collisional ionization of H0, radiative
		!     recombination, He <-> H charge exchange; with the triplet
		!     tracked, the Penning ionization source of heh_tr_rows.
		srow(1) = (ieq_cell%P_HI                                          &
		           + (ieq_cell%a_ion_HI + ieq_cell%rchiiB)*ne)*nH + cx_heh
		if (thereis_HeITR) then
			srow(1) = srow(1) + f_penning_HeI23S*ieq_cell%Q31*nHe*nH
			! (2) summed He I balance of heh_tr_rows: both photoionization
			!     channels, collisional ionization of ground and metastable
			!     He I, and both recombination paths back into He I.
			srow(2) = (ieq_cell%P_HeI + ieq_cell%P_HeITR                  &
			           + (ieq_cell%a_ion_HeI + ieq_cell%a_ion_HeITR       &
			              + ieq_cell%rcheiTR + ieq_cell%rcheiiB)*ne)*nHe  &
			        + cx_heh
			! (3) He+ <-> He++.
			srow(3) = (ieq_cell%P_HeII                                    &
			           + (ieq_cell%a_ion_HeII + ieq_cell%rcheiiiB)*ne)*nHe
			! (4) He 2^3S: populated from He+ recombination and 1^1S
			!     excitation, drained by photoionization, A31,
			!     de-excitation, electron-impact and Penning ionization.
			srow(4) = (ieq_cell%P_HeITR + ieq_cell%A31                    &
			           + (ieq_cell%rcheiTR + ieq_cell%q13 + ieq_cell%q31a &
			              + ieq_cell%q31b + ieq_cell%a_ion_HeITR)*ne)*nHe &
			        + ieq_cell%Q31*nHe*nH
		else
			! (2)(3) the standard heh_rows balances.
			srow(2) = (ieq_cell%P_HeI                                     &
			           + (ieq_cell%a_ion_HeI + ieq_cell%rcheiiB)*ne)*nHe  &
			        + cx_heh
			srow(3) = (ieq_cell%P_HeII                                    &
			           + (ieq_cell%a_ion_HeII + ieq_cell%rcheiiiB)*ne)*nHe
		endif
		if (thereis_metals) then
			! X0 <-> X+ and X+ <-> X++ of each element (the metal_rows
			! balances), plus the charge-exchange bound on every row it
			! reaches; identity rows are then forced back to their
			! dimensionless scale of 1.
			do e = 1,met_nelem
				ix = mbase + 2*(e-1)
				if (met_ntot(e) .le. 1.0d-30) cycle
				srow(ix) = met_ntot(e)*(met_g0(e)                         &
				           + (met_b0(e) + met_a1(e))*ne)
				if (met_top(e) .ge. 2)                                    &
					srow(ix+1) = met_ntot(e)*(met_g1(e)                   &
					             + (met_b1(e) + met_a2(e))*ne)
			enddo
			cx_metal_base = mbase
			call cx_add_to_turnover(srow, el_tot)
			cx_metal_base = 4
			do e = 1,met_nelem
				ix = mbase + 2*(e-1)
				if (met_ntot(e) .le. 1.0d-30) then
					srow(ix)   = 1.0d0
					srow(ix+1) = 1.0d0
				else if (met_top(e) .lt. 2) then
					srow(ix+1) = 1.0d0
				endif
			enddo
		endif
	endif

	res = 0.0d0
	do i = 1,n
		s = srow(i)
		if (s .le. 0.0d0) s = 1.0d0
		res = max(res, abs(fv(i))/s)
	enddo

	end function normalized_reaction_residual

	!----------------------------------!

	integer function residual_decade(res) result(ib)
	! Decade bin of a normalized reaction residual for the run-wide
	! histograms: bin 0 collects everything at or below 1e-16 (exact zeros
	! included), bin 15 everything at or above 1e-1, and bin k in between
	! covers [10^(k-16), 10^(k-15)).
	real*8, intent(in) :: res
	if (res .le. 1.0d-16) then
		ib = 0
	else
		ib = min(15, 16 + int(floor(log10(res))))
	endif
	end function residual_decade

	!----------------------------------!

	subroutine nonroot_streak_update(acc_class,j,step,radius,T,dens,n_e,   &
	                                 info,viol,res,x,n)
	! Relaxation amnesty for a NON-ROOT acceptance, and its limit. A cell
	! whose accepted state is a root (classes 1-3, and class 5, the root of
	! the constrained element-conserving continuation solve) resets its
	! streak; a cell
	! accepting a non-root (class 4) is reported loudly, counted, and
	! allowed to continue for at most ieq_nonroot_streak_stop consecutive
	! sweeps -- the measured signature of a cold-start transient that the
	! next sweeps repair (see ieq_res_tol / ieq_nonroot_streak_stop above).
	! One sweep beyond that the same cell resting on a non-root is a
	! solution built on a state that does not solve the reaction network:
	! continuing would produce a physically meaningless wind
	! (docs/supersonic_molecular_base.md section 11.5-A), so the run stops
	! here with the full cell diagnostics.
	integer, intent(in) :: acc_class, j, step, info, n
	real*8,  intent(in) :: radius, T, dens, n_e, viol, res
	real*8,  intent(in) :: x(n)

	! Class 4 alone is the non-root acceptance; every other class, class 5
	! included, is a root and clears the streak.
	if (acc_class .ne. 4) then
		ieq_nonroot_streak(j) = 0
		return
	endif

	ieq_nonroot_streak(j) = ieq_nonroot_streak(j) + 1
	call report_acceptance_event('NON-ROOT accepted (relaxation amnesty)', &
	                             j,step,radius,T,info,viol,res)

	if (ieq_nonroot_streak(j) .ge. ieq_nonroot_streak_stop) then
		!$omp critical (ieq_acc_report)
		write(*,'(A)') ' (ioniz_eq) STOP: a cell has rested on a '//      &
			'NON-ROOT chemical equilibrium beyond the relaxation amnesty'
		write(*,'(A,I0,A,I0,A,I0,A)') '   cell ', j, '  step ', step,     &
			'  consecutive non-root sweeps ', ieq_nonroot_streak(j), ''
		write(*,'(A,F0.6,A,ES11.4,A,ES11.4,A,ES11.4)')                    &
			'   r [R_p] ', radius, '  T [K] ', T,                         &
			'  n_tot [cm^-3] ', dens, '  n_e [cm^-3] ', n_e
		write(*,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3)')                      &
			'   solver info ', info, '  element violation ', viol,        &
			'  normalized reaction residual ', res,                       &
			'  (tolerance ', ieq_res_tol, ')'
		write(*,'(A)') '   candidate stage fractions:'
		write(*,'(6ES12.4)') x(1:n)
		write(*,'(A)') '   Bounds and element conservation are necessary'//&
			' conditions, not reaction equilibrium; continuing with a'
		write(*,'(A)') '   known non-root produces a physically'//        &
			' meaningless wind. See docs/Update_EXHALE.md section 113.'
		flush(6)
		!$omp end critical (ieq_acc_report)
		error stop 'ioniz_eq: persistent non-root chemical equilibrium'
	endif

	end subroutine nonroot_streak_update

	!----------------------------------!

	subroutine report_acceptance_event(label,j,step,radius,T,info,viol,res)
	! One line per NON-ROOT acceptance (nonroot_streak_update), with the
	! cell, step, radius, temperature, solver exit code, element-budget
	! violation and normalized reaction residual -- a non-root acceptance
	! is never silent. Capped at ieq_acc_nprint_max lines per run so a
	! pathological run cannot flood the log; the run-wide counters keep
	! the full population.
	character(len=*), intent(in) :: label
	integer, intent(in) :: j, step, info
	real*8,  intent(in) :: radius, T, viol, res

	!$omp critical (ieq_acc_report)
	if (ieq_acc_nprint .lt. ieq_acc_nprint_max) then
		ieq_acc_nprint = ieq_acc_nprint + 1
		write(*,'(A,A,A,I0,A,I0,A,F0.6,A,ES10.3,A,I0,A,ES10.3,A,ES10.3)') &
			' (ioniz_eq) ', label, ': cell ', j, ' step ', step,          &
			' r ', radius, ' T[K] ', T, ' info ', info,                   &
			' viol ', viol, ' res ', res
		if (ieq_acc_nprint .eq. ieq_acc_nprint_max)                       &
			write(*,'(A)') ' (ioniz_eq) further acceptance report '//     &
				'lines suppressed (run-wide counters keep counting)'
	endif
	!$omp end critical (ieq_acc_report)
	end subroutine report_acceptance_event

	!----------------------------------!

	subroutine clamp_fractions_to_element_budget(x,n,mbase)
	! Move a root that lies just outside the physically allowed states onto the
	! nearest allowed one: no negative populations, and no element with more
	! nuclei in its tracked stages than it has. Negative fractions are set to
	! zero and, where an element's tracked stages still sum above its nuclei,
	! they are scaled down to sum exactly to them, leaving the neutral stage
	! empty.
	!
	! This is the treatment of a root the solver has resolved to a face of the
	! allowed region -- a fully dissociated H2, a fully neutral or fully
	! ionized element -- and delivered as a small number of either sign. The
	! clamped state is that root to the accuracy the solve reached, so it
	! remains a solution of the equilibrium at THIS cell's state; what it must
	! not be replaced by is a composition carried over from an earlier
	! evaluation.

	integer, intent(in)    :: n, mbase
	real*8,  intent(inout) :: x(n)
	integer :: im, ix, iy
	real*8  :: s, hehp_he, ox_h

	do ix = 1,n
		if (x(ix) .lt. 0.0d0) x(ix) = 0.0d0
	enddo

	! Each metal element separately.  The oxygen element carries its two
	! molecular carriers in the same budget and is scaled with them, so that
	! the clamped state still has a non-negative free atomic oxygen. This
	! block runs FIRST because the oxygen carriers also hold H nuclei, and
	! the hydrogen budget below has to see them already capped. It is
	! disjoint from the hydrogen and helium unknowns, so an oxygen-chemistry
	! -free run is unaffected by the reordering.
	if (thereis_metals) then
		do im = 1,n_melem
			ix = mbase + 2*(im-1)
			s  = x(ix)
			if (melem_top(im) .ge. 2) s = s + x(ix+1)
			iy = 0
			if (thereis_oxychem .and. im .eq. iel_O) then
				iy = oxygen_row_base()
				s  = s + x(iy) + x(iy+1)
			endif
			if (s .gt. 1.0d0) then
				x(ix) = x(ix)/s
				if (melem_top(im) .ge. 2) x(ix+1) = x(ix+1)/s
				if (iy .gt. 0) then
					x(iy)   = x(iy)/s
					x(iy+1) = x(iy+1)/s
				endif
			endif
		enddo
	endif

	! Hydrogen nuclei: H II, the H bound in molecules where tracked, and the
	! H bound in the oxygen carriers. ox_h is that last part, expressed per
	! H nucleus; it has just been capped against the oxygen budget, so it is
	! taken as given here and the H species are scaled into what it leaves,
	! exactly as HeH+ is treated in the helium budget below.
	ox_h = 0.0d0
	if (thereis_oxychem .and. ieq_cell%nh .gt. 0.0d0) then
		iy   = oxygen_row_base()
		ox_h = (x(iy) + 2.0d0*x(iy+1))*ieq_cell%n_ofam/ieq_cell%nh
		if (ox_h .gt. 1.0d0) then
			x(iy)   = x(iy)/ox_h
			x(iy+1) = x(iy+1)/ox_h
			ox_h    = 1.0d0
		endif
	endif
	s = x(1)
	if (thereis_mol) s = s + x(4) + x(5) + x(6) + x(7)
	if (s .gt. 1.0d0 - ox_h) then
		x(1) = x(1)*(1.0d0 - ox_h)/s
		if (thereis_mol) then
			do ix = 4,7
				x(ix) = x(ix)*(1.0d0 - ox_h)/s
			enddo
		endif
	endif

	! Helium nuclei: He II, He III, the 2^3S metastable and, in the molecular
	! layout, the He nucleus bound in HeH+. x(7) is HeH+ per H nucleus and
	! belongs to BOTH element budgets; it has just been capped against the H
	! nuclei, so it is capped against the He nuclei here as well (which can
	! only relax the H budget), and the free-He stages are scaled into what
	! is left. hehp_he is exactly zero without molecules, where the scaling
	! reduces to the plain x/s of every non-molecular run.
	if (thereis_He) then
		hehp_he = 0.0d0
		if (thereis_mol .and. ieq_cell%nhe .gt. 0.0d0) then
			hehp_he = x(7)*ieq_cell%nh/ieq_cell%nhe
			if (hehp_he .gt. 1.0d0) then
				x(7)    = x(7)/hehp_he
				hehp_he = 1.0d0
			endif
		endif
		s = x(2) + x(3)
		ix = 0
		if (thereis_HeITR) then
			ix = 4
			if (thereis_mol) ix = 8
			s = s + x(ix)
		endif
		if (s .gt. 1.0d0 - hehp_he) then
			x(2) = x(2)*(1.0d0 - hehp_he)/s
			x(3) = x(3)*(1.0d0 - hehp_he)/s
			if (ix .gt. 0) x(ix) = x(ix)*(1.0d0 - hehp_he)/s
		endif
	endif

	end subroutine clamp_fractions_to_element_budget

	!----------------------------------!

	subroutine dissociation_ionization_balance_at_fixed_ne(x,n,mbase,n_e,   &
	                                                       p_bar,T_gas)
	! H2 dissociation equilibrium at the local (p, T), with every element in
	! its own ionization balance at a given electron density: the limit of the
	! molecular network in which the couplings between the molecular ions and
	! the ionization of H and He are dropped.
	!
	!   - The H nuclei are partitioned between H2 and atomic H by the
	!     chemical-equilibrium mixing ratio q_H2(p, T) of Koskinen et al.
	!     (2022) Eq. 11, the same fit the molecular base boundary condition
	!     uses; per H nucleus the fraction bound in H2 is
	!     2 q (1 + He/H)/(1 + q) (lower_column::mu_mixture).
	!   - The atomic remainder is ionized in the H ionization balance of
	!     ionization_balance_at_fixed_ne, which also sets He, the 2^3S
	!     metastable and the metal stages.
	!   - H2+, H3+ and HeH+ are left at zero: each is a trace intermediate
	!     whose abundance is set by the very couplings this limit drops, and
	!     each holds far fewer H nuclei than H2 wherever the gas is molecular.
	!
	! Every fraction returned is non-negative and each element's stages sum to
	! at most its nuclei, so the starting point is a physically allowed state,
	! and it depends only on the local (p, T, n_e) and on the rate
	! coefficients of the cell -- not on the cell's previous composition.
	!
	! Molecular layout only: x(4) is the H-nucleus fraction bound in H2 there,
	! and the He 2^3S metastable sits at x(8).

	integer, intent(in)  :: n, mbase
	real*8,  intent(in)  :: n_e, p_bar, T_gas
	real*8,  intent(out) :: x(n)
	integer :: iy, ix
	real*8 :: qh2, x_h2, nh2s, nhis, f_oh, f_h2o, s_free

	call ionization_balance_at_fixed_ne(x,n,mbase,n_e)

	qh2  = q_h2_equilibrium(p_bar, T_gas)
	x_h2 = 2.0d0*qh2*(1.0d0 + HeH)/(1.0d0 + qh2)
	if (x_h2 .gt. 1.0d0) x_h2 = 1.0d0
	x(4) = x_h2
	! Only the H nuclei left atomic are available to ionize.
	x(1) = x(1)*(1.0d0 - x_h2)

	! Oxygen carriers at the chemical equilibrium of the same partition:
	! O + H2 <-> OH + H and OH + H2 <-> H2O + H at this cell's H2/H ratio.
	! Zero is a valid root of the water cycle and hybrd1 is bistable from a
	! zero molecular seed, so the molecular basin has to be seeded with a
	! molecular oxygen partition as well, not only a molecular hydrogen one.
	! The free atomic oxygen then closes the oxygen budget, and the metal
	! stage fractions ionization_balance_at_fixed_ne wrote are scaled into
	! what is left of the element.
	if (thereis_oxychem) then
		iy   = oxygen_row_base()
		nh2s = 0.5d0*x_h2*ieq_cell%nh
		nhis = max(1.0d0 - x_h2, 0.0d0)*(1.0d0 - x(1))*ieq_cell%nh
		call oxygen_chemical_equilibrium_fractions(T_gas, nh2s, nhis,     &
		                                           f_oh, f_h2o)
		x(iy)   = f_oh
		x(iy+1) = f_h2o
		s_free  = max(1.0d0 - f_oh - f_h2o, 0.0d0)
		ix = mbase + 2*(iel_O-1)
		x(ix)   = x(ix)*s_free
		x(ix+1) = x(ix+1)*s_free
	endif

	end subroutine dissociation_ionization_balance_at_fixed_ne

	!----------------------------------!

	subroutine metal_ionization_balance_at_fixed_ne(x,n,mbase,n_e)
	! Metal stage fractions in each element's OWN ionization balance at a
	! given electron density, element by element in canonical order (the
	! metal part of ionization_balance_at_fixed_ne above; see there for the
	! balance itself). Written separately because the molecular retry needs
	! exactly this part: its H/He/molecular unknowns come from the
	! chemical-equilibrium H2 fit, its metal unknowns from here. Only
	! x(mbase..) is touched, so the caller's other seeds are preserved. An
	! element that is absent keeps its stages at zero, as its residual rows do.

	use System_HeH_metals, only: met_nelem, met_ntot, met_g0, met_g1,     &
	                             met_b0, met_b1, met_a1, met_a2, met_top

	integer, intent(in)    :: n, mbase
	real*8,  intent(in)    :: n_e
	real*8,  intent(inout) :: x(n)
	integer :: im, ix, top
	real*8  :: u0,u1,d1,d2,w1,w2,s

	do im = 1,met_nelem
		ix = mbase + 2*(im-1)
		x(ix)   = 0.0d0
		x(ix+1) = 0.0d0
		if (met_ntot(im) .le. 1.0d-30) cycle
		top = met_top(im)
		u0  = met_g0(im) + met_b0(im)*n_e
		d1  = met_a1(im)*n_e
		if (top .ge. 2) then
			u1 = met_g1(im) + met_b1(im)*n_e
			d2 = met_a2(im)*n_e
			w1 = u0*d2
			w2 = u0*u1
			s  = d1*d2 + w1 + w2
			if (s .gt. 0.0d0) then
				x(ix)   = w1/s
				x(ix+1) = w2/s
			endif
		else
			s = d1 + u0
			if (s .gt. 0.0d0) x(ix) = u0/s
		endif
	enddo

	end subroutine metal_ionization_balance_at_fixed_ne

	! End of module
	end module ionization_equilibrium
