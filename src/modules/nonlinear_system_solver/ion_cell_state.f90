	module ion_cell_state
	! Named-field cell state for the equilibrium ionization residuals. It
	! replaces reading the position-packed params() inside the residuals
	! (§5.2 Inc 4, docs/refactor_plan_system_composition_parser.md): params
	! stays as the MINPACK transport argument (hybrd1/fdjac1/solve_ieq), but
	! the converted equilibrium systems no longer read it -- they read the
	! named fields below instead. The field names follow the physical
	! quantities of the equilibrium params layout.
	!
	! The ionization cell sweep (ioniz_eq) is OpenMP-parallel over cells, so
	! the single module instance is threadprivate: each thread fills its own
	! copy field by field, exactly like the threadprivate met_* arrays in
	! System_HeH_metals. All components are fixed-size scalars (no allocatable
	! components), so the threadprivate copy needs no per-thread allocation.
	!
	! The advection-correction residuals (System_implicit_adv_H/HeH/HeH_TR) and
	! the post-process temperature residual (T_equation / solve_T_brent) read
	! their own named-field types (adv_rates/teq_state) below. Their loop in
	! post_process_adv is serial (master thread), but the instances are kept
	! threadprivate for uniformity with ieq_cell.

	implicit none
	public

	type ion_rates
		real*8 :: P_HI
		real*8 :: P_HeI
		real*8 :: P_HeII
		real*8 :: rchiiB
		real*8 :: rcheiiB
		real*8 :: rcheiiiB
		real*8 :: nh
		real*8 :: nhe
		real*8 :: a_ion_HI
		real*8 :: a_ion_HeI
		real*8 :: a_ion_HeII
		real*8 :: rcheiTR
		real*8 :: A31
		real*8 :: P_HeITR
		real*8 :: q13
		real*8 :: q31a
		real*8 :: q31b
		real*8 :: Q31
		real*8 :: P_H2
		real*8 :: T_K
		real*8 :: ntot
	end type ion_rates

	type(ion_rates), save :: ieq_cell
	!$omp threadprivate(ieq_cell)

	! Named-field cell state for the advection-correction ionization residuals
	! (System_implicit_adv_H/HeH/HeH_TR). One type covers BOTH advection
	! layouts: the H-only and HeH systems read the leading fields, the HeH_TR
	! system additionally reads the triplet channels (rcheiTR..xheiTR_old).
	! heh_loc is the effective He/H for the electron density (global HeH
	! normally, the local nhe/nh with He_diffusion). Filled field by field by
	! post_process_adv before each hybrd1 advection solve.
	type adv_rates
		real*8 :: c1
		real*8 :: xhi_old
		real*8 :: xhei_old
		real*8 :: xheiii_old
		real*8 :: nh
		real*8 :: P_HI
		real*8 :: P_HeI
		real*8 :: P_HeII
		real*8 :: rchiiB
		real*8 :: rcheiiB
		real*8 :: rcheiiiB
		real*8 :: a_ion_HI
		real*8 :: a_ion_HeI
		real*8 :: a_ion_HeII
		real*8 :: heh_loc
		real*8 :: rcheiTR
		real*8 :: A31
		real*8 :: P_HeITR
		real*8 :: q13
		real*8 :: q31a
		real*8 :: q31b
		real*8 :: Q31
		real*8 :: xheiTR_old
	end type adv_rates

	type(adv_rates), save :: adv_cell
	!$omp threadprivate(adv_cell)

	! Named-field cell state for the post-process temperature residual
	! (T_equation / solve_T_brent). Filled by post_process_adv before each
	! cell temperature solve.
	type teq_state
		real*8 :: nhi
		real*8 :: nhii
		real*8 :: nhei
		real*8 :: nheii
		real*8 :: nheiii
		real*8 :: mup
		real*8 :: mum
		real*8 :: rhov
		real*8 :: coeff
		real*8 :: dr
		real*8 :: Told
		real*8 :: heaold
	end type teq_state

	type(teq_state), save :: teq_cell
	!$omp threadprivate(teq_cell)

	end module ion_cell_state
