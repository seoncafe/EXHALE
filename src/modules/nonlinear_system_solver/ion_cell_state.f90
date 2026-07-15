	module ion_cell_state
	! Named-field cell state for the equilibrium ionization residuals. It
	! replaces reading the position-packed params() inside the residuals
	! (§5.2 Inc 4, docs/refactor_plan_system_composition_parser.md): params
	! stays as the MINPACK transport argument (hybrd1/fdjac1/solve_ieq), but
	! the converted equilibrium systems no longer read it -- they read the
	! named fields below instead. The field names follow the physical
	! quantities of the equilibrium params layout (params_idx IPE_*/IPH_*).
	!
	! The ionization cell sweep (ioniz_eq) is OpenMP-parallel over cells, so
	! the single module instance is threadprivate: each thread fills its own
	! copy field by field, exactly like the threadprivate met_* arrays in
	! System_HeH_metals. All components are fixed-size scalars (no allocatable
	! components), so the threadprivate copy needs no per-thread allocation.

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

	end module ion_cell_state
