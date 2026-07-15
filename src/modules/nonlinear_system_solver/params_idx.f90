	module params_idx
	! Named indices into the position-packed params() arrays used by the
	! ionization and temperature nonlinear systems. Inc 0 of
	! docs/refactor_plan_system_composition_parser.md: a pure renaming of the
	! magic params(N) slots into type-checkable named constants, applied at
	! BOTH the packing sites (ionization_equilibrium.f90, post_process_adv.f90)
	! and the System_* unpacking. The slot meanings differ between the
	! equilibrium and advection families, so each family carries its own
	! prefixed set below.

	implicit none
	public

	! --- Group 1: equilibrium HeH family (prefix IPE_) --------------------
	! Used by System_HeH, System_HeH_TR, System_HeH_metals,
	! System_HeH_TR_metals, System_HeH_mol; packed in
	! ionization_equilibrium.f90 and the metals_pp slots (7,8) of
	! post_process_adv.f90.
	integer, parameter :: IPE_PHI    = 1     ! P_HI
	integer, parameter :: IPE_PHEI   = 2     ! P_HeI
	integer, parameter :: IPE_PHEII  = 3     ! P_HeII
	integer, parameter :: IPE_AHII   = 4     ! rchiiB
	integer, parameter :: IPE_AHEII  = 5     ! rcheiiB
	integer, parameter :: IPE_AHEIII = 6     ! rcheiiiB
	integer, parameter :: IPE_NH     = 7     ! nh
	integer, parameter :: IPE_NHE    = 8     ! nhe
	integer, parameter :: IPE_BHI    = 9     ! a_ion_HI
	integer, parameter :: IPE_BHEI   = 10    ! a_ion_HeI
	integer, parameter :: IPE_BHEII  = 11    ! a_ion_HeII
	integer, parameter :: IPE_ATR    = 12    ! rcheiTR
	integer, parameter :: IPE_A31    = 13    ! A31
	integer, parameter :: IPE_PTR    = 14    ! P_HeITR
	integer, parameter :: IPE_Q13    = 15    ! q13
	integer, parameter :: IPE_Q31A   = 16    ! q31a
	integer, parameter :: IPE_Q31B   = 17    ! q31b
	integer, parameter :: IPE_Q31    = 18    ! Q31
	integer, parameter :: IPE_PH2    = 19    ! P_H2
	integer, parameter :: IPE_T      = 20    ! T_K
	integer, parameter :: IPE_NTOT   = 21    ! n_tot (M for the 3-body rates)

	! --- Group 2: equilibrium H-only (prefix IPH_) ------------------------
	! Used by System_H; packed in ionization_equilibrium.f90.
	integer, parameter :: IPH_PHI  = 1       ! P_HI
	integer, parameter :: IPH_AHII = 2       ! rchiiB
	integer, parameter :: IPH_NH   = 3       ! nh
	integer, parameter :: IPH_BHI  = 4       ! a_ion_HI

	! --- Group 3: advection H-only (prefix IPAH_) -------------------------
	! Used by System_implicit_adv_H; packed in post_process_adv.f90.
	integer, parameter :: IPAH_C1   = 1      ! dr/v
	integer, parameter :: IPAH_XHI  = 2      ! nhi/nh
	integer, parameter :: IPAH_NH   = 3      ! nh
	integer, parameter :: IPAH_PHI  = 4      ! P_HI
	integer, parameter :: IPAH_AHII = 5      ! rchiiB
	integer, parameter :: IPAH_BHI  = 6      ! a_ion_HI

	! --- Group 4: advection HeH core (prefix IPA_) ------------------------
	! Used by System_implicit_adv_HeH (all 15) and the shared 1-14 block of
	! System_implicit_adv_HeH_TR; packed in post_process_adv.f90.
	integer, parameter :: IPA_C1     = 1     ! dr/v
	integer, parameter :: IPA_XHI    = 2     ! nhi/nh
	integer, parameter :: IPA_XHEI   = 3     ! nhei/nhe
	integer, parameter :: IPA_XHEIII = 4     ! nheiii/nhe
	integer, parameter :: IPA_NH     = 5     ! nh
	integer, parameter :: IPA_PHI    = 6     ! P_HI
	integer, parameter :: IPA_PHEI   = 7     ! P_HeI
	integer, parameter :: IPA_PHEII  = 8     ! P_HeII
	integer, parameter :: IPA_AHII   = 9     ! rchiiB
	integer, parameter :: IPA_AHEII  = 10    ! rcheiiB
	integer, parameter :: IPA_AHEIII = 11    ! rcheiiiB
	integer, parameter :: IPA_BHI    = 12    ! a_ion_HI
	integer, parameter :: IPA_BHEI   = 13    ! a_ion_HeI
	integer, parameter :: IPA_BHEII  = 14    ! a_ion_HeII
	integer, parameter :: IPA_HEH    = 15    ! local He/H (non-TR adv system)

	! --- Group 5: advection HeH_TR tail (prefix IPAT_) --------------------
	! Used by System_implicit_adv_HeH_TR (its slots 15-23); packed in
	! post_process_adv.f90. NOTE: the same physical quantity, the local He/H
	! ratio, sits at IPA_HEH = 15 in the non-TR advection system but at
	! IPAT_HEH = 23 in the TR advection system, because the triplet channels
	! occupy slots 15-22 there. This displaced-slot pair is the exact
	! fragility that Inc 0 exists to name.
	integer, parameter :: IPAT_ATR   = 15    ! rcheiTR
	integer, parameter :: IPAT_A31   = 16    ! A31
	integer, parameter :: IPAT_PTR   = 17    ! P_HeITR
	integer, parameter :: IPAT_Q13   = 18    ! q13
	integer, parameter :: IPAT_Q31A  = 19    ! q31a
	integer, parameter :: IPAT_Q31B  = 20    ! q31b
	integer, parameter :: IPAT_Q31   = 21    ! Q31
	integer, parameter :: IPAT_XTR   = 22    ! xheiTR_old
	integer, parameter :: IPAT_HEH   = 23    ! local He/H (TR adv system)

	! --- Group 6: temperature equation (prefix IPT_) ----------------------
	! Used by T_equation (equation_T); packed as paramsT(...) in
	! post_process_adv.f90.
	integer, parameter :: IPT_NHI    = 1     ! nhi
	integer, parameter :: IPT_NHII   = 2     ! nhii
	integer, parameter :: IPT_NHEI   = 3     ! nhei
	integer, parameter :: IPT_NHEII  = 4     ! nheii
	integer, parameter :: IPT_NHEIII = 5     ! nheiii
	integer, parameter :: IPT_MUP    = 6     ! mup
	integer, parameter :: IPT_MUM    = 7     ! mum
	integer, parameter :: IPT_RHOV   = 8     ! rho*v
	integer, parameter :: IPT_COEFF  = 9     ! coeff
	integer, parameter :: IPT_DR     = 10    ! dr
	integer, parameter :: IPT_TOLD   = 11    ! T_old
	integer, parameter :: IPT_HEAOLD = 12    ! heat_old

	! End of module
	end module params_idx
