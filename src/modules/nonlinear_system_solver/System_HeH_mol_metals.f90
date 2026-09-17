	module System_HeH_mol_metals
	! Coupled ionization-equilibrium system for H and He with the molecular
	! network (H2, H2+, H3+, HeH+), optionally the He 2^3S metastable, AND an
	! arbitrary number of trace metals.
	!
	! This is the merge of System_HeH_mol (H/He + molecules + optional
	! triplet) and System_HeH_metals (H/He + metals). It follows the merge
	! already made for the triplet and the metals (System_HeH_TR_metals): the
	! molecular unknowns keep their rows and the metals are appended above
	! them, so both blocks are written once and shared verbatim.
	!
	!   x(1) = n_HII  /n_H(nuclei)
	!   x(2) = n_HeII /n_He        x(3) = n_HeIII/n_He
	!   x(4) = 2 n_H2 /n_H         x(5) = 2 n_H2+/n_H
	!   x(6) = 3 n_H3+/n_H         x(7) = n_HeH+ /n_H
	!   x(8) = n_HeITR/n_He        (only when thereis_HeITR)
	! with the oxygen chemistry on (thereis_oxychem), at iox = oxygen_row_base(),
	!   x(iox)   = n_OH /n_O(free family)
	!   x(iox+1) = n_H2O/n_O(free family)
	! and, for each metal element e = 1..met_nelem (canonical order),
	!   x(mbase   + 2*(e-1)) = n_(Xe)II /n_Xe
	!   x(mbase+1 + 2*(e-1)) = n_(Xe)III/n_Xe    (pinned 0 if met_top(e) < 2)
	! with mbase = metal_row_base() = 8, or 9 when the triplet occupies x(8),
	! plus 2 when the oxygen carriers occupy x(iox..iox+1).
	! N_eq = 7 (+1 with the triplet) (+2 with the oxygen carriers)
	!        + 2*met_nelem.
	!
	! Physics of the coupling. The two blocks meet in the free electron
	! density: it is the H/He/molecular sum PLUS the metal charges (X+ once,
	! X++ twice). In the shielded, molecular base the metals are the dominant
	! electron donors, so this is not a small correction there -- it is the
	! reason the two networks cannot be solved apart from each other. Every
	! recombination and electron-impact term of BOTH blocks then sees the same
	! n_e. Metal <-> H/He charge exchange is added by cx_add_to_fvec with the
	! metal rows pointed at mbase; the atomic neutral H it reacts with is the
	! free H0 (the H nuclei bound into H2/H2+/H3+/HeH+ are excluded, exactly as
	! in the molecular rows), and the He I reservoir is the free neutral He.
	!
	! With met_nelem = 0 this system reduces exactly to System_HeH_mol, and
	! with the molecular fractions at zero its H/He/metal rows reduce to
	! System_HeH_metals; both limits are shared code, not re-derived here.
	!
	! Solved with hybrd1 (numerical Jacobian), matching the molecular system;
	! no analytic Jacobian is written for the merged residual.

	use global_parameters, only: thereis_HeITR, thereis_oxychem
	use ion_cell_state,    only: ieq_cell
	use ion_residual_core, only: metal_fractions, metal_electron_sum,     &
	                             metal_rows,                              &
	                             impose_transported_ionization_fractions
	use System_HeH_mol,    only: mol_heh_rows, mol_inv_turnover,          &
	                             oxygen_carrier_rows,                     &
	                             oj3, oj4, oj5, oj7
	use System_HeH_metals, only: met_nelem, met_ntot, met_g0, met_g1,     &
	                             met_b0, met_b1, met_a1, met_a2, met_top
	use species_table,     only: iel_O
	use charge_exchange,   only: cx_add_to_fvec, he_h_cx_fvec

	implicit none

	contains

	! Row of the first metal unknown in the molecular layout: the molecular
	! network owns x(1..7), the He 2^3S metastable x(8) when tracked, and the
	! metals follow. Single definition, read by the residual here and by the
	! driver (unknown seeding, root validation, cx_metal_base).
	integer function metal_row_base()
	metal_row_base = oxygen_row_base()
	if (thereis_oxychem) metal_row_base = metal_row_base + 2
	end function metal_row_base

	! Row of the first oxygen-carrier unknown (OH; H2O follows it). It sits
	! immediately above the molecular block and the He 2^3S metastable, and
	! below the metals, so metal_row_base() above shifts by two when the
	! oxygen chemistry is on and is unchanged when it is off. Single
	! definition, read by the residual here and by the driver.
	integer function oxygen_row_base()
	oxygen_row_base = 8
	if (thereis_HeITR) oxygen_row_base = 9
	end function oxygen_row_base

	! Turnover scale of the metal rows, appended to the molecular ones so that
	! the whole system reaches hybrd1 equilibrated (see mol_inv_turnover in
	! System_HeH_mol for what the scale is and why).  A metal row runs as
	! n_Xe times a rate, and a trace element carries n_Xe ~ 1e-4 n_H, so
	! without this the metal block sits far below the helium block for the
	! same reason the molecular one does.  Called once per cell from ioniz_eq
	! AFTER set_mol_turnover_rates (which resets the array) and after
	! set_metal_coeffs (which fills met_*).
	!
	! The rows of an absent element, and the X++ row of a two-stage element,
	! are the identity rows fvec = x written by metal_rows: they are already
	! dimensionless and O(1), so their scale stays 1.  Metal <-> H/He charge
	! exchange is left out of the bound; it is one contribution among several
	! to the same row and the scale only has to be right to within a factor.
	!
	! photo_scale multiplies the two photoionization rates met_g0/met_g1 and
	! nothing else, with the same meaning and the same exactness at 1 as in
	! set_mol_turnover_rates (System_HeH_mol).
	subroutine set_mol_metal_turnover_rates(n_e_ref, photo_scale)
	real*8, intent(in) :: n_e_ref, photo_scale
	real*8  :: sc
	integer :: e, ix, mbase

	mbase = metal_row_base()
	do e = 1,met_nelem
		ix = mbase + 2*(e-1)
		if (met_ntot(e) .le. 1.0d-30) cycle      ! identity rows, scale 1
		! X0 <-> X+ : photoionization, electron-impact ionization and
		! radiative recombination of the first stage.
		sc = met_ntot(e)*(photo_scale*met_g0(e)                      &
		                  + (met_b0(e) + met_a1(e))*n_e_ref)
		if (sc .gt. 0.0d0) mol_inv_turnover(ix) = 1.0d0/sc
		if (met_top(e) .ge. 2) then
			! X+ <-> X++ : the same three channels one stage up.
			sc = met_ntot(e)*(photo_scale*met_g1(e)              &
			                  + (met_b1(e) + met_a2(e))*n_e_ref)
			if (sc .gt. 0.0d0) mol_inv_turnover(ix+1) = 1.0d0/sc
		endif
	enddo

	end subroutine set_mol_metal_turnover_rates

	subroutine ion_system_HeH_mol_metals(N_eq,x,fvec,iflag,params)

	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
	real*8  :: params(60)

	real*8  :: g_hi,g_hei,g_heii,g_heiTR,g_h2,g_h2_di  ! photoionization
	real*8  :: g_h2_dd,g_h2_nd                  ! the other two H2 channels
	real*8  :: g_lw                             ! LW photodissociation
	real*8  :: b_hi,b_hei,b_heii,b_heiTR        ! collisional ionization
	real*8  :: a_hii,a_heii,a_heiii,a_heiTR     ! recombination
	real*8  :: A31,q13,q31a,q31b,Q31            ! triplet kinetics
	real*8  :: n_h,n_he,n_e,ntot
	real*8  :: n_hi,n_hii,n_h2,n_h2p,n_h3p,n_hehp
	real*8  :: n_hei,n_heii,n_heiii,n_heiTR,n_heiSI
	! Oxygen carriers (zero and unused without the oxygen chemistry).
	real*8  :: n_oh,n_h2o,n_ofam
	! Each element's metal densities (neutral/+/++); charge exchange later.
	real*8  :: nm0(met_nelem),nm1(met_nelem),nm2(met_nelem)
	integer :: mbase,iox

	! Unpack the cell state (same named fields as System_HeH_mol)
	g_hi    = ieq_cell%P_HI
	g_hei   = ieq_cell%P_HeI
	g_heii  = ieq_cell%P_HeII
	a_hii   = ieq_cell%rchiiB
	a_heii  = ieq_cell%rcheiiB
	a_heiii = ieq_cell%rcheiiiB
	n_h     = ieq_cell%nh
	n_he    = ieq_cell%nhe
	b_hi    = ieq_cell%a_ion_HI
	b_hei   = ieq_cell%a_ion_HeI
	b_heii  = ieq_cell%a_ion_HeII
	b_heiTR = ieq_cell%a_ion_HeITR   ! He 2^3S collisional ioniz. (0 if no triplet)
	a_heiTR = ieq_cell%rcheiTR
	A31     = ieq_cell%A31
	g_heiTR = ieq_cell%P_HeITR
	q13     = ieq_cell%q13
	q31a    = ieq_cell%q31a
	q31b    = ieq_cell%q31b
	Q31     = ieq_cell%Q31
	g_h2    = ieq_cell%P_H2
	g_h2_di = ieq_cell%P_H2_di
	g_h2_dd = ieq_cell%P_H2_dd   ! double ionization (0 unless a model is on)
	g_h2_nd = ieq_cell%P_H2_nd   ! neutral dissociation (0 unless on)
	g_lw    = ieq_cell%k_LW      ! Lyman-Werner photodissociation (0 if off)
	ntot    = ieq_cell%ntot

	! Oxygen carriers, as fractions of the free oxygen family. They hold H
	! nuclei too (OH one, H2O two), so they enter the atomic-H closure below.
	if (thereis_oxychem) then
		iox    = oxygen_row_base()
		n_ofam = ieq_cell%n_ofam
		n_oh   = x(iox)*n_ofam
		n_h2o  = x(iox+1)*n_ofam
	else
		iox    = 0
		n_ofam = 0.0d0
		n_oh   = 0.0d0
		n_h2o  = 0.0d0
	endif

	! Species densities (verbatim System_HeH_mol)
	n_hi   = (1.0d0 - x(1) - x(4) - x(5) - x(6) - x(7))*n_h
	n_hii  = x(1)*n_h
	n_h2   = 0.5d0*x(4)*n_h
	n_h2p  = 0.5d0*x(5)*n_h
	n_h3p  = x(6)*n_h/3.0d0
	n_hehp = x(7)*n_h
	n_hei   = (1.0d0 - x(2) - x(3))*n_he - n_hehp   ! free neutral He
	n_heii  = x(2)*n_he
	n_heiii = x(3)*n_he
	if (thereis_HeITR) then
		n_heiTR = x(8)*n_he
	else
		n_heiTR = 0.0d0
	endif
	n_heiSI = n_hei - n_heiTR

	! The H nuclei bound in OH and H2O are not atomic H. Written as a
	! separate statement so that the expression above stays bit-for-bit the
	! one a run without the oxygen chemistry evaluates.
	if (thereis_oxychem) n_hi = n_hi - n_oh - 2.0d0*n_h2o

	! Metal ion densities (neutral/+/++) from fractions, canonical order.
	mbase = metal_row_base()
	call metal_fractions(x, mbase, met_nelem, met_ntot, nm0, nm1, nm2)

	! The oxygen row of the metal block is solved against the FREE OXYGEN
	! FAMILY (met_ntot(iel_O) = every O nucleus but the one locked in CO), so
	! its neutral stage as metal_fractions built it still contains the oxygen
	! bound in OH and H2O. Removing it here is the redefinition of O I to
	! mean FREE ATOMIC neutral oxygen (section 4.2 of the design): it is the
	! reactant of the [O I] fine-structure cooling, of the oxygen charge
	! exchange and of O2 below, none of which a bound O atom takes part in.
	! The ionization rows themselves are unchanged -- no reaction of the
	! audited set makes or destroys an oxygen ION -- so O I / O II / O III
	! above the molecular layer are untouched, which is what gate G2 tests.
	if (thereis_oxychem)                                                  &
		nm0(iel_O) = nm0(iel_O) - n_oh - n_h2o

	! Electron density: H/He/molecular first (reproducing the molecular
	! system's sum), then the metal charges.
	n_e = n_hii + n_h2p + n_h3p + n_hehp + n_heii + 2.0d0*n_heiii
	call metal_electron_sum(n_e, met_nelem, nm1, nm2)

	! --- Molecular network rows (System_HeH_mol, metal-inclusive n_e) ---
	call mol_heh_rows(fvec, n_hi, n_hii, n_h2, n_h2p, n_h3p, n_hehp,   &
	                  n_heiSI, n_heiTR, n_heii, n_heiii, n_e, ntot,     &
	                  g_hi, g_hei, g_heii, g_heiTR, g_h2, g_h2_di,      &
	                  g_h2_dd, g_h2_nd, g_lw,                           &
	                  a_hii, a_heii, a_heiii, a_heiTR,                  &
	                  b_hi, b_hei, b_heii, b_heiTR,                     &
	                  q13, q31a, q31b, Q31, A31)

	! --- Oxygen-carrier rows, and the oxygen cycle's exchange with H2 ---
	! Called after mol_heh_rows, which it adds to (fvec(4)). nm0(iel_O) is
	! the free atomic oxygen just formed above. The photolysis channels are
	! this cell's own (set_oxygen_coeffs), passed explicitly so that the
	! same rows can also be evaluated at a scaled radiation field by the
	! continuation of constrained_chemical_equilibrium.
	if (thereis_oxychem)                                                  &
		call oxygen_carrier_rows(fvec, iox, n_hi, n_h2, n_oh, n_h2o,      &
		                         nm0(iel_O), oj3, oj4, oj5, oj7)

	! --- Metal rows (System_HeH_metals, shifted to mbase..) ---
	call metal_rows(fvec, x, mbase, met_nelem, met_ntot, met_g0, met_g1, &
	                met_b0, met_b1, met_a1, met_a2, met_top,             &
	                nm0, nm1, nm2, n_e)

	! Charge exchange (Huang Table 4) on the H, He and metal rows. The
	! driver sets cx_metal_base = mbase before this solve. Rows 1 and 2 are
	! the H+ and He+ balances written production positive, the convention
	! cx_add_to_fvec assumes. Absent reactants contribute zero, preserving
	! the identity rows of absent elements.
	call cx_add_to_fvec(N_eq, fvec, nm0, nm1, nm2,                       &
	                    n_hi, n_hii, n_hei, n_heii, n_heiii)

	! He <-> H charge exchange (Huang Table 4 group B); rows 1 and 2 are
	! production positive here, so he_row_sign = +1. Group B is excluded
	! from cx_act, so it is applied only here (no double counting).
	call he_h_cx_fvec(fvec, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,    &
	                  n_hi, n_hii, n_heiSI, n_heii, 1.0d0)

	! Each row divided by its own turnover rate, molecular block and metal
	! block alike (set_mol_turnover_rates, set_mol_metal_turnover_rates).
	fvec(1:N_eq) = fvec(1:N_eq)*mol_inv_turnover(1:N_eq)

	! IMPOSED CARRIER PARTITIONS. Where a partition is not a local root, the
	! balance row that would have computed it is replaced by the value,
	! AFTER the turnover scaling so that the row is exactly x - x_fix and its
	! Jacobian is the identity; every other row keeps its balance and is
	! solved against it, which is what keeps the ionization stages consistent
	! with the imposed partition.
	!
	! The two flags are separate because the reasons are: H2 is imposed by
	! carrier transport OR by the lower-boundary reservoir composition, the
	! oxygen carriers only by carrier transport.
	if (ieq_cell%x_h2_fixed) fvec(4) = x(4) - ieq_cell%x_h2_fix
	! The transported ionization fractions, where the flow carries them
	! and not this cell's local balance (ion_residual_core). The H2
	! partition just above is imposed for its own reasons and keeps its
	! own flag.
	call impose_transported_ionization_fractions(ieq_cell, x, fvec)
	if (ieq_cell%x_ox_fixed .and. thereis_oxychem) then
		fvec(iox)   = x(iox)   - ieq_cell%x_oh_fix
		fvec(iox+1) = x(iox+1) - ieq_cell%x_h2o_fix
	endif

	return

	end subroutine ion_system_HeH_mol_metals

	! End of module
	end module System_HeH_mol_metals
