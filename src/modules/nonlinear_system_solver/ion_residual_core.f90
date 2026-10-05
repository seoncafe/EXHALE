	module ion_residual_core
	! Pure, stateless helpers for the standard H/He ionization residual rows
	! (with collisional ionization) and their analytic-Jacobian pieces, plus
	! the TR-form H/He/triplet rows (Oklopcic form; electron-impact ionization
	! of H0/He(1^1S)/He+ and of the He 2^3S metastable is included).
	! The standard blocks are verbatim-shared by System_HeH and
	! System_HeH_metals; the TR-form rows and the He 2^3S triplet
	! row are shared by System_HeH_TR, System_HeH_TR_metals and System_HeH_mol.
	! It also
	! holds impose_transported_ionization_fractions, the substitution that
	! replaces an ionization balance row by the fraction the flow carries;
	! every system that can reach those rows calls that one routine. Only
	! explicit-shape / assumed-size dummies are used here: assumed-shape (:)
	! dummies can change gfortran -O3 code generation and break the
	! byte-identical regression. The arithmetic order is
	! kept exactly as in the original inline statements.

	! The one imported quantity is the Penning : associative branching of the
	! He(2^3S) + H total ionization rate, kept in Cool_coeff.f90 next to the
	! rate coefficient itself so it has a single definition.
	use Cooling_Coefficients, only: f_penning_HeI23S
	! The cell state is READ, never written: the imposed-fraction routine
	! below takes the cell as an argument, so this module holds no state of
	! its own and stays safe to call from the OpenMP cell sweep.
	use ion_cell_state, only: ion_rates, ieq_cell

	implicit none

	! ---- THE GROSS CHANNELS OF THE THREE ION-STAGE ROWS ----------------
	! Production and loss of H II, He II and He III reaction by reaction,
	! each a volumetric rate [cm^-3 s^-1] and each a single term of the row
	! expressions below, written from the same factors.  The signed sum of
	! the channels of a stage is that stage's source.
	!
	! WHY THEY ARE HANDED OUT.  Where an ionization stage is fast its
	! production and its loss stand orders above their difference, so a
	! reader of the net source cannot recover them from it, and neither can
	! a reader of the positive and negative parts of that net (which is all
	! the helium diagnostics of the transport operator can state).  The
	! channels are written HERE, at the one place the row terms are
	! written, so that a channel is never a second spelling of a rate.
	! Diagnostic: nothing in the solution reads them.
	!
	! THE CHARGE-EXCHANGE CHANNELS ARE NOT FILLED BY THE ROWS.  Charge
	! exchange is applied by the caller after the rows return, with a
	! reservoir and an orientation that are the caller's, so the channels
	! it owns are left at zero here and filled by whoever applies it.
	! Two reaction sets own channels this way:
	!   the He <-> H charge exchange of Huang et al. (2023) Table 4 group
	!   B and He2+ + H0 -> He+ + H+ (charge_exchange::he_h_cx_fvec),
	!   channels 10, 11, 22, 23 and 33 to 35; calling that routine once
	!   with all but one of its three rate coefficients zeroed isolates
	!   each reaction from the one expression that defines them;
	!   the metal charge exchange of the same table, groups A and C and
	!   the group E electron capture
	!   (charge_exchange::charge_exchange_stage_sources), channels 27 to
	!   32, which that routine hands out as the gross production and the
	!   gross loss of each stage directly.
	! Their channels exist here so that the signed channel sum of a stage
	! is the source the transported row was assembled from, whichever of
	! those reactions the configuration activates.
	integer, parameter :: n_stage_chan = 35
	! Which stage's source each channel belongs to: 1 = H II, 2 = He II,
	! 3 = He III.
	integer, parameter :: stage_chan_row(n_stage_chan) =                  &
	     (/ 1,1,1,1,1,1,1,1,1,1,1,                                        &
	        2,2,2,2,2,2,2,2,2,2,2,2,                                      &
	        3,3,3,                                                        &
	        1,1, 2,2, 3,3,                                                &
	        1, 2, 3 /)
	! Production (+1) or loss (-1) of that stage.
	real*8, parameter :: stage_chan_sign(n_stage_chan) =                  &
	     (/  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,             &
	        -1.0d0,-1.0d0, 1.0d0,-1.0d0,                                  &
	         1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,                           &
	        -1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0, 1.0d0,-1.0d0,             &
	         1.0d0, 1.0d0,-1.0d0,                                         &
	         1.0d0,-1.0d0, 1.0d0,-1.0d0, 1.0d0,-1.0d0,                    &
	         1.0d0, 1.0d0,-1.0d0 /)
	character(len=22), parameter :: stage_chan_name(n_stage_chan) =       &
	     (/ 'HII_photo_HI          ', 'HII_collisional_HI    ',           &
	        'HII_penning_He23S     ', 'HII_H2_photo_diss_ion ',           &
	        'HII_H2_photo_double   ', 'HII_H2p_plus_HI       ',           &
	        'HII_HeII_plus_H2      ', 'HII_recombination     ',           &
	        'HII_plus_H2           ', 'HII_cx_gain_HeII_HI   ',           &
	        'HII_cx_loss_HeI_HII   ',                                     &
	        'HeII_photo_HeI_singlet', 'HeII_photo_HeI_triplet',           &
	        'HeII_coll_HeI_singlet ', 'HeII_coll_HeI_triplet ',           &
	        'HeII_rec_from_HeIII   ', 'HeII_rec_to_HeI       ',           &
	        'HeII_photo_to_HeIII   ', 'HeII_coll_to_HeIII    ',           &
	        'HeII_plus_H2          ', 'HeII_plus_CO          ',           &
	        'HeII_cx_gain_HeI_HII  ', 'HeII_cx_loss_HeII_HI  ',           &
	        'HeIII_photo_HeII      ', 'HeIII_coll_HeII       ',           &
	        'HeIII_recombination   ',                                     &
	        'HII_cx_gain_metal     ', 'HII_cx_loss_metal     ',           &
	        'HeII_cx_gain_metal    ', 'HeII_cx_loss_metal    ',           &
	        'HeIII_cx_gain_metal   ', 'HeIII_cx_loss_metal   ',           &
	        'HII_cx_gain_HeIII_HI  ', 'HeII_cx_gain_HeIII_HI ',           &
	        'HeIII_cx_loss_HeIII_HI' /)
	! The two channels of each stage the caller's He <-> H pair owns.
	integer, parameter :: ich_HII_cx_gain  = 10
	integer, parameter :: ich_HII_cx_loss  = 11
	integer, parameter :: ich_HeII_cx_gain = 22
	integer, parameter :: ich_HeII_cx_loss = 23
	! The three channels of He2+ + H0 -> He+ + H+ (the caller's as well):
	! a proton and a He+ made, a He2+ lost.
	integer, parameter :: ich_HII_cx_gain_hepp   = 33
	integer, parameter :: ich_HeII_cx_gain_hepp  = 34
	integer, parameter :: ich_HeIII_cx_loss_hepp = 35
	! The two channels of each stage the caller's metal charge exchange
	! owns.  No reaction of the present set promotes He II to He III or
	! captures an electron onto He III (Table 4 groups A, C and E, whose
	! helium reactants are He and He+ alone), so the He III pair carries
	! the zero its own gross rates return, and it exists so that the
	! ledger of that stage stays complete if the set grows.
	integer, parameter :: ich_metal_HII_gain   = 27
	integer, parameter :: ich_metal_HII_loss   = 28
	integer, parameter :: ich_metal_HeII_gain  = 29
	integer, parameter :: ich_metal_HeII_loss  = 30
	integer, parameter :: ich_metal_HeIII_gain = 31
	integer, parameter :: ich_metal_HeIII_loss = 32

	! ---- THE OPTIONAL RECORD OF THOSE CHANNELS, CELL BY CELL ----------
	! EXHALE_STAGE_CHANNELS names up to stage_chan_ncell_max grid cells, as
	! a comma-separated list of indices ("107,248").  Unset or empty, which
	! is every run by default, nothing below is reached and no file is
	! opened.  Named, every evaluation of the rows at one of those cells
	! appends one line to ./output/stage_channels.txt carrying the gas
	! branch, the cell, a call counter, the kernel's own arguments and the
	! channels, so that the rows can be re-evaluated afterwards at the very
	! state the run gave them.
	!
	! The rows run inside the OpenMP cell sweep, so the append sits in a
	! critical region; the ORDER of the lines is deterministic only with one
	! thread, which is how a record is meant to be taken.
	integer, parameter :: stage_chan_ncell_max = 8
	! The kernel's own arguments, in the order the record writes them.
	! Entries 1 to 29 are common to both gas branches (26, 27 and 29 are
	! the caller's He <-> H charge-exchange coefficients), 30 to 59 are the
	! molecular network's; the atomic branch writes 29 and stops.
	integer, parameter :: n_stage_state_common = 29
	integer, parameter :: n_stage_state_mol    = 59
	integer, save :: stage_chan_cells(stage_chan_ncell_max) = 0
	integer, save :: stage_chan_ncell = 0
	logical, save :: stage_chan_asked = .false.
	logical, save :: stage_chan_opened = .false.
	integer, save :: stage_chan_icall = 0

	contains

	! Whether this cell is one the record was asked for.  Reads the
	! environment once for the run.
	logical function stage_channels_watch(jcell)
	integer, intent(in) :: jcell
	character(len=128) :: env
	integer :: i, ios, jv, ib, ie
	stage_channels_watch = .false.
	if (.not. stage_chan_asked) then
		!$omp critical (stage_channel_record_region)
		if (.not. stage_chan_asked) then
			call get_environment_variable('EXHALE_STAGE_CHANNELS', env)
			stage_chan_ncell = 0
			ib = 1
			do while (ib .le. len_trim(env) .and.                          &
			          stage_chan_ncell .lt. stage_chan_ncell_max)
				ie = index(env(ib:), ',')
				if (ie .eq. 0) then
					ie = len_trim(env) + 1
				else
					ie = ib + ie - 1
				endif
				read(env(ib:ie-1), *, iostat = ios) jv
				if (ios .eq. 0) then
					stage_chan_ncell = stage_chan_ncell + 1
					stage_chan_cells(stage_chan_ncell) = jv
				endif
				ib = ie + 1
			enddo
			stage_chan_asked = .true.
		endif
		!$omp end critical (stage_channel_record_region)
	endif
	do i = 1, stage_chan_ncell
		if (stage_chan_cells(i) .eq. jcell) then
			stage_channels_watch = .true.
			return
		endif
	enddo
	end function stage_channels_watch

	! One line of the record: the gas branch ('A' atomic, 'M' molecular),
	! the cell, the call counter, the kernel's arguments and the channels.
	subroutine stage_channels_append(gas, jcell, st, nst, chan)
	! gas: the branch and the call site, 'MS' the molecular network as the
	! local sweep solves it, 'MC' the same network as the transport
	! operator evaluates it (the operator is the only caller that asks for
	! the production/loss split, so it is the one the kernel can tell
	! apart), 'AS' the atomic rows as the local sweep solves them and 'AC'
	! the same rows as the transport operator evaluates them, which that
	! caller states with heh_tr_rows' transport_operator argument because
	! it asks the kernel for nothing else that would tell it apart.
	character(len=2), intent(in) :: gas
	integer, intent(in) :: jcell, nst
	real*8, intent(in)  :: st(*), chan(*)
	integer :: u, k
	!$omp critical (stage_channel_record_region)
	if (stage_chan_opened) then
		open(newunit = u, file = './output/stage_channels.txt',            &
		     status = 'old', position = 'append', action = 'write')
	else
		open(newunit = u, file = './output/stage_channels.txt',            &
		     status = 'replace', action = 'write')
		write(u,'(A)') '# gross channels of the ion-stage rows '//        &
		     '(EXHALE_STAGE_CHANNELS), one line per evaluation'
		write(u,'(A)') '# gas cell icall nstate state(1:nstate) '//       &
		     'chan(1:n)   [rates cm^-3 s^-1], n listed below'
		write(u,'(A)') '# gas: MS molecular rows in the local sweep, '//  &
		     'MC the same rows in the transport operator, A- atomic '//   &
		     'rows (call site not distinguishable in the kernel)'
		write(u,'(A)') '# state 1..29 (both branches): T_K n_HI n_HII '// &
		     'n_HeI_singlet n_HeI_triplet n_HeII n_HeIII n_e '//          &
		     'g_HI g_HeI g_HeII g_HeI23S'
		write(u,'(A)') '#   a_HII a_HeII a_HeIII a_HeI23S '//             &
		     'b_HI b_HeI b_HeII b_HeI23S q13 q31a q31b Q31 A31 '//        &
		     'k_cx_He0_Hp k_cx_Hep_H0 q31g k_cx_Hepp_H0 (q13 = direct '// &
		     '1^1S -> 2^3S excitation plus the higher-triplet feed; '//   &
		     'q31g the reverse of the direct one)'
		write(u,'(A)') '# state 30..59 (molecular only): n_tot n_H2 '//   &
		     'n_H2p n_H3p n_HeHp g_H2 g_H2_di g_H2_dd g_H2_nd g_LW '//    &
		     'k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15_two_body k16 '//     &
		     'k17 k18 k19 k_H2p_He k23 k_ion_H2 k_CO_HeII n_third'
		write(u,'(A)') '# channels, in order:'
		do k = 1, n_stage_chan
			write(u,'(A,I3,1X,A,1X,A,I2,A)') '#   ', k,                    &
			     trim(stage_chan_name(k)), 'stage ', stage_chan_row(k),    &
			     merge('  production', '  loss      ',                     &
			           stage_chan_sign(k) .gt. 0.0d0)
		enddo
		write(u,'(A)') '# channels 10, 11, 22, 23, 33 to 35 (He <-> H) and '// &
		     '27 to 32 (metal charge exchange) are the caller''s, zero here'
		stage_chan_opened = .true.
	endif
	stage_chan_icall = stage_chan_icall + 1
	write(u,'(A2,1X,I6,1X,I9,1X,I4,200(1X,ES24.16E3))') gas, jcell,       &
	     stage_chan_icall, nst, (st(k), k = 1, nst),                      &
	     (chan(k), k = 1, n_stage_chan)
	close(u)
	!$omp end critical (stage_channel_record_region)
	end subroutine stage_channels_append

	! The source of one stage from its channels: the signed sum.  The
	! order is the channel order, which is the order the rows write their
	! terms in, so the sum agrees with the assembled row to rounding.
	real*8 function stage_source_from_channels(chan, irow)
	real*8, intent(in)  :: chan(*)
	integer, intent(in) :: irow
	integer :: k
	stage_source_from_channels = 0.0d0
	do k = 1, n_stage_chan
		if (stage_chan_row(k) .eq. irow)                                   &
			stage_source_from_channels = stage_source_from_channels        &
			                           + stage_chan_sign(k)*chan(k)
	enddo
	end function stage_source_from_channels

	! The gross production and the gross loss of one stage, both positive.
	subroutine stage_gross_of_channels(chan, irow, gprod, gloss)
	real*8, intent(in)   :: chan(*)
	integer, intent(in)  :: irow
	real*8, intent(out)  :: gprod, gloss
	integer :: k
	gprod = 0.0d0
	gloss = 0.0d0
	do k = 1, n_stage_chan
		if (stage_chan_row(k) .ne. irow) cycle
		if (stage_chan_sign(k) .gt. 0.0d0) then
			gprod = gprod + chan(k)
		else
			gloss = gloss + chan(k)
		endif
	enddo
	end subroutine stage_gross_of_channels

	! Standard three H/He residual rows. Writes fvec(1), fvec(2), fvec(3) in
	! this exact order. n_e is an INPUT, so the metals system can pass its own
	! metal-inclusive electron density.
	subroutine heh_rows(fvec, n_hi, n_hii, n_hei, n_heii, n_heiii, n_e,  &
	                    g_hi, g_hei, g_heii, a_hii, a_heii, a_heiii,      &
	                    b_hi, b_hei, b_heii, chan)
	real*8 :: fvec(*)
	real*8, intent(in) :: n_hi, n_hii, n_hei, n_heii, n_heiii, n_e
	real*8, intent(in) :: g_hi, g_hei, g_heii, a_hii, a_heii, a_heiii
	real*8, intent(in) :: b_hi, b_hei, b_heii
	! The gross channels of the three stages, written from the same
	! factors as the rows above (declaration of n_stage_chan).  This form
	! carries no metastable, so the triplet channels are zero; row (2) here
	! is the He I balance of a two-level helium and its channels are the
	! same singlet ones the TR form writes.
	real*8, optional, intent(out) :: chan(*)

	fvec(1) = n_hi*g_hi + (n_hi*b_hi - a_hii*n_hii)*n_e
	fvec(2) = n_hei*g_hei + (n_hei*b_hei - a_heii*n_heii)*n_e
	fvec(3) = n_heii*g_heii + (n_heii*b_heii - a_heiii*n_heiii)*n_e
	if (present(chan)) then
		chan(1:n_stage_chan) = 0.0d0
		chan(1)  = n_hi*g_hi
		chan(2)  = n_hi*b_hi*n_e
		chan(8)  = a_hii*n_hii*n_e
		chan(12) = n_hei*g_hei
		chan(14) = n_hei*b_hei*n_e
		chan(16) = a_heiii*n_heiii*n_e
		chan(17) = a_heii*n_heii*n_e
		chan(18) = n_heii*g_heii
		chan(19) = n_heii*b_heii*n_e
		chan(24) = n_heii*g_heii
		chan(25) = n_heii*b_heii*n_e
		chan(26) = a_heiii*n_heiii*n_e
	endif
	end subroutine heh_rows

	! n_e-coefficient C_i of each standard H/He row (the factor multiplying
	! n_e in fvec_i). Sets C1, C2, C3 in this exact order.
	subroutine heh_crow(C1, C2, C3, n_hi, n_hii, n_hei, n_heii, n_heiii,  &
	                    a_hii, a_heii, a_heiii, b_hi, b_hei, b_heii)
	real*8, intent(out) :: C1, C2, C3
	real*8, intent(in)  :: n_hi, n_hii, n_hei, n_heii, n_heiii
	real*8, intent(in)  :: a_hii, a_heii, a_heiii, b_hi, b_hei, b_heii

	C1 = n_hi*b_hi - a_hii*n_hii
	C2 = n_hei*b_hei - a_heii*n_heii
	C3 = n_heii*b_heii - a_heiii*n_heiii
	end subroutine heh_crow

	! Local "direct" terms (photoionization + dC_i/dx_local * n_e) added to the
	! standard H/He Jacobian rows. Applies the five updates in this exact
	! order; call AFTER the rank-1 n_e coupling has populated fjac.
	subroutine heh_jac_local(nsz, fjac, n_h, n_he, n_e, g_hi, g_hei, g_heii, &
	                         a_hii, a_heii, a_heiii, b_hi, b_hei, b_heii)
	integer, intent(in) :: nsz
	real*8 :: fjac(nsz,nsz)
	real*8, intent(in) :: n_h, n_he, n_e
	real*8, intent(in) :: g_hi, g_hei, g_heii, a_hii, a_heii, a_heiii
	real*8, intent(in) :: b_hi, b_hei, b_heii

	fjac(1,1) = fjac(1,1) - n_h*g_hi + (-n_h*b_hi - a_hii*n_h)*n_e
	fjac(2,2) = fjac(2,2) - n_he*g_hei + (-n_he*b_hei - a_heii*n_he)*n_e
	fjac(2,3) = fjac(2,3) - n_he*g_hei + (-n_he*b_hei)*n_e
	fjac(3,2) = fjac(3,2) + n_he*g_heii + (n_he*b_heii)*n_e
	fjac(3,3) = fjac(3,3) + (-a_heiii*n_he)*n_e
	end subroutine heh_jac_local

	! He 2^3S triplet balance row (Oklopcic form + electron-impact ionization
	! of the metastable). Sets ftr to the single triplet residual expression;
	! verbatim-shared as row 4 of the TR / TR_metals systems and row 8 of the
	! molecular system. b_heiTR is the He(2^3S) collisional-ionization
	! coefficient (ionization_rate_HeI_23S, threshold 4.8 eV):
	! He(2^3S)+e- -> He+ + 2e- removes the triplet, so it enters as a
	! destruction term -n_e*n_heiTR*b_heiTR.
	! q31g is the 2^3S -> 1^1S electron-impact de-excitation, the
	! detailed-balance reverse of q13 (Cool_coeff.f90), so the electron
	! collisions between the two levels drive them towards their Boltzmann
	! ratio 3 exp(-19.82 eV/kT); q31a, q31b take the metastable to 2^1S and
	! 2^1P, which decay to 1^1S.
	! Q31 is the TOTAL He(2^3S)+H ionization rate coefficient, Penning plus
	! associative: both channels quench the metastable, so the sink here takes
	! the sum. Only the terms that create a lasting proton or deposit the
	! Penning exothermicity are scaled by f_penning_HeI23S.
	subroutine tr_triplet_row(ftr, n_hi, n_heiSI, n_heiTR, n_heii, n_e,   &
	                          g_heiTR, a_heiTR, q13, q31g, q31a, q31b,     &
	                          Q31, A31, b_heiTR)
	real*8, intent(out) :: ftr
	real*8, intent(in)  :: n_hi, n_heiSI, n_heiTR, n_heii, n_e
	real*8, intent(in)  :: g_heiTR, a_heiTR, q13, q31g, q31a, q31b, Q31, A31
	real*8, intent(in)  :: b_heiTR

	ftr = - n_heiTR*g_heiTR                                         &
	      + n_e*( n_heii*a_heiTR                                    &
	            + n_heiSI*q13                                       &
	            - n_heiTR*(q31g + q31a + q31b)                      &
	            - n_heiTR*b_heiTR)                                  &
	      - n_heiTR*(A31 + n_hi*Q31)
	end subroutine tr_triplet_row

	! The same balance as its two sign-definite halves, production
	! (recombination into the triplets and collisional excitation from the
	! ground singlet) and loss (photoionization, electron-impact transfer
	! and ionization, radiative decay and the He(2^3S) + H quench), from the
	! terms of tr_triplet_row: ftr = prod - loss to rounding.  A record of
	! the channels, for the row-term record of the transported level; the
	! residual is tr_triplet_row's.
	subroutine tr_triplet_row_channels(prod, loss, n_hi, n_heiSI, n_heiTR, &
	                                   n_heii, n_e, g_heiTR, a_heiTR, q13, &
	                                   q31g, q31a, q31b, Q31, A31, b_heiTR)
	real*8, intent(out) :: prod, loss
	real*8, intent(in)  :: n_hi, n_heiSI, n_heiTR, n_heii, n_e
	real*8, intent(in)  :: g_heiTR, a_heiTR, q13, q31g, q31a, q31b, Q31, A31
	real*8, intent(in)  :: b_heiTR

	prod = n_e*(n_heii*a_heiTR + n_heiSI*q13)
	loss = n_heiTR*g_heiTR                                          &
	     + n_e*n_heiTR*(q31g + q31a + q31b + b_heiTR)              &
	     + n_heiTR*(A31 + n_hi*Q31)
	end subroutine tr_triplet_row_channels

	! TR-form four H/He/triplet residual rows (the Oklopcic & Hirata 2018
	! form, with the electron-impact ionization of H I, of He I from the
	! ground singlet and from 2^3S, and of He II added: b_hi, b_hei,
	! b_heiTR, b_heii). Writes fvec(1), fvec(2), fvec(3) in this exact order, then
	! obtains fvec(4) from tr_triplet_row so the triplet expression lives in
	! one place. Verbatim-shared by System_HeH_TR and System_HeH_TR_metals;
	! n_e is an INPUT, so the metals system can pass its metal-inclusive
	! electron density.
	subroutine heh_tr_rows(fvec, n_hi, n_hii, n_heiSI, n_heiTR, n_heii,   &
	                       n_heiii, n_e, g_hi, g_hei, g_heii, g_heiTR,     &
	                       a_hii, a_heii, a_heiii, a_heiTR,                &
	                       b_hi, b_hei, b_heii, b_heiTR,                   &
	                       q13, q31g, q31a, q31b, Q31, A31, chan,         &
	                       transport_operator)
	real*8 :: fvec(*)
	real*8, intent(in) :: n_hi, n_hii, n_heiSI, n_heiTR, n_heii, n_heiii, n_e
	real*8, intent(in) :: g_hi, g_hei, g_heii, g_heiTR
	real*8, intent(in) :: a_hii, a_heii, a_heiii, a_heiTR
	real*8, intent(in) :: b_hi, b_hei, b_heii, b_heiTR
	real*8, intent(in) :: q13, q31g, q31a, q31b, Q31, A31
	! THE STAGE ROWS CHANNEL BY CHANNEL (declaration of n_stage_chan), each
	! a single term of the three rows below and written from the same
	! factors.  The stage sources of an atomic gas are
	!     src(H II) = fvec(1),  src(He II) = -fvec(2) - fvec(3),
	!     src(He III) = fvec(3)
	! (row (2) is the SUMMED He I balance, He-I-gain positive), and the
	! signed channel sums reproduce those three.  Diagnostic only.
	real*8, optional, intent(out) :: chan(*)
	! Which evaluation of these rows the record is of: .true. from the
	! transported-stage source of diffusive_photochemistry::carrier_source,
	! absent or .false. from the local ionization sweep.  The rows do not
	! depend on it; it labels the record line, so that the two evaluations
	! of one cell can be told apart afterwards as the molecular branch's
	! already can.
	logical, optional, intent(in) :: transport_operator
	real*8  :: cloc(n_stage_chan), stt(n_stage_state_common)
	logical :: wrec, from_operator

	! Penning ionization source He(2^3S)+H0 -> He(1^1S) + H+ + e-: the same
	! collision that removes the triplet in tr_triplet_row ionizes H0, so it
	! enters here as an H+ production term. Q31 is the TOTAL ionization rate,
	! of which only the Penning branch f_penning_HeI23S leaves a proton behind.
	! The remaining 10% is associative ionization He(2^3S)+H0 -> HeH+ + e-;
	! this system carries no HeH+, and in an atomic gas HeH+ dissociatively
	! recombines back to He + H (1.3e-8 cm^3 s^-1 at 2000 K, faster than any
	! competing HeH+ reaction there), returning the H atom and consuming the
	! electron, so that branch is a metastable sink but not a lasting proton
	! or electron source. Where H2 is abundant that cycle is broken by
	! HeH+ + H2 -> H3+ + He, and System_HeH_mol therefore does carry the
	! branch explicitly into its HeH+ row.
	! n_hi*b_hi*n_e is the electron-impact ionization of H0 (Voronov b_hi),
	! restored to match the standard heh_rows.
	fvec(1) = n_hi*g_hi + f_penning_HeI23S*n_heiTR*n_hi*Q31            &
	        + n_hi*b_hi*n_e - a_hii*n_hii*n_e

	! New equation for hei - sum of the two equations of Oklopcic. The last
	! two terms are electron-impact ionization of ground-state (b_hei) and
	! metastable (b_heiTR) He I; both remove He I and produce He+, so they sit
	! on the loss side of this summed He I balance.
	fvec(2) =   n_heii*(a_heiTR + a_heii)*n_e                            &
	          - n_heiSI*g_hei                                           &
	          - n_heiTR*g_heiTR                                         &
	          - n_heiSI*b_hei*n_e                                       &
	          - n_heiTR*b_heiTR*n_e

	! n_heii*b_heii*n_e: electron-impact ionization of He+ into He++.
	fvec(3) = n_heii*g_heii + n_heii*b_heii*n_e - a_heiii*n_heiii*n_e

	call tr_triplet_row(fvec(4), n_hi, n_heiSI, n_heiTR, n_heii, n_e,     &
	                    g_heiTR, a_heiTR, q13, q31g, q31a, q31b, Q31,     &
	                    A31, b_heiTR)

	! The same terms one by one, and the record of them when this cell was
	! named.  The rows above are untouched by either.
	wrec = .false.
	if (stage_chan_ncell .gt. 0 .or. .not. stage_chan_asked)              &
		wrec = stage_channels_watch(ieq_cell%jcell)
	if (present(chan) .or. wrec) then
		cloc     = 0.0d0
		cloc(1)  = n_hi*g_hi
		cloc(2)  = n_hi*b_hi*n_e
		cloc(3)  = f_penning_HeI23S*n_heiTR*n_hi*Q31
		cloc(8)  = a_hii*n_hii*n_e
		cloc(12) = n_heiSI*g_hei
		cloc(13) = n_heiTR*g_heiTR
		cloc(14) = n_heiSI*b_hei*n_e
		cloc(15) = n_heiTR*b_heiTR*n_e
		cloc(16) = a_heiii*n_heiii*n_e
		cloc(17) = (a_heiTR + a_heii)*n_heii*n_e
		cloc(18) = n_heii*g_heii
		cloc(19) = n_heii*b_heii*n_e
		cloc(24) = n_heii*g_heii
		cloc(25) = n_heii*b_heii*n_e
		cloc(26) = a_heiii*n_heiii*n_e
		if (present(chan)) chan(1:n_stage_chan) = cloc
		if (wrec) then
			stt( 1) = ieq_cell%T_K
			stt( 2) = n_hi;     stt( 3) = n_hii
			stt( 4) = n_heiSI;  stt( 5) = n_heiTR
			stt( 6) = n_heii;   stt( 7) = n_heiii;  stt( 8) = n_e
			stt( 9) = g_hi;     stt(10) = g_hei
			stt(11) = g_heii;   stt(12) = g_heiTR
			stt(13) = a_hii;    stt(14) = a_heii
			stt(15) = a_heiii;  stt(16) = a_heiTR
			stt(17) = b_hi;     stt(18) = b_hei
			stt(19) = b_heii;   stt(20) = b_heiTR
			stt(21) = q13;      stt(22) = q31a;     stt(23) = q31b
			stt(24) = Q31;      stt(25) = A31
			stt(26) = ieq_cell%kcx_He0_Hp
			stt(27) = ieq_cell%kcx_Hep_H0
			stt(28) = q31g
			stt(29) = ieq_cell%kcx_Hepp_H0
			from_operator = .false.
			if (present(transport_operator))                               &
				from_operator = transport_operator
			if (from_operator) then
				call stage_channels_append('AC', ieq_cell%jcell, stt,       &
				                           n_stage_state_common, cloc)
			else
				call stage_channels_append('AS', ieq_cell%jcell, stt,       &
				                           n_stage_state_common, cloc)
			endif
		endif
	endif
	end subroutine heh_tr_rows

	! Metal ion densities (neutral/+/++) from fractions, in canonical element
	! order. `base` is the row of the first metal unknown (4 in the metals
	! system, 5 with the He 2^3S triplet), so ix = base + 2*(e-1) locates the
	! X+/X++ fractions of element e. Sets nm0/nm1/nm2 in this exact order.
	subroutine metal_fractions(x, base, nelem, mtot, nm0, nm1, nm2)
	integer, intent(in) :: base, nelem
	real*8, intent(in)  :: x(*)
	real*8, intent(in)  :: mtot(nelem)
	real*8, intent(out) :: nm0(nelem), nm1(nelem), nm2(nelem)
	integer :: e, ix
	real*8  :: n_X

	do e = 1,nelem
		ix     = base + 2*(e-1)
		n_X    = mtot(e)
		nm1(e) = x(ix)*n_X
		nm2(e) = x(ix+1)*n_X
		nm0(e) = (1.0 - x(ix) - x(ix+1))*n_X
	enddo
	end subroutine metal_fractions

	! Add the metal electron contribution to n_e (X+ counts once, X++ twice),
	! in canonical element order. n_e accumulates in place, reproducing the
	! original loop's add order exactly.
	subroutine metal_electron_sum(n_e, nelem, nm1, nm2)
	integer, intent(in)   :: nelem
	real*8, intent(inout) :: n_e
	real*8, intent(in)    :: nm1(nelem), nm2(nelem)
	integer :: e

	do e = 1,nelem
		n_e = n_e + nm1(e) + 2.0*nm2(e)
	enddo
	end subroutine metal_electron_sum

	! Metal ionization balance rows, one element at a time (force-zero if the
	! element is absent; otherwise the normal balance). `base` locates the
	! first metal row (4 or 5); ix = base + 2*(e-1). Three-stage elements
	! (mtop >= 2) solve both X0<->X+ and X+<->X++; two-stage elements solve
	! only X0<->X+ and pin the unused upper unknown. Expression order is
	! verbatim from the System_HeH_metals residual.
	subroutine metal_rows(fvec, x, base, nelem, mtot, mg0, mg1, mg02,  &
	                      mb0, mb1, ma1, ma2, mtop, nm0, nm1, nm2, n_e,    &
	                      gross)
	integer, intent(in) :: base, nelem
	real*8 :: fvec(*)
	real*8, intent(in)  :: x(*)
	real*8, intent(in)  :: mtot(nelem), mg0(nelem), mg1(nelem)
	! mg02: the part of the neutral's photoionization rate mg0 that ejects
	! two or more electrons (an autoionizing inner-shell vacancy) and so
	! takes the atom straight to X++. The rows are the net flows across the
	! two stage boundaries, so a direct X0 -> X++ event crosses both: it is
	! in mg0 on the X0 <-> X+ row already, and adds nm0*mg02 to the
	! X+ <-> X++ row. The stage sum is untouched (the unknowns are the
	! fractions of X+ and X++, the neutral their complement).
	real*8, intent(in)  :: mg02(nelem)
	real*8, intent(in)  :: mb0(nelem), mb1(nelem), ma1(nelem), ma2(nelem)
	integer, intent(in) :: mtop(nelem)
	real*8, intent(in)  :: nm0(nelem), nm1(nelem), nm2(nelem)
	real*8, intent(in)  :: n_e
	! The gross rate of each row, the sum of the magnitudes of its terms
	! [cm^-3 s^-1] (System_HeH_mol, mol_inv_turnover); zero for an identity
	! row, which carries no reaction.
	real*8, optional, intent(inout) :: gross(*)
	integer :: e, ix

	do e = 1,nelem
		ix = base + 2*(e-1)
		if (present(gross)) then
			gross(ix)   = 0.0d0
			gross(ix+1) = 0.0d0
			if (mtot(e) .gt. 1.0d-30) then
				gross(ix) = abs(nm0(e)*mg0(e)) + abs(nm0(e)*mb0(e)*n_e)     &
				          + abs(ma1(e)*nm1(e)*n_e)
				if (mtop(e) .ge. 2) gross(ix+1) = abs(nm1(e)*mg1(e))         &
				          + abs(nm0(e)*mg02(e)) + abs(nm1(e)*mb1(e)*n_e)    &
				          + abs(ma2(e)*nm2(e)*n_e)
			endif
		endif
		if (mtot(e) .le. 1.0d-30) then
			fvec(ix)   = x(ix)
			fvec(ix+1) = x(ix+1)
		else
			! X0 <-> X+
			fvec(ix)   = nm0(e)*mg0(e)                            &
			           + (nm0(e)*mb0(e) - ma1(e)*nm1(e))*n_e
			if (mtop(e) .ge. 2) then
				! X+ <-> X++
				fvec(ix+1) = nm1(e)*mg1(e) + nm0(e)*mg02(e)           &
				           + (nm1(e)*mb1(e) - ma2(e)*nm2(e))*n_e
			else
				! Two-stage element: no X++, pin the unused unknown.
				fvec(ix+1) = x(ix+1)
			endif
		endif
	enddo
	end subroutine metal_rows

	! Impose the ionization fractions the flow carries.
	!
	! Where an ionization stage is transported, its fraction in a cell is
	! not a root of that cell's photoionization balance: the ionization
	! time exceeds the flow time and the composition of the gas is the one
	! the flow brought. The balance row that would have computed the
	! fraction is then replaced by the carried value,
	!
	!       fvec(i) = x(i) - x_fix(i) ,
	!
	! whose Jacobian row is the identity. Every other row of the system
	! keeps its balance and is solved against it, which is what keeps the
	! remaining stages -- helium, the molecular ions, the metals, and the
	! electron density they all share -- consistent with the transported
	! composition.
	!
	! Rows 1 to 3 hold the same three stages in every system that reaches
	! them: x(1) = n(H II)/n(H nuclei), x(2) = n(He II)/n(He) and
	! x(3) = n(He III)/n(He). System_H offers row 1 alone. One flag and one
	! value per stage live in the cell state, and only a stage whose flag is
	! set is written, so a system may call this whether or not anything is
	! imposed on it.
	!
	! CALL IT LAST, after every term of the rows has been written and after
	! any scaling of them, so that the row is exactly x - x_fix.
	!
	! This is the one place the substitution is written. Seven systems reach
	! these rows, and seven transcriptions of one substitution is the drift
	! the g2s/G2s shadowing of 2026-08-12 stands as the warning about.
	subroutine impose_transported_ionization_fractions(cell, x, fvec)
	type(ion_rates), intent(in) :: cell
	real*8, intent(in) :: x(*)
	real*8 :: fvec(*)

	! H II per H nucleus, carried by the transported proton
	! (the key Ionization transport).
	if (cell%x_hp_fixed) fvec(1) = x(1) - cell%x_hp_fix
	! He II and He III per He nucleus, rows 2 and 3, carried by the same
	! transported partition. A system that offers row 1 alone (System_H) is
	! solved only in a mixture with no helium, where no helium fraction can
	! be imposed, so these two rows exist wherever the flags can be set.
	if (cell%x_heii_fixed)  fvec(2) = x(2) - cell%x_heii_fix
	if (cell%x_heiii_fixed) fvec(3) = x(3) - cell%x_heiii_fix
	! The He 2^3S level per He nucleus, carried by the transported level
	! (He 2^3S transport), in the row its balance holds in this system
	! (tr_triplet_row: row 4 of the atomic triplet systems, row 8 of the
	! molecular ones). The flag is set only where the system carries the
	! level, so the row named is that balance.
	if (cell%x_hetr_fixed .and. cell%x_hetr_row .gt. 0)                  &
		fvec(cell%x_hetr_row) = x(cell%x_hetr_row) - cell%x_hetr_fix
	end subroutine impose_transported_ionization_fractions

	! End of module
	end module ion_residual_core
