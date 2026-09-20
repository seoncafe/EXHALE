      program stage_row_balance_tests
      ! THE GROSS CHANNELS OF THE ION-STAGE ROWS, AND WHAT THEY SAY ABOUT
      ! THE BALANCE A TRANSPORTED STAGE IS ASKED TO SATISFY
      ! (PLAN_20260918_rev2 item D3a, steps 1, 2 and 4).
      !
      ! The channels are an optional output of the two production kernels
      ! the stage sources come from -- mol_heh_rows in molecular gas and
      ! heh_tr_rows (heh_rows without the metastable) in atomic gas -- so
      ! nothing here writes a rate a second time.  The He <-> H charge
      ! exchange is applied by the CALLER of those kernels, with a reservoir
      ! and an orientation that are the caller's, and the two directions are
      ! isolated here by calling charge_exchange::he_h_cx_fvec once with
      ! each of its two rate coefficients set to zero, which is again the
      ! one expression that defines them.
      !
      ! ROWS.
      !   the_channel_sums_are_the_assembled_source_molecular
      !   the_channel_sums_are_the_assembled_source_atomic
      !   the_channel_sums_are_the_assembled_source_atomic_no_metastable
      !       the ledger identity in both gas branches and with the He 2^3S
      !       level tracked and untracked: the signed sum of the channels of
      !       a stage is that stage's source AFTER the charge-exchange pair,
      !       to rounding.
      !   the_recorded_channels_are_the_rows_of_the_recorded_state
      !   the_recorded_channel_sums_are_the_recorded_source
      !       the same identity on the record a run left behind
      !       (EXHALE_STAGE_CHANNELS), with the kernel re-called from the
      !       state the record carries, so that the record is shown to be
      !       complete: the channels come back bitwise.
      !   the_stage_sum_identity_is_within_its_rounding_bound
      !       the floating-point bound of the stage-nucleus-sum entry,
      !       derived from the number and the size of the summed terms, on
      !       manufactured columns whose eddy term spans four decades.
      !
      ! MEASUREMENT.  With EXHALE_STAGE_ROW_RECORD naming a
      ! stage_channels.txt and EXHALE_STAGE_ROW_TERMS a carrier_row_terms.txt
      ! of the same state, the tables of D3a step 2 are printed: the gross
      ! production and loss of each stage beside the transport terms and the
      ! row scale, the sensitivity of the source to the stage fraction and
      ! to the electron density, the abundance correction that would close
      ! the transport, the cancellation floor of the net source, and the
      ! smoothness of the row under budget-respecting perturbations.
      !
      ! Every row prints one
      !     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
      ! line; lines beginning with two spaces or with DIAGNOSTIC are
      ! context.  Exit status is nonzero if any row fails.

      use global_parameters
      use ion_cell_state,     only: ieq_cell
      use ion_residual_core,  only: heh_rows, heh_tr_rows,                &
                                    n_stage_chan, stage_chan_name,        &
                                    stage_chan_row, stage_chan_sign,      &
                                    stage_source_from_channels,           &
                                    stage_gross_of_channels,              &
                                    ich_HII_cx_gain, ich_HII_cx_loss,     &
                                    ich_HeII_cx_gain, ich_HeII_cx_loss,   &
                                    ich_metal_HII_gain, ich_metal_HII_loss, &
                                    ich_metal_HeII_gain,                  &
                                    ich_metal_HeII_loss,                  &
                                    ich_metal_HeIII_gain,                 &
                                    ich_metal_HeIII_loss,                 &
                                    n_stage_state_common, n_stage_state_mol
      use System_HeH_mol,     only: mol_heh_rows, set_mol_coeffs,         &
                                    mk5, mk6, mk7, mk8, mk9, mk10, mk11,  &
                                    mk12, mk13, mk14, mk15, mk16, mk17,   &
                                    mk18, mk19, mk_h2p_he, mk23,          &
                                    mk_ion_H2
      use charge_exchange,    only: he_h_cx_fvec, he_h_charge_exchange,   &
                                    cx_full, cx_init, cx_set_cell,        &
                                    charge_exchange_stage_sources
      use species_table,      only: n_melem, iel_C, iel_O, iel_Si
      use mol_rates,          only: h2_thermochemistry_init
      use ionization_stage_transport, only: ionization_stage_face_flux
      use Cooling_Coefficients, only: f_penning_HeI23S

      implicit none

      integer, parameter :: dp = kind(1.0d0)
      integer, parameter :: qp = selected_real_kind(30)
      real(dp), parameter :: eps_d = epsilon(1.0d0)

      ! The largest record a fixture measurement reads.
      integer, parameter :: nrec_max = 40000
      integer, parameter :: nst_max  = 64

      integer :: nfail

      nfail = 0
      ! The H2 equilibrium-constant table is built serially first, as the
      ! sweep does before it enters its parallel region.
      call h2_thermochemistry_init()

      call ledger_of_a_molecular_cell()
      call ledger_of_an_atomic_cell(.true.)
      call ledger_of_an_atomic_cell(.false.)
      call ledger_of_a_metal_bearing_cell(.true.)
      call ledger_of_a_metal_bearing_cell(.false.)
      call stage_sum_rounding_bound()
      call ledger_and_measurement_of_a_record()

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'FAILED ', nfail, ' row(s)'
         call exit(1)
      endif
      call exit(0)

      contains

      ! ----------------------------------------------------------------- !

      subroutine verdict(name, ok, val, ref, tol)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: ok
      real(dp),         intent(in) :: val, ref, tol
      if (ok) then
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES9.2)') 'PASS ', name,        &
              ' measured=', val, ' reference=', ref, ' tol=', tol
      else
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES9.2)') 'FAIL ', name,        &
              ' measured=', val, ' reference=', ref, ' tol=', tol
         nfail = nfail + 1
      endif
      end subroutine verdict

      ! ----------------------------------------------------------------- !

      ! The charge-exchange pair as the caller applies it, direction by
      ! direction.  he_h_cx_fvec returns the two directions' DIFFERENCE in
      ! each row; calling it once with each rate coefficient zeroed leaves
      ! one direction alone, so the two come from the one expression that
      ! defines them and never from a copy of it.  he_row_sign is the
      ! caller's orientation: +1 where row 2 is written He-II-gain positive
      ! (the molecular kernel), -1 where it is the He I balance written
      ! He-I-gain positive (the atomic kernel).
      subroutine charge_exchange_channels(k1, k2, n_hi, n_hii, n_hei_cx,  &
                                          n_heii, he_row_sign, chan)
      real(dp), intent(in)    :: k1, k2, n_hi, n_hii, n_hei_cx, n_heii
      real(dp), intent(in)    :: he_row_sign
      real(dp), intent(inout) :: chan(*)
      real(dp) :: f1(4), f2(4)
      f1 = 0.0d0
      f2 = 0.0d0
      ! Direction B1 alone, He0 + H+ -> He+ + H0: a proton lost, a He+ made.
      call he_h_cx_fvec(f1, k1, 0.0d0, n_hi, n_hii, n_hei_cx, n_heii,     &
                        he_row_sign)
      ! Direction B2 alone, He+ + H0 -> He0 + H+: a proton made, a He+ lost.
      call he_h_cx_fvec(f2, 0.0d0, k2, n_hi, n_hii, n_hei_cx, n_heii,     &
                        he_row_sign)
      chan(ich_HII_cx_loss)  = -f1(1)
      chan(ich_HII_cx_gain)  =  f2(1)
      chan(ich_HeII_cx_gain) =  he_row_sign*f1(2)
      chan(ich_HeII_cx_loss) = -he_row_sign*f2(2)
      end subroutine charge_exchange_channels

      ! ----------------------------------------------------------------- !

      ! The three stage sources of a molecular cell from the rows the
      ! transport operator takes whole: src(H II) = fvec(1),
      ! src(He II) = fvec(2), src(He III) = fvec(3), with the pair applied
      ! at he_row_sign = +1 on the ground singlet, which is the reservoir
      ! System_HeH_mol passes.
      subroutine molecular_stage_sources(st, src, chan)
      real(dp), intent(in)  :: st(*)
      real(dp), intent(out) :: src(3), chan(*)
      real(dp) :: fv(10)
      fv = 0.0d0
      chan(1:n_stage_chan) = 0.0d0
      call set_mol_coeffs(st(1), st(28))
      call mol_heh_rows(fv, st(2), st(3), st(29), st(30), st(31), st(32), &
                        st(4), st(5), st(6), st(7), st(8), st(28),        &
                        st(9), st(10), st(11), st(12), st(33),            &
                        st(34), st(35), st(36), st(37),                   &
                        st(13), st(14), st(15), st(16),                   &
                        st(17), st(18), st(19), st(20),                   &
                        st(21), st(22), st(23), st(24), st(25),           &
                        chan = chan)
      call he_h_cx_fvec(fv, st(26), st(27), st(2), st(3), st(4), st(6),   &
                        1.0d0)
      call charge_exchange_channels(st(26), st(27), st(2), st(3), st(4),  &
                                    st(6), 1.0d0, chan)
      src(1) = fv(1)
      src(2) = fv(2)
      src(3) = fv(3)
      end subroutine molecular_stage_sources

      ! The same for an atomic cell.  Row (2) of heh_tr_rows is the SUMMED
      ! He I balance written He-I-gain positive, so
      ! src(He II) = -fvec(2) - fvec(3), and the pair is applied at
      ! he_row_sign = -1 on the GROUND SINGLET He(1^1S), which is the
      ! reservoir System_HeH_TR passes: the Table 4 rate of He + H+ carries
      ! the barrier exp(-12.75/T4), the ionization-potential difference of
      ! ground-state helium against hydrogen, and the metastable's own
      ! charge exchange is a different reaction this code carries nowhere.
      subroutine atomic_stage_sources(st, src, chan)
      real(dp), intent(in)  :: st(*)
      real(dp), intent(out) :: src(3), chan(*)
      real(dp) :: fv(4)
      fv = 0.0d0
      chan(1:n_stage_chan) = 0.0d0
      call heh_tr_rows(fv, st(2), st(3), st(4), st(5), st(6), st(7),      &
                       st(8), st(9), st(10), st(11), st(12),              &
                       st(13), st(14), st(15), st(16),                    &
                       st(17), st(18), st(19), st(20),                    &
                       st(21), st(22), st(23), st(24), st(25),            &
                       chan = chan)
      call he_h_cx_fvec(fv, st(26), st(27), st(2), st(3), st(4), st(6),   &
                        -1.0d0)
      call charge_exchange_channels(st(26), st(27), st(2), st(3), st(4),  &
                                    st(6), -1.0d0, chan)
      src(1) =  fv(1)
      src(2) = -fv(2) - fv(3)
      src(3) =  fv(3)
      end subroutine atomic_stage_sources

      ! ----------------------------------------------------------------- !

      ! THE METAL CHARGE EXCHANGE AS THE TRANSPORT OPERATOR APPLIES IT.
      ! charge_exchange_stage_sources returns the stoichiometric source of
      ! each hydrogen and helium stage and the gross production and loss the
      ! same reactions give it, so the channels come from the one expression
      ! that defines them; the stage sources are added to the rows in the
      ! basis of the unknowns, which is the stage density itself, exactly as
      ! diffusive_photochemistry::carrier_source does.
      subroutine metal_charge_exchange_channels(nm0, nm1, nm2, n_hi,      &
                          n_hii, n_heiSI, n_heii, n_heiii, T, chan, dsrc)
      real(dp), intent(in)    :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real(dp), intent(in)    :: n_hi, n_hii, n_heiSI, n_heii, n_heiii, T
      real(dp), intent(inout) :: chan(*)
      real(dp), intent(out)   :: dsrc(3)
      real(dp) :: sH(0:1), sHe(0:2)
      real(dp) :: pH(0:1), lH(0:1), pHe(0:2), lHe(0:2)
      call cx_set_cell(T)
      call charge_exchange_stage_sources(nm0, nm1, nm2, n_hi, n_hii,      &
                       n_heiSI, n_heii, n_heiii, T, sH, sHe,              &
                       p_H = pH, l_H = lH, p_He = pHe, l_He = lHe)
      chan(ich_metal_HII_gain)   = pH(1)
      chan(ich_metal_HII_loss)   = lH(1)
      chan(ich_metal_HeII_gain)  = pHe(1)
      chan(ich_metal_HeII_loss)  = lHe(1)
      chan(ich_metal_HeIII_gain) = pHe(2)
      chan(ich_metal_HeIII_loss) = lHe(2)
      dsrc(1) = sH(1)
      dsrc(2) = sHe(1)
      dsrc(3) = sHe(2)
      end subroutine metal_charge_exchange_channels

      ! ----------------------------------------------------------------- !

      ! |signed channel sum - assembled source| over the largest channel of
      ! the stage: the ledger identity's own relative measure.
      real(dp) function ledger_departure(chan, irow, src)
      real(dp), intent(in) :: chan(*), src
      integer,  intent(in) :: irow
      real(dp) :: s, gp, gl
      s = stage_source_from_channels(chan, irow)
      call stage_gross_of_channels(chan, irow, gp, gl)
      ledger_departure = abs(s - src)/max(gp, gl, abs(src), 1.0d-300)
      end function ledger_departure

      ! ----------------------------------------------------------------- !

      subroutine ledger_of_a_molecular_cell()
      ! A MOLECULAR CELL OF THE HOT-URANUS KIND: a partly dissociated,
      ! weakly ionized gas at 1800 K with helium at 8 per cent of hydrogen.
      ! The numbers are a state, not a measurement; what the row asserts is
      ! the identity, which holds at any state.
      real(dp) :: st(nst_max), src(3), chan(n_stage_chan)
      real(dp) :: d1, d2, d3
      logical  :: mol_save, tr_save, oxy_save

      mol_save = thereis_mol
      tr_save  = thereis_HeITR
      oxy_save = thereis_oxychem
      thereis_mol      = .true.
      thereis_HeITR    = .true.
      thereis_oxychem  = .false.
      he_h_charge_exchange = .true.

      st = 0.0d0
      st( 1) = 1.8d3                 ! T [K]
      st( 2) = 2.0d9                 ! n(H I)
      st( 3) = 1.2d7                 ! n(H II)
      st( 4) = 2.4d9                 ! n(He I) ground singlet
      st( 5) = 1.0d2                 ! n(He 2^3S)
      st( 6) = 1.3d3                 ! n(He II)
      st( 7) = 3.6d0                 ! n(He III)
      st( 8) = 1.3d7                 ! n_e
      st( 9) = 2.0d-6                ! g(H I)
      st(10) = 3.0d-7                ! g(He I)
      st(11) = 4.0d-8                ! g(He II)
      st(12) = 1.0d-4                ! g(He 2^3S)
      st(13) = 2.5d-13               ! a(H II)
      st(14) = 1.5d-13               ! a(He II)
      st(15) = 8.0d-13               ! a(He III)
      st(16) = 5.0d-14               ! a(He 2^3S)
      st(17) = 1.0d-18               ! b(H I)
      st(18) = 2.0d-20               ! b(He I)
      st(19) = 1.0d-22               ! b(He II)
      st(20) = 3.0d-11               ! b(He 2^3S)
      st(21) = 1.0d-14               ! q13
      st(22) = 2.0d-12               ! q31a
      st(23) = 1.0d-12               ! q31b
      st(24) = 5.0d-10               ! Q31
      st(25) = 1.272d-4              ! A31
      st(26) = 1.0d-15               ! k(He0 + H+)
      st(27) = 1.0d-9                ! k(He+ + H0)
      st(28) = 5.0d9                 ! n_tot
      st(29) = 1.2d9                 ! n(H2)
      st(30) = 3.0d1                 ! n(H2+)
      st(31) = 2.0d3                 ! n(H3+)
      st(32) = 5.0d-3                ! n(HeH+)
      st(33) = 1.0d-7                ! g(H2), total destruction
      st(34) = 2.0d-8                ! g(H2) dissociative ionization
      st(35) = 0.0d0                 ! g(H2) double ionization
      st(36) = 0.0d0                 ! g(H2) neutral dissociation
      st(37) = 0.0d0                 ! Lyman-Werner
      ieq_cell%T_K  = st(1)
      ieq_cell%ntot = st(28)
      ieq_cell%n_co = 0.0d0
      ieq_cell%jcell = 0

      call molecular_stage_sources(st, src, chan)
      d1 = ledger_departure(chan, 1, src(1))
      d2 = ledger_departure(chan, 2, src(2))
      d3 = ledger_departure(chan, 3, src(3))
      call report_channels('molecular cell', st, src, chan)
      call verdict('the_channel_sums_are_the_assembled_source_molecular', &
                   max(d1, d2, d3) .le. 8.0d0*eps_d,                      &
                   max(d1, d2, d3), 0.0d0, 8.0d0*eps_d)

      thereis_mol     = mol_save
      thereis_HeITR   = tr_save
      thereis_oxychem = oxy_save
      end subroutine ledger_of_a_molecular_cell

      ! ----------------------------------------------------------------- !

      subroutine ledger_of_an_atomic_cell(with_metastable)
      ! AN ATOMIC CELL OF THE WIND: hydrogen half ionized at 6000 K, helium
      ! carrying all three stages.  With the metastable untracked its four
      ! coefficients are zero and the rows reduce to the standard H/He
      ! balances, which the second call checks against heh_rows itself.
      logical, intent(in) :: with_metastable
      real(dp) :: st(nst_max), src(3), chan(n_stage_chan)
      real(dp) :: d1, d2, d3, fv(4), dstd
      logical  :: mol_save, tr_save
      character(len=80) :: nm

      mol_save = thereis_mol
      tr_save  = thereis_HeITR
      thereis_mol   = .false.
      thereis_HeITR = with_metastable
      he_h_charge_exchange = .true.

      st = 0.0d0
      st( 1) = 6.0d3
      st( 2) = 1.0d5
      st( 3) = 1.0d5
      st( 4) = 7.0d3
      st( 5) = 0.0d0
      st( 6) = 1.0d3
      st( 7) = 1.0d1
      st( 8) = 1.01d5
      st( 9) = 1.0d-4
      st(10) = 2.0d-5
      st(11) = 1.0d-6
      st(12) = 0.0d0
      st(13) = 1.6d-13
      st(14) = 1.0d-13
      st(15) = 5.0d-13
      st(16) = 0.0d0
      st(17) = 2.0d-11
      st(18) = 1.0d-12
      st(19) = 1.0d-14
      st(20) = 0.0d0
      st(26) = 1.0d-15
      st(27) = 1.0d-9
      if (with_metastable) then
         st( 5) = 3.0d0
         st(12) = 1.0d-4
         st(16) = 5.0d-14
         st(20) = 3.0d-11
         st(21) = 1.0d-14
         st(22) = 2.0d-12
         st(23) = 1.0d-12
         st(24) = 5.0d-10
         st(25) = 1.272d-4
      endif
      ieq_cell%T_K   = st(1)
      ieq_cell%jcell = 0

      call atomic_stage_sources(st, src, chan)
      d1 = ledger_departure(chan, 1, src(1))
      d2 = ledger_departure(chan, 2, src(2))
      d3 = ledger_departure(chan, 3, src(3))
      if (with_metastable) then
         nm = 'the_channel_sums_are_the_assembled_source_atomic'
      else
         nm = 'the_channel_sums_are_the_assembled_source_atomic_'//       &
              'no_metastable'
      endif
      call report_channels('atomic cell', st, src, chan)
      call verdict(trim(nm), max(d1, d2, d3) .le. 8.0d0*eps_d,            &
                   max(d1, d2, d3), 0.0d0, 8.0d0*eps_d)

      ! Without the metastable the TR form is the standard form, so the
      ! standard kernel's own channels have to give the same three sources.
      if (.not. with_metastable) then
         fv = 0.0d0
         chan(1:n_stage_chan) = 0.0d0
         call heh_rows(fv, st(2), st(3), st(4), st(6), st(7), st(8),      &
                       st(9), st(10), st(11), st(13), st(14), st(15),     &
                       st(17), st(18), st(19), chan)
         dstd = abs((fv(2) - fv(3)) -                                     &
                    stage_source_from_channels(chan, 2))                  &
              / max(abs(fv(2)), abs(fv(3)), 1.0d-300)
         call verdict('the_standard_form_channels_are_its_own_He_II_row', &
                      dstd .le. 8.0d0*eps_d, dstd, 0.0d0, 8.0d0*eps_d)
      endif

      thereis_mol   = mol_save
      thereis_HeITR = tr_save
      end subroutine ledger_of_an_atomic_cell

      ! ----------------------------------------------------------------- !

      subroutine ledger_of_a_metal_bearing_cell(molecular)
      ! THE LEDGER WITH THE METAL CHARGE EXCHANGE PRESENT (item D7b).
      !
      ! A cell of a metal-bearing run, in each gas branch, assembled the way
      ! the transport operator assembles it: the kernel's rows, then the
      ! He <-> H pair the caller applies, then the metal charge exchange of
      ! Huang et al. (2023) Table 4 as stage sources.  cx_full is on, so
      ! group C (metal + He and He+) is active beside group A (metal + H and
      ! H+) and the group E electron capture, and all three stages of every
      ! element carry atoms so that no reaction of the set is silently zero.
      !
      ! What the row asserts is the identity: the signed sum of the
      ! channels of a stage is the source the transported row was assembled
      ! from.  Without the metal channels in the layout that sum would be
      ! short by exactly the metal terms, which is the RED this row replaces.
      logical, intent(in) :: molecular
      real(dp) :: st(nst_max), src(3), chan(n_stage_chan), dsrc(3)
      real(dp) :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
      real(dp) :: d1, d2, d3, mag
      logical  :: mol_save, tr_save, oxy_save, full_save
      character(len=80) :: nm

      mol_save  = thereis_mol
      tr_save   = thereis_HeITR
      oxy_save  = thereis_oxychem
      full_save = cx_full
      thereis_mol     = molecular
      thereis_HeITR   = .true.
      thereis_oxychem = .false.
      he_h_charge_exchange = .true.
      cx_full = .true.
      call cx_init()

      ! A partly ionized cell of the launch region at 8000 K, with the
      ! metals at their canonical order of magnitude relative to hydrogen.
      st = 0.0d0
      st( 1) = 8.0d3                 ! T [K]
      st( 2) = 6.0d8                 ! n(H I)
      st( 3) = 4.0d8                 ! n(H II)
      st( 4) = 7.0d7                 ! n(He I) ground singlet
      st( 5) = 5.0d2                 ! n(He 2^3S)
      st( 6) = 2.0d7                 ! n(He II)
      st( 7) = 1.0d6                 ! n(He III)
      st( 8) = 4.2d8                 ! n_e
      st( 9) = 1.0d-4                ! g(H I)
      st(10) = 2.0d-5                ! g(He I)
      st(11) = 1.0d-6                ! g(He II)
      st(12) = 1.0d-4                ! g(He 2^3S)
      st(13) = 1.6d-13               ! a(H II)
      st(14) = 1.0d-13               ! a(He II)
      st(15) = 5.0d-13               ! a(He III)
      st(16) = 5.0d-14               ! a(He 2^3S)
      st(17) = 2.0d-11               ! b(H I)
      st(18) = 1.0d-12               ! b(He I)
      st(19) = 1.0d-14               ! b(He II)
      st(20) = 3.0d-11               ! b(He 2^3S)
      st(21) = 1.0d-14               ! q13
      st(22) = 2.0d-12               ! q31a
      st(23) = 1.0d-12               ! q31b
      st(24) = 5.0d-10               ! Q31
      st(25) = 1.272d-4              ! A31
      st(26) = 1.0d-15               ! k(He0 + H+)
      st(27) = 1.0d-9                ! k(He+ + H0)
      st(28) = 1.1d9                 ! n_tot
      if (molecular) then
         st(29) = 1.0d8              ! n(H2)
         st(30) = 2.0d2              ! n(H2+)
         st(31) = 1.0d3              ! n(H3+)
         st(32) = 1.0d-2             ! n(HeH+)
         st(33) = 1.0d-7             ! g(H2), total destruction
         st(34) = 2.0d-8             ! g(H2) dissociative ionization
      endif
      ieq_cell%T_K   = st(1)
      ieq_cell%ntot  = st(28)
      ieq_cell%n_co  = 0.0d0
      ieq_cell%jcell = 0

      ! Every metal present and spread over its three stages, so that both
      ! directions of every group A, C and D row and the group E capture
      ! have their reactants.
      nm0 = 0.55d0*metal_totals()
      nm1 = 0.35d0*metal_totals()
      nm2 = 0.10d0*metal_totals()

      if (molecular) then
         call molecular_stage_sources(st, src, chan)
      else
         call atomic_stage_sources(st, src, chan)
      endif
      call metal_charge_exchange_channels(nm0, nm1, nm2, st(2), st(3),    &
                          st(4), st(6), st(7), st(1), chan, dsrc)
      src(1) = src(1) + dsrc(1)
      src(2) = src(2) + dsrc(2)
      src(3) = src(3) + dsrc(3)

      mag = max(abs(dsrc(1)), abs(dsrc(2)))
      call verdict('the_metal_charge_exchange_is_nonzero_here',           &
                   mag .gt. 0.0d0, mag, 0.0d0, 0.0d0)
      d1 = ledger_departure(chan, 1, src(1))
      d2 = ledger_departure(chan, 2, src(2))
      d3 = ledger_departure(chan, 3, src(3))
      if (molecular) then
         nm = 'the_channel_sums_are_the_assembled_source_metal_molecular'
      else
         nm = 'the_channel_sums_are_the_assembled_source_metal_atomic'
      endif
      call report_channels('metal-bearing cell', st, src, chan)
      call verdict(trim(nm), max(d1, d2, d3) .le. 8.0d0*eps_d,            &
                   max(d1, d2, d3), 0.0d0, 8.0d0*eps_d)

      thereis_mol     = mol_save
      thereis_HeITR   = tr_save
      thereis_oxychem = oxy_save
      cx_full         = full_save
      call cx_init()
      end subroutine ledger_of_a_metal_bearing_cell

      ! Element totals [cm^-3] of the canonical metal list at solar-order
      ! abundance relative to the 1e9 cm^-3 of hydrogen nuclei above, in the
      ! canonical order C, O, N, Mg, Si, Ca, Na, K, S, Fe.  They set the
      ! scale of the metal terms and nothing else asserted here depends on
      ! them.
      function metal_totals() result(nt)
      real(dp) :: nt(n_melem)
      nt = 1.0d9*[ 2.7d-4, 4.9d-4, 6.8d-5, 4.0d-5, 3.2d-5,                &
                   2.2d-6, 1.7d-6, 1.2d-7, 1.3d-5, 3.2d-5 ]
      end function metal_totals

      ! ----------------------------------------------------------------- !

      subroutine report_channels(label, st, src, chan)
      character(len=*), intent(in) :: label
      real(dp),         intent(in) :: st(*), src(3), chan(*)
      integer  :: k, irow
      real(dp) :: gp, gl
      write(*,'(A,A)') '  DIAGNOSTIC gross channels of the ', label
      do irow = 1, 3
         call stage_gross_of_channels(chan, irow, gp, gl)
         write(*,'(A,I2,A,ES13.6,A,ES13.6,A,ES13.6,A,ES9.2)')             &
              '    stage ', irow, ' production=', gp, ' loss=', gl,       &
              ' net=', src(irow), ' net/gross=',                          &
              abs(src(irow))/max(gp, gl, 1.0d-300)
      enddo
      do k = 1, n_stage_chan
         if (chan(k) .eq. 0.0d0) cycle
         write(*,'(A,I3,1X,A,A,I2,A,ES13.6)') '    channel ', k,          &
              stage_chan_name(k), ' stage ', stage_chan_row(k),           &
              ' rate=', stage_chan_sign(k)*chan(k)
      enddo
      end subroutine report_channels

      ! ----------------------------------------------------------------- !

      subroutine stage_sum_rounding_bound()
      ! THE FLOATING-POINT BOUND OF THE STAGE-NUCLEUS-SUM ENTRY (D3a step 4).
      !
      ! sum_k F_k(f) + F_close(f) = N_el(f) is an algebraic identity: it has
      ! no truncation error, so what a tolerance on it has to clear is the
      ! rounding of the sums and nothing else.  Counting the rounding events
      ! on the path from the face fractions to the measured difference
      ! (ionization_stage_face_fractions and ionization_stage_face_flux):
      !
      !   nk   forming xclose = 1 - sum_k x_k
      !    1   xclose * N
      !   nk   subtracting each eddy term from the closing flux
      !   2nk  each carried stage's product and its sum with the eddy term
      !   nk   summing the nk + 1 fluxes in the measure
      !    1   the final difference
      !
      ! that is 5 nk + 2 events, each bounded by eps times the magnitude of
      ! the quantity it rounds.  Every such quantity is bounded by
      !
      !   S = |N| + |F_close| + sum_k |F_k| + sum_k |E_k| ,
      !
      ! the eddy terms appearing because they are subtracted from the
      ! closing flux and added to the carried ones, so
      !
      !   |sum F - N|  <=  (5 nk + 2) eps S .
      !
      ! The entry divides by sc = max(|N|, |F_close| + sum_k|F_k|), which is
      ! at least half of S' = |N| + |F_close| + sum_k|F_k|, so with
      ! g = S/S' >= 1 the entry itself is bounded by
      !
      !   d  <=  2 (5 nk + 2) eps g .
      !
      ! With g = 1 that is 3.109e-15 for hydrogen (nk = 1) and 5.329e-15 for
      ! helium (nk = 2).  The rows below run the production face-flux
      ! routine on columns whose eddy term spans four decades, so that g is
      ! MEASURED and not assumed, and check the bound at each.
      integer, parameter :: ngrid = 5
      real(dp) :: kscale(ngrid)
      real(dp), allocatable :: Nel(:), xion(:,:), n_elf(:), Kf(:), drf(:)
      real(dp), allocatable :: xf(:,:), Fk(:,:), xclose(:), Fclose(:)
      real(dp) :: dr, Fsum, sc, d, S1, S2, g, ed, bnd
      real(dp) :: dworst(2), gworst(2), ratio(2)
      integer  :: ig, j, k, nk, ie

      dworst = 0.0d0
      gworst = 1.0d0
      ratio  = 0.0d0
      kscale = (/ 0.0d0, 1.0d9, 1.0d11, 1.0d13, 1.0d15 /)

      do ie = 1, 2
         nk = ie                       ! 1 = hydrogen, 2 = helium
         do ig = 1, ngrid
            N  = 400
            R0 = 1.0d10
            spherical_domain = .true.
            grid_type = 'Uniform'
            r_max     = 2.0d0
            use_plm = .false.;  use_weno3 = .true.;  rec_method = 'WENO3'
            recon_lambda_on = .false.
            j_min = 1
            if (allocated(r))     deallocate(r)
            if (allocated(r_edg)) deallocate(r_edg)
            if (allocated(dr_j))  deallocate(dr_j)
            allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
            dr = (r_max - 1.0d0)/dble(N)
            do j = 1-Ng, N+Ng
               r(j)     = 1.0d0 + (dble(j) - 0.5d0)*dr
               r_edg(j) = 1.0d0 +  dble(j)         *dr
               dr_j(j)  = dr
            enddo
            allocate(Nel(1-Ng:N+Ng), xion(nk,1-Ng:N+Ng))
            allocate(n_elf(0:N), Kf(0:N), drf(0:N))
            allocate(xf(nk,0:N), Fk(nk,0:N), xclose(0:N), Fclose(0:N))
            do j = 1-Ng, N+Ng
               xion(1,j) = 0.05d0 + 0.80d0*(r(j) - 1.0d0)
               if (nk .eq. 2) xion(2,j) = 1.0d-3                          &
                                        + 0.10d0*(r(j) - 1.0d0)**2
               Nel(j) = 1.0d13/(r_edg(j)*r_edg(j))
            enddo
            do j = 0, N
               n_elf(j) = 1.0d10*exp(-4.0d0*(r(j) - 1.0d0))
               Kf(j)    = kscale(ig)*r(j)*r(j)
               drf(j)   = (r(j+1) - r(j))*R0
            enddo
            call ionization_stage_face_flux(xion, Nel, n_elf, Kf, drf,    &
                                            xf, xclose, Fk, Fclose)
            do j = 0, N
               Fsum = Fclose(j)
               sc   = abs(Fclose(j))
               S1   = abs(Fclose(j))
               ed   = 0.0d0
               do k = 1, nk
                  Fsum = Fsum + Fk(k,j)
                  sc   = sc + abs(Fk(k,j))
                  S1   = S1 + abs(Fk(k,j))
                  ! The eddy term of this stage, re-formed from the face
                  ! quantities the routine was given; it is the one term
                  ! shared between a carried stage and the closing one.
                  if (j .gt. 0 .and. j .lt. N)                            &
                     ed = ed + abs(n_elf(j)*Kf(j)/drf(j)                  &
                                  *(xion(k,j+1) - xion(k,j)))
               enddo
               sc = max(abs(Nel(j)), sc)
               S1 = S1 + abs(Nel(j))
               if (sc .le. 0.0d0) cycle
               S2 = S1 + ed
               g  = S2/max(S1, 1.0d-300)
               d  = abs(Fsum - Nel(j))/sc
               ! The bound is a statement about EVERY face, so the row
               ! below is the largest ratio of an entry to its own face's
               ! bound, and not the largest entry against the bound of the
               ! face that happened to attain it.
               ratio(ie)  = max(ratio(ie),                                &
                                d/(2.0d0*dble(5*nk + 2)*eps_d*g))
               gworst(ie) = max(gworst(ie), g)
               dworst(ie) = max(dworst(ie), d)
            enddo
            write(*,'(A,I1,A,ES9.2,A,ES11.4,A,ES9.2,A,F8.4)')             &
                 '  DIAGNOSTIC stage sum, nk=', nk, '  K_0=', kscale(ig), &
                 '  worst entry=', dworst(ie), '  largest g=',            &
                 gworst(ie), '  worst entry/its own bound=', ratio(ie)
            deallocate(Nel, xion, n_elf, Kf, drf, xf, Fk, xclose, Fclose)
         enddo
      enddo

      call verdict('the_stage_sum_identity_is_within_its_rounding_'//     &
                   'bound_hydrogen', ratio(1) .le. 1.0d0, ratio(1),       &
                   1.0d0, 1.0d0)
      call verdict('the_stage_sum_identity_is_within_its_rounding_'//     &
                   'bound_helium', ratio(2) .le. 1.0d0, ratio(2),         &
                   1.0d0, 1.0d0)
      bnd = 0.0d0
      write(*,'(A,ES11.4,A,ES11.4)')                                      &
           '  DIAGNOSTIC stage-sum bound at g = 1: hydrogen ',            &
           2.0d0*7.0d0*eps_d, '  helium ', 2.0d0*12.0d0*eps_d
      end subroutine stage_sum_rounding_bound

      ! ----------------------------------------------------------------- !

      subroutine ledger_and_measurement_of_a_record()
      ! THE RECORD A RUN LEFT BEHIND (EXHALE_STAGE_CHANNELS), re-evaluated
      ! through the production kernels at the state it carries.
      character(len=512) :: recfile, trmfile
      character(len=8192):: line
      real(dp) :: st(nst_max), chan(n_stage_chan), chanr(n_stage_chan)
      real(dp) :: src(3), dmax_chan, dmax_led, dled
      real(dp), allocatable :: strec(:,:), chrec(:,:)
      integer,  allocatable :: cellrec(:), nstrec(:), callrec(:)
      character(len=2), allocatable :: gasrec(:)
      integer  :: u, ios, nrec, i, k, nst, jc, icall, ibase
      character(len=2) :: gas
      logical  :: have

      call get_environment_variable('EXHALE_STAGE_ROW_RECORD', recfile)
      if (len_trim(recfile) .eq. 0) then
         write(*,'(A)') '  DIAGNOSTIC no EXHALE_STAGE_ROW_RECORD: the '// &
              'fixture measurement is not taken in this run'
         return
      endif
      call get_environment_variable('EXHALE_STAGE_ROW_TERMS', trmfile)

      allocate(strec(nst_max, nrec_max), chrec(n_stage_chan, nrec_max),  &
               cellrec(nrec_max), nstrec(nrec_max), callrec(nrec_max),   &
               gasrec(nrec_max))
      strec = 0.0d0
      chrec = 0.0d0
      open(newunit=u, file=trim(recfile), status='old', action='read',    &
           iostat=ios)
      if (ios .ne. 0) then
         call verdict('the_record_named_by_the_environment_is_readable',  &
                      .false., 1.0d0, 0.0d0, 0.0d0)
         return
      endif
      nrec = 0
      do
         read(u,'(A)', iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#') cycle
         if (len_trim(line) .eq. 0) cycle
         if (nrec .ge. nrec_max) exit
         nrec = nrec + 1
         read(line,*,iostat=ios) gas, jc, icall, nst
         if (ios .ne. 0) then
            nrec = nrec - 1
            cycle
         endif
         read(line,*,iostat=ios) gas, jc, icall, nst,                     &
              (strec(k,nrec), k = 1, nst),                                &
              (chrec(k,nrec), k = 1, n_stage_chan)
         if (ios .ne. 0) then
            nrec = nrec - 1
            cycle
         endif
         gasrec(nrec)  = gas
         cellrec(nrec) = jc
         callrec(nrec) = icall
         nstrec(nrec)  = nst
      enddo
      close(u)
      write(*,'(A,I0,A,A)') '  DIAGNOSTIC ', nrec,                        &
           ' record(s) read from ', trim(recfile)
      if (nrec .eq. 0) then
         call verdict('the_record_named_by_the_environment_is_readable',  &
                      .false., 0.0d0, 1.0d0, 0.0d0)
         return
      endif

      ! The rows re-called at the recorded state give the recorded
      ! channels; the signed channel sums give the assembled source.
      dmax_chan = 0.0d0
      dmax_led  = 0.0d0
      do i = 1, nrec
         st = 0.0d0
         st(1:nstrec(i)) = strec(1:nstrec(i), i)
         ieq_cell%T_K   = st(1)
         ieq_cell%ntot  = st(28)
         ieq_cell%n_co  = 0.0d0
         ieq_cell%jcell = 0
         he_h_charge_exchange = .true.
         if (gasrec(i)(1:1) .eq. 'M') then
            thereis_mol   = .true.
            thereis_HeITR = (st(5) .gt. 0.0d0 .or. st(12) .gt. 0.0d0)
            thereis_oxychem = .false.
            call molecular_stage_sources(st, src, chan)
         else
            thereis_mol   = .false.
            thereis_HeITR = (st(5) .gt. 0.0d0 .or. st(12) .gt. 0.0d0)
            call atomic_stage_sources(st, src, chan)
         endif
         do k = 1, n_stage_chan
            if (k .eq. ich_HII_cx_gain  .or. k .eq. ich_HII_cx_loss .or.  &
                k .eq. ich_HeII_cx_gain .or. k .eq. ich_HeII_cx_loss) cycle
            dmax_chan = max(dmax_chan, abs(chan(k) - chrec(k,i))          &
                            /max(abs(chrec(k,i)), 1.0d-300))
         enddo
         do k = 1, 3
            dled = ledger_departure(chan, k, src(k))
            dmax_led = max(dmax_led, dled)
         enddo
      enddo
      call verdict('the_recorded_channels_are_the_rows_of_the_'//         &
                   'recorded_state', dmax_chan .eq. 0.0d0, dmax_chan,     &
                   0.0d0, 0.0d0)
      call verdict('the_recorded_channel_sums_are_the_recorded_source',   &
                   dmax_led .le. 8.0d0*eps_d, dmax_led, 0.0d0,            &
                   8.0d0*eps_d)

      ! THE SAME COMPOSITION AND THE SAME RADIATION FIELD.  Every
      ! evaluation of the rows at one cell -- the local sweep's own and the
      ! transport operator's -- reads its rate coefficients from that
      ! cell's stored state, so the two must agree BITWISE in the
      ! temperature, the four photoionization rates, the recombination and
      ! collisional coefficients, the He 2^3S rates, the two
      ! charge-exchange coefficients and the molecular photo rates.  The
      ! densities are NOT expected to agree: the operator evaluates the row
      ! at the transported trial, which is the difference this item
      ! measures.
      call fields_across_the_records(nrec, gasrec, cellrec, nstrec, strec)

      ! The measurement itself, on the FIRST record of each cell, which is
      ! the unperturbed evaluation of the assembly (the ones that follow it
      ! are the difference probes of the Jacobian).
      do i = 1, nrec
         have = .true.
         do k = 1, i-1
            if (cellrec(k) .eq. cellrec(i)) have = .false.
         enddo
         if (.not. have) cycle
         ! The TRANSPORT OPERATOR's own evaluation of this cell, which is
         ! the one the certification measured; where the record cannot tell
         ! the call sites apart (the atomic rows) the first evaluation is
         ! the operator's, because nothing else evaluates the rows in a
         ! stationary evaluation.
         ibase = i
         do k = i, nrec
            if (cellrec(k) .ne. cellrec(i)) cycle
            if (gasrec(k)(2:2) .eq. 'C') then
               ibase = k
               exit
            endif
         enddo
         call measure_one_cell(gasrec(ibase), cellrec(ibase),             &
                               strec(:,ibase), trmfile)
      enddo
      end subroutine ledger_and_measurement_of_a_record


      ! ----------------------------------------------------------------- !

      subroutine fields_across_the_records(nrec, gasrec, cellrec, nstrec, &
                                           strec)
      integer,          intent(in) :: nrec
      character(len=2), intent(in) :: gasrec(*)
      integer,          intent(in) :: cellrec(*), nstrec(*)
      real(dp),         intent(in) :: strec(nst_max,*)
      ! The entries that are rate coefficients or the temperature: the
      ! cell's stored state, which the transport operator reads from the
      ! sweep's own cell state and may not change.
      integer, parameter :: nfix = 25
      integer, parameter :: ifix(nfix) = (/ 1, 9,10,11,12,13,14,15,16,    &
                                           17,18,19,20,21,22,23,24,25,    &
                                           26,27, 33,34,35,36,37 /)
      real(dp) :: dfix, dden, dv
      integer  :: i, k, m, isw, npair
      dfix  = 0.0d0
      dden  = 0.0d0
      npair = 0
      do i = 1, nrec
         if (gasrec(i)(2:2) .ne. 'C') cycle
         ! The last evaluation the LOCAL SWEEP made at this cell before the
         ! operator's; the two must read one radiation field.
         isw = 0
         do k = i-1, 1, -1
            if (cellrec(k) .ne. cellrec(i)) cycle
            if (gasrec(k)(2:2) .eq. 'S') then
               isw = k
               exit
            endif
         enddo
         if (isw .eq. 0) cycle
         npair = npair + 1
         do m = 1, nfix
            k = ifix(m)
            if (k .gt. nstrec(i)) cycle
            dv = abs(strec(k,i) - strec(k,isw))                           &
                /max(abs(strec(k,isw)), 1.0d-300)
            dfix = max(dfix, dv)
         enddo
         do k = 2, 8
            dv = abs(strec(k,i) - strec(k,isw))                           &
                /max(abs(strec(k,isw)), 1.0d-300)
            dden = max(dden, dv)
         enddo
      enddo
      if (npair .eq. 0) then
         write(*,'(A)') '  DIAGNOSTIC the record holds no sweep and '//   &
              'operator evaluation of one cell to compare (a stationary'//&
              ' evaluation runs the operator alone)'
         return
      endif
      call verdict('the_sweep_and_the_transported_source_read_one_'//     &
                   'radiation_field', dfix .eq. 0.0d0, dfix, 0.0d0, 0.0d0)
      write(*,'(A,I0,A,ES12.5)') '  DIAGNOSTIC over ', npair,             &
           ' sweep/operator pairs the largest relative difference of '//  &
           'the densities is ', dden
      end subroutine fields_across_the_records

      ! ----------------------------------------------------------------- !

      subroutine measure_one_cell(gas, jc, st_in, trmfile)
      ! D3a STEP 2 AT ONE CELL.
      character(len=2), intent(in) :: gas
      integer,          intent(in) :: jc
      real(dp),         intent(in) :: st_in(nst_max)
      character(len=*), intent(in) :: trmfile
      real(dp) :: st(nst_max), src(3), chan(n_stage_chan)
      real(dp) :: stp(nst_max), srcp(3), chanp(n_stage_chan)
      real(dp) :: gp(3), gl(3), nHe, nHnuc, xstage(3)
      real(dp) :: dif(3), adv(3), res(3), scl(3), prod(3), loss(3), flo(3)
      real(dp) :: nc(3), dsdx, dsdne, dlt, sref, sfloor, qsum
      real(dp) :: quadsum, corr, D, delta(6)
      real(dp) :: dtab(6), stab(6)
      integer  :: irow, k, id
      logical  :: haveterms

      st = 0.0d0
      st(1:nst_max) = st_in(1:nst_max)
      ieq_cell%T_K   = st(1)
      ieq_cell%ntot  = st(28)
      ieq_cell%n_co  = 0.0d0
      ieq_cell%jcell = 0
      he_h_charge_exchange = .true.
      if (gas(1:1) .eq. 'M') then
         thereis_mol = .true.
         thereis_oxychem = .false.
      else
         thereis_mol = .false.
      endif
      thereis_HeITR = (st(5) .gt. 0.0d0 .or. st(12) .gt. 0.0d0)
      if (gas(1:1) .eq. 'M') then
         call molecular_stage_sources(st, src, chan)
      else
         call atomic_stage_sources(st, src, chan)
      endif
      do irow = 1, 3
         call stage_gross_of_channels(chan, irow, gp(irow), gl(irow))
      enddo
      nHnuc = st(2) + st(3) + 2.0d0*st(29) + 2.0d0*st(30) + 3.0d0*st(31)  &
            + st(32)
      nHe   = st(4) + st(5) + st(6) + st(7) + st(32)
      xstage(1) = st(3)/max(nHnuc, 1.0d-300)
      xstage(2) = st(6)/max(nHe,   1.0d-300)
      xstage(3) = st(7)/max(nHe,   1.0d-300)

      call read_row_terms(trmfile, jc, nc, dif, adv, prod, loss, res,     &
                          flo, scl, haveterms)

      write(*,'(A)') ' '
      write(*,'(A,A2,A,I0,A,ES12.5,A)') '  DIAGNOSTIC cell (gas ', gas,   &
           ') ', jc, ', T = ', st(1), ' K'
      write(*,'(A)') '    stage   x            P            L       '//   &
           '     P-L          |P-L|/max(P,L)'
      do irow = 1, 3
         write(*,'(A,I1,5(1X,ES13.6))') '      ', irow, xstage(irow),     &
              gp(irow), gl(irow), src(irow),                              &
              abs(src(irow))/max(gp(irow), gl(irow), 1.0d-300)
      enddo
      if (haveterms) then
         write(*,'(A)') '    stage   diffusive     advective    '//       &
              'net source   residual     floor        scale        measure'
         do irow = 1, 3
            write(*,'(A,I1,7(1X,ES13.6))') '      ', irow, dif(irow),     &
                 adv(irow), prod(irow) - loss(irow), res(irow),           &
                 flo(irow), scl(irow),                                    &
                 abs(res(irow))/max(scl(irow), 1.0d-300)
         enddo
         write(*,'(A)') '    the row scale on the gross channels '//      &
              'instead of the net source, stage by stage:'
         do irow = 1, 3
            write(*,'(A,I1,3(1X,ES13.6))') '      ', irow,                &
                 scl(irow) - abs(src(irow)) + gp(irow) + gl(irow),        &
                 abs(res(irow))                                           &
                 /max(scl(irow) - abs(src(irow)) + gp(irow) + gl(irow),   &
                      1.0d-300),                                          &
                 (gp(irow) + gl(irow))/max(abs(src(irow)), 1.0d-300)
         enddo
      endif

      ! THE CANCELLATION FLOOR OF THE NET SOURCE.  The kernels are written
      ! in double precision, so a quadruple-precision evaluation of the
      ! rows is not available without a second spelling of every rate; what
      ! IS available is the quadruple-precision SUM of the channels the
      ! kernel returned, which isolates the rounding of the assembly, and
      ! the term-magnitude estimate of the rounding inside each channel.
      ! Each channel is a product of two or three doubles, so its own
      ! relative error is at most 2 eps, and the sum of m of them carries at
      ! most (m-1) eps of their magnitudes.
      do irow = 1, 3
         quadsum = 0.0d0
         qsum    = 0.0d0
         k = 0
         do id = 1, n_stage_chan
            if (stage_chan_row(id) .ne. irow) cycle
            if (chan(id) .eq. 0.0d0) cycle
            k = k + 1
            quadsum = quadsum + real(stage_chan_sign(id), qp)             &
                                *real(chan(id), qp)
            qsum    = qsum + abs(chan(id))
         enddo
         sfloor = dble(k + 2)*eps_d*qsum
         write(*,'(A,I1,A,ES12.5,A,ES12.5,A,ES12.5,A,ES9.2)')             &
              '    stage ', irow,                                         &
              ' cancellation: |net|=', abs(src(irow)),                    &
              ' floor=', sfloor,                                          &
              ' quad-sum departure=', abs(dble(quadsum) - src(irow)),     &
              ' |net|/floor=', abs(src(irow))/max(sfloor, 1.0d-300)
      enddo

      ! THE SENSITIVITY OF EACH STAGE'S SOURCE TO ITS OWN FRACTION AND TO
      ! THE ELECTRON DENSITY, and the abundance correction the transport
      ! would need.  The perturbation respects the budgets: a helium
      ! nucleus moved into He II comes out of the ground singlet and takes
      ! its electron with it.
      delta = (/ 1.0d-2, 1.0d-3, 1.0d-4, 1.0d-5, 1.0d-6, 1.0d-8 /)
      do irow = 1, 3
         do id = 1, 6
            dlt = delta(id)
            stp = st
            ! A BUDGET-RESPECTING PERTURBATION of the stage this row owns:
            ! the nucleus moved into the stage comes out of the stage below
            ! it in the same element, and its electron is added to n_e.
            if (irow .eq. 1) then
               stp(3) = st(3)*(1.0d0 + dlt)
               stp(2) = st(2) - st(3)*dlt
               stp(8) = st(8) + st(3)*dlt
            else if (irow .eq. 2) then
               stp(6) = st(6)*(1.0d0 + dlt)
               stp(4) = st(4) - st(6)*dlt
               stp(8) = st(8) + st(6)*dlt
            else
               stp(7) = st(7)*(1.0d0 + dlt)
               stp(6) = st(6) - st(7)*dlt
               stp(8) = st(8) + st(7)*dlt
            endif
            if (gas(1:1) .eq. 'M') then
               call molecular_stage_sources(stp, srcp, chanp)
            else
               call atomic_stage_sources(stp, srcp, chanp)
            endif
            dtab(id) = dlt
            stab(id) = (srcp(irow) - src(irow))                           &
                     /max(xstage(irow)*dlt, 1.0d-300)
         enddo
         dsdx = stab(3)
         write(*,'(A,I1,A)') '    stage ', irow,                          &
              ' d(P-L)/dx by relative step (budget-respecting):'
         do id = 1, 6
            write(*,'(A,ES9.2,A,ES13.6)') '      step=', dtab(id),        &
                 '  d(P-L)/dx=', stab(id)
         enddo
         stp = st
         stp(8) = st(8)*(1.0d0 + 1.0d-4)
         if (gas(1:1) .eq. 'M') then
            call molecular_stage_sources(stp, srcp, chanp)
         else
            call atomic_stage_sources(stp, srcp, chanp)
         endif
         dsdne = (srcp(irow) - src(irow))/1.0d-4
         write(*,'(A,ES13.6,A,ES13.6)')                                   &
              '      d(P-L)/dln(n_e)=', dsdne, '   d(P-L)/dx=', dsdx
         if (haveterms) then
            D    = dif(irow) + adv(irow)
            corr = D/sign(max(abs(dsdx), 1.0d-300), dsdx)
            write(*,'(A,ES13.6,A,ES13.6,A,ES13.6)')                       &
                 '      transport divergence D=', D,                      &
                 '   abundance correction D/(d(P-L)/dx)=', corr,          &
                 '   residual correction R/(d(P-L)/dx)=',                 &
                 res(irow)/sign(max(abs(dsdx), 1.0d-300), dsdx)
            write(*,'(A,ES13.6,A,ES13.6)')                                &
                 '      D correction as a fraction of x: ',              &
                 corr/max(xstage(irow), 1.0d-300),                        &
                 '   R correction as a fraction of x: ',                  &
                 res(irow)/sign(max(abs(dsdx), 1.0d-300), dsdx)           &
                 /max(xstage(irow), 1.0d-300)
         endif
      enddo
      end subroutine measure_one_cell

      ! ----------------------------------------------------------------- !

      subroutine read_row_terms(fname, jc, nc, dif, adv, prod, loss, res, &
                                flo, scl, have)
      ! The carrier row record of the same state, for the three stage rows
      ! of the named cell.  The file is the run's own output, written by
      ! carrier_row_terms_write.
      character(len=*), intent(in)  :: fname
      integer,          intent(in)  :: jc
      real(dp),         intent(out) :: nc(3), dif(3), adv(3), prod(3)
      real(dp),         intent(out) :: loss(3), res(3), flo(3), scl(3)
      logical,          intent(out) :: have
      character(len=1024) :: line
      character(len=8)    :: cname
      real(dp) :: rj, tj, ncj, xcj, dj, aj, pj, lj, phj, nsj, rrj, flj
      real(dp) :: scj, msj
      integer  :: u, ios, jj, irow
      have = .false.
      nc = 0.0d0;  dif = 0.0d0;  adv = 0.0d0;  prod = 0.0d0
      loss = 0.0d0;  res = 0.0d0;  flo = 0.0d0;  scl = 0.0d0
      if (len_trim(fname) .eq. 0) return
      open(newunit=u, file=trim(fname), status='old', action='read',      &
           iostat=ios)
      if (ios .ne. 0) return
      do
         read(u,'(A)', iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#') cycle
         ! The face and channel lines of the same record, which carry a
         ! label where a cell index stands on a row line.
         if (line(1:6) .eq. '  face') cycle
         if (line(1:6) .eq. '  chan') cycle
         read(line,'(I6,1X,ES14.7,1X,ES12.5,1X,A8,12(1X,ES14.7))',        &
              iostat=ios) jj, rj, tj, cname, ncj, xcj, dj, aj, pj, lj,    &
              phj, nsj, rrj, flj, scj, msj
         if (ios .ne. 0) cycle
         if (jj .ne. jc) cycle
         irow = 0
         if (trim(adjustl(cname)) .eq. 'H+')   irow = 1
         if (trim(adjustl(cname)) .eq. 'He+')  irow = 2
         if (trim(adjustl(cname)) .eq. 'He++') irow = 3
         if (irow .eq. 0) cycle
         nc(irow)   = ncj
         dif(irow)  = dj
         adv(irow)  = aj
         prod(irow) = pj
         loss(irow) = lj
         res(irow)  = rrj
         flo(irow)  = flj
         scl(irow)  = scj
         have = .true.
      enddo
      close(u)
      end subroutine read_row_terms

      end program stage_row_balance_tests
