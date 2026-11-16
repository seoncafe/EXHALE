      module field_spectral_basis_reuse
      ! REUSE OF THE DENSITY-FREE SPECTRAL BASIS OF THE ATTENUATED STELLAR
      ! FIELD ACROSS THE PASSES OF ONE COUPLED-SOURCE ENTRY.
      !
      ! In the coupled composition-temperature loop of a marching step
      ! (EXHALE_main, coupled_source) every pass forms the attenuated XUV
      ! field of every cell from the composition the pass starts from.
      ! Between passes the columns and densities change little. For an
      ! interior cell, pass k >= 2 may take the density-free spectral
      ! integrals (the basis, spectral_basis_t of utils_ion_eq:
      ! spectral_basis_of_cell) from the pass that formed them; the
      ! photoelectron partition, the density assembly, the heating channels
      ! and the efficiency are formed from the current arguments in every
      ! call (photoionization_rates_of_spectral_basis, the same routine the
      ! recomputation uses), so only the spectral weights carry the lag.
      !
      ! THE BOUND. The integrated field of the cell is F_XUV opa_pf A, with
      ! the cell mean of the attenuation A = exp(-t_out) (1 - exp(-t))/t,
      ! t_out the star-ward and t the cell's own optical depth of the bin.
      ! d ln A / d t_out = -1, and d ln[(1 - exp(-t))/t]/dt = 1/(exp(t) - 1)
      ! - 1/t, whose magnitude is at most 1/2 for t >= 0 (its limit at t = 0
      ! is -1/2). With t = sigma n w, w = dr_j R0 opa_pf the opacity-weighted
      ! width, at fixed grid geometry a change of opa_pf is the change of w,
      ! and |d ln w| bounds the change of the prefactor opa_pf. So, for every
      ! bin,
      !   |d ln(opa_pf A)| <= |d t_out| + |d t|/2 + |d ln w| <= M_hat,
      !   M_hat = 1e-18 [ sum_a sigma_max,a |dN_a| + (1/2) sum_a sigma_max,a
      !           |d(n_a w)| ] + |d ln w|,
      ! with N_a the star-ward column and n_a the own density of absorber a
      ! and sigma_max,a its largest cross section (scalar_attenuation_bound).
      ! A cell reuses only when M_hat <= ln(1 + eps_F): then every spectral
      ! weight of the stored basis lies within a factor 1 + eps_F of the
      ! current one. This is an exact-arithmetic argument for nonnegative
      ! depths and a positive width, with the grid, the stellar spectrum and
      ! the tables fixed within the entry. The pointwise bound concerns the
      ! spectral weights; it does not imply a relative eps_F accuracy of a
      ! signed combination (the H2 heating, a difference of channels), of a
      ! quantity with a vanishing denominator, or of the accepted state of
      ! the step. A generating spectrum with an underflowed bin (int_f <
      ! tiny where F_XUV > 0) is never reused; a bin that was normal when the
      ! basis was formed can fall below tiny only from within a factor
      ! 1 + eps_F of it, so its contribution is of order tiny (absolute
      ! smallness, not a relative floating-point guarantee).
      !
      ! VALIDITY. Reuse only when all hold; everything else recomputes:
      !   - a coupled pass >= 2 of the same coupled-loop entry; a new entry
      !     (a new step, a refused, restored or sampled attempt) starts from
      !     a recomputed pass 1; the basis must have been formed by an
      !     earlier pass of this entry and is never rewritten by a reuse;
      !   - an interior cell (ghost cells recompute) in the first ghost
      !     fixed-point application of the sweep; a sweep with more than one
      !     application recomputes in its later applications and invalidates
      !     every basis at its end;
      !   - neither the fixed-point probe (csm_err_probe) nor the geometric
      !     extrapolation (csm_extrap_on) is enabled for the entry, nor is the
      !     pass a probe pass;
      !   - a_tau = 0 (the bound is derived for that attenuation) and no
      !     photon-fraction heating;
      !   - the same branch as the generating call (photoelectron partition
      !     on or off, H2 present, Pm2 wanted); an active low-energy
      !     partition is never stored and always recomputes;
      !   - no underflowed bin and no nonfinite value in the generating
      !     basis; finite current inputs;
      !   - M_hat <= ln(1 + eps_F), eps_F = spectral_basis_reuse_eps.
      ! ASSUMED FIXED WITHIN AN ENTRY (not part of the comparison): the
      ! stellar spectrum and energy grid, the cross-section and spectral
      ! tables, the grid geometry other than the width, the species set and
      ! the switches other than those the basis records; the marching loop
      ! changes none of them inside an entry.
      !
      ! VERIFIED SCOPE (update log stage 3, section 189). eps_F = 1e-9
      ! against the recomputation of every pass, for the four LHS 1140 b
      ! states He/H 2.09, 1.2, 14 and 0.1 (ionization transport, molecular
      ! carriers, metals; CFL 0.15, the C2P estimator) to 2100 s: the final
      ! conditional verdicts of the frozen comparison are PASS in all four
      ! cases (He/H 1.2 by the frozen addendum with the h0/16 level), with
      ! the matched-time field lag at most 3e-5 of the allocation 0.25 tau
      ! and the inner allowance widened by the fixed-step difference of the
      ! two maps; held-state, short-run and event tests cover the probe, the
      ! extrapolation, a_tau > 0, photon-fraction heating, the active
      ! low-energy partition, repeated ghost applications, refused entries
      ! and an atomic and a molecular configuration without ionization
      ! transport. Not exercised: the self-field path (xuv_self_field_passes
      ! > 1), the H-only field, and the imposed-stage rounding repair at
      ! field level.
      !
      ! spectral_basis_reuse_on = .false. builds the recomputation control
      ! (the same dynamic schedule of the field loop, no reuse).
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use global_parameters
      use species_table, only: n_mion, n_mphot, mion_isphot, mion_iphot
      use utils_ion_eq, only: photoionization_field_at_cell_HHe,            &
                              photoionization_field_from_spectral_basis,    &
                              spectral_basis_t
      implicit none
      private
      public :: spectral_basis_reuse_on, spectral_basis_reuse_eps
      public :: spectral_basis_reuse_begin_pass, spectral_basis_reuse_end_pass
      public :: spectral_basis_reuse_application, spectral_basis_reuse_application_end
      public :: spectral_basis_reuse_applications
      public :: field_with_spectral_basis_reuse, spectral_basis_reuse_report
      public :: spectral_basis_reuse_counts, sbr_ncode, sbr_code_name

      logical, parameter :: spectral_basis_reuse_on = .true.
      ! eps_F: the largest admitted change of a spectral weight, as a factor
      ! 1 + eps_F (the verified value).
      real*8, parameter :: spectral_basis_reuse_eps = 1.0d-9

      ! The comparison inputs of a field call: columns (5 + n_mphot), own
      ! densities (5 + n_mion), xion, f_vib_quench, e_vib_bound, width.
      integer, parameter :: sbr_iNcol = 5, sbr_in = 5 + n_mphot,          &
                            sbr_inm = sbr_in + 5,                        &
                            sbr_ix = sbr_inm + n_mion
      integer, parameter :: sbr_nin = sbr_ix + 4

      ! The basis of each interior cell, the identity of the call that formed
      ! it (coupled-loop entry, pass) and that call's inputs.
      type(spectral_basis_t), allocatable :: sbr_basis(:)
      integer, allocatable :: sbr_basis_entry(:), sbr_basis_pass(:)
      real*8,  allocatable :: sbr_basis_inp(:,:)
      ! Largest cross section of each absorber [1e-18 cm^2]: H I, He I,
      ! He II, He 2^3S, H2, and every metal photoionization column.
      real*8 :: sbr_shat(5) = 0.0d0
      real*8, allocatable :: sbr_shat_m(:)
      logical :: sbr_ready = .false.
      logical :: sbr_in_pass = .false.
      logical :: sbr_entry_reuse = .true.
      integer :: sbr_entry = 0, sbr_pass = 0
      integer :: sbr_cur_app = 0, sbr_apps = 0

      ! The decision of each field request of the HHe field loop, for the
      ! interior cell j of the application being run (cell-owned in the
      ! parallel loop; counted serially after the loop): 1 reused,
      ! 2 reused with numerically equal inputs; recomputed: 3 first pass of
      ! the entry, 4 no basis of this entry, 5 later ghost application,
      ! 6 option (a_tau > 0 or photon-fraction heating), 7 branch (active
      ! low-energy partition or basis of another branch), 8 underflow in the
      ! generating basis, 9 M_hat above ln(1 + eps_F), 10 nonfinite input,
      ! 11 ghost cell, 12 entry with the probe or the extrapolation, 13 a
      ! request outside a coupled pass or in a probe pass.
      integer, parameter :: sbr_ncode = 13
      character(len=30), parameter :: sbr_code_name(sbr_ncode) = [character(len=30) :: &
         'reused', 'reused_exact_inputs', 'recomputed_first_pass', 'recomputed_no_basis', &
         'recomputed_ghost_application', 'recomputed_option', 'recomputed_branch',         &
         'recomputed_underflow', 'recomputed_indicator', 'recomputed_nonfinite',           &
         'recomputed_ghost_cell', 'recomputed_probe_or_extrap',                            &
         'recomputed_outside_pass']
      integer, allocatable :: sbr_code(:)
      integer*8 :: sbr_count(sbr_ncode) = 0

      contains

      subroutine initialize_spectral_basis_reuse()
      integer :: k
      allocate(sbr_basis(N), sbr_basis_entry(N), sbr_basis_pass(N),          &
               sbr_basis_inp(sbr_nin, N), sbr_code(N))
      sbr_basis_entry = -1
      sbr_basis_pass  = -1
      sbr_basis_inp   = 0.0d0
      sbr_code        = 0
      sbr_shat(1) = maxval(s_hi)
      sbr_shat(2) = maxval(s_hei)
      sbr_shat(3) = maxval(s_heii)
      if (allocated(s_heiTR)) sbr_shat(4) = maxval(s_heiTR)
      if (allocated(s_h2))    sbr_shat(5) = maxval(s_h2)
      allocate(sbr_shat_m(n_mphot))
      sbr_shat_m = 0.0d0
      if (allocated(sigma_tab)) then
         do k = 1, min(n_mphot, size(sigma_tab, 2))
            sbr_shat_m(k) = maxval(sigma_tab(:, k))
         enddo
      endif
      sbr_ready = .true.
      end subroutine initialize_spectral_basis_reuse


      subroutine spectral_basis_reuse_begin_pass(i_csm, probe, reuse_allowed)
      ! Before the sweep of coupled pass i_csm (EXHALE_main). Pass 1 opens a
      ! new entry; reuse_allowed (false when the fixed-point probe or the
      ! geometric extrapolation is enabled) holds for the whole entry; a
      ! probe pass recomputes every cell.
      integer, intent(in) :: i_csm
      logical, intent(in) :: probe, reuse_allowed
      if (.not. sbr_ready) call initialize_spectral_basis_reuse()
      if (i_csm .eq. 1) then
         sbr_entry_reuse = reuse_allowed
         sbr_entry = sbr_entry + 1
      endif
      sbr_pass = i_csm
      sbr_apps = 0
      sbr_in_pass = .not. probe
      end subroutine spectral_basis_reuse_begin_pass


      subroutine spectral_basis_reuse_application(i_app)
      ! Before an HHe field loop of the sweep (ioniz_eq, serial): the ghost
      ! fixed-point application it belongs to.
      integer, intent(in) :: i_app
      if (.not. sbr_ready) call initialize_spectral_basis_reuse()
      sbr_cur_app = i_app
      sbr_code = 0
      end subroutine spectral_basis_reuse_application


      subroutine spectral_basis_reuse_application_end(jlo, jhi)
      ! After an HHe field loop over the cells jlo..jhi (ioniz_eq, serial):
      ! every request of that loop is counted, by its decision; the ghost
      ! cells of the range by code 11.
      integer, intent(in) :: jlo, jhi
      integer :: j
      do j = max(jlo, 1), min(jhi, N)
         if (sbr_code(j) .ge. 1) sbr_count(sbr_code(j)) = sbr_count(sbr_code(j)) + 1
      enddo
      sbr_count(11) = sbr_count(11) + max(0, min(jhi, 0) - jlo + 1)           &
                                    + max(0, jhi - max(jlo, N + 1) + 1)
      end subroutine spectral_basis_reuse_application_end


      subroutine spectral_basis_reuse_applications(n_app)
      ! The number of ghost fixed-point applications the sweep made.
      integer, intent(in) :: n_app
      if (sbr_in_pass) sbr_apps = n_app
      end subroutine spectral_basis_reuse_applications


      subroutine spectral_basis_reuse_end_pass()
      ! After the sweep of the pass (serial): a pass with more than one
      ! ghost application leaves no basis to be used again.
      if (.not. sbr_in_pass) return
      sbr_in_pass = .false.
      if (sbr_apps .ne. 1) sbr_basis(:)%valid = .false.
      end subroutine spectral_basis_reuse_end_pass


      pure real*8 function scalar_attenuation_bound(inp, inp0, has_h2) result(mh)
      ! M_hat between the inputs inp and the generating inputs inp0 of a
      ! basis (module header).
      real*8, intent(in) :: inp(sbr_nin), inp0(sbr_nin)
      logical, intent(in) :: has_h2
      real*8 :: s, d, wn, wo, dlnw
      integer :: k, i
      wn = inp(sbr_ix+4)
      wo = inp0(sbr_ix+4)
      dlnw = 0.0d0
      if (wn .gt. 0.0d0 .and. wo .gt. 0.0d0) dlnw = abs(log(wn/wo))
      s = sbr_shat(1)*abs(inp(1) - inp0(1)) + sbr_shat(2)*abs(inp(2) - inp0(2)) &
        + sbr_shat(3)*abs(inp(3) - inp0(3))
      if (thereis_He .and. thereis_HeITR) s = s + sbr_shat(4)*abs(inp(4) - inp0(4))
      if (has_h2) s = s + sbr_shat(5)*abs(inp(5) - inp0(5))
      do k = 1, n_mphot
         s = s + sbr_shat_m(k)*abs(inp(sbr_iNcol+k) - inp0(sbr_iNcol+k))
      enddo
      d = sbr_shat(1)*abs(inp(sbr_in+1)*wn - inp0(sbr_in+1)*wo)              &
        + sbr_shat(2)*abs(inp(sbr_in+2)*wn - inp0(sbr_in+2)*wo)              &
        + sbr_shat(3)*abs(inp(sbr_in+3)*wn - inp0(sbr_in+3)*wo)
      if (thereis_He .and. thereis_HeITR)                                     &
         d = d + sbr_shat(4)*abs(inp(sbr_in+4)*wn - inp0(sbr_in+4)*wo)
      if (has_h2) d = d + sbr_shat(5)*abs(inp(sbr_in+5)*wn - inp0(sbr_in+5)*wo)
      do i = 1, n_mion
         if (.not. mion_isphot(i)) cycle
         k = mion_iphot(i)
         d = d + sbr_shat_m(k)*abs(inp(sbr_inm+i)*wn - inp0(sbr_inm+i)*wo)
      enddo
      mh = 1.0d-18*(s + 0.5d0*d) + dlnw
      end function scalar_attenuation_bound


      subroutine field_with_spectral_basis_reuse(j, N1o, N15o, N2o, NTRo, NH2o, Nmo,  &
                   nhi_j, nheiS_j, nheii_j, nheiTR_j, nh2_j, nm_j, xion_j, &
                   fvq_j, evq_j, has_h2, sec_on, mol_sec, D0e, Ekere,      &
                   P1, P15, P2, PTR, Pm_r, PH2, PH2di, PH2dd, PH2nd,       &
                   h1HI, h1HeI, h1HeII, h1HeTR, h1H2, h1m_r, heat_r,       &
                   chan_r, q_r, qabs_r, Pm2_r)
      ! The field of cell j in the HHe field loop of ioniz_eq, with the
      ! arguments and results of photoionization_field_at_cell_HHe: the
      ! basis of the cell is reused or formed (validity: module header).
      ! Writes only the elements j of the module arrays (the loop over cells
      ! is parallel).
      integer, intent(in) :: j
      real*8, intent(in) :: N1o, N15o, N2o, NTRo, NH2o, Nmo(n_mphot)
      real*8, intent(in) :: nhi_j, nheiS_j, nheii_j, nheiTR_j, nh2_j
      real*8, intent(in) :: nm_j(n_mion), xion_j, fvq_j, evq_j
      logical, intent(in) :: has_h2, sec_on, mol_sec
      real*8, intent(in) :: D0e, Ekere
      real*8, intent(out) :: P1, P15, P2, PTR, Pm_r(n_mion), PH2, PH2di
      real*8, intent(out) :: PH2dd, PH2nd, h1HI, h1HeI, h1HeII, h1HeTR, h1H2
      real*8, intent(out) :: h1m_r(n_mion), heat_r, chan_r(6), q_r, qabs_r
      real*8, intent(out) :: Pm2_r(n_mion)
      real*8 :: inp(sbr_nin), mh, ln_eps
      integer :: code
      logical :: want, used_b, reused_b

      if (.not. sbr_in_pass .or. j .lt. 1 .or. j .gt. N .or.               &
          (a_tau .ne. 0.0d0) .or. (photoheat_photon_fraction .gt. 0.0d0)) then
         call photoionization_field_at_cell_HHe(j, N1o, N15o, N2o, NTRo,    &
                  NH2o, Nmo, nhi_j, nheiS_j, nheii_j, nheiTR_j, nh2_j, nm_j,   &
                  xion_j, fvq_j, evq_j, has_h2, sec_on, mol_sec, D0e, Ekere,   &
                  P1, P15, P2, PTR, Pm_r, PH2, PH2di, PH2dd, PH2nd,            &
                  h1HI, h1HeI, h1HeII, h1HeTR, h1H2, h1m_r, heat_r, chan_r,    &
                  q_r, qabs_r, Pm2_r)
         if (j .ge. 1 .and. j .le. N) then
            if (sbr_in_pass) then
               sbr_basis(j)%valid = .false.
               sbr_code(j) = 6
            else
               sbr_code(j) = 13
            endif
         endif
         return
      endif
      inp(1) = N1o; inp(2) = N15o; inp(3) = N2o; inp(4) = NTRo; inp(5) = NH2o
      inp(sbr_iNcol+1:sbr_iNcol+n_mphot) = Nmo
      inp(sbr_in+1) = nhi_j; inp(sbr_in+2) = nheiS_j; inp(sbr_in+3) = nheii_j
      inp(sbr_in+4) = nheiTR_j; inp(sbr_in+5) = nh2_j
      inp(sbr_inm+1:sbr_inm+n_mion) = nm_j
      inp(sbr_ix+1) = xion_j; inp(sbr_ix+2) = fvq_j; inp(sbr_ix+3) = evq_j
      inp(sbr_ix+4) = dr_j(j)*R0*opa_pf(j)
      ! ln(1 + eps_F) by its series (eps_F = 1e-9: the cubic term is 3e-28)
      ln_eps = spectral_basis_reuse_eps*(1.0d0 - spectral_basis_reuse_eps*(0.5d0 - spectral_basis_reuse_eps/3.0d0))
      want = .false.
      ! finite inputs: the sum of their magnitudes is finite (one test; a sum
      ! overflowing from finite inputs of order huge also recomputes)
      if (.not. ieee_is_finite(sum(abs(inp)))) then
         code = 10
      else if (sbr_pass .le. 1) then
         code = 3
      else if (.not. sbr_entry_reuse) then
         code = 12
      else if (sbr_cur_app .ne. 1) then
         code = 5
      else if (.not. sbr_basis(j)%valid .or. sbr_basis_entry(j) .ne. sbr_entry &
               .or. sbr_basis_pass(j) .ge. sbr_pass) then
         code = 4
      else if (sbr_basis(j)%underflow) then
         code = 8
      else
         mh = scalar_attenuation_bound(inp, sbr_basis_inp(:, j), has_h2)
         if (ieee_is_finite(mh) .and. mh .le. ln_eps) then
            want = .true.
            code = 1
            if (all(inp .eq. sbr_basis_inp(:, j))) code = 2
         else
            code = 9
         endif
      endif
      call photoionization_field_from_spectral_basis(j, N1o, N15o, N2o, NTRo, &
               NH2o, Nmo, nhi_j, nheiS_j, nheii_j, nheiTR_j, nh2_j, nm_j,       &
               xion_j, fvq_j, evq_j, has_h2, sec_on, mol_sec, D0e, Ekere,       &
               P1, P15, P2, PTR, Pm_r, PH2, PH2di, PH2dd, PH2nd,                &
               h1HI, h1HeI, h1HeII, h1HeTR, h1H2, h1m_r, heat_r, chan_r,        &
               q_r, qabs_r, Pm2_r, sbr_basis(j), want, used_b, reused_b)
      if (.not. used_b) then
         code = 7
      else if (want .and. .not. reused_b) then
         code = 7
      endif
      if (used_b .and. .not. reused_b) then
         ! a basis formed by this call: it carries this call's identity
         sbr_basis_entry(j) = sbr_entry
         sbr_basis_pass(j)  = sbr_pass
         sbr_basis_inp(:, j) = inp
         if (code .eq. 10) sbr_basis(j)%valid = .false.
      endif
      sbr_code(j) = code
      end subroutine field_with_spectral_basis_reuse


      subroutine spectral_basis_reuse_counts(c)
      ! The counts of the run, by decision (sbr_code_name order).
      integer*8, intent(out) :: c(sbr_ncode)
      c = sbr_count
      end subroutine spectral_basis_reuse_counts


      subroutine spectral_basis_reuse_report()
      ! The footer of run.log (nothing is printed when no HHe field loop ran).
      integer :: k
      if (all(sbr_count .eq. 0)) return
      write(*,'(A,ES9.2,A)') '   spectral-basis reuse (eps_F =', spectral_basis_reuse_eps, &
           '): field requests of the HHe field loop by decision,'
      write(*,'(A)') '     counted once for every cell of every field loop run (every ghost application)'
      do k = 1, sbr_ncode
         if (sbr_count(k) .gt. 0) write(*,'(A,A30,1X,I0)') '     ', sbr_code_name(k), sbr_count(k)
      enddo
      end subroutine spectral_basis_reuse_report

      end module field_spectral_basis_reuse
