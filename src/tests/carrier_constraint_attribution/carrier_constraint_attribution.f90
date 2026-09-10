      program carrier_constraint_attribution
      ! WHICH RESIDUALS BELONG TO THE SOLVER AND WHICH TO A CONSTRAINT.
      !
      ! photochemical_transport_step stops the run when the carrier Newton
      ! ends short in a carrier that "is not on an element constraint".  The
      ! state it measures that on is not the Newton's last iterate:
      ! limit_to_element_budget runs in between, and where an element clamp
      ! acts it rescales carriers of that cell.  A row measured after such a
      ! clamp is not the row the Newton solved, so the residual standing in
      ! it is the constraint's.
      !
      ! THE CONSTRAINT USED HERE IS THE CARBON CLAMP: CO is the only carbon
      ! carrier, so a cell whose carried CO exceeds the free carbon of the
      ! cell has its CO and nothing else rescaled.  It is conservation, and
      ! it is the only constraint that acts on CO alone.
      !
      ! THE PROPERTY UNDER TEST IS THAT THE ATTRIBUTION IS PER CELL.  The
      ! clamps act cell by cell.  Two statements follow, and both are
      ! measured here:
      !
      !  (1) a cell the carbon clamp acted in is marked, and
      !  (2) a cell no constraint touched is NOT marked, however many other
      !      cells of the same column were clamped.
      !
      ! (2) is the one an attribution indexed by CARRIER cannot make: a
      ! carrier-indexed record is a maximum over the grid, so one clamped
      ! cell sets the CO entry for the whole domain and excuses a genuine
      ! solver failure in any cell of it.
      !
      !  (3) THE PREMISE, measured rather than assumed: a clamp on CO alone
      !      moves the OH row of the same cell.  carrier_source closes the
      !      free atomic oxygen as n_o0 = nO_free - n(O ions) - n(OH) -
      !      n(H2O) - n(CO), so removing CO hands its oxygen to O I and
      !      changes the O + H2 -> OH + H production the OH row is built
      !      from.  The carriers the clamp did not itself rescale are
      !      therefore not free of it either, which is why the question is
      !      asked of the cell and not of the carrier.
      !
      ! The gas is the metal-free hydrogen and helium column of
      ! carrier_reference_scales -- the frozen background of a cell is what
      ! carrier_source reads, and this suite needs it defined, not
      ! realistic.  The free oxygen and carbon densities the limiter is
      ! given are arguments of it, so they are set here directly.

      use global_parameters
      use species_table
      use ion_cell_state, only: ion_rates
      use ionization_equilibrium, only: bg_cell, ioniz_eq_allocate_arrays
      use diffusive_photochemistry, only: carrier_set_init, carrier_state, &
                                          carrier_source,                  &
                                          limit_to_element_budget,         &
                                          carrier_cell_is_constrained,     &
                                          n_carrier_max,                   &
                                          ic_H2, ic_OH, ic_H2O, ic_CO
      use assertion_report
      implicit none

      ! The cell whose CO exceeds the free carbon, and one that stays inside
      ! it.  They are not neighbours: the limiter is cell-local, so a mark
      ! that leaked would leak to a neighbour first.
      integer, parameter :: j_over = 2
      integer, parameter :: j_free = 5
      ! A third cell holding CO a relative 1e-12 above the free carbon: the
      ! clamp moves it, and by less than the limit_report = 1e-10 threshold
      ! the limiter's counters are reported above.
      integer, parameter :: j_tiny = 4
      real*8,  parameter :: over_tiny = 1.0d-12
      ! How far above the free carbon the carried CO of j_over sits.
      real*8,  parameter :: over_big  = 0.5d0
      ! Free element densities handed to the limiter [cm^-3].  Carbon is the
      ! scarcer of the two, so it and not oxygen bounds CO.
      real*8, parameter :: nO_use = 4.8d10
      real*8, parameter :: nC_use = 2.6d10
      ! The gas temperature of the column [K].  Nothing asserted here
      ! depends on it: the carbon clamp is conservation and carries no
      ! temperature.
      real*8, parameter :: T_gas  = 1.2d3

      real*8, allocatable :: rho(:), f_sp(:,:), fc(:,:)
      real*8, allocatable :: ntot(:), TK(:), mbar(:), nrho(:), wfac(:)
      real*8, allocatable :: nH_free(:), nO_free(:), nC_free(:)
      real*8, allocatable :: TK_use(:), nO_arr(:), nC_arr(:)
      real*8 :: nc(n_carrier_max), src_before(n_carrier_max)
      real*8 :: src_after(n_carrier_max)
      real*8 :: nco_before, nco_after, nco_free_before
      real*8 :: doh
      integer :: j

      call setup_globals()
      call build_oxygen_bearing_column()

      ! ---- the clamp really binds in one cell and not in the other ----- !
      nco_before      = fc(j_over,ic_CO)*nrho(j_over)
      nco_free_before = fc(j_free,ic_CO)*nrho(j_free)
      call check_positive('carried_CO_above_the_free_carbon',             &
                          nco_before - nC_use)
      call check_positive('the_other_cell_stays_inside_its_carbon',       &
                          nC_use - nco_free_before)

      ! A clamp of that size is far below the reporting threshold, and it
      ! is still a change to the state the Newton measured its rows at.
      fc(j_tiny,ic_CO) = nC_use*(1.0d0 + over_tiny)/nrho(j_tiny)

      ! ---- the OH row of the clamped cell, before and after ------------ !
      nc = fc(j_over,:)*nrho(j_over)
      call carrier_source(j_over, nc, nH_free(j_over), nO_use, src_before)

      call limit_to_element_budget(fc, nrho, TK_use, nH_free, nO_arr,     &
                                   nC_arr)

      nco_after = fc(j_over,ic_CO)*nrho(j_over)
      call check_relative('the_cell_CO_is_clamped_to_its_free_carbon',    &
                          nco_after, nC_use, 1.0d-12)
      call check_relative('the_other_cell_CO_is_untouched',               &
                          fc(j_free,ic_CO)*nrho(j_free),                  &
                          nco_free_before, 1.0d-12)

      ! ---- (1) and (2): the attribution ------------------------------- !
      call check_absolute('clamped_cell_is_marked_constrained',           &
                          logical_as_double(                              &
                             carrier_cell_is_constrained(j_over)),        &
                          1.0d0, 0.0d0)
      call check_absolute('untouched_cell_is_not_marked_constrained',     &
                          logical_as_double(                              &
                             carrier_cell_is_constrained(j_free)),        &
                          0.0d0, 0.0d0)

      ! ---- (2b) a clamp below the reporting threshold marks its cell --- !
      ! The mark and the counter answer different questions. The counter
      ! says the constraint did something worth reporting; the mark says
      ! this cell no longer holds the state the Newton measured, and that is
      ! true of a rescaling of any size, because a carrier row is a stiff
      ! difference of large terms and a clamp of 1e-12 in the density is not
      ! a change of 1e-12 in the row.
      call check_relative('tiny_clamp_moved_the_cell',                    &
                          fc(j_tiny,ic_CO)*nrho(j_tiny), nC_use,          &
                          1.0d-12)
      call check_positive('tiny_clamp_is_below_the_reporting_threshold',  &
                          1.0d-10 - over_tiny)
      call check_absolute('a_clamp_below_the_threshold_marks_its_cell',   &
                          logical_as_double(                              &
                             carrier_cell_is_constrained(j_tiny)),        &
                          1.0d0, 0.0d0)

      ! ---- (3) the premise: the CO clamp moved the OH row -------------- !
      nc = fc(j_over,:)*nrho(j_over)
      call carrier_source(j_over, nc, nH_free(j_over), nO_use, src_after)
      doh = abs(src_after(ic_OH) - src_before(ic_OH))                     &
            /max(abs(src_before(ic_OH)), abs(src_after(ic_OH)), 1.0d-300)
      ! A tenth of the row is far more than the round-off a clamp of no
      ! consequence could produce, and far less than the change measured
      ! here, so the assertion separates the two without pinning a number
      ! the oxygen network fixes.
      call check_positive('CO_clamp_moves_the_OH_row_of_its_cell',        &
                          doh - 0.1d0)

      if (assertion_failures .gt. 0) stop 1

      contains

      !--------------!

      double precision function logical_as_double(l) result(x)
      logical, intent(in) :: l
      x = 0.0d0
      if (l) x = 1.0d0
      end function logical_as_double

      !--------------!

      subroutine setup_globals()
      ! A five-cell hydrogen and helium column with the molecular chemistry
      ! and the oxygen cycle on, so that OH, H2O and CO are carriers.
      integer :: j
      N   = 5
      n0  = 1.0d10
      R0  = 1.0d10
      v0  = 1.0d6
      T0  = 1000.0d0
      HeH = 0.0793d0
      p_base_bar         = 1.0d-6
      thereis_He         = .true.
      thereis_HeITR      = .false.
      thereis_metals     = .true.
      thereis_mol        = .true.
      thereis_oxychem    = .true.
      eos_include_metals = .false.
      he_diffusion       = .false.
      carrier_transport  = .true.
      ionization_transport = .false.
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         r(j)     = 1.0d0 + 0.01d0*dble(j)
         r_edg(j) = 1.0d0 + 0.01d0*(dble(j) + 0.5d0)
         dr_j(j)  = 0.01d0
      enddo
      allocate(melem_ab(n_melem))
      melem_ab = 0.0d0
      call carrier_set_init()
      call ioniz_eq_allocate_arrays()
      end subroutine setup_globals

      !--------------!

      subroutine build_oxygen_bearing_column()
      ! Molecular base gas carrying CO at 90 per cent of the free carbon and
      ! a trace of the water family, with cell j_over carrying more CO than
      ! the cell has carbon for.
      integer :: j
      allocate(rho(1-Ng:N+Ng), f_sp(1-Ng:N+Ng,n_species))
      allocate(fc(1-Ng:N+Ng,n_carrier_max))
      allocate(ntot(1-Ng:N+Ng), TK(1-Ng:N+Ng), mbar(1-Ng:N+Ng))
      allocate(nrho(1-Ng:N+Ng), wfac(1-Ng:N+Ng))
      allocate(nH_free(1-Ng:N+Ng), nO_free(1-Ng:N+Ng), nC_free(1-Ng:N+Ng))
      allocate(TK_use(1-Ng:N+Ng), nO_arr(1-Ng:N+Ng), nC_arr(1-Ng:N+Ng))
      rho  = 1.0d0
      f_sp = 0.0d0
      f_sp(:,isp_HI)   = 0.20d0
      f_sp(:,isp_HII)  = 1.0d-6
      f_sp(:,isp_H2)   = 0.40d0
      f_sp(:,isp_H2p)  = 1.0d-10
      f_sp(:,isp_H3p)  = 1.0d-10
      f_sp(:,isp_HeI)  = 0.0793d0
      f_sp(:,isp_HeII) = 1.0d-8
      do j = 1-Ng, N+Ng
         call set_frozen_cell_rates(j)
      enddo
      call carrier_state(rho, f_sp, fc, ntot, nrho, wfac, TK, mbar,       &
                         nH_free, nO_free, nC_free)
      ! The element headroom the limiter is tested against.  The column is
      ! metal-free, so carrier_state returns zero oxygen and carbon; the
      ! limiter takes both as arguments and they are set here.
      nO_arr = nO_use
      nC_arr = nC_use
      TK_use = T_gas
      ! The carried state: CO at 90 per cent of the carbon everywhere, the
      ! water family at a trace of the oxygen, and one cell carrying half
      ! again as much CO as its carbon allows.  The oxygen budget is still
      ! satisfied there, so the carbon clamp is the only constraint acting.
      fc(:,ic_CO)  = 0.90d0*nC_use/nrho
      fc(:,ic_OH)  = 1.0d-4*nO_use/nrho
      fc(:,ic_H2O) = 1.0d-3*nO_use/nrho
      fc(j_over,ic_CO) = (1.0d0 + over_big)*nC_use/nrho(j_over)
      end subroutine build_oxygen_bearing_column

      !--------------!

      subroutine set_frozen_cell_rates(j)
      ! Every field of the frozen cell state, so that none is read
      ! undefined.  The values are the ones carrier_reference_scales uses;
      ! nothing asserted here depends on their size.
      integer, intent(in) :: j
      bg_cell(j)%P_HI        = 1.0d-6      ! [1/s]
      bg_cell(j)%P_HeI       = 1.0d-7
      bg_cell(j)%P_HeII      = 0.0d0
      bg_cell(j)%P_HeITR     = 0.0d0
      bg_cell(j)%P_H2        = 1.0d-8
      bg_cell(j)%P_H2_di     = 0.0d0
      bg_cell(j)%P_H2_dd     = 0.0d0
      bg_cell(j)%P_H2_nd     = 0.0d0
      bg_cell(j)%k_LW        = 0.0d0
      bg_cell(j)%rchiiB      = 2.6d-13     ! [cm^3/s] at ~1e4 K
      bg_cell(j)%rcheiiB     = 4.3d-13
      bg_cell(j)%rcheiiiB    = 2.2d-12
      bg_cell(j)%rcheiTR     = 0.0d0
      bg_cell(j)%a_ion_HI    = 0.0d0
      bg_cell(j)%a_ion_HeI   = 0.0d0
      bg_cell(j)%a_ion_HeII  = 0.0d0
      bg_cell(j)%a_ion_HeITR = 0.0d0
      bg_cell(j)%q13         = 0.0d0
      bg_cell(j)%q31a        = 0.0d0
      bg_cell(j)%q31b        = 0.0d0
      bg_cell(j)%Q31         = 0.0d0
      bg_cell(j)%A31         = 0.0d0
      bg_cell(j)%kcx_He0_Hp  = 0.0d0
      bg_cell(j)%kcx_Hep_H0  = 0.0d0
      bg_cell(j)%nh          = 0.0d0
      bg_cell(j)%nhe         = 0.0d0
      bg_cell(j)%n_ofam      = 0.0d0
      bg_cell(j)%n_co        = 0.0d0
      bg_cell(j)%x_h2_fixed  = .false.
      bg_cell(j)%x_ox_fixed  = .false.
      bg_cell(j)%x_hp_fixed  = .false.
      bg_cell(j)%x_h2_fix    = 0.0d0
      bg_cell(j)%x_oh_fix    = 0.0d0
      bg_cell(j)%x_h2o_fix   = 0.0d0
      bg_cell(j)%x_hp_fix    = 0.0d0
      bg_cell(j)%T_K         = T_gas
      bg_cell(j)%ntot        = 0.5d0*rho(j)*n0
      end subroutine set_frozen_cell_rates

      end program carrier_constraint_attribution
