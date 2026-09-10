      program heating_channel_closure
      ! One assembly of the volumetric heating rate, and its channels close.
      !
      ! WHAT IS BEING ASSERTED.  heating_of_composition
      ! (src/modules/radiation/util_ion_eq.f90) is the one place the heating
      ! of a cell is put together: the ionization sweep takes the total it
      ! returns for the energy equation, output/Heating_breakdown.txt writes
      ! the channel array it fills, and the advection-corrected post-process
      ! calls it again on its own composition.  The property that makes one
      ! assembly worth having is that the total is the running sum of the
      ! channels, so a deposit cannot reach the energy equation and be absent
      ! from the breakdown.  That is assertion (1).
      !
      ! Assertion (2) is what makes assertion (1) mean anything: the
      ! synthetic cell below has to run EVERY channel.  A closure test on a
      ! cell whose molecular, oxygen and metastable deposits are all zero
      ! closes trivially, and it was exactly two silent zeros in a second
      ! copy of this sum -- the associative He(2^3S) branch and the
      ! collisional oxygen channels -- that the breakdown file was missing.
      !
      ! THE CELL.  A molecular base: 1400 K, 10^13 cm^-3 of H2 over a partly
      ! ionized H/He gas with the He 2^3S metastable, the oxygen carriers OH,
      ! H2O and CO, trace metals, an attenuated XUV field, a Lyman-Werner
      ! band flux and FUV photolysis rates.  The numbers are order-of-
      ! magnitude values of the hot-Uranus gate, not a solution of anything:
      ! nothing here is asserted against a physical reference, only the
      ! closure of the sum and the fact that every channel fires.
      use global_parameters
      use species_table,   only: n_mion, im_OI, mion_fsp
      use water_photolysis, only: n_fuv_band, water_photolysis_init
      use utils_ion_eq,    only: heating_of_composition, n_heat_channel,   &
                                 heat_channel_name
      use assertion_report, only: check_absolute, assertion_failures

      implicit none

      integer, parameter :: ncell = 3
      real*8,  parameter :: t_cell = 1.4d3      ! K
      real*8,  parameter :: n_h2   = 1.0d13     ! cm^-3
      real*8,  parameter :: n_hi   = 2.0d12
      real*8,  parameter :: n_hii  = 1.0d10
      real*8,  parameter :: n_hei  = 1.0d12
      real*8,  parameter :: n_heii = 1.0d9
      real*8,  parameter :: n_hetr = 1.0d5
      real*8,  parameter :: n_e    = 1.1d10
      ! Relative closure the sum of 17 channels can be held to.  The total
      ! and the check below add the same 17 numbers in different orders, so
      ! what separates them is the reassociation of a floating sum, a few
      ! units in the last place of the largest term.  1e-14 is about 45 of
      ! them on a double.
      real*8,  parameter :: tol_close = 1.0d-14

      integer :: jlo, jhi, j, ic, n_zero
      real*8, allocatable :: T_K(:), nhi(:), nhii(:), nhei(:), nheii(:),   &
                             nheiii(:), nheiTR(:), ne(:), n_tot(:)
      real*8, allocatable :: nm(:,:), nmol(:,:), nox(:,:)
      real*8, allocatable :: h1_HI(:), h1_HeI(:), h1_HeII(:), h1_HeTR(:),  &
                             h1_H2(:), h1_m(:,:)
      real*8, allocatable :: q31a(:), q31b(:), Q31(:)
      real*8, allocatable :: k_lw(:), p_lw(:), j_h2o(:,:), j_oh(:,:)
      ! CO photodissociation rate, the fourth absorber of the same beam.
      real*8, allocatable :: k_co(:)
      real*8, allocatable :: heat(:), heat_chan(:,:)
      real*8 :: A31, csum, worst, dev, scale

      ! ---- Run configuration: every optional channel on.
      N  = ncell
      ! Planet radius: the grid spacing is carried in units of it, and the
      ! He recombination photons are absorbed over one cell width.
      R0 = 1.0d10
      thereis_He      = .true.
      thereis_HeITR   = .true.
      thereis_metals  = .true.
      thereis_mol     = .true.
      thereis_oxychem = .true.
      mol_reaction_heat   = .true.
      use_excited_H       = .true.
      use_he_rec_coupling = .true.
      F_LW_star  = 3.43d2
      F_Lya_star = 0.0d0
      ! The FUV band thresholds, without which one photolysis event deposits
      ! nothing (the guard in heat_per_water_dissociation).
      call water_photolysis_init

      jlo = 1 - Ng
      jhi = N + Ng
      allocate(T_K(jlo:jhi), nhi(jlo:jhi), nhii(jlo:jhi), nhei(jlo:jhi),   &
               nheii(jlo:jhi), nheiii(jlo:jhi), nheiTR(jlo:jhi),           &
               ne(jlo:jhi), n_tot(jlo:jhi))
      allocate(nm(jlo:jhi,n_mion), nmol(jlo:jhi,4), nox(jlo:jhi,3))
      allocate(h1_HI(jlo:jhi), h1_HeI(jlo:jhi), h1_HeII(jlo:jhi),          &
               h1_HeTR(jlo:jhi), h1_H2(jlo:jhi), h1_m(jlo:jhi,n_mion))
      allocate(q31a(jlo:jhi), q31b(jlo:jhi), Q31(jlo:jhi))
      allocate(k_lw(jlo:jhi), p_lw(jlo:jhi), k_co(jlo:jhi))
      allocate(j_h2o(jlo:jhi,n_fuv_band), j_oh(jlo:jhi,n_fuv_band))
      allocate(heat(jlo:jhi), heat_chan(jlo:jhi,n_heat_channel))
      ! The two H(n=2) deposits are already-contracted volumetric rates that
      ! excited_H_update fills; the assembly reads them from the module.
      allocate(Hpe_arr(jlo:jhi), Hdx_arr(jlo:jhi))
      ! The He recombination coupling escapes its photons over one cell
      ! width, so it reads the grid spacing.
      allocate(dr_j(jlo:jhi), r(jlo:jhi))
      dr_j = 1.0d7
      do j = jlo, jhi
         r(j) = 1.0d0 + dble(j - 1)*1.0d-3
      enddo

      T_K    = t_cell
      nhi    = n_hi
      nhii   = n_hii
      nhei   = n_hei
      nheii  = n_heii
      nheiii = 1.0d6
      nheiTR = n_hetr
      ne     = n_e
      n_tot  = n_h2 + n_hi + n_hii + n_hei

      nmol = 0.0d0
      nmol(:,1) = n_h2         ! H2
      nmol(:,2) = 1.0d6        ! H2+
      nmol(:,3) = 1.0d7        ! H3+
      nmol(:,4) = 1.0d4        ! HeH+
      nox  = 0.0d0
      nox(:,1) = 1.0d7         ! OH
      nox(:,2) = 1.0d8         ! H2O
      nox(:,3) = 1.0d8         ! CO
      nm = 0.0d0
      nm(:,im_OI) = 1.0d8      ! free atomic oxygen, the O2 channel runs on it
      nm(:,1)     = 1.0d7

      ! Photoheating of one particle of each absorber [erg s^-1]: values of
      ! an attenuated XUV field at the top of such a layer.
      h1_HI   = 1.0d-24
      h1_HeI  = 3.0d-24
      h1_HeII = 1.0d-25
      h1_HeTR = 5.0d-23
      h1_H2   = 8.0d-25
      h1_m    = 0.0d0
      h1_m(:,1) = 2.0d-24

      ! He 2^3S coefficients and the Lyman-Werner and FUV rates.
      A31  = 1.272d-4
      q31a = 1.0d-11
      q31b = 1.0d-11
      Q31  = 5.0d-10
      k_lw = 1.0d-9
      p_lw = 0.135d0
      ! Shielded but not extinguished, so both CO channels deposit on this
      ! cell and the every-channel assertion covers them.
      k_co = 3.0d-9
      j_h2o = 1.0d-8
      j_oh  = 1.0d-8

      Hpe_arr = 1.0d-14
      Hdx_arr = 4.0d-15

      call heating_of_composition(T_K,                                     &
               nhi,nhii,nhei,nheii,nheiii,nheiTR, nm, nmol, nox,           &
               ne, n_tot,                                                  &
               h1_HI,h1_HeI,h1_HeII,h1_HeTR,h1_H2,h1_m,                    &
               A31,q31a,q31b,Q31, k_lw, p_lw, k_co, j_h2o, j_oh,          &
               .true., .true., heat, heat_chan)

      ! ---- (1) The total is the sum of the channels.
      worst = 0.0d0
      do j = jlo, jhi
         csum = 0.0d0
         do ic = 1, n_heat_channel
            csum = csum + heat_chan(j,ic)
         enddo
         scale = max(abs(heat(j)), abs(csum))
         if (scale .gt. 0.0d0) then
            dev = abs(csum - heat(j))/scale
            worst = max(worst, dev)
         endif
      enddo
      call check_absolute('heat_total_is_channel_sum', worst, 0.0d0,     &
                          tol_close)

      ! ---- (2) Every channel of the list deposits on this cell, so the
      ! closure above is a statement about all seventeen of them.
      n_zero = 0
      do ic = 1, n_heat_channel
         if (heat_chan(1,ic) .eq. 0.0d0) then
            n_zero = n_zero + 1
            write(*,'(a,i3,1x,a)') '  channel with no deposit: ', ic,      &
                                   trim(heat_channel_name(ic))
         endif
      enddo
      call check_absolute('every_heat_channel_exercised', dble(n_zero),  &
                          0.0d0, 0.0d0)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'heating_channel_closure: ',                  &
              assertion_failures, ' assertion(s) failed'
         stop 1
      endif
      write(*,'(a)') 'heating_channel_closure: every assertion passed'

      end program heating_channel_closure
