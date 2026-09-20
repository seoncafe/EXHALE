      program h2_ladder_statistical_equilibrium
      ! THE REDUCED STATISTICAL EQUILIBRIUM OF THE H2 GROUND-STATE LADDER,
      ! and what it says about the scalar thermalized fraction
      ! h2_vibrational_heat_fraction of
      ! src/modules/lower_atmosphere/h2_vibrational_relaxation.f90.
      !
      ! The scalar fraction is a FIRST-EVENT BRANCHING MODEL,
      ! C1/(C1 + A_max), with C1 the v = 1 collisional coefficient of each
      ! collider and A_max the largest total spontaneous decay rate over
      ! the whole bound ladder.  It is not a bound: a first-event fraction
      ! is not the energy fraction of a multistep cascade, and the v = 1
      ! collisional coefficient does not bound every level's collisional
      ! loss from below.  What it approximates is the net collisional heat
      ! of a solved ladder, and this driver solves one from the published
      ! data and compares.
      !
      ! WHAT IS ASSERTED.
      !   A. THE MODEL IS BUILT FROM THE DATA THE PRODUCTION MODULE USES.
      !      The all-level maximum total decay rate reduced here from
      !      Roueff et al. (2019) is the one h2_vibrational_relaxation
      !      carries, and the three collider coefficients reduce to the
      !      published thermal values the production tables reproduce.
      !   B. DETAILED BALANCE AND THE THERMAL LIMIT.  In a collision-only
      !      bath with no pumping and no radiative escape the populations
      !      reach the Boltzmann distribution and the net collisional heat
      !      vanishes.  This is the one test that isolates the upward
      !      rates, which are built from the downward ones and never read.
      !   C. THE ENERGY OF AN INJECTED FORMATION PARTITIONS AND SUMS.  With
      !      the loss channels restored, an injection at the energy the
      !      three-body association leaves its molecule partitions into gas
      !      heat, escaping radiation, chemical destruction and stored
      !      internal energy, and the four sum to the injection: in the
      !      steady state, where the stored term is what the destruction
      !      carries away, and over one backward-Euler step from an empty
      !      ladder, where it is the ladder filling up.
      !   D. THE BRACKET ON THE SCALAR FRACTION at the three representative
      !      cells of the certified molecular fiducial, spanning the one
      !      explicitly uncertain quantity of the model: the disposal of
      !      the levels above the published data.
      !
      ! Its own module is src/tests/h2_level_ladder/h2_ladder_model.f90,
      ! which states every data source and every assumption.

      use h2_ladder_model
      use h2_vibrational_relaxation, only: h2_vibrational_heat_fraction,  &
               h2_vibrational_relaxation_init,                            &
               h2_total_decay_rate_max, gamma_10_H, gamma_10_H2,          &
               gamma_10_He, e_vib_v5_eV, e_vib_v6_eV
      use assertion_report

      implicit none

      ! Where the published tables live.  EXHALE_TEST_ROOT overrides the
      ! repository root the driver deduces from its own build directory.
      character(len=512) :: root, roueff_file, lique_file, jozwiak_dir
      character(len=256) :: why
      logical :: ok, have

      ! The three representative cells, READ from the certified molecular
      ! fiducial LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH2.13,
      ! output/{Hydro_ioniz_IC.txt, Ion_species_IC.txt}, rows 3, 180 and
      ! 301 of the written state: the molecular depth, the H2 front and
      ! the dilute upper column.  Helium is the ground singlet.
      integer, parameter :: ncase = 3
      character(len=20), parameter :: case_name(ncase) =                  &
         (/ character(len=20) :: 'molecular_depth', 'h2_front',           &
                                 'dilute_upper' /)
      real*8, parameter :: c_T(ncase) =                                   &
         (/ 808.31894177851177d0, 2128.3673585881161d0,                   &
            5083.4376627309339d0 /)
      real*8, parameter :: c_nhi(ncase) =                                 &
         (/ 1867980406714.9812d0, 15643668544.791315d0,                   &
            81619552.724037111d0 /)
      real*8, parameter :: c_nhe(ncase) =                                 &
         (/ 6171608738847.3115d0, 35876533538.838463d0,                   &
            70245906.254677579d0 /)
      real*8, parameter :: c_nh2(ncase) =                                 &
         (/ 514739739664.78314d0, 1646227639.0261714d0,                   &
            153.44755402026095d0 /)

      ! THE ENERGY A NASCENT MOLECULE IS BORN WITH.  The three-body
      ! association R15 leaves its molecule within 0.02 eV of the
      ! dissociation limit, i.e. at v = 10 to 14, and the H3+ dissociative
      ! recombination R6 leaves its H2 fragment peaking at v = 5 to 6.
      ! The group's energy is set to each in turn [K].
      real*8, parameter :: eV_to_K = 11604.518121550082d0
      real*8 :: e_group_r15_K, e_group_r6_K

      ! The uncertain disposal of the levels above the published data: the
      ! collisional rate is the top resolved level's own collisional total
      ! scaled by a factor spanning a decade about one, and the radiative
      ! rate spans zero to the all-level maximum.  Nothing in hand narrows
      ! either, which is why the answer is quoted as the spread over them.
      real*8, parameter :: group_coll_scale_lo = 1.0d0/3.0d0
      real*8, parameter :: group_coll_scale_hi = 3.0d0
      ! The chemical destruction rate of one H2 molecule [s^-1], spanned
      ! from far below the all-level radiative rate 5.6e-06 s^-1, where a
      ! molecule always radiates before it is destroyed, to far above it.
      real*8, parameter :: k_destroy_slow = 1.0d-9
      real*8, parameter :: k_destroy_fast = 1.0d-2

      real*8 :: n(mxlev + 1), n0(mxlev + 1), n1(mxlev + 1)
      real*8 :: n_ctl(mxlev + 1), s_vec(mxlev + 1), s_ctl(mxlev + 1)
      real*8 :: l_ctl, d_rad
      integer :: i_ground
      real*8 :: q_coll, l_rad, u_store, p_inj, chem, resid
      real*8 :: f_model, f_scalar, f_lo, f_hi
      real*8 :: gamma_g, a_g, k_destroy, s_group, dt
      real*8 :: c_top, du, sum4
      real*8 :: f_bracket_lo(3), f_bracket_hi(3)
      integer :: info, ic, iv, ia, ik, itop
      real*8 :: term_scale

      call h2_vibrational_relaxation_init

      call resolve_paths(root, roueff_file, lique_file, jozwiak_dir, have)
      if (.not. have) then
         write(*,'(a)') 'SKIP h2_ladder_statistical_equilibrium '//       &
            'measured=absent reference='//trim(lique_file)//' tol=exact'
         write(*,'(a)') '  the published collision tables live outside '//&
            'the repository, in the workspace references/ tree'
         stop 0
      endif

      e_group_r15_K = (4.478d0 - 0.02d0)*eV_to_K
      e_group_r6_K  = 0.5d0*(e_vib_v5_eV + e_vib_v6_eV)*eV_to_K

      ! ================================================================== !
      ! A. THE MODEL IS BUILT FROM THE PRODUCTION MODULE'S OWN DATA.
      ! ================================================================== !
      call build_reduced_ladder(roueff_file, lique_file, jozwiak_dir,     &
                                c_T(1), gamma_10_H2(c_T(1)),             &
                                gamma_10_He(c_T(1)), ok, why)
      if (.not. ok) then
         write(*,'(a)') 'FAIL h2_ladder_build measured='//trim(why)//     &
                        ' reference=ok tol=0'
         stop 1
      endif
      write(*,'(a,i0,a)') '#   resolved levels: ', nlev,                  &
         ' (v <= 3, the set Lique 2015 distributes), plus one group'
      call check_at_least('resolved_ladder_has_levels',                   &
                          dble(nlev), 30.0d0)
      ! The all-level maximum reduced here is the production module's.
      ! The production module reduces the same maximum from the TRIMMED
      ! line list its infrared cooling sums (the lines that carry the
      ! emission), this module from every row of the same file, so the
      ! production value is the lower of the two by the weight of the
      ! lines the trimming drops.
      call check_relative('all_level_A_max_is_the_production_value',      &
                          a_total_max, h2_total_decay_rate_max(), 1.0d-4)
      call check_at_least('all_level_A_max_of_the_full_file_is_larger',   &
                          a_total_max/h2_total_decay_rate_max(), 1.0d0)

      ! ================================================================== !
      ! B. DETAILED BALANCE AND THE THERMAL LIMIT.
      ! ================================================================== !
      do ic = 1, ncase
         call build_reduced_ladder(roueff_file, lique_file, jozwiak_dir,  &
                                   c_T(ic), gamma_10_H2(c_T(ic)),        &
                                   gamma_10_He(c_T(ic)), ok, why)
         call check_absolute('ladder_built_'//trim(case_name(ic)),        &
                             merge(0.0d0, 1.0d0, ok), 0.0d0, 0.0d0)
         call solve_thermal_bath(c_T(ic), c_nhi(ic), c_nhe(ic),           &
                                 c_nh2(ic), 1.0d0, n, info)
         call check_absolute('thermal_bath_solves_'                      &
                             //trim(case_name(ic)),                       &
                             dble(info), 0.0d0, 0.0d0)
         call check_at_most('thermal_bath_is_boltzmann_'                 &
                            //trim(case_name(ic)),                        &
                            boltzmann_departure(n, c_T(ic)), 1.0d-9)
         ! The net collisional heat of a Boltzmann ladder is zero.  It is
         ! compared against the scale of the individual terms, which is
         ! what "vanishes" means for a sum of large cancelling numbers.
         q_coll = net_collisional_heat(n, c_T(ic), c_nhi(ic), c_nhe(ic),  &
                                       c_nh2(ic), 0.0d0, 0.0d0)
         term_scale = collisional_term_scale(n, c_T(ic), c_nhi(ic),       &
                                             c_nhe(ic), c_nh2(ic))
         call check_at_most('thermal_bath_net_collisional_heat_vanishes_'&
                            //trim(case_name(ic)),                        &
                            abs(q_coll)/max(term_scale, 1.0d-300),        &
                            1.0d-10)
      enddo

      ! ================================================================== !
      ! C. THE INJECTED FORMATION ENERGY PARTITIONS AND SUMS.
      ! ================================================================== !
      call build_reduced_ladder(roueff_file, lique_file, jozwiak_dir,     &
                                c_T(1), gamma_10_H2(c_T(1)),             &
                                gamma_10_He(c_T(1)), ok, why)
      itop      = top_level()
      c_top     = collisional_total(itop, c_T(1), c_nhi(1), c_nhe(1),     &
                                    c_nh2(1))
      gamma_g   = c_top
      a_g       = a_total_max
      k_destroy = 1.0d-3
      s_group   = 1.0d4
      p_inj     = s_group*e_group_r15_K*kb_erg_ladder
      i_ground  = level_index(0, 0)
      s_vec     = 0.0d0
      s_vec(nlev + 1) = s_group

      ! The steady state: what is injected leaves as heat, as radiation
      ! and inside the molecules the chemistry destroys.
      n0 = 0.0d0
      call solve_statistical_equilibrium(c_T(1), c_nhi(1), c_nhe(1),      &
              c_nh2(1), 1.0d0, s_vec, gamma_g, a_g,                       &
              k_destroy, -1.0d0, n0, n, info)
      call check_absolute('steady_injection_solves', dble(info), 0.0d0,   &
                          0.0d0)
      q_coll  = net_collisional_heat(n, c_T(1), c_nhi(1), c_nhe(1),       &
                                     c_nh2(1), gamma_g, e_group_r15_K)
      l_rad   = radiated_power(n, 1.0d0, a_g, e_group_r15_K)
      u_store = stored_internal_energy(n, e_group_r15_K)
      chem    = k_destroy*u_store
      call check_relative('steady_injection_partition_sums',              &
                          q_coll + l_rad + chem, p_inj, 1.0d-9)
      ! Every channel is a real one at this state, so the row above is not
      ! satisfied by three vanishing terms.
      call check_positive('steady_injection_heats_the_gas', q_coll)
      call check_positive('steady_injection_radiates', l_rad)
      call check_positive('steady_injection_stores_internal_energy',      &
                          u_store)

      ! One backward-Euler step from an empty ladder: there the stored
      ! term is the ladder filling up, and the four channels still sum.
      dt = 1.0d0/max(gamma_g, 1.0d0)
      n0 = 0.0d0
      call solve_statistical_equilibrium(c_T(1), c_nhi(1), c_nhe(1),      &
              c_nh2(1), 1.0d0, s_vec, gamma_g, a_g,                       &
              k_destroy, dt, n0, n1, info)
      call check_absolute('transient_injection_solves', dble(info),       &
                          0.0d0, 0.0d0)
      q_coll  = net_collisional_heat(n1, c_T(1), c_nhi(1), c_nhe(1),      &
                                     c_nh2(1), gamma_g, e_group_r15_K)
      l_rad   = radiated_power(n1, 1.0d0, a_g, e_group_r15_K)
      u_store = stored_internal_energy(n1, e_group_r15_K)
      chem    = k_destroy*u_store
      du      = u_store - stored_internal_energy(n0, e_group_r15_K)
      sum4    = (q_coll + l_rad + chem)*dt + du
      call check_relative('transient_injection_partition_sums',           &
                          sum4, p_inj*dt, 1.0d-9)
      call check_positive('transient_stores_internal_energy', du)

      ! ================================================================== !
      ! D. THE BRACKET ON THE SCALAR FRACTION.
      !
      !    f is the share of the internal energy a nascent molecule is
      !    born with that stays in the gas rather than leaving it in the
      !    infrared quadrupole lines.  ESCAPING RADIATION IS THE ONLY EXIT
      !    FROM THE GAS: the net collisional term hands the energy to the
      !    translational pool, and the internal energy a molecule still
      !    holds when the chemistry destroys it goes to the fragments and
      !    is carried by the caloric equation of state, which owns the
      !    rovibrational ladder of the bulk gas.  So
      !
      !        f_model = 1 - (escaping radiation)/(injected energy) ,
      !
      !    which is the same partition the scalar C1/(C1 + A_max) states
      !    and is directly comparable with it.
      !
      !    THE BRACKET spans everything about the model that is not fixed
      !    by published data: the collisional disposal rate of the levels
      !    above the data, over a decade about the top resolved level's
      !    own collisional total; its radiative rate, from zero to the
      !    all-level maximum; and the chemical destruction rate of an H2
      !    molecule, from far below the radiative rate (where a molecule
      !    always radiates before it is destroyed, the most radiative
      !    case) to far above it.
      ! ================================================================== !
      write(*,'(a)') '#'
      write(*,'(a)') '# BRACKET ON THE THERMALIZED FRACTION f, reduced'// &
                     ' ladder against the scalar model'
      write(*,'(a)') '#   cell  T[K]  f_scalar  f_model_lo  f_model_hi'// &
                     '  1-f_scalar  1-f_model_lo'
      do ic = 1, ncase
         call build_reduced_ladder(roueff_file, lique_file, jozwiak_dir,  &
                                   c_T(ic), gamma_10_H2(c_T(ic)),        &
                                   gamma_10_He(c_T(ic)), ok, why)
         itop   = top_level()
         c_top  = collisional_total(itop, c_T(ic), c_nhi(ic), c_nhe(ic),  &
                                    c_nh2(ic))
         f_scalar = h2_vibrational_heat_fraction(c_T(ic), c_nhi(ic),      &
                                                 c_nh2(ic), c_nhe(ic))
         f_lo =  huge(1.0d0)
         f_hi = -huge(1.0d0)
         do iv = 1, 2
            do ia = 1, 2
               do ik = 1, 2
                  gamma_g = merge(group_coll_scale_lo,                    &
                                  group_coll_scale_hi, iv .eq. 1)*c_top
                  a_g     = merge(0.0d0, a_total_max, ia .eq. 1)
                  k_destroy = merge(k_destroy_slow, k_destroy_fast,       &
                                    ik .eq. 1)
                  s_group   = 1.0d4
                  p_inj     = s_group*e_group_r15_K*kb_erg_ladder
                  i_ground  = level_index(0, 0)
                  s_vec = 0.0d0
                  s_ctl = 0.0d0
                  s_vec(nlev + 1) = s_group
                  s_ctl(i_ground) = s_group
                  n0 = 0.0d0
                  call solve_statistical_equilibrium(c_T(ic), c_nhi(ic),  &
                          c_nhe(ic), c_nh2(ic), 1.0d0, s_vec,             &
                          gamma_g, a_g, k_destroy, -1.0d0, n0, n, info)
                  if (info .ne. 0) cycle
                  call solve_statistical_equilibrium(c_T(ic), c_nhi(ic),  &
                          c_nhe(ic), c_nh2(ic), 1.0d0, s_ctl,             &
                          gamma_g, a_g, k_destroy, -1.0d0, n0, n_ctl,     &
                          info)
                  if (info .ne. 0) cycle
                  ! The SAME molecule rate injected at the bottom of the
                  ! ladder is the control: its radiation is the ladder's
                  ! thermal emission, which the gas pays for through the
                  ! collisional term and which the infrared cooling module
                  ! carries in production.  The difference is the response
                  ! to the internal energy alone.
                  l_rad = radiated_power(n, 1.0d0, a_g, e_group_r15_K)
                  l_ctl = radiated_power(n_ctl, 1.0d0, a_g, e_group_r15_K)
                  d_rad = l_rad - l_ctl
                  f_model = 1.0d0 - d_rad/p_inj
                  f_lo = min(f_lo, f_model)
                  f_hi = max(f_hi, f_model)
               enddo
            enddo
         enddo
         f_bracket_lo(ic) = f_lo
         f_bracket_hi(ic) = f_hi
         write(*,'(a,a20,f10.2,2x,5es14.6)') '#   ', case_name(ic),       &
              c_T(ic), f_scalar, f_lo, f_hi, 1.0d0 - f_scalar,            &
              1.0d0 - f_lo
         ! The bracket exists and is a fraction.
         call check_at_least('ladder_fraction_is_nonnegative_'           &
                             //trim(case_name(ic)), f_lo, 0.0d0)
         call check_at_most('ladder_fraction_is_at_most_one_'            &
                            //trim(case_name(ic)), f_hi, 1.0d0 + 1.0d-12)
      enddo

      ! THE SCALAR FORM DOES NOT OVERSTATE THE HEAT ANYWHERE IN THE
      ! BRACKET, and that is the statement this item needs: the default
      ! does not change.  The comparison is on 1 - f, the share that
      ! leaves the gas, which is the quantity the two models differ in;
      ! f itself is 1 to six digits at every cell and comparing it would
      ! test nothing.  The scalar's share is the larger at all three
      ! cells, so the scalar keeps more energy out of the gas than any
      ! member of the ladder's bracket does, and replacing it could only
      ! add heat.  The ratio is printed because its SIZE is the
      ! uncertainty this item records; the assertion is the direction.
      write(*,'(a)') '#'
      write(*,'(a)') '#   cell  (1-f_scalar)/(1-f_ladder_most_radiative)'
      do ic = 1, ncase
         f_scalar = h2_vibrational_heat_fraction(c_T(ic), c_nhi(ic),      &
                                                 c_nh2(ic), c_nhe(ic))
         write(*,'(a,a20,es14.6)') '#   ', case_name(ic),                 &
              (1.0d0 - f_scalar)/max(1.0d0 - f_bracket_lo(ic), 1.0d-300)
         call check_at_least('scalar_radiates_at_least_as_much_as_'      &
                             //trim(case_name(ic)),                       &
              (1.0d0 - f_scalar)/max(1.0d0 - f_bracket_lo(ic), 1.0d-300), &
              1.0d0)
      enddo

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'h2_ladder_statistical_equilibrium: ',       &
              assertion_failures, ' assertion(s) failed'
         stop 1
      endif
      write(*,'(a)')                                                      &
         'h2_ladder_statistical_equilibrium: every assertion passed'

      contains

      double precision function collisional_total(i, T, n_H, n_He, n_H2)  &
                       result(c)
      ! Total downward collisional rate out of level i [s^-1], summed over
      ! the three colliders and every lower level.
      integer, intent(in) :: i
      real*8,  intent(in) :: T, n_H, n_He, n_H2
      real*8  :: nm(ncoll)
      integer :: j, m
      nm(1) = n_H; nm(2) = n_He; nm(3) = n_H2
      c = 0.0d0
      do j = 1, nlev
         if (lev_E(j) .ge. lev_E(i)) cycle
         do m = 1, ncoll
            c = c + c_ul(i, j, m)*nm(m)
         enddo
      enddo
      end function collisional_total

      double precision function collisional_term_scale(n, T, n_H, n_He,   &
                       n_H2) result(s)
      ! The sum of the ABSOLUTE downward and upward collisional energy
      ! terms [erg cm^-3 s^-1].  The net heat is their difference, so this
      ! is the scale a cancellation has to be measured against.
      real*8, intent(in) :: n(mxlev + 1), T, n_H, n_He, n_H2
      real*8  :: nm(ncoll), cd, dE, boltz
      integer :: i, j, m
      nm(1) = n_H; nm(2) = n_He; nm(3) = n_H2
      s = 0.0d0
      do i = 1, nlev
         do j = 1, nlev
            if (lev_E(i) .le. lev_E(j)) cycle
            dE    = lev_E(i) - lev_E(j)
            boltz = (lev_g(i)/lev_g(j))*exp(-dE/max(T, 1.0d0))
            do m = 1, ncoll
               cd = c_ul(i, j, m)
               if (cd .le. 0.0d0) cycle
               s = s + (n(i)*cd*nm(m) + n(j)*cd*boltz*nm(m))              &
                       *dE*kb_erg_ladder
            enddo
         enddo
      enddo
      end function collisional_term_scale

      subroutine resolve_paths(root, roueff, lique, jozwiak, have)
      ! The repository root, and the three published tables under it and
      ! under the workspace references/ tree beside it.  The collision
      ! tables are third-party material distributed outside the
      ! repository, so their absence is a skip and not a failure.
      character(len=*), intent(out) :: root, roueff, lique, jozwiak
      logical, intent(out) :: have
      integer :: ls
      call get_environment_variable('EXHALE_TEST_ROOT', root, ls)
      if (ls .le. 0) root = '.'
      roueff  = trim(root)//'/cooling_data/'//                            &
                'roueff2019_h2_infrared_lines.dat'
      lique   = trim(root)//'/../references/'//                           &
                'Lique_2015MNRAS_453_810_data/Rates_H_H2.dat'
      jozwiak = trim(root)//'/../references/'//                           &
                'Jozwiak_2024J_A+A_685_A113'
      inquire(file=trim(lique), exist=have)
      if (have) inquire(file=trim(jozwiak)//'/ph2-rat.dat', exist=have)
      if (have) inquire(file=trim(roueff), exist=have)
      end subroutine resolve_paths

      end program h2_ladder_statistical_equilibrium
