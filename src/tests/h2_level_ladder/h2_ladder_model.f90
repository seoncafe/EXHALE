      module h2_ladder_model
      ! A REDUCED STATISTICAL EQUILIBRIUM OF THE H2 GROUND-STATE LADDER,
      ! built only from published data, to bracket the scalar thermalized
      ! fraction h2_vibrational_heat_fraction of
      ! src/modules/lower_atmosphere/h2_vibrational_relaxation.f90.
      !
      ! WHAT THE SCALAR FRACTION APPROXIMATES.  The exact statement of the
      ! same physics is the net collisional heat of the ladder,
      !
      !     Q = sum over level pairs (u, l) and colliders M of
      !         [ n_u C_ul(M) - n_l C_lu(M) ] n_M (E_u - E_l) ,
      !
      ! with n_u a statistical equilibrium of the whole X ladder and the
      ! upward coefficients from the downward ones by detailed balance.  It
      ! carries no critical density, no single level standing in for the
      ! cascade and no branching fraction at all.  This module evaluates
      ! that expression on the levels the published collision data reach.
      !
      ! WHY IT IS REDUCED, AND WHERE THE MISSING DATA GO.  The three
      ! quantum calculations in hand stop well below the dissociation
      ! limit: Lique (2015) distributes H2-H state to state for v <= 3,
      ! Jozwiak et al. (2024) H2-He for the 53 levels below 15000 cm^-1,
      ! and the molecules whose energy this fraction disposes of are born
      ! ABOVE both, at v = 10 to 14 (the three-body association R15) and at
      ! v = 5 to 6 (the H3+ dissociative recombination R6).  Rather than
      ! invent rates for the levels with none, everything above the
      ! resolved set is carried as ONE group with an explicitly uncertain
      ! disposal, and the uncertainty of the answer is the spread of the
      ! answer over that group's stated range.  This is the substitute the
      ! level-resolved cascade is blocked on, not the cascade itself.
      !
      ! THE DATA, each read from the file its own publication distributes:
      !   levels, weights and A values
      !       cooling_data/roueff2019_h2_infrared_lines.dat, Roueff,
      !       Abgrall, Czachorowski, Pachucki, Puchalski & Komasa (2019),
      !       A&A 630, A58, table 2.  The same file the code's own ladder
      !       and infrared line list are reduced from, so the model and the
      !       production module describe one molecule.  Its all-level
      !       maximum total decay rate agrees with the complete quadrupole
      !       set of Wolniewicz, Simbotin & Dalgarno (1998), ApJS 115, 293,
      !       to 1.1e-04 relative (READ, the validity block of
      !       h2_vibrational_relaxation).
      !   H collider
      !       ../references/Lique_2015MNRAS_453_810_data/Rates_H_H2.dat,
      !       Lique (2015), MNRAS 453, 810: k(v j -> v' j') on
      !       T = 100 to 5000 K in steps of 100 K.
      !   He collider
      !       ../references/Jozwiak_2024J_A+A_685_A113/{ph2,oh2}-{lev,rat}.dat,
      !       Jozwiak, Thibault, Viel, Wcislo & Lique (2024), A&A 685,
      !       A113: 43 temperatures from 20 to 8000 K, both directions.
      !   H2 collider
      !       Le Bourlot, Pineau des Forets & Flower (1999), MNRAS 305,
      !       802.  Their state-to-state table is not in this tree; what is
      !       in the tree is its reduction to the thermal v = 1 -> v' = 0
      !       coefficient, gamma_10_H2 of the production module.  The
      !       helium state-to-state matrix is therefore used for the SHAPE
      !       of the H2 collider and scaled to that published thermal
      !       value.  THIS IS AN ASSUMPTION OF THE MODEL and is stated
      !       wherever the model's answer is quoted; at the cells this
      !       model is run on, H2 carries 6e-05 of the collider sum
      !       (MEASURED), so it cannot carry the answer either way.
      !
      ! DEGENERACIES.  g = g_I (2J+1) with g_I = 3 for odd J and 1 for
      ! even J, the weights of the code's own ladder.  Within the helium
      ! files the nuclear-spin factor is one constant per file and cancels
      ! out of every detailed-balance ratio; it enters only the
      ! ortho-to-para transitions, which only the reactive H exchange has,
      ! and there it is the thermodynamically correct weight.
      !
      ! UNITS.  Level energies in K above (v = 0, J = 0); rate
      ! coefficients in cm^3 s^-1; densities in cm^-3; energy rates in
      ! erg cm^-3 s^-1.

      implicit none
      public

      integer, parameter :: mxlev  = 96
      integer, parameter :: ncoll  = 3      ! 1 = H, 2 = He, 3 = H2
      real*8,  parameter :: kb_erg_ladder = 1.380649d-16
      real*8,  parameter :: cm_to_K       = 1.4387768775039337d0

      ! The resolved ladder.
      integer :: nlev = 0
      integer :: lev_v(mxlev), lev_J(mxlev)
      real*8  :: lev_E(mxlev)              ! [K] above (0,0)
      real*8  :: lev_g(mxlev)
      ! Spontaneous decay, upper to lower [s^-1], and the downward
      ! collisional rate coefficient of each collider [cm^3 s^-1].
      real*8  :: a_ul(mxlev, mxlev)
      real*8  :: c_ul(mxlev, mxlev, ncoll)
      ! The all-level maximum total decay rate of the full Roueff ladder
      ! [s^-1], reduced from the same file over every level it carries,
      ! not only the resolved ones.
      real*8  :: a_total_max = 0.0d0
      logical :: ladder_ready = .false.

      contains

      ! ------------------------------------------------------------------ !

      subroutine build_reduced_ladder(roueff_file, lique_file,            &
                                      jozwiak_dir, T, gamma10_H2_pub,     &
                                      gamma10_He_pub, ok, why)
      ! Fill the resolved ladder at one temperature.  The level set is the
      ! one Lique's H2-H table names, which is the collider that carries
      ! three orders more of the collider sum than the other two together;
      ! a level the other files do not reach simply gets no contribution
      ! from them, which leaves every collider's own detailed balance
      ! intact.
      character(len=*), intent(in)  :: roueff_file, lique_file, jozwiak_dir
      real*8,  intent(in)  :: T
      real*8,  intent(in)  :: gamma10_H2_pub, gamma10_He_pub
      logical, intent(out) :: ok
      character(len=*), intent(out) :: why
      real*8  :: scale_h2
      integer :: i, k

      nlev  = 0
      a_ul  = 0.0d0
      c_ul  = 0.0d0
      ladder_ready = .false.
      ok  = .false.
      why = ''

      call collect_lique_levels(lique_file, ok, why)
      if (.not. ok) return
      call read_roueff(roueff_file, ok, why)
      if (.not. ok) return
      call read_lique_rates(lique_file, T, ok, why)
      if (.not. ok) return
      call read_jozwiak_rates(jozwiak_dir, T, ok, why)
      if (.not. ok) return

      ! THE H2 COLLIDER, shaped like the helium one and normalized to the
      ! published Le Bourlot thermal v = 1 -> v' = 0 coefficient.  The
      ! scale is the ratio of the two published thermal values the
      ! production module carries, so the model's H2 collider reproduces
      ! that value by construction.
      scale_h2 = 0.0d0
      if (gamma10_He_pub .gt. 0.0d0)                                      &
         scale_h2 = gamma10_H2_pub/gamma10_He_pub
      do i = 1, nlev
         do k = 1, nlev
            c_ul(i, k, 3) = scale_h2*c_ul(i, k, 2)
         enddo
      enddo

      ladder_ready = (nlev .gt. 1)
      ok = ladder_ready
      if (.not. ok) why = 'the resolved level set is empty'

      end subroutine build_reduced_ladder

      ! ------------------------------------------------------------------ !

      subroutine collect_lique_levels(path, ok, why)
      ! The (v, J) pairs Lique's table names, in the order they first
      ! appear.  Energies and weights are filled from the Roueff table
      ! afterwards.
      character(len=*), intent(in)  :: path
      logical, intent(out) :: ok
      character(len=*), intent(out) :: why
      character(len=4096) :: line
      integer :: iu, ios, v1, j1, v2, j2
      ok = .false.; why = ''
      open(newunit=iu, file=path, status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
         why = 'cannot open '//trim(path); return
      endif
      do
         read(iu,'(a)',iostat=ios) line
         if (ios .ne. 0) exit
         if (len_trim(line) .lt. 20) cycle
         read(line,*,iostat=ios) v1, j1, v2, j2
         if (ios .ne. 0) cycle
         call add_level(v1, j1)
         call add_level(v2, j2)
      enddo
      close(iu)
      ok = (nlev .gt. 1)
      if (.not. ok) why = 'no (v, J) pairs in '//trim(path)
      end subroutine collect_lique_levels

      ! ------------------------------------------------------------------ !

      subroutine add_level(v, J)
      integer, intent(in) :: v, J
      integer :: i
      do i = 1, nlev
         if (lev_v(i) .eq. v .and. lev_J(i) .eq. J) return
      enddo
      if (nlev .ge. mxlev) return
      nlev = nlev + 1
      lev_v(nlev) = v
      lev_J(nlev) = J
      lev_E(nlev) = -1.0d0
      lev_g(nlev) = merge(3.0d0, 1.0d0, mod(J, 2) .eq. 1)*(2*J + 1)
      end subroutine add_level

      integer function level_index(v, J) result(i)
      integer, intent(in) :: v, J
      integer :: k
      i = 0
      do k = 1, nlev
         if (lev_v(k) .eq. v .and. lev_J(k) .eq. J) then
            i = k; return
         endif
      enddo
      end function level_index

      ! ------------------------------------------------------------------ !

      subroutine read_roueff(path, ok, why)
      ! Term energies, weights and spontaneous decay rates of the resolved
      ! levels, from Roueff et al. (2019) table 2 in the fixed-column form
      ! the VizieR catalog distributes.  The upper level of each row
      ! carries its own term energy and weight; (0, 0) is the origin and
      ! (0, 1) is placed by the (0, 3) -> (0, 1) row, exactly as
      ! cooling_data/molecular_infrared_bands.py places them for the
      ! production ladder.  The all-level maximum total decay rate is
      ! reduced over EVERY level of the file and not only the resolved
      ! ones, so it is the same quantity the production module carries.
      character(len=*), intent(in)  :: path
      logical, intent(out) :: ok
      character(len=*), intent(out) :: why
      character(len=400) :: rec
      integer :: iu, ios, vu, ju, vl, jl, iu_lev, il_lev, i
      real*8  :: sigma, aul, t_up, g_up, e01
      ! Total decay rate of every level the file carries, keyed by the
      ! level's term energy through a running sum per (v, J).
      integer, parameter :: mxall = 400
      integer :: nall, all_v(mxall), all_J(mxall), k
      real*8  :: all_A(mxall)

      ok = .false.; why = ''
      e01 = -1.0d0
      nall = 0
      all_A = 0.0d0
      do i = 1, nlev
         lev_E(i) = -1.0d0
      enddo
      open(newunit=iu, file=path, status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
         why = 'cannot open '//trim(path); return
      endif
      do
         read(iu,'(a)',iostat=ios) rec
         if (ios .ne. 0) exit
         if (len_trim(rec) .lt. 164) cycle
         read(rec(2:3),  *, iostat=ios) vu
         if (ios .ne. 0) cycle
         read(rec(5:6),  *, iostat=ios) ju
         if (ios .ne. 0) cycle
         read(rec(8:9),  *, iostat=ios) vl
         if (ios .ne. 0) cycle
         read(rec(11:12),*, iostat=ios) jl
         if (ios .ne. 0) cycle
         read(rec(14:29),*, iostat=ios) sigma
         if (ios .ne. 0) cycle
         read(rec(97:105),*, iostat=ios) aul
         if (ios .ne. 0) cycle
         read(rec(149:159),*, iostat=ios) t_up
         if (ios .ne. 0) cycle
         read(rec(162:164),*, iostat=ios) g_up
         if (ios .ne. 0) cycle

         ! The term energy and weight of the upper level.
         i = level_index(vu, ju)
         if (i .gt. 0) then
            lev_E(i) = t_up
            lev_g(i) = g_up
         endif
         if (vu .eq. 0 .and. ju .eq. 3 .and. vl .eq. 0 .and. jl .eq. 1)   &
            e01 = t_up - sigma*cm_to_K
         ! The total decay rate of the upper level, over the whole file.
         k = 0
         do i = 1, nall
            if (all_v(i) .eq. vu .and. all_J(i) .eq. ju) then
               k = i; exit
            endif
         enddo
         if (k .eq. 0 .and. nall .lt. mxall) then
            nall = nall + 1
            k = nall
            all_v(k) = vu; all_J(k) = ju; all_A(k) = 0.0d0
         endif
         if (k .gt. 0) all_A(k) = all_A(k) + aul
         ! The resolved transition itself.
         iu_lev = level_index(vu, ju)
         il_lev = level_index(vl, jl)
         if (iu_lev .gt. 0 .and. il_lev .gt. 0)                           &
            a_ul(iu_lev, il_lev) = a_ul(iu_lev, il_lev) + aul
      enddo
      close(iu)

      a_total_max = 0.0d0
      do i = 1, nall
         a_total_max = max(a_total_max, all_A(i))
      enddo

      ! The two levels no row has as an upper level.
      i = level_index(0, 0)
      if (i .gt. 0) then
         lev_E(i) = 0.0d0
         lev_g(i) = 1.0d0
      endif
      i = level_index(0, 1)
      if (i .gt. 0) then
         if (e01 .lt. 0.0d0) then
            why = 'the (0,3) -> (0,1) row is missing from '//trim(path)
            return
         endif
         lev_E(i) = e01
         lev_g(i) = 9.0d0
      endif

      do i = 1, nlev
         if (lev_E(i) .lt. 0.0d0) then
            why = 'a Lique level has no Roueff term energy'
            return
         endif
      enddo
      ok = .true.
      end subroutine read_roueff

      ! ------------------------------------------------------------------ !

      subroutine read_lique_rates(path, T, ok, why)
      ! k(v j -> v' j') of the H collider at T, log-interpolated on the
      ! file's own 100 to 5000 K grid and held at the end values outside
      ! it, which is the treatment the production module gives the same
      ! table.  Only the DOWNWARD direction is taken; the upward rates are
      ! built by detailed balance where the equilibrium is assembled.
      character(len=*), intent(in)  :: path
      real*8,  intent(in)  :: T
      logical, intent(out) :: ok
      character(len=*), intent(out) :: why
      integer, parameter :: nT = 50
      real*8  :: tgrid(nT), k(nT), kT
      character(len=4096) :: line
      integer :: iu, ios, i, v1, j1, v2, j2, i1, i2
      ok = .false.; why = ''
      do i = 1, nT
         tgrid(i) = 100.0d0*i
      enddo
      open(newunit=iu, file=path, status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
         why = 'cannot open '//trim(path); return
      endif
      do
         read(iu,'(a)',iostat=ios) line
         if (ios .ne. 0) exit
         if (len_trim(line) .lt. 20) cycle
         read(line,*,iostat=ios) v1, j1, v2, j2, k
         if (ios .ne. 0) cycle
         i1 = level_index(v1, j1)
         i2 = level_index(v2, j2)
         if (i1 .le. 0 .or. i2 .le. 0) cycle
         kT = log_interp(T, tgrid, k, nT)
         if (lev_E(i1) .gt. lev_E(i2)) then
            c_ul(i1, i2, 1) = c_ul(i1, i2, 1) + kT
         else if (lev_E(i2) .gt. lev_E(i1)) then
            ! An upward row; detailed balance rebuilds it from the
            ! downward one, so it is not taken twice.
            continue
         endif
      enddo
      close(iu)
      ok = .true.
      end subroutine read_lique_rates

      ! ------------------------------------------------------------------ !

      subroutine read_jozwiak_rates(dir, T, ok, why)
      ! k(i -> f) of the He collider at T, from the ortho and para files
      ! of the VizieR catalog, on their own 43-point grid.  Both
      ! directions are distributed; only the downward rows are taken, for
      ! the same reason as above.
      character(len=*), intent(in)  :: dir
      real*8,  intent(in)  :: T
      logical, intent(out) :: ok
      character(len=*), intent(out) :: why
      integer, parameter :: nT = 43
      real*8, parameter :: tgrid(nT) = (/                                 &
         20.d0, 30.d0, 40.d0, 50.d0, 60.d0, 70.d0, 80.d0, 90.d0, 100.d0,  &
         120.d0, 140.d0, 160.d0, 180.d0, 200.d0, 250.d0, 300.d0, 350.d0,  &
         400.d0, 450.d0, 500.d0, 550.d0, 600.d0, 650.d0, 700.d0, 750.d0,  &
         800.d0, 850.d0, 900.d0, 950.d0, 1000.d0, 1100.d0, 1200.d0,       &
         1300.d0, 1400.d0, 1500.d0, 1750.d0, 2000.d0, 3000.d0, 4000.d0,   &
         5000.d0, 6000.d0, 7000.d0, 8000.d0 /)
      ok = .false.; why = ''
      call one_jozwiak_set(trim(dir)//'/ph2-lev.dat',                     &
                           trim(dir)//'/ph2-rat.dat', T, tgrid, nT,       &
                           ok, why)
      if (.not. ok) return
      call one_jozwiak_set(trim(dir)//'/oh2-lev.dat',                     &
                           trim(dir)//'/oh2-rat.dat', T, tgrid, nT,       &
                           ok, why)
      end subroutine read_jozwiak_rates

      subroutine one_jozwiak_set(levfile, ratfile, T, tgrid, nT, ok, why)
      character(len=*), intent(in)  :: levfile, ratfile
      real*8,  intent(in)  :: T, tgrid(nT)
      integer, intent(in)  :: nT
      logical, intent(out) :: ok
      character(len=*), intent(out) :: why
      integer, parameter :: mxset = 64
      integer :: set_v(mxset), set_J(mxset), nset
      character(len=4096) :: line
      integer :: iu, ios, idx, i, f, i1, i2, n
      real*8  :: k(nT), kT
      ok = .false.; why = ''
      nset = 0
      open(newunit=iu, file=levfile, status='old', action='read',         &
           iostat=ios)
      if (ios .ne. 0) then
         why = 'cannot open '//trim(levfile); return
      endif
      do
         read(iu,'(a)',iostat=ios) line
         if (ios .ne. 0) exit
         if (len_trim(line) .lt. 5) cycle
         read(line,*,iostat=ios) idx, i, f     ! label, J, v
         if (ios .ne. 0) cycle
         if (idx .lt. 1 .or. idx .gt. mxset) cycle
         nset = max(nset, idx)
         set_J(idx) = i
         set_v(idx) = f
      enddo
      close(iu)

      open(newunit=iu, file=ratfile, status='old', action='read',         &
           iostat=ios)
      if (ios .ne. 0) then
         why = 'cannot open '//trim(ratfile); return
      endif
      do
         read(iu,'(a)',iostat=ios) line
         if (ios .ne. 0) exit
         if (len_trim(line) .lt. 30) cycle
         read(line,*,iostat=ios) n, i, f, k
         if (ios .ne. 0) cycle
         if (i .lt. 1 .or. i .gt. nset) cycle
         if (f .lt. 1 .or. f .gt. nset) cycle
         if (i .eq. f) cycle
         i1 = level_index(set_v(i), set_J(i))
         i2 = level_index(set_v(f), set_J(f))
         if (i1 .le. 0 .or. i2 .le. 0) cycle
         if (lev_E(i1) .le. lev_E(i2)) cycle
         kT = log_interp(T, tgrid, k, nT)
         c_ul(i1, i2, 2) = c_ul(i1, i2, 2) + kT
      enddo
      close(iu)
      ok = .true.
      end subroutine one_jozwiak_set

      ! ------------------------------------------------------------------ !

      double precision function log_interp(T, tg, y, n) result(v)
      ! log y against T on a non-uniform grid, held at the end values
      ! outside it.  A zero or negative entry is returned as zero, which
      ! is what an unconverged transition of either table is.
      real*8,  intent(in) :: T, tg(n), y(n)
      integer, intent(in) :: n
      integer :: i, kk
      real*8  :: tc, frac
      tc = min(max(T, tg(1)), tg(n))
      i  = 1
      do kk = 1, n - 1
         if (tc .ge. tg(kk)) i = kk
      enddo
      if (y(i) .le. 0.0d0 .or. y(i + 1) .le. 0.0d0) then
         v = 0.0d0
         return
      endif
      frac = (tc - tg(i))/(tg(i + 1) - tg(i))
      v = exp(log(y(i)) + frac*(log(y(i + 1)) - log(y(i))))
      end function log_interp

      ! ------------------------------------------------------------------ !

      subroutine solve_statistical_equilibrium(T, n_H, n_He, n_H2,        &
                     escape, s_vec, gamma_group,                         &
                     a_group, k_destroy, dt, n_in, n_out, info)
      ! Populations of the resolved ladder and of the one high-level group
      ! [cm^-3], as the steady state (dt <= 0) or as one backward-Euler
      ! step of length dt from the state n_in.  n_out(1:nlev) are the
      ! resolved levels and n_out(nlev+1) the group.
      !
      ! THE SYSTEM.  For each resolved level i,
      !
      !   dn_i/dt = sum_j [ n_j R_ji - n_i R_ij ]
      !             + (the group's disposal, if i is the top level)
      !             - k_destroy n_i ,
      !
      ! with R_ij the sum over colliders of n_M C_ij(M) in both
      ! directions, the upward coefficient from the downward one by
      ! detailed balance, plus the spontaneous decay n_i A_ij times the
      ! escape probability `escape`.  `escape` = 1 is the optically thin
      ! ladder and 0 the completely trapped one, which is the limit the
      ! thermal-bath test needs.
      !
      ! s_vec [cm^-3 s^-1] is where molecules enter: entry nlev+1 is the
      ! high-level group, entries 1 to nlev the resolved levels.  Injecting
      ! the same molecule rate into the GROUND level instead of the group
      ! is the control the fraction is read against: the system is linear,
      ! so the difference of the two solutions is the response to the
      ! internal energy alone and the thermal background cancels exactly.
      !
      ! THE GROUP is every bound level above the resolved set: the levels
      ! the molecules of R15 (v = 10 to 14) and R6 (v = 5 to 6) are born
      ! in.  It carries its own energy, which enters the heat and the
      ! radiated power rather than the populations, and it empties into the
      ! TOP resolved level at the
      ! collisional rate gamma_group [s^-1] and the radiative rate a_group
      ! [s^-1].  Both are explicitly uncertain: no published set reaches
      ! those levels, and the spread of the answer over their stated range
      ! IS this model's uncertainty.  Nothing climbs back into the group,
      ! which is the statement that its levels stand far above kT.
      !
      ! k_destroy [s^-1] removes a molecule from the ladder with whatever
      ! internal energy it carries: the chemical destruction channel.  It
      ! is also what makes the steady system nonsingular, which is the
      ! physical statement that a molecule does not stay forever.
      real*8,  intent(in)  :: T, n_H, n_He, n_H2
      real*8,  intent(in)  :: escape, s_vec(mxlev + 1)
      real*8,  intent(in)  :: gamma_group, a_group, k_destroy, dt
      real*8,  intent(in)  :: n_in(mxlev + 1)
      real*8,  intent(out) :: n_out(mxlev + 1)
      integer, intent(out) :: info
      real*8, allocatable :: amat(:,:), rhs(:)
      integer, allocatable :: ipiv(:)
      integer :: i, j, nn, itop

      n_out = 0.0d0
      info  = -1
      if (.not. ladder_ready) return

      nn   = nlev + 1
      itop = top_level()
      allocate(amat(nn, nn), rhs(nn), ipiv(nn))
      call build_rate_matrix(T, n_H, n_He, n_H2, escape, gamma_group,     &
                             a_group, k_destroy, itop, amat, nn)

      if (dt .le. 0.0d0) then
         ! Steady state: M n = -s.
         do i = 1, nn
            rhs(i) = -s_vec(i)
         enddo
      else
         ! One backward-Euler step: (I/dt - M) n = n_in/dt + s.
         do i = 1, nn
            do j = 1, nn
               amat(i, j) = -amat(i, j)
            enddo
            amat(i, i) = amat(i, i) + 1.0d0/dt
            rhs(i) = n_in(i)/dt + s_vec(i)
         enddo
      endif

      call dgesv(nn, 1, amat, nn, ipiv, rhs, nn, info)
      if (info .eq. 0) then
         do i = 1, nn
            n_out(i) = rhs(i)
         enddo
      endif
      deallocate(amat, rhs, ipiv)

      end subroutine solve_statistical_equilibrium

      ! ------------------------------------------------------------------ !

      subroutine build_rate_matrix(T, n_H, n_He, n_H2, escape,            &
                     gamma_group, a_group, k_destroy, itop, amat, nn)
      ! The rate matrix M of dn/dt = M n + s, resolved levels first and
      ! the high-level group last.
      real*8,  intent(in)  :: T, n_H, n_He, n_H2, escape
      real*8,  intent(in)  :: gamma_group, a_group, k_destroy
      integer, intent(in)  :: itop, nn
      real*8,  intent(out) :: amat(nn, nn)
      real*8  :: r_ij, r_ji
      integer :: i, j
      amat = 0.0d0
      do i = 1, nlev
         do j = 1, nlev
            if (i .eq. j) cycle
            call pair_rates(i, j, T, n_H, n_He, n_H2, escape, r_ij, r_ji)
            amat(i, i) = amat(i, i) - r_ij
            amat(i, j) = amat(i, j) + r_ji
         enddo
         amat(i, i) = amat(i, i) - k_destroy
      enddo
      amat(itop, nlev + 1) = gamma_group + a_group
      amat(nlev + 1, nlev + 1) = -(gamma_group + a_group + k_destroy)
      end subroutine build_rate_matrix

      ! ------------------------------------------------------------------ !

      subroutine solve_thermal_bath(T, n_H, n_He, n_H2, n_total, n_out,   &
                                    info)
      ! The populations a collision-only bath settles at: no spontaneous
      ! decay, no injection, no chemical destruction, and the ladder's own
      ! total density fixed.  The one physical content of the solve is the
      ! detailed balance of the upward rates against the downward ones, so
      ! this is the test that isolates it.
      !
      ! SOLVED IN THE BOLTZMANN-SCALED VARIABLE n_i = b_i x_i, with
      ! b_i = g_i exp(-E_i/kT).  The exact answer is then x_i independent
      ! of i, and the test reads how far from constant the computed x is.
      ! In the raw variable the populations of a ladder at 800 K span
      ! forty orders, the linear solve resolves the smallest of them only
      ! to its own conditioning, and the departure a test would read would
      ! be that conditioning and not the physics.
      real*8,  intent(in)  :: T, n_H, n_He, n_H2, n_total
      real*8,  intent(out) :: n_out(mxlev + 1)
      integer, intent(out) :: info
      real*8, allocatable :: amat(:,:), rhs(:)
      integer, allocatable :: ipiv(:)
      real*8  :: b(mxlev), zsum, scal
      integer :: i, j, nn, inorm
      n_out = 0.0d0
      info  = -1
      if (.not. ladder_ready) return
      nn = nlev
      allocate(amat(nn, nn), rhs(nn), ipiv(nn))
      call build_bath_matrix(T, n_H, n_He, n_H2, amat, nn)
      zsum = 0.0d0
      do i = 1, nlev
         b(i) = lev_g(i)*exp(-lev_E(i)/max(T, 1.0d0))
         zsum = zsum + b(i)
      enddo
      ! Columns scaled by b_j, and then each balance row divided by its
      ! own b_i.  Detailed balance makes R_ij b_i = R_ji b_j, so the
      ! doubly scaled matrix is symmetric and its entries are all of order
      ! the rates themselves; without the row scaling the high levels'
      ! rows are uniformly forty orders below the low levels' and the
      ! solve resolves them only to its own conditioning.
      do i = 1, nn
         do j = 1, nn
            amat(i, j) = amat(i, j)*b(j)/b(i)
         enddo
      enddo
      ! A collision-only bath conserves molecules, so one row of the
      ! balance is redundant; the total density replaces it.  The row
      ! dropped is the MOST POPULATED level's, whose population the
      ! normalization then fixes almost by itself; dropping any other
      ! level's leaves that level with no equation of its own and its
      ! population read as the residual of the normalization, which is a
      ! subtraction of the ladder's whole population from itself.
      inorm = 1
      do i = 2, nlev
         if (b(i) .gt. b(inorm)) inorm = i
      enddo
      rhs = 0.0d0
      do j = 1, nn
         amat(inorm, j) = b(j)
      enddo
      rhs(inorm) = n_total
      ! Row equilibration.  Gaussian elimination with partial pivoting is
      ! not invariant under a row rescaling, and the rows of this system
      ! span the whole range of the rates: the rows of the high levels
      ! carry downward rates of order unity and those of the low levels
      ! upward rates forty orders below.  Each row is therefore divided by
      ! its own largest entry before the solve.
      do i = 1, nn
         scal = 0.0d0
         do j = 1, nn
            scal = max(scal, abs(amat(i, j)))
         enddo
         if (scal .le. 0.0d0) cycle
         do j = 1, nn
            amat(i, j) = amat(i, j)/scal
         enddo
         rhs(i) = rhs(i)/scal
      enddo
      call dgesv(nn, 1, amat, nn, ipiv, rhs, nn, info)
      if (info .eq. 0) then
         do i = 1, nn
            n_out(i) = rhs(i)*b(i)
         enddo
      endif
      deallocate(amat, rhs, ipiv)
      end subroutine solve_thermal_bath

      subroutine build_bath_matrix(T, n_H, n_He, n_H2, amat, nn)
      real*8,  intent(in)  :: T, n_H, n_He, n_H2
      integer, intent(in)  :: nn
      real*8,  intent(out) :: amat(nn, nn)
      real*8  :: r_ij, r_ji
      integer :: i, j
      amat = 0.0d0
      do i = 1, nlev
         do j = 1, nlev
            if (i .eq. j) cycle
            call pair_rates(i, j, T, n_H, n_He, n_H2, 0.0d0, r_ij, r_ji)
            amat(i, i) = amat(i, i) - r_ij
            amat(i, j) = amat(i, j) + r_ji
         enddo
      enddo
      end subroutine build_bath_matrix

      ! ------------------------------------------------------------------ !

      subroutine pair_rates(i, j, T, n_H, n_He, n_H2, escape, r_ij, r_ji)
      ! The total rate i -> j and j -> i [s^-1] between two resolved
      ! levels: the collisional terms of the three colliders in both
      ! directions, and the spontaneous decay of the upper one times the
      ! escape probability.
      integer, intent(in)  :: i, j
      real*8,  intent(in)  :: T, n_H, n_He, n_H2, escape
      real*8,  intent(out) :: r_ij, r_ji
      real*8 :: nm(ncoll), cd, dE, boltz
      integer :: m
      nm(1) = n_H; nm(2) = n_He; nm(3) = n_H2
      r_ij = 0.0d0
      r_ji = 0.0d0
      if (lev_E(i) .gt. lev_E(j)) then
         dE    = lev_E(i) - lev_E(j)
         boltz = (lev_g(i)/lev_g(j))*exp(-dE/max(T, 1.0d0))
         do m = 1, ncoll
            cd = c_ul(i, j, m)
            if (cd .le. 0.0d0) cycle
            r_ij = r_ij + cd*nm(m)                 ! down, i -> j
            r_ji = r_ji + cd*boltz*nm(m)           ! up, by detailed balance
         enddo
         r_ij = r_ij + a_ul(i, j)*escape
      else
         dE    = lev_E(j) - lev_E(i)
         boltz = (lev_g(j)/lev_g(i))*exp(-dE/max(T, 1.0d0))
         do m = 1, ncoll
            cd = c_ul(j, i, m)
            if (cd .le. 0.0d0) cycle
            r_ji = r_ji + cd*nm(m)                 ! down, j -> i
            r_ij = r_ij + cd*boltz*nm(m)           ! up, by detailed balance
         enddo
         r_ji = r_ji + a_ul(j, i)*escape
      endif
      end subroutine pair_rates

      ! ------------------------------------------------------------------ !

      integer function top_level() result(itop)
      ! The highest resolved level, where the group lands its molecules.
      integer :: i
      itop = 1
      do i = 2, nlev
         if (lev_E(i) .gt. lev_E(itop)) itop = i
      enddo
      end function top_level

      ! ------------------------------------------------------------------ !

      double precision function net_collisional_heat(n, T,                &
                       n_H, n_He, n_H2, gamma_group, e_group_K)           &
                       result(q)
      ! The exact form the scalar fraction approximates
      ! [erg cm^-3 s^-1]: downward minus upward collisional rate, summed
      ! over level pairs and colliders, times the energy gap.  Positive
      ! means the ladder is heating the gas.  The group's collisional
      ! disposal into the top resolved level carries its own gap.
      real*8, intent(in) :: n(mxlev + 1), T, n_H, n_He, n_H2
      real*8, intent(in) :: gamma_group, e_group_K
      real*8 :: nm(ncoll), cd, dE, boltz, up, dn
      integer :: i, j, m, itop
      q = 0.0d0
      nm(1) = n_H; nm(2) = n_He; nm(3) = n_H2
      do i = 1, nlev
         do j = 1, nlev
            if (lev_E(i) .le. lev_E(j)) cycle
            dE    = lev_E(i) - lev_E(j)
            boltz = (lev_g(i)/lev_g(j))*exp(-dE/max(T, 1.0d0))
            do m = 1, ncoll
               cd = c_ul(i, j, m)
               if (cd .le. 0.0d0) cycle
               dn = n(i)*cd*nm(m)
               up = n(j)*cd*boltz*nm(m)
               q  = q + (dn - up)*dE*kb_erg_ladder
            enddo
         enddo
      enddo
      itop = top_level()
      q = q + n(nlev + 1)*gamma_group*(e_group_K - lev_E(itop))           &
              *kb_erg_ladder
      end function net_collisional_heat

      ! ------------------------------------------------------------------ !

      double precision function radiated_power(n, escape,                 &
                       a_group, e_group_K) result(l)
      ! The energy the ladder sends away in its quadrupole lines and the
      ! group sends away in its own decay [erg cm^-3 s^-1].
      real*8, intent(in) :: n(mxlev + 1), escape, a_group
      real*8, intent(in) :: e_group_K
      integer :: i, j, itop
      l = 0.0d0
      do i = 1, nlev
         do j = 1, nlev
            if (a_ul(i, j) .le. 0.0d0) cycle
            l = l + n(i)*a_ul(i, j)*escape                                &
                    *(lev_E(i) - lev_E(j))*kb_erg_ladder
         enddo
      enddo
      itop = top_level()
      l = l + n(nlev + 1)*a_group*(e_group_K - lev_E(itop))*kb_erg_ladder
      end function radiated_power

      ! ------------------------------------------------------------------ !

      double precision function stored_internal_energy(n, e_group_K)      &
                       result(u)
      ! Internal energy held by the ladder and the group [erg cm^-3].
      real*8, intent(in) :: n(mxlev + 1), e_group_K
      integer :: i
      u = n(nlev + 1)*e_group_K*kb_erg_ladder
      do i = 1, nlev
         u = u + n(i)*lev_E(i)*kb_erg_ladder
      enddo
      end function stored_internal_energy

      ! ------------------------------------------------------------------ !

      double precision function boltzmann_departure(n, T) result(d)
      ! The largest relative departure of the resolved populations from
      ! the Boltzmann distribution at T, read in the scaled variable
      ! n_i/(g_i exp(-E_i/kT)), which is constant in the thermal limit
      ! whatever the ladder's dynamic range.  Zero is that limit.
      real*8, intent(in) :: n(mxlev + 1), T
      real*8 :: x(mxlev), xmean
      integer :: i
      d = 0.0d0
      xmean = 0.0d0
      do i = 1, nlev
         x(i) = n(i)/(lev_g(i)*exp(-lev_E(i)/max(T, 1.0d0)))
         xmean = xmean + x(i)
      enddo
      if (nlev .lt. 1 .or. xmean .le. 0.0d0) then
         d = huge(1.0d0)
         return
      endif
      xmean = xmean/dble(nlev)
      do i = 1, nlev
         d = max(d, abs(x(i) - xmean)/xmean)
      enddo
      end function boltzmann_departure

      end module h2_ladder_model
