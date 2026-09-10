      program lya_escape_probability
      ! The Ly-alpha escape probability of the two-level closure is the one
      ! of the published static plane-parallel damping-wing slab solution.
      !
      ! Production routines exercised:
      ! src/modules/radiation/lya_rt.f90 --
      !   lya_slab_source_position_factor (the source-position factor of the
      !   slab solution), lya_slab_scatterings_before_escape (its mean
      !   number of scatterings), lya_wing_escape_probability (the escape
      !   chance per emission that closes the 2p budget),
      !   lya_starward_escape_fraction (the split of the escaping population
      !   between the two faces), lya_wing_domain_record (the record of the
      !   cells that left the solution's stated domain), and
      !   jlya_escape_prob, which is what the H(n=2) pumping reads.
      !
      ! PHYSICS AND REFERENCES.
      !  * Ly-alpha scatters with partial redistribution (R_II-A, Hummer
      !    1962): a photon shifts by about one Doppler width per scattering
      !    and cannot be re-emitted directly in the far wing, so escape from
      !    a thick static slab is a random walk in frequency.  The mean
      !    number of scatterings before escape is then proportional to the
      !    line-centre optical depth and INDEPENDENT of the Voigt parameter:
      !      <N> = (4 sqrt(6)/pi^2) u_2 B = 0.909316 B
      !    for a mid-plane source in a slab of line-centre optical half
      !    thickness B (Harrington 1973, MNRAS 162, 43, eq. 40, with
      !    u_2 = 1 - 3^-2 + 5^-2 - ... = 0.9159656, Catalan's constant), and
      !      <N> = (4 sqrt(6)/pi^2) [(tau_up+tau_dn)/2] Phi(xi),
      !      Phi(xi) = sum_{n odd} sin(n pi xi)/n^2, xi = tau_dn/(tau_up+tau_dn)
      !    for a source plane anywhere between the two faces (Neufeld 1990,
      !    ApJ 350, 216, eq. 3.27 with his continuum destruction opacity set
      !    to zero).  Both are asserted below against sums written out here.
      !  * <N> counts ABSORPTIONS (Harrington 1973 eq. 38), so a generated
      !    photon is EMITTED 1 + <N> times before it leaves, and the escape
      !    chance per emission that closes the two-level budget
      !    n2p A beta = P is beta = 1/(1 + <N>).  That is exact in the thin
      !    limit as well, which is asserted here, so the closure needs no
      !    interpolation between the thin and thick regimes.
      !  * The two faces of the slab share ONE escaping population, split by
      !    Neufeld eq. (2.25) in the ratio (tau_0 -+ tau_s)/(2 tau_0); the
      !    two shares sum to one, which is asserted here, and they are never
      !    two independent escape probabilities to be multiplied.
      !  * The solution is a damping-wing one and holds for
      !    (a tau)^(1/3) > 10 (Neufeld section Va).  The count of cell
      !    visits that left that domain is asserted against a column built
      !    to straddle it.
      !
      ! Tolerances: the slab identities are exact relations between the
      ! production routine and a sum written out here, and are tested at
      ! 1e-12 relative; the source-position factor is a quadrature and is
      ! tested at 1e-9 against a direct 4e6-term summation of its series;
      ! the closure identity recovered from the production field carries the
      ! round-off of that reconstruction and is tested at 1e-10.
      use global_parameters, only: N, Ng, R0, r, dr_j, v0, pi, appx_mth,   &
                                   F_Lya_star, dv_star_lya,               &
                                   lya_star_boost, lya_bottom_absorber,   &
                                   gamma2_bal, c_light, kb_erg, mu, hp_erg
      use hydrogen_n2_rates, only: nu_lya, C_lya, A_2p1s, g1s, g2p,       &
                                   alpha_B_hydrogen, c1s2p_rate,          &
                                   n2p_destruction_rate
      use lya_rt, only: lya_rt_allocate_arrays, jint_arr,                 &
                        lya_line_center_optical_depth,                    &
                        lya_slab_source_position_factor,                  &
                        lya_slab_scatterings_before_escape,               &
                        lya_wing_escape_probability,                      &
                        lya_starward_escape_fraction,                     &
                        lya_wing_domain_reset, lya_wing_domain_record,    &
                        jlya_escape_prob
      use assertion_report
      implicit none

      ! 4 sqrt(6)/pi^2 and Catalan's constant, written out here rather than
      ! taken from the module, so that the assertions are an independent
      ! statement of Harrington eq. (40).
      real*8, parameter :: k_slab = 4.0d0*sqrt(6.0d0)/(4.0d0*atan(1.0d0))**2
      real*8, parameter :: u2     = 0.9159655941772190d0

      integer, parameter :: ncell   = 60
      real*8,  parameter :: dr_grid = 0.02d0     ! cell width [R0]
      real*8,  parameter :: r0_cm   = 1.0d10     ! planet radius [cm]
      real*8,  parameter :: t_layer = 8.0d3      ! layer temperature [K]
      real*8,  parameter :: h_scale = 0.06d0     ! H I scale height [R0]
      real*8,  parameter :: tau_base_target = 1.0d8

      integer, parameter :: ndepth = 3
      real*8, parameter :: tau_wing(ndepth) = (/ 1.0d8, 1.0d6, 1.0d4 /)

      real*8, allocatable :: nHI(:), T_K(:), nHII(:), ne(:), v_in(:)
      real*8, allocatable :: tau(:), dtau(:), tau_out(:), Jlya(:)
      real*8 :: DnuD, kap_base, avoigt, tau_c, nbar, beta, want, worst
      real*8 :: Jpref, P_src, D2p, budget, atau13, thr
      integer :: j, jlo, jhi, id, nout, nseen, nout_want
      character(len=64) :: label

      ! ------------------------------------------------------------------ !
      ! (1) The source-position factor of the slab solution.
      ! ------------------------------------------------------------------ !
      ! Phi(1/2) is Catalan's constant, the mid-plane case of Harrington
      ! eq. (40); the factor is symmetric under exchanging the two faces;
      ! and a source lying on a face is not trapped at all.
      call check_relative('lya_slab_midplane_factor_is_catalan',           &
           lya_slab_source_position_factor(0.5d0), u2, 1.0d-12)
      call check_absolute('lya_slab_factor_face_symmetry',                 &
           lya_slab_source_position_factor(0.3d0)                          &
         - lya_slab_source_position_factor(0.7d0), 0.0d0, 1.0d-14)
      call check_absolute('lya_slab_factor_on_a_face_vanishes',            &
           lya_slab_source_position_factor(0.0d0), 0.0d0, 1.0d-99)
      ! Against the series it sums, term by term.
      worst = 0.0d0
      do id = 1, 5
         want  = phi_series(0.1d0*dble(id))
         worst = max(worst, abs(lya_slab_source_position_factor            &
                                (0.1d0*dble(id)) - want)/want)
      enddo
      call check_absolute('lya_slab_factor_matches_its_series', worst,     &
                          0.0d0, 1.0d-9)

      ! ------------------------------------------------------------------ !
      ! (2) The mean number of scatterings, Harrington eq. (40).
      ! ------------------------------------------------------------------ !
      ! A mid-plane source (the reflecting lower boundary of the default
      ! closure, mirrored) in a slab of half thickness B = tau: <N> is
      ! 0.909316 tau, proportional to the depth and with NO dependence on
      ! the Voigt parameter.
      do id = 1, ndepth
         nbar = lya_slab_scatterings_before_escape(tau_wing(id),           &
                                                   tau_wing(id))
         write(label,'(a,i0)') 'lya_slab_scatterings_harrington_eq40_', id
         call check_relative(label, nbar, k_slab*u2*tau_wing(id), 1.0d-12)
         write(label,'(a,i0)') 'lya_wing_escape_is_one_over_one_plus_N_',id
         call check_relative(label,                                        &
              lya_wing_escape_probability(tau_wing(id), tau_wing(id)),     &
              1.0d0/(1.0d0 + k_slab*u2*tau_wing(id)), 1.0d-12)
      enddo
      call check_absolute('lya_slab_scatterings_face_symmetry',            &
           lya_slab_scatterings_before_escape(1.0d7, 3.0d6)                &
         - lya_slab_scatterings_before_escape(3.0d6, 1.0d7),               &
           0.0d0, 1.0d-8)

      ! ------------------------------------------------------------------ !
      ! (3) The thin limit.
      ! ------------------------------------------------------------------ !
      ! A photon created in a slab of no depth escapes on its first flight.
      call check_relative('lya_wing_escape_thin_limit',                    &
           lya_wing_escape_probability(0.0d0, 0.0d0), 1.0d0, 1.0d-99)
      call check_relative('lya_wing_escape_thin_limit_small_depth',        &
           lya_wing_escape_probability(1.0d-9, 1.0d-9),                    &
           1.0d0/(1.0d0 + k_slab*u2*1.0d-9), 1.0d-12)
      ! And it never exceeds one at any depth.
      worst = 0.0d0
      do id = 1, 14
         worst = max(worst, lya_wing_escape_probability                    &
                            (10.0d0**(dble(id) - 4.0d0),                   &
                             10.0d0**(dble(id) - 4.0d0)))
      enddo
      call check_at_most('lya_wing_escape_never_exceeds_one', worst, 1.0d0)

      ! ------------------------------------------------------------------ !
      ! (4) The two faces share one escaping population.
      ! ------------------------------------------------------------------ !
      call check_relative('lya_face_partition_sums_to_one',                &
           lya_starward_escape_fraction(2.0d7, 5.0d6)                      &
         + lya_starward_escape_fraction(5.0d6, 2.0d7), 1.0d0, 1.0d-14)
      call check_relative('lya_face_partition_midplane_is_half',           &
           lya_starward_escape_fraction(4.0d6, 4.0d6), 0.5d0, 1.0d-14)
      ! A source plane lying on the absorbing face loses everything there.
      call check_absolute('lya_face_partition_source_on_the_face',         &
           lya_starward_escape_fraction(3.0d6, 0.0d0), 0.0d0, 1.0d-99)

      ! ------------------------------------------------------------------ !
      ! (5) The field the H(n=2) pumping reads carries that probability.
      ! ------------------------------------------------------------------ !
      N  = ncell
      R0 = r0_cm
      v0 = 1.0d5
      appx_mth   = 'None'
      dv_star_lya = 70.0d0
      F_Lya_star  = 0.0d0          ! internal sources alone
      lya_star_boost      = 1.0d0
      lya_bottom_absorber = .false.
      gamma2_bal          = 0.0d0

      jlo = 1 - Ng
      jhi = N + Ng
      allocate(r(jlo:jhi), dr_j(jlo:jhi))
      allocate(nHI(jlo:jhi), T_K(jlo:jhi), nHII(jlo:jhi), ne(jlo:jhi))
      allocate(v_in(jlo:jhi), tau(jlo:jhi), dtau(jlo:jhi))
      allocate(tau_out(jlo:jhi), Jlya(jlo:jhi))
      call lya_rt_allocate_arrays

      do j = jlo, jhi
         r(j)    = 1.0d0 + dble(j - 1)*dr_grid
         dr_j(j) = dr_grid
      enddo
      T_K  = t_layer
      DnuD     = nu_lya*sqrt(2.0d0*kb_erg*t_layer/mu)/c_light
      kap_base = tau_base_target/(h_scale*R0)
      do j = jlo, jhi
         nHI(j)  = kap_base*DnuD/C_lya*exp(-(r(j) - r(jlo))/h_scale)
         nHII(j) = 1.0d-3*nHI(j)
         ne(j)   = nHII(j)
      enddo
      ! No velocity gradient, so the Sobolev channel is closed and the
      ! escape probability of the field is the slab one alone.
      v_in = 0.0d0

      call lya_line_center_optical_depth(T_K, nHI, tau, dtau, tau_out)
      call lya_wing_domain_reset
      call jlya_escape_prob(T_K, nHI, nHII, ne, v_in, Jlya, tau)

      ! The internal field is the trapped source function of the two-level
      ! budget n2p = P/(A beta + D) closed with THIS escape probability:
      !   Jint n1s (A beta + D) = (2 h nu^3/c^2)(g1s/g2p) P (1 - beta) .
      ! Written as the equality of the two sides, so that no step of the
      ! check divides by 1 - beta, which vanishes in the thin outer cells.
      ! The published escape probability is the only one that satisfies it.
      Jpref = 2.0d0*hp_erg*nu_lya**3.0d0/c_light**2.0d0
      worst = 0.0d0
      do j = jlo, jhi
         tau_c = max(tau(j) - 0.5d0*dtau(j), 0.0d0)
         beta  = lya_wing_escape_probability(tau_c, tau_c)
         P_src = alpha_B_hydrogen(t_layer)*ne(j)*nHII(j)                   &
               + c1s2p_rate(t_layer)*ne(j)*nHI(j)
         D2p   = n2p_destruction_rate(t_layer, ne(j), 0.0d0, 0.0d0)
         budget = jint_arr(j)*nHI(j)*(A_2p1s*beta + D2p)
         want   = Jpref*(g1s/g2p)*P_src*(1.0d0 - beta)
         worst  = max(worst, abs(budget - want)/want)
      enddo
      call check_absolute('lya_closure_carries_the_published_escape',      &
                          worst, 0.0d0, 1.0d-13)

      ! ------------------------------------------------------------------ !
      ! (6) The domain record of the damping-wing solution.
      ! ------------------------------------------------------------------ !
      ! (a tau)^(1/3) > 10 (Neufeld section Va).  The exponential column
      ! above straddles that limit; the record must count exactly the cells
      ! that fail it, on the full slab depth 2 tau_c of the mirrored source.
      avoigt = (A_2p1s/(4.0d0*pi))/DnuD
      nout_want = 0
      do j = jlo, jhi
         tau_c  = max(tau(j) - 0.5d0*dtau(j), 0.0d0)
         atau13 = (avoigt*2.0d0*tau_c)**(1.0d0/3.0d0)
         if (atau13 .lt. 10.0d0) nout_want = nout_want + 1
      enddo
      call lya_wing_domain_record(nout, nseen, want, thr)
      call check_relative('lya_wing_domain_cells_seen', dble(nseen),       &
                          dble(jhi - jlo + 1), 1.0d-14)
      call check_relative('lya_wing_domain_cells_out_of_domain',           &
                          dble(nout), dble(nout_want), 1.0d-14)
      call check_relative('lya_wing_domain_threshold', thr, 10.0d0,        &
                          1.0d-14)
      ! Reference-only: the worst (a tau)^(1/3) of this column, and that the
      ! record found cells on both sides of the limit.
      call check_positive('lya_wing_domain_worst_atau_one_third', want)
      call check_positive('lya_wing_domain_straddles_the_limit',           &
                          dble(nout)*dble(nseen - nout))

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'lya_escape_probability: ',                   &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'lya_escape_probability: all assertions passed'

      contains

      ! ---------------------------------------------------------------- !

      double precision function phi_series(xi) result(s)
      ! Phi(xi) = sum_{n odd} sin(n pi xi)/n^2 summed term by term, the
      ! series the production routine replaces by its integrated form.  Four
      ! million terms leave a truncation below 1e-12 at the arguments used
      ! here; an independent statement of the same function.
      real*8, intent(in) :: xi
      integer, parameter :: nmax = 4000001
      integer :: n
      s = 0.0d0
      do n = nmax, 1, -2
         s = s + sin(dble(n)*pi*xi)/dble(n)**2
      enddo
      end function phi_series

      end program lya_escape_probability
