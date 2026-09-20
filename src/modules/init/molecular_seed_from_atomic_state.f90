      module molecular_seed
      ! THE MOLECULAR STATE IMPLIED BY A SOLVED ATOMIC STATE OF THE SAME
      ! WIND, written as an initialization file and nothing else.
      !
      ! WHY IT EXISTS. A molecular run carries four unknowns an atomic run
      ! does not (H2, H2+, H3+, HeH+), so the restart contract calls an
      ! atomic state file "a state of another equation set" and refuses it
      ! (load_IC, tokens mol/carrier/carrier_newton/oxychem/iontrans). The
      ! refusal is right: the file does not carry the molecular unknowns.
      ! The opening is therefore to WRITE them -- to take the atomic
      ! solution's rho, v, p and its 34 species columns, state a hydrogen
      ! partition explicitly, close the equation of state on the result, and
      ! write the 38-column state that the molecular run can load.
      !
      ! WHAT IT IS NOT. The written state is an initialization seed and says
      ! so in its own header (certified=F, cert_reason=molecular_seed,
      ! mode=init). It is not stationary under the molecular equations: the
      ! partition is imposed, the molecular ions start at zero, and the
      ! first chemical sweep moves all of them. Nothing here claims that
      ! only the chemistry will move.
      !
      ! THE PARTITION. The fraction of hydrogen NUCLEI bound into H2 is the
      ! one the lower-atmosphere handoff states,
      !
      !     x2 = 2 q (1 + He/H)/(1 + q),    q = q_H2_base,
      !
      ! the same conversion from a mixture mixing ratio to a nucleus
      ! fraction that composition.f90 (base_h2_nuclei_fraction) and
      ! lower_column.f90 (mu_mixture) use, so there is one definition of
      ! "how molecular is the base". It is applied to each cell's own
      ! hydrogen-nucleus budget, taken from the production census
      ! (element_nuclei_and_charge), and capped by the neutral hydrogen that
      ! cell actually has:
      !
      !     delta = min(x2 n_H,nuclei , n_HI),  n_H2 += delta/2, n_HI -= delta
      !
      ! which is the form set_IC uses for a cold molecular start. A fraction
      ! of the NEUTRAL budget would be a different statement and would agree
      ! with this one only in fully neutral gas; in the ionized wind above
      ! the front the two differ by the ionization fraction.
      !
      ! Carrying a base partition through the whole column is an
      ! initialization choice and not a boundary statement: the cap is what
      ! keeps the hot, ionized cells atomic, and the molecular ions are left
      ! at zero for the first sweep to find their roots.
      !
      ! THE THERMODYNAMIC INVARIANT, STATED. Two H atoms into one H2
      ! molecule keep the nuclei and the mass but remove half a particle per
      ! converted pair, so with the electrons and the other species
      ! unchanged n_new = n_old - delta/2 and p = n k_B T: rho, v, p and T
      ! cannot all be kept. Which one is dropped is a choice and is written
      ! into the file:
      !
      !   invariant=p  keep rho, v, p; T = p/((n_tot + n_e) k_B) follows from
      !                the new particle census. Keeps the pressure-force
      !                profile of the source wind, which is what a nearly
      !                hydrostatic column is held up by.
      !   invariant=T  keep rho, v, T; p = (n_tot + n_e) k_B T follows. Keeps
      !                the rate coefficients of the source state and moves
      !                the pressure gradient instead.
      !
      ! Neither keeps the chemical or the energy equilibrium of the source.
      ! The energy is closed with the production caloric equation of state
      ! in both cases (W_to_U / U_to_W of Conversion), and the round trip
      ! through it is measured and reported.
      !
      ! WHERE THE REQUESTED FRACTION COMES FROM, and why there are three
      ! answers. EXHALE_MOLECULAR_SEED_X2 takes
      !
      !   (absent)  the handoff's x2, one number for the whole column: the
      !             base partition extended to the outer boundary. It is a
      !             boundary statement carried where the boundary does not
      !             reach, and at x2 near 1 it takes every neutral hydrogen
      !             atom the wind has.
      !   local     the H2 content that HOLDS in each cell: the smaller of
      !             the thermochemical fit q_H2(p, T) (Visscher/Koskinen,
      !             converted with x2 = 2 q (1 + He/H)/(1 + q) and clipped
      !             at the element-ratio ceiling, the branch set_IC takes
      !             when no handoff states the base) and the root of that
      !             cell's own H2 carrier row, production over loss rate,
      !             with the photodissociation, the photoionization and the
      !             ion channels in it (carrier_h2_chemical_root).
      !
      !             WHY THE SMALLER OF THE TWO, and it is a physical
      !             statement and not a safeguard. The fit is the balance of
      !             the three-body association against the thermal
      !             dissociation; it is what holds at the base, where the
      !             gas is dense, neutral and shielded. In an irradiated
      !             outer wind it is not the balance that holds -- the gas
      !             is there at a low temperature because it expanded, the
      !             star dissociates its H2, and half its hydrogen is
      !             ionized -- and the fit knows none of that and returns a
      !             large H2 because the temperature alone is low. MEASURED
      !             on LHS 1140 b (item L7e, docs/lhs1140b_stationary_L7e_
      !             20260915.md section 13): the fit alone put 13 to 21
      !             percent of the gas into H2 from 5.9 R_p out to 29 R_p,
      !             80 to 360 times the root of its own row, with atomic H
      !             four decades BELOW the molecule it is made from, and no
      !             substep of the carrier operator could integrate that
      !             state, so the relaxation of the stationary route kept no
      !             step at all. A seed may not assert an equilibrium that
      !             does not hold in the cell it asserts it for; where the
      !             two disagree the one that holds is the smaller, because
      !             each is a balance the gas would reach if it were the
      !             only one, and the faster of the two sets the content.
      !             The crossover radius and the number of cells the root
      !             governs are reported and written into the state files.
      !
      !             THE ROOT NEEDS THE MOLECULAR IONS. The production of H2
      !             runs through H2+ and HeH+, which an atomic state does
      !             not carry, so the root is taken after the fit has been
      !             applied and one equilibrium sweep has solved the ions
      !             against it; that sweep's composition is then discarded
      !             and the partition is re-applied to the ATOMIC state, so
      !             the ionization state the seed writes is the one it was
      !             handed and not one this module solved.
      !   a number in [0,1]  that fraction, for every cell. 0 is the
      !             conversion identity (nothing is transferred, so every
      !             species column, the nuclei and the chosen primitive
      !             invariants come back unchanged).
      !
      ! A zero H2 column is a legitimate test and a poor seed: with the
      ! carrier transported, a LOADED state's H2 column is the authority and
      ! the sweep is handed it rather than solving for it
      ! (ionization_equilibrium, "with the carriers transported, their
      ! partition is not a local root any more"), so H2 = 0 stays 0. The value used, and where it came from, is written
      ! into both state files.
      !
      ! HOW IT IS ASKED FOR. Environment variable EXHALE_MOLECULAR_SEED set
      ! to the directory holding the atomic state (its Hydro_ioniz_IC.txt
      ! and Ion_species_IC.txt). The run is configured by the TARGET
      ! molecular input.inp with "Load IC? True"; the loader reads the
      ! atomic pair from that directory, this module builds the molecular
      ! state, the run writes output/Hydro_ioniz_IC.txt and
      ! output/Ion_species_IC.txt and stops. EXHALE_MOLECULAR_SEED_INVARIANT
      ! selects p (default) or T. Documented in docs/input_schema.md
      ! appendix D.
      !
      ! Item L7 of docs/PLAN_20260913_lhs_stationary.md.

      use global_parameters
      use species_table, only: isp_HI, isp_HII, isp_H2, n_mion
      use element_census, only: element_nuclei_and_charge, n_element,     &
                                ie_H, ie_He
      use composition, only: get_species_densities, comp_T_from_p,        &
                             comp_p_from_T, h2_mixing_ratio_base,         &
                             h2_mixing_ratio_ceiling,                     &
                             base_h2_nuclei_fraction,                     &
                             base_h2_composition_imposed
      ! The Visscher/Koskinen chemical-equilibrium H2 mixing ratio, the same
      ! fit set_IC uses when no handoff states the base partition.
      use lower_column, only: q_h2_equilibrium
      ! The root of each cell's own H2 carrier row, and the sweep that fills
      ! the molecular ions the row makes its H2 from.  The seed adopts the
      ! smaller of the fit and that root, see "WHERE THE REQUESTED FRACTION
      ! COMES FROM" above.
      use diffusive_photochemistry, only: carrier_h2_chemical_root
      use ionization_equilibrium, only: ioniz_eq
      use Conversion, only: W_to_U, U_to_W
      use BC_Apply, only: Apply_BC

      implicit none
      private

      ! Is a molecular seed being built, and from where.
      logical            :: seed_on  = .false.
      character(len=512) :: seed_dir = ''
      ! 'p' or 'T': the primitive variable kept across the conversion.
      character(len=1)   :: seed_invariant = 'p'
      ! Where the requested fraction comes from: 'handoff', 'local' (each
      ! cell's own chemical equilibrium) or 'stated'.
      character(len=8)   :: seed_x2_mode = 'handoff'
      ! The stated fraction, when the mode is 'stated'.
      real*8             :: seed_x2_stated = -1.0d0
      ! Was the state now in memory built by this module? Only then do the
      ! two record lines belong in a state file.
      logical            :: seed_built = .false.

      ! WHAT THE PARTITION DID, for the record lines of the written files.
      real*8 :: rec_q = 0.0d0, rec_qmax = 0.0d0, rec_x2 = 0.0d0
      ! Where the row's own root became the smaller of the two, and over how
      ! many cells it governs the partition ('local' only).
      integer :: rec_n_root_governs = 0
      real*8  :: rec_r_crossover = 0.0d0
      ! The base layer the handoff states, and what the network's own root
      ! would have said there: the top cell of the layer, its radius, and
      ! the smallest and largest ratio handoff/root inside it.  The layer is
      ! NOT revised by the root (see the header); the ratio is a measurement
      ! reported beside it.
      integer :: rec_j_base_top = 0
      real*8  :: rec_r_base_top = 0.0d0
      real*8  :: rec_base_ratio_min = 0.0d0, rec_base_ratio_max = 0.0d0
      character(len=16) :: rec_x2_source = 'handoff'
      integer :: rec_n_capped = 0
      ! Cells whose local chemical-equilibrium fit was clipped at the ceiling.
      integer :: rec_n_qmax = 0
      real*8 :: rec_x2_eff_max = 0.0d0
      real*8 :: rec_dT_rel = 0.0d0, rec_dp_rel = 0.0d0
      real*8 :: rec_eos_close = 0.0d0
      real*8 :: rec_dnH = 0.0d0, rec_dnHe = 0.0d0, rec_drho = 0.0d0

      ! Every conservation identity of the conversion is exact in the
      ! fraction space it is carried out in (2 x 0.5 delta = delta for the
      ! nuclei and, with bsp_mass(H2) = 2, for the mass), so what is
      ! admitted here is the round-off of the census sum and nothing more.
      real*8, parameter, public :: molecular_seed_closure_tol = 1.0d-12

      public :: molecular_seed_configure, molecular_seed_on
      public :: molecular_seed_state_file
      public :: molecular_seed_option_may_differ
      public :: molecular_seed_from_atomic_state
      public :: write_molecular_seed_header

      contains

      ! ------------------------------------------------------!

      subroutine molecular_seed_configure()
      ! Read the two environment variables once, at startup, and refuse an
      ! invariant this module does not implement rather than falling back
      ! to a default the caller did not ask for.
      character(len=512) :: env
      integer :: ios

      call get_environment_variable('EXHALE_MOLECULAR_SEED', env)
      if (len_trim(env) .eq. 0) then
         seed_on  = .false.
         seed_dir = ''
         return
      endif
      seed_on  = .true.
      seed_dir = trim(adjustl(env))

      call get_environment_variable('EXHALE_MOLECULAR_SEED_X2', env)
      if (len_trim(env) .eq. 0) then
         seed_x2_mode   = 'handoff'
         seed_x2_stated = -1.0d0
      else if (trim(adjustl(env)) .eq. 'local') then
         seed_x2_mode   = 'local'
         seed_x2_stated = -1.0d0
      else
         seed_x2_mode = 'stated'
         read(env,*,iostat=ios) seed_x2_stated
         if (ios .ne. 0 .or. seed_x2_stated .lt. 0.0d0 .or.               &
             seed_x2_stated .gt. 1.0d0) then
            write(*,'(A)') ' (molecular_seed) ERROR:'//                   &
                 ' EXHALE_MOLECULAR_SEED_X2 is "'//trim(adjustl(env))//   &
                 '"; it is the fraction of the'
            write(*,'(A)') '   hydrogen NUCLEI to bind into H2 and lies'//&
                 ' between 0 and 1, or the word "local" for each cell''s'//&
                 ' own chemical equilibrium.'
            error stop 1
         endif
      endif

      call get_environment_variable('EXHALE_MOLECULAR_SEED_INVARIANT', env)
      if (len_trim(env) .eq. 0) then
         seed_invariant = 'p'
      else if (trim(adjustl(env)) .eq. 'p') then
         seed_invariant = 'p'
      else if (trim(adjustl(env)) .eq. 'T') then
         seed_invariant = 'T'
      else
         write(*,'(A)') ' (molecular_seed) ERROR:'//                      &
              ' EXHALE_MOLECULAR_SEED_INVARIANT is "'//trim(adjustl(env))//&
              '"; the two primitive variables this'
         write(*,'(A)') '   conversion can keep are p (pressure, the'//   &
              ' default) and T (temperature). Two H atoms into one H2'
         write(*,'(A)') '   remove half a particle per pair, so rho, v,'//&
              ' p and T cannot all be kept and one has to be named.'
         error stop 1
      endif

      end subroutine molecular_seed_configure

      ! ------------------------------------------------------!

      ! One line per 25 cells of the fit against the row's own root, under
      ! the switch the carrier operator already uses.  Default off.
      logical function molecular_seed_debug_on()
      character(len=8) :: env
      call get_environment_variable('EXHALE_CARRIER_DEBUG', env)
      molecular_seed_debug_on = (trim(env) .eq. '1')
      end function molecular_seed_debug_on

      ! ------------------------------------------------------------- !

      logical function molecular_seed_on()
      molecular_seed_on = seed_on
      end function molecular_seed_on

      ! ------------------------------------------------------!

      function molecular_seed_state_file(base) result(path)
      ! Where the restart pair is read from. The seed run reads the ATOMIC
      ! state of another case directory and writes its own molecular one
      ! into ./output, so the two never collide.
      character(len=*), intent(in) :: base
      character(len=640) :: path
      if (seed_on) then
         path = trim(seed_dir)//'/'//trim(base)//'_IC.txt'
      else
         path = 'output/'//trim(base)//'_IC.txt'
      endif
      end function molecular_seed_state_file

      ! ------------------------------------------------------!

      logical function molecular_seed_option_may_differ(name)
      ! WHICH OPTION TOKENS A SEED CONVERSION MAY FIND DIFFERENT between the
      ! atomic file and this run, and no others.
      !
      ! These six decide how many unknowns a state has, which is exactly why
      ! an ordinary restart may never name them (load_IC, decision 21): a
      ! state whose rows are not this run's rows is not this run's state.
      ! A seed conversion is the one operation whose whole purpose is to
      ! ADD those rows, and it does not integrate the state it built: it
      ! writes an uncertified initialization file and stops. Every other
      ! token -- the elements carried, the constants, the grid, the
      ! reservoir, the boundary, the reconstruction -- is compared and
      ! refused as it always is, so the seed cannot come from a wind of
      ! another planet, another composition or another discretization.
      character(len=*), intent(in) :: name
      molecular_seed_option_may_differ = seed_on .and.                    &
           (  trim(name) .eq. 'mol'                                       &
         .or. trim(name) .eq. 'molbase'                                   &
         .or. trim(name) .eq. 'carrier'                                   &
         .or. trim(name) .eq. 'carrier_newton'                            &
         .or. trim(name) .eq. 'oxychem'                                   &
         .or. trim(name) .eq. 'iontrans' )
      end function molecular_seed_option_may_differ

      ! ------------------------------------------------------!

      subroutine molecular_seed_from_atomic_state(W, f_sp, T)
      ! Build the molecular state implied by the loaded atomic one, in
      ! place, and report every check of section "THE PARTITION" above.
      !
      ! On entry (W, f_sp) is what load_IC restored from the atomic pair:
      ! rho rebuilt from the species, v and p from the hydro file, the four
      ! molecular columns zero because the file has no such columns. On exit
      ! the molecular columns carry the partition, T is the temperature of
      ! the chosen invariant, and W is closed under the production equation
      ! of state.
      real*8, dimension(3,1-Ng:N+Ng),         intent(inout) :: W
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out)   :: T

      real*8, dimension(1-Ng:N+Ng) :: rho, v, p
      real*8, dimension(1-Ng:N+Ng) :: nhi, nhii, nhei, nheii, nheiii, nheiTR
      real*8, dimension(1-Ng:N+Ng) :: ne, n_tot, T_src, p_src
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
      real*8, dimension(1-Ng:N+Ng,n_element) :: nnuc0, nnuc1
      real*8, dimension(1-Ng:N+Ng) :: nchg0, nchg1, rcomp0, rcomp1
      real*8, dimension(3,1-Ng:N+Ng) :: u, W_rt
      real*8  :: q, qmax, x2, fH, df, x2_eff, d, q_cell
      integer :: j, jworst, n_capped, n_neg, n_qmax
      ! THE PARTITION OF EVERY CELL, formed before it is applied so that the
      ! transfer is one operation on one column and the 'local' mode can
      ! revise it without a second copy of the transfer.
      real*8, dimension(1-Ng:N+Ng) :: x2_col
      ! The fit's own request, kept beside the adopted one for the report.
      real*8, dimension(1-Ng:N+Ng) :: x2_fit_col
      ! The fit's mixture mixing ratio per cell, which is what the top of
      ! the base layer is found with.
      real*8, dimension(1-Ng:N+Ng) :: q_fit_col
      real*8  :: q_base, x2_base, ratio
      integer :: j_base_top
      ! The atomic composition the transfer always starts from, so that
      ! applying the partition twice gives the same state as applying the
      ! second one once.
      real*8, dimension(:,:), allocatable :: f_atomic
      ! The root of each cell's own H2 row, the sweep's outputs that are
      ! discarded with it, and what the comparison of the two found.
      real*8, dimension(1-Ng:N+Ng) :: n_root, heat_s, cool_s, eta_s
      ! |P - L|/(P + L) at the density carrier_h2_chemical_root hands back,
      ! so the report states that the number IS the root of the row.
      real*8, dimension(1-Ng:N+Ng) :: root_imbalance
      real*8  :: x2_root
      integer :: n_root_governs, j_cross

      if (.not. seed_on) return

      if (.not. thereis_mol) then
         write(*,'(A)') ' (molecular_seed) ERROR: EXHALE_MOLECULAR_SEED'// &
              ' is set, but this run has "Molecular chemistry: False",'
         write(*,'(A)') '   so there are no molecular unknowns to write.'//&
              ' The seed is built by the TARGET molecular input.'
         error stop 1
      endif
      if (.not. do_load_IC) then
         write(*,'(A)') ' (molecular_seed) ERROR: EXHALE_MOLECULAR_SEED'// &
              ' is set, but this run has "Load IC? False", so no atomic'
         write(*,'(A)') '   state was read to convert. State'//           &
              ' "Load IC? True" in the seed run''s input.'
         error stop 1
      endif
      if (seed_x2_mode .eq. 'handoff' .and.                               &
          .not. base_h2_composition_imposed()) then
         write(*,'(A)') ' (molecular_seed) ERROR: the seed partition is'// &
              ' the hydrogen-nucleus fraction the lower-atmosphere'
         write(*,'(A)') '   handoff states, and this run has no'//        &
              ' q_H2_base (base.inp key, or a lower-atmosphere profile).'
         write(*,'(A)') '   Without it the partition would be a local'//  &
              ' chemical-equilibrium estimate of the seed itself, which'
         write(*,'(A)') '   is not upstream information and is not'//     &
              ' written into a state file as though it were.'
         error stop 1
      endif

      rho = W(1,:)
      v   = W(2,:)
      p   = W(3,:)

      ! ---- the partition, and its two admissibility statements ----
      qmax = h2_mixing_ratio_ceiling()
      if (seed_x2_mode .eq. 'stated') then
         ! A stated fraction is the request itself; the mixture mixing
         ! ratio it corresponds to is reported beside it by inverting
         ! x2 = 2 q (1+He/H)/(1+q).
         x2 = seed_x2_stated
         q  = x2/(2.0d0*(1.0d0 + HeH) - x2)
         rec_x2_source = 'stated'
      else if (seed_x2_mode .eq. 'local') then
         ! Cell by cell below; the reported q and x2 are the column's
         ! largest, so that the file states the strongest request made.
         q  = 0.0d0
         x2 = 0.0d0
         rec_x2_source = 'local'
      else
         q  = h2_mixing_ratio_base()
         x2 = base_h2_nuclei_fraction()
         rec_x2_source = 'handoff'
      endif
      if (q .gt. qmax) then
         write(*,'(A,ES23.16,A,ES23.16)') ' (molecular_seed) ERROR:'//    &
              ' q_H2_base ', q, ' is above the element-ratio ceiling ', qmax
         write(*,'(A)') '   0.5/(0.5 + He/H): it asks for more hydrogen'//&
              ' than the element ratio contains.'
         error stop 1
      endif
      if (x2 .gt. 1.0d0) then
         write(*,'(A,ES23.16)') ' (molecular_seed) ERROR: the hydrogen'// &
              '-nucleus fraction bound into H2 is ', x2
         write(*,'(A)') '   and a fraction cannot exceed one.'
         error stop 1
      endif

      ! ---- the elemental census of the state as it was loaded ----
      call element_nuclei_and_charge(rho, f_sp, nnuc0, nchg0, rcomp0)

      ! ---- the temperature of the source state, by this code's EOS ----
      ! Formed before the partition so that the "keep T" invariant keeps the
      ! temperature of the ATOMIC state, and so that both invariants are
      ! stated against one number rather than against the file's column.
      nhei = 0.0d0;  nheii = 0.0d0;  nheiii = 0.0d0;  nheiTR = 0.0d0
      call get_species_densities(rho, f_sp, nhi, nhii, nhei, nheii,       &
                                 nheiii, nheiTR, nm, ne, n_tot)
      call comp_T_from_p(p, n_tot, ne, T_src)
      p_src = p

      ! ---- 2 H -> H2, cell by cell, on the cell's own budget ----
      ! Carried out on the mass fractions f = n/(rho n0): the nuclei
      ! (2 x 0.5 delta = delta) and the mass (bsp_mass(H2) = 2) are then
      ! conserved to the last bit rather than to the round-off of a
      ! multiply and divide by rho n0.
      allocate(f_atomic(1-Ng:N+Ng,n_species))
      f_atomic = f_sp
      n_qmax   = 0
      n_root_governs = 0
      j_cross        = 0
      rec_n_root_governs = 0
      rec_r_crossover    = 0.0d0

      ! ---- the partition asked for, cell by cell ----
      x2_col     = x2
      j_base_top = 0
      if (seed_x2_mode .eq. 'local') then
         do j = 1-Ng, N+Ng
            ! The fit is calibrated for a solar-like mixture and its
            ! asymptote lies above the element-ratio ceiling for any
            ! He/H >= 0.167, so it is clipped AT the ceiling, which is the
            ! physical statement "every H nucleus is in H2" (set_IC).
            q_cell = q_h2_equilibrium(p(j)*p0/1.0d6, T_src(j)*T0)
            if (q_cell .gt. qmax) then
               q_cell = qmax
               n_qmax = n_qmax + 1
            endif
            q_fit_col(j)  = q_cell
            x2_col(j)     = 2.0d0*q_cell*(1.0d0 + HeH)/(1.0d0 + q_cell)
         enddo
         ! The fit's own request, kept before the base layer overwrites it,
         ! so the report can show what the fit said where the handoff is
         ! what is carried.
         x2_fit_col = x2_col
         ! ---- WHERE THE HANDOFF STILL DESCRIBES THE GAS, MEASURED ----
         ! The base layer runs from the base to the cell BELOW the first at
         ! which the thermochemical fit falls under the handoff's own q_H2,
         ! which is where thermal dissociation starts to take the molecule
         ! apart.  It is located and REPORTED and the partition inside it is
         ! NOT the handoff's, for a reason that is about the solution and not
         ! about the seed: the handoff partition is imposed during a run on
         ! the GHOST cells alone (ionization_equilibrium, the
         ! base_h2_composition_imposed block, j <= 0), so cells 1 to the top
         ! of the layer are free, and the fixed point of the solve puts
         ! their x2 at 0.33 on LHS 1140 b against the handoff's 0.99998 --
         ! not because the network's chemistry wants it there (its own root
         ! is 0.98 at the base, item L7f) but because the layer is eddy-mixed
         ! with the dissociation region above it, and the base face carries
         ! no diffusive flux to feed it from below.  A seed that imposes the
         ! handoff there is therefore not seeding the fixed point; it is
         ! adding a transient the relaxation must then undo, and MEASURED
         ! (item L7e) that transient is a 50 percent change of the particle
         ! count of those cells, taken at one percent a pass by the movement
         ! bound, which holds the whole column while it runs.
         !
         ! WHAT IS LEFT IS A COMPARISON, and item L7f has settled it.  The
         ! ratio this block reports was 7 to 21 while the "root" was
         ! production over loss RATE at the trial; with the root SOLVED for
         ! (carrier_h2_chemical_root, which the row's quadratic dependence on
         ! n(H2) requires) it is 1.02 at the first cell and 5.87 at the top
         ! of the layer.  The lower-atmosphere model and the wind network
         ! AGREE about the bottom of the layer, and the composition step the
         ! base boundary condition makes there is a transport question and
         ! not a chemistry one (docs/lhs1140b_stationary_L7f_20260915.md).
         if (base_h2_composition_imposed()) then
            q_base  = h2_mixing_ratio_base()
            x2_base = base_h2_nuclei_fraction()
            ! A column whose fit never falls under the handoff never gets
            ! hot enough to take the molecule apart: the layer is then the
            ! whole column, which is why the search starts at N and not 0.
            j_base_top = N
            do j = 1, N
               if (q_fit_col(j) .lt. q_base) then
                  j_base_top = j - 1
                  exit
               endif
            enddo
            if (j_base_top .lt. 0) j_base_top = 0
         endif
      endif
      if (seed_x2_mode .ne. 'local') x2_fit_col = x2_col
      call transfer_h2(x2_col)

      ! ---- and, in the 'local' mode, the root of each cell's own row ----
      ! The state the fit just built carries the molecular ions only after a
      ! sweep, so one is taken here and its composition is then thrown away:
      ! f_sp goes back to the atomic state and the revised partition is
      ! applied to that, so nothing this sweep decided about the ionization
      ! reaches the file.
      if (seed_x2_mode .eq. 'local') then
         heat_s = 0.0d0;  cool_s = 0.0d0;  eta_s = 0.0d0
         call ioniz_eq(T_src, rho, f_sp, heat_s, cool_s, eta_s)
         call carrier_h2_chemical_root(rho, f_sp, n_root,               &
                                       row_imbalance = root_imbalance)
         f_sp = f_atomic
         rec_base_ratio_min = 0.0d0
         rec_base_ratio_max = 0.0d0
         do j = 1-Ng, N+Ng
            if (nnuc0(j,ie_H) .le. 0.0d0) cycle
            if (n_root(j) .le. 0.0d0)     cycle
            x2_root = min(2.0d0*n_root(j)/nnuc0(j,ie_H), 1.0d0)
            ! INSIDE THE BASE LAYER the handoff is compared with the root
            ! and the comparison is reported; the partition adopted there is
            ! the same min as everywhere else (see the block above).
            if (j .ge. 1 .and. j .le. j_base_top .and. x2_root .gt. 0.0d0) &
               then
               ratio = x2_base/x2_root
               if (rec_base_ratio_min .le. 0.0d0 .or.                    &
                   ratio .lt. rec_base_ratio_min)                        &
                  rec_base_ratio_min = ratio
               if (ratio .gt. rec_base_ratio_max)                        &
                  rec_base_ratio_max = ratio
            endif
            if (x2_root .lt. x2_col(j)) then
               x2_col(j) = x2_root
               if (j .ge. 1 .and. j .le. N) then
                  n_root_governs = n_root_governs + 1
                  if (j_cross .eq. 0) j_cross = j
               endif
            endif
         enddo
         ! The two statements side by side at a few radii, so that a
         ! reader can see WHERE they cross and by how much they differ
         ! (EXHALE_CARRIER_DEBUG=1, default off).
         if (molecular_seed_debug_on()) then
            write(*,'(A)') '   (molecular_seed) fit against the row''s'// &
                 ' own root:'
            write(*,'(A)') '      cell    r[R_p]     T[K]      x2 fit'//  &
                 '      x2 root    adopted   |P-L|/(P+L)'
            do j = 1, N
               if (mod(j, 25) .ne. 0 .and. j .ne. 1) cycle
               x2_root = 0.0d0
               if (nnuc0(j,ie_H) .gt. 0.0d0 .and. n_root(j) .gt. 0.0d0)   &
                  x2_root = min(2.0d0*n_root(j)/nnuc0(j,ie_H), 1.0d0)
               write(*,'(A,I6,F10.4,F10.1,4ES12.4)') '     ', j, r(j),    &
                    T_src(j)*T0, x2_fit_col(j), x2_root, x2_col(j),       &
                    root_imbalance(j)
            enddo
         endif
         call transfer_h2(x2_col)
         rec_n_root_governs = n_root_governs
         if (j_cross .gt. 0) rec_r_crossover = r(j_cross)
         rec_j_base_top = j_base_top
         if (j_base_top .ge. 1) rec_r_base_top = r(j_base_top)
      endif
      deallocate(f_atomic)

      ! ---- what the transfer did to the invariants of the state ----
      call element_nuclei_and_charge(rho, f_sp, nnuc1, nchg1, rcomp1)
      rec_dnH  = max_rel_change(nnuc0(:,ie_H),  nnuc1(:,ie_H))
      rec_dnHe = max_rel_change(nnuc0(:,ie_He), nnuc1(:,ie_He))
      rec_drho = max_rel_change(rcomp0, rcomp1)
      n_neg = 0
      do j = 1-Ng, N+Ng
         if (f_sp(j,isp_HI) .lt. 0.0d0 .or. f_sp(j,isp_H2) .lt. 0.0d0)    &
            n_neg = n_neg + 1
      enddo
      if (n_neg .gt. 0) then
         write(*,'(A,I0,A)') ' (molecular_seed) ERROR: ', n_neg,          &
              ' cell(s) carry a negative H I or H2 fraction after the'
         write(*,'(A)') '   transfer, so the cap did not hold.'
         error stop 1
      endif
      if (rec_dnH .gt. molecular_seed_closure_tol .or.                    &
          rec_dnHe .gt. molecular_seed_closure_tol .or.                   &
          rec_drho .gt. molecular_seed_closure_tol) then
         write(*,'(A)') ' (molecular_seed) ERROR: the transfer did not'// &
              ' conserve what 2 H -> H2 conserves.'
         write(*,'(A,3ES12.4)') '   max relative change of the H nuclei,'//&
              ' the He nuclei and the mass: ', rec_dnH, rec_dnHe, rec_drho
         error stop 1
      endif

      ! ---- the thermodynamic invariant ----
      call get_species_densities(rho, f_sp, nhi, nhii, nhei, nheii,       &
                                 nheiii, nheiTR, nm, ne, n_tot)
      if (seed_invariant .eq. 'p') then
         call comp_T_from_p(p, n_tot, ne, T)
      else
         T = T_src
         call comp_p_from_T(T, n_tot, ne, p)
         W(3,:) = p
      endif
      ! The two numbers that say how far the chosen invariant moved the
      ! other variable, so that the two runs of the pilot are comparable
      ! from the files alone.
      rec_dT_rel = max_rel_change(T_src, T)
      rec_dp_rel = max_rel_change(p_src, p)

      ! ---- the boundary, and the closure of the equation of state ----
      ! Apply_BC fills the ghost rows with the TARGET run's own lower and
      ! outer boundary rather than leaving the atomic run's in the file,
      ! and the round trip back through U_to_W is the closure test: the
      ! pressure that comes back out of the caloric EOS against the one
      ! that went in.
      call W_to_U(W, u)
      call Apply_BC(u)
      call U_to_W(u, W_rt)
      rec_eos_close = 0.0d0
      jworst = 1
      do j = 1, N
         d = abs(W_rt(3,j) - W(3,j))/max(abs(W(3,j)), tiny(1.0d0))
         if (d .gt. rec_eos_close) then
            rec_eos_close = d
            jworst = j
         endif
      enddo
      if (rec_eos_close .gt. molecular_seed_closure_tol) then
         write(*,'(A,ES12.4,A,I0)') ' (molecular_seed) ERROR: the'//      &
              ' written state does not close under the production'//      &
              ' equation of state: ', rec_eos_close, ' at cell ', jworst
         error stop 1
      endif
      ! The ghosts the boundary just wrote are part of the state file, so
      ! the primitive and the composition have to be the ones they imply.
      W = W_rt
      rho = W(1,:)
      v   = W(2,:)
      p   = W(3,:)
      call get_species_densities(rho, f_sp, nhi, nhii, nhei, nheii,       &
                                 nheiii, nheiTR, nm, ne, n_tot)
      call comp_T_from_p(p, n_tot, ne, T)

      rec_q        = q
      rec_qmax     = qmax
      ! rec_x2 and q are what transfer_h2 recorded of the partition it
      ! applied, in every mode, so nothing is restated here.
      rec_n_qmax   = n_qmax
      rec_n_capped = n_capped
      seed_built   = .true.

      call report_molecular_seed()

      contains

      ! 2 H -> H2 ON EACH CELL'S OWN BUDGET, from the ATOMIC composition.
      ! Carried out on the mass fractions f = n/(rho n0): the nuclei
      ! (2 x 0.5 delta = delta) and the mass (bsp_mass(H2) = 2) are then
      ! conserved to the last bit rather than to the round-off of a multiply
      ! and divide by rho n0.  Starting from f_atomic every time is what
      ! lets the 'local' mode revise its partition and apply it once more
      ! without the two applications compounding.
      subroutine transfer_h2(x2_of_cell)
      real*8, dimension(1-Ng:N+Ng), intent(in) :: x2_of_cell
      integer :: jc
      real*8  :: x2c, fHc, dfc, x2e, f_pool
      f_sp = f_atomic
      n_capped       = 0
      rec_x2         = 0.0d0
      rec_x2_eff_max = 0.0d0
      do jc = 1-Ng, N+Ng
         x2c = x2_of_cell(jc)
         if (x2c .gt. rec_x2) then
            rec_x2 = x2c
            ! The mixture mixing ratio the adopted fraction corresponds to,
            ! by inverting x2 = 2 q (1+He/H)/(1+q).  The handoff and the
            ! stated modes already carry their own q, which is the number
            ! the run was configured with, so only the cell-by-cell mode
            ! forms one here.
            if (seed_x2_mode .eq. 'local')                               &
               q = x2c/(2.0d0*(1.0d0 + HeH) - x2c)
         endif
         ! THE NEUTRAL HYDROGEN THE PARTITION ACTS ON is H I plus the
         ! nuclei already bound in H2. In a cell of the atomic state the
         ! second term is zero and this is H I alone; the lower ghosts are
         ! not atomic, because load_IC does not read the ghost rows of the
         ! pair and fills them with this run's molecular reservoir, which
         ! already carries H2 (f_H2 = 5.3e-2 against f_HI = 1.8e-6 on the
         ! LHS 1140 b ghost at He/H = 2.13). Partitioning H I alone there
         ! and overwriting f_H2 would destroy those nuclei. The molecular
         ! ions are left as they are for the same reason: zero in every
         ! atomic cell, and part of the reservoir's budget in a ghost.
         fHc    = nnuc0(jc,ie_H)/(rho(jc)*n0)
         f_pool = f_sp(jc,isp_HI) + 2.0d0*f_sp(jc,isp_H2)
         dfc    = x2c*fHc
         if (dfc .gt. f_pool) then
            dfc = f_pool
            n_capped = n_capped + 1
         endif
         f_sp(jc,isp_H2)   = 0.5d0*dfc
         f_sp(jc,isp_HI)   = f_pool - dfc
         if (fHc .gt. 0.0d0) then
            x2e = dfc/fHc
            if (x2e .gt. rec_x2_eff_max) rec_x2_eff_max = x2e
         endif
      enddo
      end subroutine transfer_h2

      end subroutine molecular_seed_from_atomic_state

      ! ------------------------------------------------------!

      real*8 function max_rel_change(a, b)
      ! Largest relative change between two columns over the whole grid,
      ! ghosts included; a zero entry contributes only where the other is
      ! zero too, and then nothing.
      real*8, dimension(1-Ng:N+Ng), intent(in) :: a, b
      integer :: j
      real*8  :: den
      max_rel_change = 0.0d0
      do j = 1-Ng, N+Ng
         den = max(abs(a(j)), abs(b(j)))
         if (den .le. 0.0d0) cycle
         max_rel_change = max(max_rel_change, abs(a(j) - b(j))/den)
      enddo
      end function max_rel_change

      ! ------------------------------------------------------!

      subroutine report_molecular_seed()
      ! Everything the conversion did, on the run log, in the order the
      ! checks were made.
      write(*,'(A)') ' (molecular_seed) the molecular state implied by'//  &
           ' the atomic state in'
      write(*,'(A)') '   '//trim(seed_dir)
      write(*,'(A,ES23.16,A,ES23.16)') '   q_H2_base  ', rec_q,           &
           '   ceiling 0.5/(0.5+He/H) ', rec_qmax
      write(*,'(A,A)')                 '   fraction requested by: ',      &
           trim(rec_x2_source)
      write(*,'(A,ES23.16,A,ES23.16)') '   x2 = 2q(1+He/H)/(1+q) ',       &
           rec_x2, '   largest x2 actually transferred ', rec_x2_eff_max
      write(*,'(A,I0,A,I0,A)') '   cells capped by their own neutral'//   &
           ' hydrogen: ', rec_n_capped, ' of ', N + 2*Ng, ''
      if (rec_x2_source .eq. 'local') then
         write(*,'(A,I0,A)') '   cells whose equilibrium fit was clipped'//&
              ' at the element-ratio ceiling: ', rec_n_qmax, ''
         ! WHERE THE TWO STATEMENTS CROSS.  Inside the crossover the
         ! thermochemical fit is the smaller and is what the seed carries;
         ! outside it the row's own root is, and the difference there is
         ! what the fit alone would have asserted about a gas the star is
         ! dissociating.
         if (rec_j_base_top .ge. 1) then
            write(*,'(A,I0,A,F8.4,A,ES23.16)') '   layer the handoff'//  &
                 ' still describes: cells 1 to ', rec_j_base_top,        &
                 ', up to r = ', rec_r_base_top, ' R_p, its x2 = ',       &
                 base_h2_nuclei_fraction()
            write(*,'(A)') '     (the thermochemical fit falls under'//   &
                 ' q_H2_base above it, which is where thermal'//          &
                 ' dissociation takes over; the partition carried there'// &
                 ' is the same min as everywhere else, because the run'//  &
                 ' imposes the handoff on the GHOST cells alone and the'// &
                 ' fixed point of the solve puts that layer at the'//      &
                 ' network''s own root)'
            ! WHAT THE TWO CHEMISTRIES SAY ABOUT THIS LAYER, reported and
            ! not acted on: the handoff against the root of the wind
            ! network's own H2 row, cell by cell, as a ratio.  Item L7f
            ! measured it at 1.02 at the first cell of LHS 1140 b and 5.87
            ! at the top of the layer, i.e. the two agree at the bottom and
            ! part company as the temperature rises through it; a ratio far
            ! from 1 at the FIRST cell would be the disagreement that item
            ! looked for and did not find.
            if (rec_base_ratio_max .gt. 0.0d0)                            &
               write(*,'(A,F7.2,A,F7.2,A)') '   the handoff stands ',     &
                 rec_base_ratio_min, ' to ', rec_base_ratio_max,          &
                 ' times the root of the wind network''s own H2 row'//    &
                 ' through that layer (item L7f)'
         else
            write(*,'(A)') '   no handoff states a base layer; the'//     &
                 ' partition is the fit against the root everywhere'
         endif
         if (rec_n_root_governs .gt. 0) then
            write(*,'(A,I0,A,F8.4,A)') '   the row''s own root is the'//  &
                 ' smaller in ', rec_n_root_governs, ' cell(s), from'//   &
                 ' r = ', rec_r_crossover, ' R_p outward'
         else
            write(*,'(A)') '   the thermochemical fit is the smaller in'//&
                 ' every cell above the base layer; the root revises'//   &
                 ' nothing'
         endif
      endif
      write(*,'(A,3ES12.4)') '   max relative change of H nuclei, He'//   &
           ' nuclei, mass: ', rec_dnH, rec_dnHe, rec_drho
      write(*,'(A,A1)') '   thermodynamic invariant kept: ', seed_invariant
      write(*,'(A,ES12.4,A,ES12.4)') '   it moved T by (max, relative) ', &
           rec_dT_rel, '   and p by ', rec_dp_rel
      write(*,'(A,ES12.4)') '   equation-of-state closure of the'//       &
           ' written state (p round trip): ', rec_eos_close
      end subroutine report_molecular_seed

      ! ------------------------------------------------------!

      subroutine write_molecular_seed_header(unit)
      ! The two lines that say a state file is a molecular seed and what
      ! was done to build it. Written into BOTH state files, so the pair
      ! cannot describe different constructions. '#' comments, read by no
      ! parser and compared by nothing.
      integer, intent(in) :: unit
      if (.not. seed_built) return
      write(unit,'(A)') '# molecular-seed-from: '//trim(seed_dir)
      write(unit,'(A,ES23.16,A,ES23.16,A,ES23.16,A,I0,A,A1,A,A)')         &
           '# molecular_partition: q_H2 ', rec_q,                         &
           ' q_max ', rec_qmax,                                           &
           ' x2 ', rec_x2,                                                &
           ' capped_cells ', rec_n_capped,                                &
           ' invariant ', seed_invariant,                                 &
           ' x2_source ', trim(rec_x2_source)
      if (rec_x2_source .eq. 'local')                                     &
         write(unit,'(A,I0,A,F10.5,A,I0,A,F10.5,A,F8.2,A,F8.2)')          &
           '# molecular_partition_local: cells where the row''s own'//    &
           ' root is the smaller ', rec_n_root_governs,                   &
           '  crossover r[R_p] ', rec_r_crossover,                        &
           '  layer the handoff describes cells 1..', rec_j_base_top,     &
           ' up to r[R_p] ', rec_r_base_top,                              &
           '  handoff/root there ', rec_base_ratio_min, ' to ',           &
           rec_base_ratio_max
      write(unit,'(A)') '# molecular_partition_rule: H2 from'//           &
           ' min(x2 n_H,nuclei, n_HI) of each cell, its own element'//    &
           ' census; H2+ H3+ HeH+ start at zero'//                        &
           '; x2_source local means min(thermochemical fit, root of'//  &
           ' the cell''s own H2 row) in every cell'
      end subroutine write_molecular_seed_header

      end module molecular_seed
