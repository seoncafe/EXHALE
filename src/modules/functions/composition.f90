      module composition
      ! Single point that turns (rho, f_sp) into all species number
      ! densities plus the free-electron (ne) and total-particle (n_tot)
      ! densities, and converts between pressure and temperature.
      !
      ! This collapses the three near-identical extraction blocks that were
      ! repeated in EXHALE_main (and the inline copies in energy_semi_implicit
      ! / post_process), so that the electron-density and total-density
      ! POLICY lives in exactly one place.
      !
      ! The eos_include_metals policy (global_parameters, default .true.;
      ! runtime key 'eos_metals 0|1' in metals.inp) is applied here by
      ! passing nm through to calc_ne / calc_ntot, which add the metal
      ! electrons and metal nuclei. With eos_metals 0 (or metals off) the
      ! legacy H/He-only behavior is reproduced exactly.

      use global_parameters
      use species_table, only: n_mion, mion_fsp,                       &
                               isp_H2, isp_H2p, isp_H3p, isp_HeHp,      &
                               isp_HI, isp_HII, isp_HeI, isp_HeII,      &
                               isp_HeIII, isp_HeTR, bsp_mass, melem_A,  &
                               isp_OH, isp_H2O, isp_CO,                 &
                               n_bsp, bsp_fsp, bsp_nH, bsp_nHe,         &
                               bsp_is_excited_level
      use utils, only: calc_ne, calc_ntot
      use caloric_eos, only: caloric_state_from_composition
      use lower_column, only: q_h2_equilibrium

      implicit none
      private
      public :: get_species_densities, comp_T_from_p, comp_p_from_T
      public :: comp_mass_per_H, comp_ntot_bc, comp_rho_bc
      public :: mass_per_H_nucleus_without_He
      public :: element_ratio_HeH
      public :: h2_mixing_ratio_base, h2_bound_fraction
      public :: h2_mixing_ratio_ceiling
      public :: base_h2_nuclei_fraction, base_h2_composition_imposed
      public :: he_ground_singlet_density, n_cells_he_singlet_clamped

      ! Cells at which the ground singlet came out negative and was set to
      ! its physical floor, zero, summed over the whole run (marching,
      ! equilibrium diagnostics and the advection post-process alike).
      ! Reported at the end of a run; zero means the two helium columns never
      ! crossed and the run is the one an unguarded build would have produced.
      integer :: n_cells_he_singlet_clamped = 0

      ! n(1^1S) from the summed neutral helium and the metastable. The state
      ! vector carries n(He I) with the 2^3S population inside it
      ! (bsp_is_excited_level), so the ground singlet is a difference, and
      ! this is the only place the difference is taken.
      interface he_ground_singlet_density
         module procedure he_ground_singlet_density_cell
         module procedure he_ground_singlet_density_column
      end interface he_ground_singlet_density

      contains

      ! ------------------------------------------------------!

      subroutine get_species_densities(rho, f_sp, nhi, nhii, nhei,     &
                                       nheii, nheiii, nheiTR, nm,       &
                                       ne, n_tot)
      ! (rho, f_sp) -> all number densities + ne + n_tot, reproducing the
      ! legacy EXHALE_main extraction blocks exactly. The He arrays are
      ! intent(inout): like the original module-scope locals they keep
      ! their previous values when helium is absent (they are then never
      ! read by calc_ne/calc_ntot, which guard on thereis_He).

      real*8, dimension(1-Ng:N+Ng),          intent(in)    :: rho
      real*8, dimension(1-Ng:N+Ng,n_species),intent(in)    :: f_sp
      real*8, dimension(1-Ng:N+Ng),          intent(out)   :: nhi, nhii
      real*8, dimension(1-Ng:N+Ng),          intent(inout) :: nhei, nheii
      real*8, dimension(1-Ng:N+Ng),          intent(inout) :: nheiii, nheiTR
      real*8, dimension(1-Ng:N+Ng,n_mion),   intent(out)   :: nm
      real*8, dimension(1-Ng:N+Ng),          intent(out)   :: ne, n_tot
      real*8, dimension(1-Ng:N+Ng,4) :: nmol_l   ! molecules
      real*8, dimension(1-Ng:N+Ng,3) :: nox_l    ! OH H2O CO
      integer :: im

      nhi  = rho*f_sp(:,isp_HI)
      nhii = rho*f_sp(:,isp_HII)
      if (thereis_He) then
         nhei   = rho*f_sp(:,isp_HeI)
         nheii  = rho*f_sp(:,isp_HeII)
         nheiii = rho*f_sp(:,isp_HeIII)
         if (thereis_HeITR) nheiTR = rho*f_sp(:,isp_HeTR)
      endif
      do im = 1, n_mion
         nm(:,im) = rho*f_sp(:,mion_fsp(im))
      enddo

      ! molecular species: include their electrons and their (one
      ! particle each) contribution to the EOS particle count.  nmol_l is
      ! zero when thereis_mol is off, so the legacy path is unchanged.
      if (thereis_mol) then
         nmol_l(:,1) = rho*f_sp(:,isp_H2)
         nmol_l(:,2) = rho*f_sp(:,isp_H2p)
         nmol_l(:,3) = rho*f_sp(:,isp_H3p)
         nmol_l(:,4) = rho*f_sp(:,isp_HeHp)
         call calc_ne(nhii, nheii, nheiii, ne, nm, nmol_l)
         ! The oxygen-chemistry carriers are gas particles too, and the
         ! oxygen and carbon nuclei they hold have been taken OUT of the
         ! metal ion columns by the ionization solve -- so leaving them out
         ! here does not merely lose 5e-4 of the particle count, it makes
         ! this routine and ioniz_eq (which does pass them) disagree about
         ! what a particle is for the same state, and the temperature
         ! T = p/((n_tot + n_e) k) then depends on which of the two last
         ! wrote it. They are neutral, so calc_ne is unaffected.
         if (thereis_oxychem) then
            nox_l(:,1) = rho*f_sp(:,isp_OH)
            nox_l(:,2) = rho*f_sp(:,isp_H2O)
            nox_l(:,3) = rho*f_sp(:,isp_CO)
            call calc_ntot(nhi, nhii, nhei, nheii, nheiii, n_tot,       &
                           nm, nmol_l, nox_l)
         else
            call calc_ntot(nhi, nhii, nhei, nheii, nheiii, n_tot,       &
                           nm, nmol_l)
         endif
      else
         call calc_ne(nhii, nheii, nheiii, ne, nm)
         call calc_ntot(nhi, nhii, nhei, nheii, nheiii, n_tot, nm)
      endif

      ! Particle count of the first interior cell, kept where n_tot and n_e are
      ! defined so the lower boundary cannot disagree with the EOS about what a
      ! particle is. Read by the continuous-temperature base ghost (Apply_BC),
      ! which needs T(1) = p(1)/n_part_cell1; unused otherwise.
      n_part_cell1 = n_tot(1) + ne(1)

      ! Caloric half of the equation of state.  It reads the same
      ! (rho, f_sp, n_e, n_tot) as the thermal half above, from here, so the
      ! two cannot describe different gas -- and, being refreshed only where
      ! the composition is refreshed, the map from energy to pressure never
      ! depends on a previous state.
      call caloric_state_from_composition(rho, f_sp, ne, n_tot)

      end subroutine get_species_densities

      ! ------------------------------------------------------!

      real*8 function he_ground_singlet_density_cell(nhei, nheiTR)
      ! Number density of neutral helium in its ground singlet, n(1^1S), in
      ! one cell: the summed neutral helium minus the metastable level it
      ! contains, n(He I) - n(2^3S), floored at zero.
      !
      ! Why the floor. Zero is the physical lower bound of a population, and
      ! the difference can fall below it: where the metastable holds most of
      ! the neutral helium the two columns share their leading digits, so what
      ! is left is round-off, and a restart whose two columns were written by
      ! different solves can put the difference on the wrong side of zero
      ! outright. Every consumer here is LINEAR in the singlet -- the He I
      ! opacity, the He I photoelectric heating, the He I secondary-ionization
      ! source, the absorbed energy in the denominator of q -- so a negative
      ! value is not an instability but a negative opacity and a negative
      ! heating rate, both unphysical. Clamping is the enforcement of the
      ! bound, not an approximation.
      !
      ! The test is "less than zero", not the negation of "greater than or
      ! equal to zero": a NaN handed in must pass through as a NaN so the
      ! run's NaN detector reports it, rather than be silently turned into a
      ! zero singlet.
      !
      ! The counter is incremented atomically, so the function may be called
      ! from inside a parallel region (the self-consistent field loop of
      ! ionization_equilibrium does) as well as from the serial column
      ! formations, and the run-wide total stays exact without a lock around
      ! the call.
      real*8, intent(in) :: nhei, nheiTR

      he_ground_singlet_density_cell = nhei - nheiTR
      if (he_ground_singlet_density_cell .lt. 0.0d0) then
         he_ground_singlet_density_cell = 0.0d0
         !$omp atomic update
         n_cells_he_singlet_clamped = n_cells_he_singlet_clamped + 1
      endif

      end function he_ground_singlet_density_cell

      ! ------------------------------------------------------!

      function he_ground_singlet_density_column(nhei, nheiTR) result(nheiS)
      ! n(1^1S) over the whole grid, ghost cells included. Element by element
      ! this is the cell function above, so it is bitwise the array
      ! expression nhei - nheiTR wherever that expression is already
      ! non-negative.
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei, nheiTR
      real*8, dimension(1-Ng:N+Ng) :: nheiS
      integer :: j

      do j = 1-Ng, N+Ng
         nheiS(j) = he_ground_singlet_density_cell(nhei(j), nheiTR(j))
      enddo

      end function he_ground_singlet_density_column

      ! ------------------------------------------------------!

      subroutine comp_T_from_p(p, n_tot, ne, T)
      ! Adimensional ideal-gas inverse: T = p / (n_tot + ne).
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: p, n_tot, ne
      real*8, dimension(1-Ng:N+Ng), intent(out) :: T
      T = p/(n_tot + ne)
      end subroutine comp_T_from_p

      ! ------------------------------------------------------!

      subroutine comp_p_from_T(T, n_tot, ne, p)
      ! Adimensional ideal-gas law: p = (n_tot + ne) * T.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: T, n_tot, ne
      real*8, dimension(1-Ng:N+Ng), intent(out) :: p
      p = (n_tot + ne)*T
      end subroutine comp_p_from_T

      ! ------------------------------------------------------!
      ! Base composition scalars. These functions are the SINGLE
      ! SOURCE of the base composition policy: mass_per_H, ntot_bc, rho_bc
      ! (set in input_read) all flow from here, so the policy cannot
      ! disagree between code paths (the §3.4 root cause).  The molecular
      ! base belongs here for the same reason: it is part of the base
      ! particle count, not a separate correction applied afterwards.
      !
      ! The H and He masses are the species table's and nothing else:
      ! bsp_mass(isp_HI) = 1 and bsp_mass(isp_HeI) = 3.9715259, the
      ! helium-4 atom in units of the hydrogen atom, so this base
      ! composition and the rho that calc_rho forms from the same table
      ! cannot disagree.  The metal weights melem_A the sum below adds are
      ! in the same unit, the hydrogen atom, not in u.
      ! ------------------------------------------------------!

      real*8 function comp_mass_per_H()
      ! Gas mass per H nucleus [m_H]. H + He (bsp metadata) plus the trace
      ! metals when they are in the EOS budget (eos_include_metals and
      ! metals present); H/He-only otherwise.
      ! Indexing note: bsp_mass is indexed by bsp POSITION, and for the six
      ! base atomic species bsp_fsp(1:6) = 1..6, so isp_HI/isp_HeI coincide
      ! with the bsp positions here. Do NOT extend this shortcut to the
      ! molecular species (isp_H2 = 34 but bsp position 7).
      comp_mass_per_H = bsp_mass(isp_HI) + bsp_mass(isp_HeI)*HeH
      if (eos_include_metals .and. thereis_metals) then
         comp_mass_per_H = comp_mass_per_H + sum(melem_ab*melem_A)
      endif
      end function comp_mass_per_H

      ! ------------------------------------------------------!

      real*8 function mass_per_H_nucleus_without_He()
      ! Mass carried by one hydrogen nucleus together with the trace metals
      ! slaved to it [m_H], i.e. comp_mass_per_H() minus its helium term. This
      ! is the m_1 of the two-component (H+metals vs He) split used by the
      ! binary element diffusion: with it, m_1 n_H + m_He n_He = rho exactly,
      ! under either eos_metals setting, because calc_rho drops the metal mass
      ! from rho by the same policy that drops it from here.
      ! comp_mass_per_H keeps its own literal expression rather than calling
      ! this function: adding the helium term last instead of first would
      ! reorder the floating-point sum and move the goldens.
      mass_per_H_nucleus_without_He = bsp_mass(isp_HI)
      if (eos_include_metals .and. thereis_metals) then
         mass_per_H_nucleus_without_He =                                 &
              mass_per_H_nucleus_without_He + sum(melem_ab*melem_A)
      endif
      end function mass_per_H_nucleus_without_He

      ! ------------------------------------------------------!

      real*8 function comp_ntot_bc()
      ! Total nuclei density at the base in units of n0 (n0 = H+He nuclei).
      ! Legacy H/He-only value is 1; the trace metals add their nuclei when
      ! in the EOS budget.  With a molecular base the H nuclei bound into H2
      ! no longer count as separate particles and are removed here, so the
      ! base particle count has a single definition (the metal terms first,
      ! then the H2 binding, as in the original input_read sequence).
      comp_ntot_bc = 1.0d0
      if (eos_include_metals .and. thereis_metals) then
         comp_ntot_bc = (1.0d0 + HeH + sum(melem_ab))/(1.0d0 + HeH)
      endif
      if (molecular_base) comp_ntot_bc = comp_ntot_bc - h2_bound_fraction()
      end function comp_ntot_bc

      ! ------------------------------------------------------!

      real*8 function h2_mixing_ratio_base()
      ! H2 volume mixing ratio q_H2 = n_H2/(n_H2+n_H+n_He) at the base level.
      ! Taken from the lower-atmosphere photochemistry when base.inp supplied
      ! one (q_H2_base > 0), and from the Visscher/Koskinen chemical-
      ! equilibrium fit at (p_base_bar, T0) otherwise.  Chemical equilibrium
      ! underestimates H2 dissociation at T_eq ~ 1000-2000 K, which is why
      ! the photochemical value takes precedence when it exists.
      if (q_h2_base .gt. 0.0d0) then
         h2_mixing_ratio_base = q_h2_base
      else
         h2_mixing_ratio_base = q_h2_equilibrium(p_base_bar, T0)
      endif
      end function h2_mixing_ratio_base

      ! ------------------------------------------------------!

      real*8 function h2_mixing_ratio_ceiling()
      ! Largest H2 volume mixing ratio q_H2 = n_H2/(n_H2+n_H+n_He) a mixture
      ! can carry at this helium-to-hydrogen ratio, reached when every H
      ! nucleus is bound into H2.  Two H nuclei then make 0.5 molecules per H
      ! nucleus against HeH helium atoms, so
      !
      !     q_H2,max = 0.5/(0.5 + He/H),
      !
      ! which is 0.863 at the solar-like He/H = 0.0793, 1/3 at He/H = 1 and
      ! 0.048 at He/H = 10.  A requested q_H2 above this is not a large value
      ! but an impossible one: it asks for more hydrogen than the element
      ! ratio contains.  input_read refuses such a value at startup, which is
      ! why h2_bound_fraction below can evaluate its expression unguarded.
      h2_mixing_ratio_ceiling = 0.5d0/(0.5d0 + HeH)
      end function h2_mixing_ratio_ceiling

      ! ------------------------------------------------------!

      real*8 function h2_bound_fraction()
      ! Particles removed from the base budget, per (H+He) nucleus, by the H
      ! nuclei bound into H2: two H nuclei make one molecule, so a fraction
      ! x2 of the H nuclei costs x2/2 particles.  The fit returns the MIXTURE
      ! mixing ratio, hence x2 = 2 q (1+HeH)/(1+q) per H nucleus (see
      ! mu_mixture in lower_column.f90), and the result is expressed per
      ! (H+He) nucleus.
      !
      ! x2 <= 1 needs no clamp here: q <= h2_mixing_ratio_ceiling() is
      ! established at startup (input_read, "q_H2_base above the mixture
      ! ceiling"), and x2 = 1 is exactly that ceiling.  The clamp this
      ! function used to carry silently turned an impossible request into a
      ! fully molecular base and let the run continue, so the input error it
      ! concealed reached the wind as a base state nobody had asked for.
      h2_bound_fraction = 0.5d0*base_h2_nuclei_fraction()/(1.0d0 + HeH)
      end function h2_bound_fraction

      ! ------------------------------------------------------!

      real*8 function base_h2_nuclei_fraction()
      ! Fraction x2 of the base HYDROGEN NUCLEI that are bound into H2,
      !
      !     x2 = 2 q (1+He/H)/(1+q),   q = h2_mixing_ratio_base(),
      !
      ! the same conversion from a mixture mixing ratio to a hydrogen-nucleus
      ! fraction that mu_mixture uses in lower_column.f90.  x2 = 1 is fully
      ! molecular hydrogen and is reached exactly at q = q_H2,max, which
      ! startup has already established this q does not exceed.
      !
      ! This is the single definition of "how molecular is the base": the
      ! equation of state removes x2/2 particles per H nucleus through
      ! h2_bound_fraction above, and the species state imposes the same x2 on
      ! the inflowing ghosts (ionization_equilibrium). One number, so the two
      ! cannot describe different gas.
      real*8 :: q
      q = h2_mixing_ratio_base()
      base_h2_nuclei_fraction = 2.0d0*q*(1.0d0 + HeH)/(1.0d0 + q)
      end function base_h2_nuclei_fraction

      ! ------------------------------------------------------!

      logical function base_h2_composition_imposed()
      ! Is the base molecular partition incoming data rather than something
      ! the base cell derives for itself?
      !
      ! True only when a lower-atmosphere handoff stated it (base.inp key
      ! q_H2_base, or the profile's value at the matching level) AND the
      ! molecular network is actually solved, so there are H2 species to
      ! impose it on. Without a handoff the partition comes from the
      ! chemical-equilibrium fit, which is a local estimate and not upstream
      ! information; imposing that on the solver would only hand the solver
      ! back its own answer, so the fit branch leaves the species free and
      ! reaches consistency the other way, through the particle count.
      base_h2_composition_imposed = (q_h2_base .gt. 0.0d0) .and. thereis_mol
      end function base_h2_composition_imposed

      ! ------------------------------------------------------!

      function element_ratio_HeH(f_sp) result(heh_cell)
      ! Helium-to-hydrogen ELEMENT ratio n_He/n_H per cell: nuclei counted
      ! over every species that carries them, with the bsp_nH / bsp_nHe
      ! weights of species_table (H2 and H2+ carry two H nuclei, H3+ three,
      ! HeH+ one of each).  The molecular species are therefore inside the
      ! count, so the ratio is the same physical quantity in the atomic and
      ! molecular regions; the columns of the species that a run does not
      ! carry are zero, so nothing has to be switched on the flags.  The
      ! metastable He 2^3S column is skipped, not because its nucleus does
      ! not count but because it is an excited level of He I and its nucleus
      ! is already counted there (bsp_is_excited_level).  This is the same
      ! element count load_IC applies to a restart file.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng) :: heh_cell
      real*8, dimension(1-Ng:N+Ng) :: nuc_H, nuc_He
      integer :: ib

      nuc_H  = 0.0d0
      nuc_He = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         if (bsp_nH(ib)  .gt. 0)                                        &
            nuc_H  = nuc_H  + dble(bsp_nH(ib)) *f_sp(:,bsp_fsp(ib))
         if (bsp_nHe(ib) .gt. 0)                                        &
            nuc_He = nuc_He + dble(bsp_nHe(ib))*f_sp(:,bsp_fsp(ib))
      enddo
      heh_cell = nuc_He/max(nuc_H, 1.0d-30)

      end function element_ratio_HeH

      ! ------------------------------------------------------!

      real*8 function comp_rho_bc()
      ! Base mass density [rho normalization]: gas mass per (H+He) nucleus.
      comp_rho_bc = comp_mass_per_H()/(1.0d0 + HeH)
      end function comp_rho_bc

      ! End of module
      end module composition
