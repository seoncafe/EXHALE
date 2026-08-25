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
                               n_bsp, bsp_fsp, bsp_nH, bsp_nHe
      use utils, only: calc_ne, calc_ntot
      use lower_column, only: q_h2_equilibrium

      implicit none
      private
      public :: get_species_densities, comp_T_from_p, comp_p_from_T
      public :: comp_mass_per_H, comp_ntot_bc, comp_rho_bc
      public :: mass_per_H_nucleus_without_He
      public :: element_ratio_HeH
      public :: h2_mixing_ratio_base, h2_bound_fraction

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
         call calc_ntot(nhi, nhii, nhei, nheii, nheiii, nheiTR, n_tot,  &
                        nm, nmol_l)
      else
      call calc_ne(nhii, nheii, nheiii, ne, nm)
      call calc_ntot(nhi, nhii, nhei, nheii, nheiii, nheiTR, n_tot, nm)
      endif

      ! Particle count of the first interior cell, kept where n_tot and n_e are
      ! defined so the lower boundary cannot disagree with the EOS about what a
      ! particle is. Read by the continuous-temperature base ghost (Apply_BC),
      ! which needs T(1) = p(1)/n_part_cell1; unused otherwise.
      n_part_cell1 = n_tot(1) + ne(1)

      end subroutine get_species_densities

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
      ! Byte-identity: bsp_mass(isp_HI) is exactly 1.0d0 and
      ! bsp_mass(isp_HeI) exactly 4.0d0, and 1.0 / 4.0 are exact in double,
      ! so bsp_mass(isp_HI) + bsp_mass(isp_HeI)*HeH reproduces the legacy
      ! literal "1.0 + 4.0*HeH" bitwise. The metal terms and the division
      ! order match the original expressions exactly.
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
      ! the photochemical value takes precedence when it exists
      ! (docs/base_composition_handoff_plan.md).
      if (q_h2_base .gt. 0.0d0) then
         h2_mixing_ratio_base = q_h2_base
      else
         h2_mixing_ratio_base = q_h2_equilibrium(p_base_bar, T0)
      endif
      end function h2_mixing_ratio_base

      ! ------------------------------------------------------!

      real*8 function h2_bound_fraction()
      ! Particles removed from the base budget, per (H+He) nucleus, by the H
      ! nuclei bound into H2: two H nuclei make one molecule, so a fraction
      ! x2 of the H nuclei costs x2/2 particles.  The fit returns the MIXTURE
      ! mixing ratio, hence x2 = 2 q (1+HeH)/(1+q) per H nucleus (see
      ! mu_mixture in lower_column.f90), capped at full molecular hydrogen,
      ! and the result is expressed per (H+He) nucleus.
      real*8 :: q, x2
      q  = h2_mixing_ratio_base()
      x2 = 2.0d0*q*(1.0d0 + HeH)/(1.0d0 + q)
      if (x2 .gt. 1.0d0) x2 = 1.0d0
      h2_bound_fraction = 0.5d0*x2/(1.0d0 + HeH)
      end function h2_bound_fraction

      ! ------------------------------------------------------!

      function element_ratio_HeH(f_sp) result(heh_cell)
      ! Helium-to-hydrogen ELEMENT ratio n_He/n_H per cell: nuclei counted
      ! over every species that carries them, with the bsp_nH / bsp_nHe
      ! weights of species_table (H2 and H2+ carry two H nuclei, H3+ three,
      ! HeH+ one of each).  The metastable He 2^3S triplet and the molecular
      ! species are therefore inside the count, so the ratio is the same
      ! physical quantity in the atomic, triplet-helium and molecular
      ! regions; the columns of the species that a run does not carry are
      ! zero, so nothing has to be switched on the flags.  This is the same
      ! element count load_IC applies to a restart file.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng) :: heh_cell
      real*8, dimension(1-Ng:N+Ng) :: nuc_H, nuc_He
      integer :: ib

      nuc_H  = 0.0d0
      nuc_He = 0.0d0
      do ib = 1, n_bsp
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
