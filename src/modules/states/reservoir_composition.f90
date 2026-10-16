      module reservoir_composition
      ! Elemental reservoir abundances with the first physical cell's partition.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use global_parameters
      use species_table, only: isp_HI, isp_HII, isp_HeI, isp_HeII, isp_HeIII, &
           isp_HeTR, isp_H2, isp_H2p, isp_H3p, isp_HeHp, isp_OH, isp_H2O, &
           isp_CO, n_melem, n_mion, mion_fsp, mion_elem, melem_A,        &
           melem_i0, melem_top, iel_O, iel_C, n_bsp, bsp_fsp, bsp_nH,    &
           bsp_nHe, bsp_nO, bsp_nC, bsp_mass, bsp_is_excited_level
      use composition, only: base_h2_nuclei_fraction, base_h2_composition_imposed
      implicit none
      private
      public :: base_reservoir_composition_row
      contains
      subroutine base_reservoir_composition_row(f_cell1, f_row)
      ! THE COMPOSITION OF THE GAS THE RESERVOIR HOLDS AT THE BASE LEVEL,
      ! per unit mass, in the layout of f_sp: what the lower ghost cells of a
      ! restart carry instead of the rows the file wrote there.
      !
      ! TWO STATEMENTS MAKE IT, and both are declared inputs of the boundary.
      !
      !   THE ELEMENTAL ABUNDANCES ARE THE RESERVOIR'S: He/H and each trace
      !   element's El/H, the numbers the input states and the restart
      !   metadata block carries as its reservoir field. That is the rule
      !   this loader already applies to the base cells of a diffused state,
      !   where the column keeps its own separation and cells 1-Ng to 1 are
      !   carried onto the input's He/H (the He_diffusion branch above): the
      !   gas below the base is the inflow, and the lower atmosphere states
      !   what it is made of.
      !
      !   THE PARTITION WITHIN EACH ELEMENT IS THE FIRST PHYSICAL CELL'S:
      !   the ionization stages, the molecular ions and the oxygen carriers,
      !   as fractions of their own element's nuclei. The lower-atmosphere
      !   model never saw the wind's field, so it states no ionization; the
      !   ghost's own balance does, and the sweep solves it. What is needed
      !   here is a starting point in the right basin, and the cell half a
      !   grid spacing above the ghost is one. MEASURED without it, with the
      !   gas seeded neutral instead: on the hot WASP-121b base (T = 2358 K
      !   at the ghost) the ghost's ionization solve returned a non-root
      !   above the amnesty cap at both lower cells and the cells kept the
      !   neutral seed, which is not a state of the network.
      !
      !   The one partition the first cell does NOT state is the molecular
      !   one: the handoff states x2 for the inflowing gas
      !   (base_h2_nuclei_fraction, the single definition the equation of
      !   state shares), and the sweep closes it against the ghost's own
      !   ionization. Without a handoff the first cell's own H2 share is the
      !   seed, since there is then nothing upstream to impose.
      !
      ! MASS. Every species below is a density per unit mass whose nuclei are
      ! counted exactly once, and each element's neutral ground species takes
      ! the remainder of that element's nuclei, so the row weighs
      ! (1 + m_He He/H + sum El/H m_El)/mass_per_H = 1 exactly, under either
      ! metal policy, as the cold start's own base row does (set_IC builds the
      ! same composition for a whole column).
      !
      ! HeH+ carries one nucleus of each element, so no single factor
      ! rescales it: it takes the smaller of the two element factors and the
      ! nuclei that leaves short stay in the neutral ground species, which is
      ! how the He/H rescale of this loader treats it too.
      real*8, intent(in)  :: f_cell1(n_species)
      real*8, intent(out) :: f_row(n_species)
      real*8  :: nH_1, nHe_1, nH_r, nHe_r, s_he, x2
      real*8  :: n_hehp, n_h2, n_h2p, n_h3p, n_hii, n_oh, n_h2o, n_co
      real*8  :: nE_1, s_E, sO, sC
      real*8  :: n_h_available, n_h_molecular, n_h_oxygen, oxygen_scale
      integer :: e, k, i0

      f_row = 0.0d0

      ! Element nuclei of the first physical cell, over every species that
      ! carries them, and the reservoir's own counts per unit mass.
      nH_1 = f_cell1(isp_HI) + f_cell1(isp_HII)                           &
           + 2.0d0*(f_cell1(isp_H2) + f_cell1(isp_H2p))                   &
           + 3.0d0*f_cell1(isp_H3p) + f_cell1(isp_HeHp)                   &
           + f_cell1(isp_OH) + 2.0d0*f_cell1(isp_H2O)
      nHe_1 = f_cell1(isp_HeI) + f_cell1(isp_HeII)                        &
            + f_cell1(isp_HeIII) + f_cell1(isp_HeHp)
      nH_r  = 1.0d0/mass_per_H
      nHe_r = HeH/mass_per_H

      ! A first cell with no hydrogen states no partition; the reservoir is
      ! then neutral and the sweep starts from that.
      if (.not. (nH_1 .gt. 0.0d0)) then
         f_row(isp_HI) = nH_r
         if (thereis_He) f_row(isp_HeI) = nHe_r
         do e = 1, n_melem
            f_row(mion_fsp(melem_i0(e))) = melem_ab(e)/mass_per_H
         enddo
         call verify_reservoir_nuclei(f_row)
         return
      endif

      s_he = 0.0d0
      if (thereis_He .and. nHe_1 .gt. 0.0d0) s_he = nHe_r/nHe_1

      ! ---- helium, and the molecule that carries a nucleus of each -------
      n_hehp = 0.0d0
      if (thereis_mol .and. f_cell1(isp_HeHp) .gt. 0.0d0)                 &
         n_hehp = min(f_cell1(isp_HeHp)*nH_r/nH_1,                        &
                      f_cell1(isp_HeHp)*s_he)
      if (thereis_He) then
         f_row(isp_HeII)  = f_cell1(isp_HeII) *s_he
         f_row(isp_HeIII) = f_cell1(isp_HeIII)*s_he
         f_row(isp_HeTR)  = f_cell1(isp_HeTR) *s_he
         ! He I is the TOTAL neutral helium, the 2^3S level included, so the
         ! triplet is a share of it and not a further nucleus.
         f_row(isp_HeI)   = max(nHe_r - f_row(isp_HeII)                   &
                                - f_row(isp_HeIII) - n_hehp, 0.0d0)
      endif
      if (thereis_mol) f_row(isp_HeHp) = n_hehp

      ! ---- the trace elements, stage by stage ----------------------------
      do e = 1, n_melem
         i0   = melem_i0(e)
         nE_1 = 0.0d0
         do k = 0, melem_top(e)
            nE_1 = nE_1 + f_cell1(mion_fsp(i0+k))
         enddo
         if (thereis_oxychem .and. e .eq. iel_O)                          &
            nE_1 = nE_1 + f_cell1(isp_OH) + f_cell1(isp_H2O)              &
                 + f_cell1(isp_CO)
         if (thereis_oxychem .and. e .eq. iel_C)                          &
            nE_1 = nE_1 + f_cell1(isp_CO)
         if (nE_1 .gt. 0.0d0) then
            s_E = (melem_ab(e)/mass_per_H)/nE_1
            do k = 1, melem_top(e)
               f_row(mion_fsp(i0+k)) = f_cell1(mion_fsp(i0+k))*s_E
            enddo
         else
            s_E = 0.0d0
         endif
         ! The neutral stage takes the remainder of this element's nuclei;
         ! the oxygen and carbon carriers are removed from it below.
         f_row(mion_fsp(i0)) = melem_ab(e)/mass_per_H
         do k = 1, melem_top(e)
            f_row(mion_fsp(i0)) = f_row(mion_fsp(i0))                     &
                                - f_row(mion_fsp(i0+k))
         enddo
         f_row(mion_fsp(i0)) = max(f_row(mion_fsp(i0)), 0.0d0)
      enddo

      ! ---- the oxygen carriers, on the same rule --------------------------
      n_oh = 0.0d0;  n_h2o = 0.0d0;  n_co = 0.0d0
      if (thereis_oxychem) then
         sO = 0.0d0;  sC = 0.0d0
         nE_1 = f_cell1(isp_OH) + f_cell1(isp_H2O) + f_cell1(isp_CO)
         do k = 0, melem_top(iel_O)
            nE_1 = nE_1 + f_cell1(mion_fsp(melem_i0(iel_O)+k))
         enddo
         if (nE_1 .gt. 0.0d0) sO = (melem_ab(iel_O)/mass_per_H)/nE_1
         nE_1 = f_cell1(isp_CO)
         do k = 0, melem_top(iel_C)
            nE_1 = nE_1 + f_cell1(mion_fsp(melem_i0(iel_C)+k))
         enddo
         if (nE_1 .gt. 0.0d0) sC = (melem_ab(iel_C)/mass_per_H)/nE_1
         n_oh  = f_cell1(isp_OH) *sO
         n_h2o = f_cell1(isp_H2O)*sO
         ! CO holds a nucleus of each, and takes the smaller share.
         n_co  = f_cell1(isp_CO)*min(sO, sC)
         f_row(isp_OH)  = n_oh
         f_row(isp_H2O) = n_h2o
         f_row(isp_CO)  = n_co
         f_row(mion_fsp(melem_i0(iel_O))) =                               &
              max(f_row(mion_fsp(melem_i0(iel_O)))                        &
                  - n_oh - n_h2o - n_co, 0.0d0)
         f_row(mion_fsp(melem_i0(iel_C))) =                               &
              max(f_row(mion_fsp(melem_i0(iel_C))) - n_co, 0.0d0)
      endif

      ! ---- hydrogen: the handoff's molecular partition, the first cell's
      !      ionization, and the neutral remainder -------------------------
      !
      ! WHAT THE HANDOFF STATES IS THE PARTITION OF THE NON-IONIZED
      ! HYDROGEN. The handoff's x2 is the non-atomic share of neutral H; it
      ! includes H bound into oxygen carriers. Those carriers take their H
      ! from x2 of the non-ionized budget before H2 is assigned. The
      ! remaining (1-x2) share stays atomic.
      n_hii = f_cell1(isp_HII)*nH_r/nH_1
      n_h2  = 0.0d0;  n_h2p = 0.0d0;  n_h3p = 0.0d0
      if (thereis_mol) then
         n_h2p = f_cell1(isp_H2p)*nH_r/nH_1
         n_h3p = f_cell1(isp_H3p)*nH_r/nH_1
      endif
      n_h_available = max(nH_r - n_hii - 2.0d0*n_h2p - 3.0d0*n_h3p &
                          - n_hehp, 0.0d0)
      x2 = 1.0d0
      if (thereis_mol) then
         if (base_h2_composition_imposed()) then
            x2 = base_h2_nuclei_fraction()
         else
            ! No handoff fixes x2: the first cell's H2 and oxygen carriers
            ! together state its non-atomic share of non-ionized hydrogen.
            x2 = (2.0d0*f_cell1(isp_H2) + f_cell1(isp_OH) +            &
                  2.0d0*f_cell1(isp_H2O))/                         &
                  max(nH_1-f_cell1(isp_HII)-2.0d0*f_cell1(isp_H2p) &
                      -3.0d0*f_cell1(isp_H3p)-f_cell1(isp_HeHp),     &
                      tiny(1.0d0))
         endif
      endif
      n_h_molecular = x2*n_h_available
      n_h_oxygen = n_oh + 2.0d0*n_h2o
      if (n_h_oxygen .gt. n_h_molecular) then
         ! The oxygen carriers draw from the non-atomic budget. Return
         ! excess O to its neutral stage, preserving the atomic H share.
         oxygen_scale = n_h_molecular/n_h_oxygen
         f_row(mion_fsp(melem_i0(iel_O))) =                            &
              f_row(mion_fsp(melem_i0(iel_O))) +                        &
              (1.0d0-oxygen_scale)*(n_oh+n_h2o)
         n_oh = n_oh*oxygen_scale
         n_h2o = n_h2o*oxygen_scale
         f_row(isp_OH) = n_oh
         f_row(isp_H2O) = n_h2o
         n_h_oxygen = n_oh + 2.0d0*n_h2o
      endif
      if (thereis_mol) then
         ! H2 takes the part of the declared non-atomic budget not carried
         ! by OH and H2O.
         n_h2 = 0.5d0*max(n_h_molecular-n_h_oxygen,0.0d0)
      endif
      f_row(isp_HII) = n_hii
      f_row(isp_H2)  = n_h2
      f_row(isp_H2p) = n_h2p
      f_row(isp_H3p) = n_h3p
      f_row(isp_HI)  = max(nH_r - n_hii - 2.0d0*(n_h2 + n_h2p)            &
                           - 3.0d0*n_h3p - n_hehp - n_oh - 2.0d0*n_h2o,   &
                           0.0d0)

      call verify_reservoir_nuclei(f_row)

      end subroutine base_reservoir_composition_row

      subroutine verify_reservoir_nuclei(f_row)
      ! The species table is the same stoichiometric source as the element
      ! census. Check this one row directly, without a grid-sized census.
      real*8, intent(in) :: f_row(n_species)
      real*8 :: nuclei(2+n_melem), expected(2+n_melem), rho_comp, tol
      integer :: ib, im, e, k
      nuclei = 0.0d0
      rho_comp = 0.0d0
      do ib=1,n_bsp
         if (bsp_is_excited_level(ib)) cycle
         k=bsp_fsp(ib)
         nuclei(1)=nuclei(1)+dble(bsp_nH(ib))*f_row(k)
         nuclei(2)=nuclei(2)+dble(bsp_nHe(ib))*f_row(k)
         nuclei(2+iel_O)=nuclei(2+iel_O)+dble(bsp_nO(ib))*f_row(k)
         nuclei(2+iel_C)=nuclei(2+iel_C)+dble(bsp_nC(ib))*f_row(k)
         rho_comp=rho_comp+bsp_mass(ib)*f_row(k)
      enddo
      do im=1,n_mion
         nuclei(2+mion_elem(im))=nuclei(2+mion_elem(im))              &
                               +f_row(mion_fsp(im))
         if (eos_include_metals .and. thereis_metals)                  &
            rho_comp=rho_comp+melem_A(mion_elem(im))*f_row(mion_fsp(im))
      enddo
      expected(1)=1.0d0/mass_per_H
      expected(2)=0.0d0
      if (thereis_He) expected(2)=HeH/mass_per_H
      do e=1,n_melem
         expected(2+e)=melem_ab(e)/mass_per_H
      enddo
      tol=256.0d0*epsilon(1.0d0)
      if (.not.all(ieee_is_finite(f_row)) .or. any(f_row.lt.0.0d0) .or. &
          any(abs(nuclei-expected).gt.tol*max(expected,1.0d-30)) .or. &
          abs(rho_comp-1.0d0).gt.tol) then
         write(*,'(A)') ' ERROR: reservoir composition violates nuclei or mass'
         stop 1
      endif
      end subroutine verify_reservoir_nuclei
      end module reservoir_composition
