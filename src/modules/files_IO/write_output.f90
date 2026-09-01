      module output_write
      ! Write the output to the standard output files

      use global_parameters
      use species_table, only: n_mion, mion_name, im_OI, melem_i0,       &
                               iel_O, iel_C
      use ionization_equilibrium, only: nmol_eq,   &  ! molecular columns
                                       NH2_col_lw, f_shield_lw, k_lw_diss, &
                                       tr_lines_lw, a_lines_lw, P_H2_eq,   &
                                       NH2_db96_max, NH2_richings_max,     &
                                       nox_eq, n_o1d_eq,                   &
                                       NH2O_col, NOH_col,                  &
                                       j_h2o_fuv, j_oh_fuv, tau_fuv,       &
                                       heat_fuv
      use utils_ion_eq, only: fuv_band_flux
      use Cooling_Coefficients, only: base_sky_fraction
      use molecular_infrared_cooling, only: h2_line_emission_lte,          &
                            h2o_band_emission_lte, co_band_emission_lte,   &
                            h2_line_net_cooling_rate,                     &
                            h2o_band_net_cooling_rate,                    &
                            co_band_net_cooling_rate
      use lyman_werner_photodissociation, only: e_lw_fragment_erg,         &
                                       e_lw_photon_erg, p_diss_lw
      use water_photolysis, only: n_fuv_band, fuv_band_name,               &
                                       sigma_H2O_band, sigma_OH_band,      &
                                       e_photon_H2O_band, e_photon_OH_band,&
                                       e_photon_flat_band,                 &
                                       fuv_band_photon_flux,               &
                                       heat_per_water_dissociation,        &
                                       heat_per_hydroxyl_dissociation,     &
                                       ib_LW
      use oxygen_rates, only: rk_O1_OH_H2_water, rk_O2_O_H2_hydroxyl,     &
                              rate_from_detailed_balance,                 &
                              ith_H, ith_H2, ith_O, ith_OH, ith_H2O
      use mol_rates, only: rk_R12_H2_thdis, rk_R10_Hp_H2v4,               &
                           rk_R13_Hp_H2_M, rk_R14_H2_edis, rk_R8_H2p_H2,  &
                           rk_R17_Hep_H2_diss, rk_R20_Hep_H2_HeHp,        &
                           rk_R23_H2_Hep_cx, rk_R18_HeHp_H2,              &
                           rk_R15_3body_H2, rk_R6_H3p_dr_H2,              &
                           rk_R9_H2p_H, rk_R11_H3p_H
      use Cooling_Coefficients, only: ioniz_HeI23S_H2
      use diffusive_photochemistry, only: carrier_diffusion_coefficient,   &
                                          carrier_transport_diagnostics,   &
                                          carrier_co_ceiling_cells, ic_H2
      use utils, only: calc_ne, calc_ntot
      ! O I ground-term statistical equilibrium: the same solution the
      ! [O I] fine-structure cooling is built on (Cool_coeff.f90).
      use Cooling_Coefficients, only: n_fsline,                        &
                                      fine_structure_line_transfer,    &
                                      oxygen_ground_term_levels

      contains

      subroutine write_output(rho,v,p,T,heat,cool,eta,                &
                              nhi,nhii,nhei,nheii,nheiii,nheiTR,      &
                              nm,flag)
      ! Metal ion densities are passed as the 2D array nm(:, 1:n_mion),
      ! one column per metal ion stage in the canonical species_table
      ! order (CI, CII, CIII, OI, ..., MgIII). This keeps the argument
      ! list fixed as metals are added.

      character(len = 2) :: flag
      integer :: j,i
      real*8  :: nh_lw     ! H nuclei density [cm^-3], LW diagnostic
      ! O I ground-term level output: dimensional state, the line transfer
      ! the populations depend on, and the three fractional populations.
      real*8, dimension(1-Ng:N+Ng) :: T_K_oi,ne_oi,f_3P2,f_3P1,f_3P0
      real*8, dimension(1-Ng:N+Ng,n_mion)   :: nm_oi
      real*8, dimension(1-Ng:N+Ng,4)        :: nmol_oi
      real*8, dimension(1-Ng:N+Ng,n_fsline) :: beta_fs_oi, nbar_fs_oi
      real*8 :: noi_sum, oi_close
      real*8, dimension(1-Ng:N+Ng), intent(in) :: rho,v,p,T
      real*8, dimension(1-Ng:N+Ng), intent(in) :: heat,cool
      real*8, dimension(1-Ng:N+Ng), intent(in) :: eta
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi,nhii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei,nheii,nheiii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nheiTR
      real*8, dimension(1-Ng:N+Ng,n_mion), intent(in) :: nm

      
      !---- Write thermodynamic profiles ----! 
          
      if (flag.eq.'eq') then
      	open(unit = 2, file = './output/Hydro_ioniz.txt')
	   else	! Change output file after postprocessing
		   open(unit = 2, file = './output/Hydro_ioniz_adv.txt')
	   endif

      ! Schema header ('#' comment lines; readers that predate the header
      ! can skip them, numeric content is unchanged)
      ! Column 2 is rho*n0: the MASS density in units of m_H per cm^3 (metals
      ! included under the eos_metals policy), not a number density. Multiply by
      ! m_H to get g/cm^3.
      write(2,'(A)') '# EXHALE schema 2'
      write(2,'(A)') '# columns r[Rp] rho[mH/cm3] v[cm/s] p[cgs] T[K] '//   &
                     'heat[erg/cm3/s] cool[erg/cm3/s]'

         do j = 1-Ng,N+Ng
            write(2,*) r(j),        &     ! Rad. dist.
                     rho(j)*n0,     &     ! Density
                     v(j)*v0,       &     ! Velocity
                     p(j)*p0,       &     ! Pressure
                     T(j)*T0,       &     ! Temperature
                     heat(j)*q0,    &     ! Rad. heat.
                     cool(j)*q0           ! Rad. cool.
         enddo
      close(2)
      
      !---- Write ionization profiles ----!  
      if (flag.eq.'eq') then
      	open(unit = 3, file = './output/Ion_species.txt')
	   else	! Change output file after postprocessing
		   open(unit = 3, file = './output/Ion_species_adv.txt')
	   endif

      ! Schema header: species labels in column order, generated from the
      ! species table so they stay correct when species are added.
      write(3,'(A)') '# EXHALE schema 2'
      ! The advection post-process reconstructs an H/He + trace-metal gas and
      ! does not carry the molecular or oxygen species: their columns below
      ! are the EQUILIBRIUM values, written beside advection-corrected H, He
      ! and metal columns. Say so in the file rather than only in the module
      ! that produces it, so that a reader of Ion_species_adv.txt cannot take
      ! them for corrected profiles. Lifting the limitation is separate work.
      if (flag .ne. 'eq' .and. (thereis_mol .or. thereis_oxychem)) then
         write(3,'(A)') '# NOTE the molecular columns (H2 H2p H3p HeHp)'// &
                        ' and, when present, the oxygen columns'
         write(3,'(A)') '#   (OH H2O CO) are NOT advection-corrected:'//   &
                        ' the post-process reconstructs an'
         write(3,'(A)') '#   H/He + trace-metal gas, so those columns'//   &
                        ' are the equilibrium solution and the'
         write(3,'(A)') '#   metal columns inside the molecular layer'//   &
                        ' inherit that approximation.'
      endif
      write(3,'(A)', advance='no') '# columns r[Rp] HI HII HeI HeII '// &
                                   'HeIII HeITR'
      do i = 1,n_mion
         write(3,'(A)', advance='no') ' '//trim(mion_name(i))
      enddo
      ! molecular columns (present only when thereis_mol)
      if (thereis_mol) write(3,'(A)', advance='no') ' H2 H2p H3p HeHp'
      ! oxygen-carrier columns (present only when thereis_oxychem)
      if (thereis_oxychem) write(3,'(A)', advance='no') ' OH H2O CO'
      write(3,'(A)') ''

      do j = 1-Ng,N+Ng

         if (thereis_oxychem) then
            write(3,*) r(j), nhi(j)*n0, nhii(j)*n0, nhei(j)*n0,        &
                     nheii(j)*n0, nheiii(j)*n0, nheiTR(j)*n0,          &
                     (nm(j,i)*n0, i = 1,n_mion),                       &
                     (nmol_eq(j,i), i = 1,4),                          &
                     (nox_eq(j,i), i = 1,3)   ! OH H2O CO (already cm^-3)
         else if (thereis_mol) then
            write(3,*) r(j), nhi(j)*n0, nhii(j)*n0, nhei(j)*n0,        &
                     nheii(j)*n0, nheiii(j)*n0, nheiTR(j)*n0,          &
                     (nm(j,i)*n0, i = 1,n_mion),                       &
                     (nmol_eq(j,i), i = 1,4)  ! H2 H2+ H3+ HeH+ (already cm^-3)
         else
         write(3,*) r(j),       & ! Rad. dist.
                  nhi(j)*n0,    & ! HI
                  nhii(j)*n0,   & ! HII
                  nhei(j)*n0,   & ! HeI
                  nheii(j)*n0,  & ! HeII
                  nheiii(j)*n0, & ! HeIII
                  nheiTR(j)*n0, & ! HeITR
                  (nm(j,i)*n0, i = 1,n_mion)  ! metal ions (canonical order)
         endif
      enddo
      close(3)

      !---- Lyman-Werner photodissociation diagnostic ----!
      ! Written only for a molecular run that carries a Lyman-Werner band
      ! flux, and only for the equilibrium state (the advection
      ! post-process is atomic and does not re-solve the molecules).
      ! It is the record of how deep the band penetrates: f_shield -> 1 in
      ! the thin wind above the H2 -> H front and collapses in the
      ! self-shielded molecular base.  The f_shield column is the factor
      ! the RATE carries, i.e. the Richings, Schaye & Oppenheimer (2014)
      ! temperature-dependent fit; the band share the H2 lines take from
      ! the shared FUV beam is a separate quantity on a separate fit
      ! (lyman_werner.f90 sec. 2); it matters only when the oxygen
      ! chemistry shares the band, and is reported in output/FUV_bands.txt
      ! when it does.
      if (thereis_mol .and. F_LW_star .gt. 0.0d0 .and. flag .eq. 'eq') then
         open(unit = 4, file = './output/Lyman_Werner.txt')
         write(4,'(A)') '# EXHALE schema 2'
         write(4,'(A,ES12.5,A)') '# Lyman-Werner band flux at the planet: ', &
                        F_LW_star, ' erg cm^-2 s^-1 (912-1110 A)'
         write(4,'(A)') '# columns r[Rp] T[K] x_H2[2nH2/nH] nH2[cm^-3] '//  &
                        'NH2_star[cm^-2] f_shield k_LW[1/s] '//             &
                        'heat_LW[erg/cm3/s]'
         do j = 1-Ng,N+Ng
            nh_lw = (nhi(j) + nhii(j))*n0                                  &
                  + 2.0d0*(nmol_eq(j,1) + nmol_eq(j,2))                    &
                  + 3.0d0*nmol_eq(j,3) + nmol_eq(j,4)
            write(4,*) r(j), T(j)*T0,                                      &
                       2.0d0*nmol_eq(j,1)/max(nh_lw,1.0d-99),              &
                       nmol_eq(j,1), NH2_col_lw(j), f_shield_lw(j),        &
                       k_lw_diss(j),                                       &
                       k_lw_diss(j)*nmol_eq(j,1)*e_lw_fragment_erg
         enddo
         close(4)
      endif

      !---- Oxygen chemistry and the FUV bands ----!
      ! Written only for an oxygen-chemistry run and only for the
      ! equilibrium state (the advection post-process is molecule-free and
      ! does not re-solve the oxygen carriers).
      if (thereis_oxychem .and. flag .eq. 'eq') call write_oxygen_chemistry(&
                              v, T, nhi, nhii, nhei, nheii, nheiii,       &
                              nheiTR, nm)

      !---- O I ground-term level populations ----!
      ! The lower levels of the O I 1302.168 / 1304.858 / 1306.029 A
      ! resonance triplet are the three fine-structure levels of the 2p4 3P
      ! ground term, so a transit forward model needs them resolved: the
      ! total O I density applied to all three components would count the
      ! same atoms three times. The populations written here are not a
      ! second calculation of that equilibrium -- oxygen_ground_term_levels
      ! is the same statistical-equilibrium solver, with the same
      ! collisional rates and the same escape probabilities / incident
      ! field, that the [O I] 63/145/44um cooling is assembled from
      ! (Cool_coeff.f90). Written for metal-bearing runs only, in both the
      ! equilibrium and the advection-corrected state, so the transit tool
      ! reads the '_adv' file next to Ion_species_adv.txt. A metals-off run
      ! writes no file and the transit tool skips the line, exactly as it
      ! already does for the other metal lines.
      if (thereis_metals) then
         T_K_oi = T*T0
         do i = 1,n_mion
            nm_oi(:,i) = nm(:,i)*n0
         enddo
         nmol_oi = 0.0d0
         if (thereis_mol) nmol_oi = nmol_eq
         call calc_ne(nhii*n0,nheii*n0,nheiii*n0,ne_oi,nm_oi,nmol_oi)
         call fine_structure_line_transfer(T_K_oi,nm_oi,               &
                                           beta_fs_oi,nbar_fs_oi)
         call oxygen_ground_term_levels(T_K_oi,ne_oi,nhi*n0,           &
                                        beta_fs_oi,nbar_fs_oi,         &
                                        f_3P2,f_3P1,f_3P0)

         ! Gate: the three level densities must sum to the total O I.
         oi_close = 0.0d0
         do j = 1-Ng,N+Ng
            noi_sum = (f_3P2(j) + f_3P1(j) + f_3P0(j))*nm_oi(j,im_OI)
            oi_close = max(oi_close, abs(noi_sum - nm_oi(j,im_OI))     &
                                     /max(nm_oi(j,im_OI),1.0d-99))
         enddo
         write(*,'(a,a,a,es9.2)') ' (write_output/',flag,              &
              ') max |sum(O I levels)/n(O I) - 1| = ', oi_close

         if (flag.eq.'eq') then
            open(unit = 73, file = './output/OI_levels.txt')
         else
            open(unit = 73, file = './output/OI_levels_adv.txt')
         endif
         write(73,'(A)') '# EXHALE schema 2'
         write(73,'(A)') '# O I 2p4 3P ground-term level populations, from '// &
                        'the same three-level statistical'
         write(73,'(A)') '# equilibrium as the [O I] 63.2/145.5/44.1um '//     &
                        'cooling (Cool_coeff.f90).'
         write(73,'(A)') '# The 3P2 / 3P1 / 3P0 levels are the lower levels '//&
                        'of O I 1302.168 / 1304.858 / 1306.029 A.'
         write(73,'(A)') '# columns r[Rp] T[K] ne[cm^-3] nHI[cm^-3] '//        &
                        'nOI[cm^-3] f_3P2 f_3P1 f_3P0 '//                     &
                        'n_3P2[cm^-3] n_3P1[cm^-3] n_3P0[cm^-3]'
         do j = 1-Ng,N+Ng
            write(73,*) r(j), T_K_oi(j), ne_oi(j), nhi(j)*n0,           &
                       nm_oi(j,im_OI),                                 &
                       f_3P2(j), f_3P1(j), f_3P0(j),                   &
                       f_3P2(j)*nm_oi(j,im_OI),                        &
                       f_3P1(j)*nm_oi(j,im_OI),                        &
                       f_3P0(j)*nm_oi(j,im_OI)
         enddo
         close(73)
      endif

      ! End of subroutine
      end subroutine write_output

      !-------------------------------------------------------------!

      subroutine write_oxygen_chemistry(v, T, nhi, nhii, nhei, nheii,     &
                                        nheiii, nheiTR, nm)
      ! The two records of an oxygen-chemistry run:
      !
      !   output/Oxygen_chemistry.txt  the solved oxygen partition, the
      !       O(1D) steady state, and the chemical against the advection
      !       time scale, cell by cell;
      !   output/FUV_bands.txt         how deep each FUV band penetrates,
      !       and the band-by-band energy ledger of gate G4.
      !
      ! WHY THE TIME SCALES ARE A REQUIRED OUTPUT AND NOT A DIAGNOSTIC.
      ! With "Oxygen transport: True" -- the default of the option -- the
      ! carriers are transported, so a cell in which the chemistry is slow
      ! is SOLVED rather than assumed away. The time scales stay a required
      ! output for the other direction: they say where the answer is
      ! chemistry, where it is flow, and where the two are comparable, and
      ! that is what tells a reader which part of the profile the reaction
      ! set is responsible for. With "Oxygen transport: False" they say
      ! something stronger -- where that limit is simply wrong.
      !
      !   tau_chem  the H2 lifetime against the oxygen cycle,
      !             1/(k_O1 n_OH + k_O2 n_O), i.e. the time in which the
      !             two H2-consuming reactions of the network would remove
      !             the H2 of the cell;
      !   tau_adv   the radial flow time r/|v|;
      !   tau_diff  the diffusive time of the cell for H2,
      !             dr^2/(D_H2 + K_zz);
      !   Da        tau_adv/tau_chem. Da >> 1 is the local-equilibrium
      !             limit; Da <~ 1 means transport sets the partition, and
      !             it is the measurement that made the transport operator
      !             necessary (design sec. 3.1).
      !
      ! All are written for every cell, including where the wind is
      ! subsonic and |v| is small; tau_adv is then large and Da large, which
      ! is the correct statement (a stagnant cell is in local equilibrium).

      integer :: j, ib
      real*8, dimension(1-Ng:N+Ng), intent(in) :: v, T
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi, nhii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei, nheii, nheiii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nheiTR
      real*8, dimension(1-Ng:N+Ng,n_mion), intent(in) :: nm
      ! Gas-particle and electron densities of the cell, rebuilt here for
      ! the H2 loss budget: the third body M of the thermal channel and the
      ! electron of the impact-dissociation one.
      real*8, dimension(1-Ng:N+Ng) :: n_tot_ox, ne_ox
      real*8 :: nh_ox, tau_chem, tau_adv, da_num, k1, k2, n_o0
      real*8 :: d_h2_cell, tau_dif, ct_resid, ct_worst
      integer :: ct_steps, ct_nlim
      ! H2 loss budget at the base cell: which channel runs the partition.
      integer :: jb
      character(len=400) :: colhdr
      real*8  :: nH2b, Tb, nOb, nHIb, h2loss(11), h2tot
      real*8  :: ok1b, ok2b, ok1rb, ok2rb
      character(len=22), parameter :: h2ch_name(11) = (/                  &
           'thermal H2 + M net R15', 'Lyman-Werner photodiss', &
           'H2 photoionization    ', 'H+ + H2               ', &
           'e- impact dissociation', 'H2+ + H2              ', &
           'He+ + H2              ', 'HeH+ + H2             ', &
           'He(2^3S) + H2         ', 'oxygen O1+O2 net rev  ', &
           'mol-ion return R6R9R11' /)
      real*8 :: dr_cm, absph(n_fuv_band), absen(n_fuv_band)
      real*8 :: heat_col(n_fuv_band), bond_col(n_fuv_band)
      real*8 :: nph_b, expected, relerr, worst, e_lw_abs, e_lw_in, ph_lw
      ! How much of the oxygen and the carbon the option has moved out of
      ! the atomic coolants, and where.
      ! Column-integrated infrared exchange of the molecular bands, per unit
      ! area: spontaneous emission, absorption of the field from below, and
      ! their difference. Index 1 H2, 2 H2O, 3 CO.
      real*8 :: ir_emit(3), ir_abs(3), ir_e, ir_n, w_ir, ir_bound
      real*8 :: nO_el, nC_el, fO_mol, fC_mol, fO_worst, fC_worst
      real*8 :: rO_worst, rC_worst, TO_worst, TC_worst

      open(unit = 74, file = './output/Oxygen_chemistry.txt')
      write(74,'(A)') '# EXHALE schema 2'
      write(74,'(A)') '# In-code oxygen chemistry (docs/a2_oxygen_option'// &
                     '_design.md). Densities in cm^-3.'
      write(74,'(A)') '# n_O is FREE ATOMIC neutral oxygen; the element'//  &
                     ' total is n_O + n_OII + n_OIII + n_OH + n_H2O + n_CO.'
      write(74,'(A)') '# tau_chem = 1/(k_O1 n_OH + k_O2 n_O), the H2'//     &
                     ' lifetime against the oxygen cycle; tau_adv = r/|v|;'
      write(74,'(A)') '# Da = tau_adv/tau_chem. Da >> 1 is the local'//     &
                     '-equilibrium limit; Da <~ 1 means transport sets'
      write(74,'(A)') '# the partition. tau_diff = dr^2/(D_H2 + K_zz) is'//&
                     ' the diffusive time of the cell for H2.'
      if (oxygen_transport) then
         write(74,'(A)') '# The carriers ARE transported (implicit'//      &
                     ' diffusion-advection with the chemistry).'
      else
         write(74,'(A)') '# The carriers are a LOCAL steady state'//       &
                     ' ("Oxygen transport: False"): the chemistry alone.'
      endif
      write(74,'(A)') '# columns r[Rp] T[K] n_O n_OII n_OIII n_OH n_H2O'//  &
                     ' n_CO n_O1D x_H2[2nH2/nH] tau_chem[s] tau_adv[s]'//   &
                     ' Da j_H2O[1/s] j_OH[1/s] D_H2[cm2/s] Kzz[cm2/s]'//    &
                     ' tau_diff[s]'
      do j = 1-Ng,N+Ng
         nh_ox = (nhi(j) + nhii(j))*n0                                     &
               + 2.0d0*(nmol_eq(j,1) + nmol_eq(j,2))                       &
               + 3.0d0*nmol_eq(j,3) + nmol_eq(j,4)                         &
               + nox_eq(j,1) + 2.0d0*nox_eq(j,2)
         k1   = rk_O1_OH_H2_water(T(j)*T0)
         k2   = rk_O2_O_H2_hydroxyl(T(j)*T0)
         n_o0 = nm(j,melem_i0(iel_O))*n0
         da_num = k1*nox_eq(j,1) + k2*n_o0
         if (da_num .gt. 0.0d0) then
            tau_chem = 1.0d0/da_num
         else
            tau_chem = huge(1.0d0)
         endif
         if (abs(v(j)) .gt. 0.0d0) then
            tau_adv = r(j)*R0/(abs(v(j))*v0)
         else
            tau_adv = huge(1.0d0)
         endif
         if (tau_chem .gt. 0.0d0 .and. tau_chem .lt. huge(1.0d0)) then
            da_num = tau_adv/tau_chem
         else
            da_num = 0.0d0
         endif
         d_h2_cell = carrier_diffusion_coefficient(j, ic_H2)
         dr_cm     = dr_j(j)*R0
         tau_dif   = dr_cm*dr_cm/max(d_h2_cell + kzz_cell(j), 1.0d-30)
         write(74,*) r(j), T(j)*T0, n_o0, nm(j,melem_i0(iel_O)+1)*n0,      &
                    nm(j,melem_i0(iel_O)+2)*n0,                            &
                    nox_eq(j,1), nox_eq(j,2), nox_eq(j,3), n_o1d_eq(j),    &
                    2.0d0*nmol_eq(j,1)/max(nh_ox,1.0d-99),                 &
                    tau_chem, tau_adv, da_num,                             &
                    sum(j_h2o_fuv(j,:)), sum(j_oh_fuv(j,:)),               &
                    d_h2_cell, kzz_cell(j), tau_dif
      enddo
      if (oxygen_transport) then
         call carrier_transport_diagnostics(ct_steps, ct_resid, ct_nlim,  &
                                            ct_worst)
         write(74,'(A)') '#'
         write(74,'(A)') '# Carrier transport, last step: Newton'//        &
                     ' iterations, the relative residual it stopped at,'
         write(74,'(A)') '# and how many cells the element limiter'//      &
                     ' touched (a carrier that would have taken more'
         write(74,'(A)') '# nuclei than its element has). A limiter'//     &
                     ' count that is not small is a statement about'
         write(74,'(A)') '# the step, not a rounding: it is reported'//    &
                     ' rather than absorbed.'
         write(74,'(A,I0)')     '# newton_steps   ', ct_steps
         write(74,'(A,ES12.4)') '# newton_resid   ', ct_resid
         write(74,'(A,I0)')     '# limited_cells  ', ct_nlim
         write(74,'(A,ES12.4)') '# worst_overshoot', ct_worst
         write(74,'(A)') '#'
         write(74,'(A)') '# co_ceiling_cells: cells in which the'//       &
                     ' transported CO was cut back to the CO <-> C + O'
         write(74,'(A)') '# chemical equilibrium of their own (n, T).'//  &
                     ' Decision D4 makes CO chemically inert, and'
         write(74,'(A)') '# transported that would mean'//                &
                     ' INDESTRUCTIBLE: the wind would carry it to'
         write(74,'(A)') '# 2e4 K and hold the whole oxygen and carbon'//&
                     ' inventory there. Where the equilibrium'
         write(74,'(A)') '# forbids CO it is removed; where the'//        &
                     ' equilibrium allows it the transported value'
         write(74,'(A)') '# stands, so this is a one-sided constraint'// &
                     ' and not a return to equilibrium.'
         write(74,'(A,I0)') '# co_ceiling_cells ',                       &
                     carrier_co_ceiling_cells()
      endif

      ! ---- what actually runs the H2 partition at the base ----
      ! The option exists to make the OXYGEN cycle set the base H2/H
      ! partition, and whether it does is a property of the run rather than
      ! of the code: on a 2331 K base the thermal channel H2 + M -> H + H + M
      ! carries 70-97% of the net and the oxygen cycle a few percent, while
      ! on a 864 K base the oxygen family carries 96-99.6%
      ! (docs/a2_oxygen_option_design.md sec. 2.2). A run whose base is
      ! hotter than it should be -- and every molecular run is, until the
      ! H2O and CO infrared bands of item (G) exist -- is answering a
      ! different question, so the decomposition is written out instead of
      ! being left to be inferred from the partition.
      !
      ! Each entry is the NET rate at which that channel removes H2 from the
      ! base cell [cm^-3 s^-1] -- net of its own reverse where it has one,
      ! which is the convention section 2.2 of the design measures its
      ! budget in, because the reversible pairs run three decades above
      ! their net on a hot base. A negative entry is a net H2 SOURCE, and a
      ! share can therefore exceed 1 or be negative; they sum to 1 by
      ! construction.
      call calc_ne(nhii*n0, nheii*n0, nheiii*n0, ne_ox, nm*n0, nmol_eq)
      call calc_ntot(nhi*n0, nhii*n0, nhei*n0, nheii*n0, nheiii*n0,       &
                     n_tot_ox, nm*n0, nmol_eq, nox_eq)
      jb = 1
      nH2b = nmol_eq(jb,1)
      Tb   = T(jb)*T0
      nOb  = nm(jb,melem_i0(iel_O))*n0
      nHIb = nhi(jb)*n0
      ok1b = rk_O1_OH_H2_water(Tb)
      ok2b = rk_O2_O_H2_hydroxyl(Tb)
      ok1rb = rate_from_detailed_balance(ok1b, (/ ith_OH, ith_H2 /),      &
                                                (/ ith_H2O, ith_H /), Tb)
      ok2rb = rate_from_detailed_balance(ok2b, (/ ith_O, ith_H2 /),       &
                                                (/ ith_OH, ith_H /), Tb)
      h2loss    = 0.0d0
      ! (1) thermal dissociation net of the three-body association R15
      h2loss(1) = rk_R12_H2_thdis(Tb)*n_tot_ox(jb)*nH2b                   &
                - rk_R15_3body_H2(Tb, n_tot_ox(jb))*nHIb*nHIb
      h2loss(2) = k_lw_diss(jb)*nH2b
      h2loss(3) = P_H2_eq(jb)*nH2b
      h2loss(4) = (rk_R10_Hp_H2v4(Tb)                                      &
                 + rk_R13_Hp_H2_M(n_tot_ox(jb)))*nhii(jb)*n0*nH2b
      h2loss(5) = rk_R14_H2_edis(Tb)*ne_ox(jb)*nH2b
      h2loss(6) = rk_R8_H2p_H2()*nmol_eq(jb,2)*nH2b
      h2loss(7) = (rk_R17_Hep_H2_diss(Tb) + rk_R20_Hep_H2_HeHp()          &
                 + rk_R23_H2_Hep_cx())*nheii(jb)*n0*nH2b
      h2loss(8) = rk_R18_HeHp_H2()*nmol_eq(jb,4)*nH2b
      h2loss(9) = ioniz_HeI23S_H2(Tb)*nheiTR(jb)*n0*nH2b
      ! (10) the oxygen cycle, O1 and O2 net of their thermodynamic reverses
      h2loss(10)= (ok1b*nox_eq(jb,1) + ok2b*nOb)*nH2b                     &
                - (ok1rb*nox_eq(jb,2) + ok2rb*nox_eq(jb,1))*nHIb
      ! (11) the molecular ions handing H2 back (R6, R9, R11): a source
      h2loss(11)= -( rk_R6_H3p_dr_H2(Tb)*ne_ox(jb)*nmol_eq(jb,3)          &
                   + rk_R9_H2p_H()*nmol_eq(jb,2)*nHIb                      &
                   + rk_R11_H3p_H(Tb)*nmol_eq(jb,3)*nHIb )
      h2tot = sum(h2loss)
      write(74,'(A)') '#'
      write(74,'(A,F9.2,A,ES12.4)') '# NET H2 loss budget at the base'//  &
                     ' cell, T = ', Tb, ' K, net total [cm^-3 s^-1] ',    &
                     h2tot
      write(74,'(A)') '# A negative rate is a net H2 source; the shares'//&
                     ' sum to 1 and may exceed it individually.'
      write(74,'(A)') '# channel                  net rate      share'
      do ib = 1,11
         write(74,'(A,A,ES14.6,F10.4)') '# ', h2ch_name(ib), h2loss(ib),  &
                    h2loss(ib)/sign(max(abs(h2tot), 1.0d-99), h2tot)
      enddo
      close(74)
      write(*,'(a,f8.4,a,f8.4,a,f7.1,a)') ' (write_output/eq) base-cell'//&
           ' NET H2 budget: the oxygen cycle carries ',                   &
           h2loss(10)/sign(max(abs(h2tot), 1.0d-99), h2tot),              &
           ' of it and the thermal channel H2 + M ',                      &
           h2loss(1)/sign(max(abs(h2tot), 1.0d-99), h2tot),               &
           ' (T = ', Tb, ' K)'

      ! ---- what the composition does to the atomic coolants ----
      ! The oxygen and the carbon of a molecular layer are in H2O and CO,
      ! and that is the answer this option exists to give. But the code's
      ! metal cooling is atomic: the [O I] 63/145/44 um fine-structure
      ! lines, the O I and C I/C II lines. Moving those nuclei into
      ! molecules therefore SWITCHES THOSE COOLANTS OFF -- correctly, the
      ! atoms are not there -- and the code has no H2O or CO infrared bands
      ! to put in their place. H3+ is the only molecular coolant it carries.
      !
      ! That is item (G) of TO_BE_DONE.md arriving through the composition
      ! instead of through the temperature, and section 9 of
      ! docs/a2_oxygen_option_design.md says so in advance: A2 is the
      ! composition of the molecular layer and (G) is its energy. Measured
      ! on the hot-Uranus molecular gate, the layer just above the base
      ! cools 15x more slowly with the option on than with it off, at
      ! comparable heating.
      !
      ! It is reported here, at the end of every oxygen-chemistry run,
      ! because a user who does not know it will read the resulting
      ! temperature as a prediction.
      fO_worst = 0.0d0; fC_worst = 0.0d0
      rO_worst = 0.0d0; rC_worst = 0.0d0
      TO_worst = 0.0d0; TC_worst = 0.0d0
      do j = 1,N
         nO_el = nm(j,melem_i0(iel_O))*n0 + nm(j,melem_i0(iel_O)+1)*n0     &
               + nm(j,melem_i0(iel_O)+2)*n0                                &
               + nox_eq(j,1) + nox_eq(j,2) + nox_eq(j,3)
         nC_el = nm(j,melem_i0(iel_C))*n0 + nm(j,melem_i0(iel_C)+1)*n0     &
               + nm(j,melem_i0(iel_C)+2)*n0 + nox_eq(j,3)
         if (nO_el .gt. 0.0d0) then
            fO_mol = (nox_eq(j,1) + nox_eq(j,2) + nox_eq(j,3))/nO_el
            if (fO_mol .gt. fO_worst) then
               fO_worst = fO_mol
               rO_worst = r(j)
               TO_worst = T(j)*T0
            endif
         endif
         if (nC_el .gt. 0.0d0) then
            fC_mol = nox_eq(j,3)/nC_el
            if (fC_mol .gt. fC_worst) then
               fC_worst = fC_mol
               rC_worst = r(j)
               TC_worst = T(j)*T0
            endif
         endif
      enddo
      write(*,'(a,f6.3,a,f8.4,a,f8.1,a)')                                  &
         ' (write_output/eq) oxygen chemistry: at most ', fO_worst,        &
         ' of the oxygen is molecular (r = ', rO_worst, ' Rp, T = ',       &
         TO_worst, ' K)'
      write(*,'(a,f6.3,a,f8.4,a,f8.1,a)')                                  &
         ' (write_output/eq)                   at most ', fC_worst,        &
         ' of the carbon  is molecular (r = ', rC_worst, ' Rp, T = ',      &
         TC_worst, ' K)'
      ! The top-of-domain value, because that is the one a transit forward
      ! model integrates: EXHALE_transit.py reads the O I column, and with
      ! this option on that column is FREE ATOMIC oxygen. A cold outer wind
      ! re-forms CO at the local chemical equilibrium -- equilibrium is not
      ! kinetics, and CO up there would in reality be quenched, not
      ! re-formed -- so a run whose outer wind is cool hands the transit
      ! tool an O I column short by whatever this fraction is.
      nO_el = nm(N,melem_i0(iel_O))*n0 + nm(N,melem_i0(iel_O)+1)*n0        &
            + nm(N,melem_i0(iel_O)+2)*n0                                   &
            + nox_eq(N,1) + nox_eq(N,2) + nox_eq(N,3)
      fO_mol = 0.0d0
      if (nO_el .gt. 0.0d0)                                                &
         fO_mol = (nox_eq(N,1) + nox_eq(N,2) + nox_eq(N,3))/nO_el
      write(*,'(a,f6.3,a,f8.1,a)')                                         &
         ' (write_output/eq)                   at the top of the domain, ',&
         fO_mol, ' of the oxygen is molecular (T = ', T(N)*T0, ' K)'
      if (max(fO_worst, fC_worst) .gt. 0.5d0 .and. .not. mol_ir_bands)      &
         write(*,'(a)') ' (write_output/eq) NOTE those nuclei no longer'// &
            ' cool: the metal cooling is atomic and this run has'//        &
            ' "Molecular IR bands" off, so nothing replaces the [O I]'//   &
            ' and C I/C II lines the molecules switched off and the'//     &
            ' layer they are in is warmer than the same composition'//     &
            ' would really be (TO_BE_DONE.md item (G)).'

      ! ---- FUV band penetration and the G4 energy ledger ----
      ! For a beam attenuated as exp(-tau) the photons absorbed per unit area
      ! over the whole column are exactly N_b (1 - exp(-tau_b(r_min))),
      ! whatever the density profile. The ledger below measures that identity
      ! on the discrete grid: the sum of the absorbed rate over the cells,
      ! with the same opa_pf weighting the columns were built with, against
      ! the closed form. Because each cell's rate is the MEAN over the cell
      ! rather than its face value (water_photolysis_rate), the two agree to
      ! round-off at any grid spacing; a residual here means the rates and
      ! the columns have come apart, not that the grid is coarse.
      !
      ! The energy ledger is then exact by construction: every absorbed
      ! photon carries <E>_b, of which the threshold goes into the bond and
      ! the rest into the gas, so absorbed = heat + bond cell by cell. It is
      ! computed anyway, because a mismatch would mean the heating and the
      ! rates had drifted apart.
      !
      ! THE FIRST BAND HAS THREE ABSORBERS AND ONE BEAM. Over 912-1110 A
      ! the same photons are taken by H2 in the Lyman-Werner lines and by
      ! H2O and OH in a continuum. They share the beam rather than each
      ! attenuating a private copy of it (water_photolysis.f90 sec. 3), so
      ! the closed form the sum is checked against is
      ! N_b (1 - (1 - A) exp(-tau)) with A the line-removed fraction, and
      ! the H2 share belongs INSIDE that sum. The H2 photon count is the
      ! DISSOCIATION rate divided by the DB96 dissociation probability per
      ! pump, because the other 86.5% of the pumps also take a photon out of
      ! the band; of the energy they carry, only the fragment kinetic energy
      ! of the dissociating fraction reaches the gas, so for this band the
      ! 'returned' column is the fluorescence as well as the bond energy.
      absph    = 0.0d0
      absen    = 0.0d0
      heat_col = 0.0d0
      bond_col = 0.0d0
      e_lw_abs = 0.0d0
      do j = 1-Ng,N+Ng
         dr_cm = dr_j(j)*R0*opa_pf(j)
         do ib = 1,n_fuv_band
            absph(ib) = absph(ib)                                          &
                      + (j_h2o_fuv(j,ib)*nox_eq(j,2)                       &
                       + j_oh_fuv(j,ib) *nox_eq(j,1))*dr_cm
            ! Each absorbed photon carries the band's own mean energy
            ! <hv>_b, so the total can never exceed the incident flux; see
            ! the conservation argument in water_photolysis.f90 sec. 2.
            absen(ib) = absen(ib)                                          &
                      + (j_h2o_fuv(j,ib)*nox_eq(j,2)                       &
                       + j_oh_fuv(j,ib)*nox_eq(j,1))                       &
                        *e_photon_flat_band(ib)*dr_cm
            heat_col(ib) = heat_col(ib)                                    &
                      + (j_h2o_fuv(j,ib)*nox_eq(j,2)                       &
                         *heat_per_water_dissociation(ib)                  &
                       + j_oh_fuv(j,ib)*nox_eq(j,1)                        &
                         *heat_per_hydroxyl_dissociation(ib))*dr_cm
         enddo
         if (F_LW_star .gt. 0.0d0) then
            ph_lw    = k_lw_diss(j)/p_diss_lw*nmol_eq(j,1)*dr_cm
            e_lw_abs = e_lw_abs + ph_lw*e_lw_photon_erg
            ! The H2 share of the shared LW beam, on the same ledger as the
            ! continuum absorbers: every pumped photon leaves the beam, and
            ! p_diss of them end in a dissociation that gives the fragment
            ! pair e_lw_fragment_erg.
            absph(ib_LW)    = absph(ib_LW)    + ph_lw
            absen(ib_LW)    = absen(ib_LW)    + ph_lw*e_lw_photon_erg
            heat_col(ib_LW) = heat_col(ib_LW)                              &
                            + ph_lw*p_diss_lw*e_lw_fragment_erg
         endif
      enddo
      bond_col = absen - heat_col

      ! ---- Infrared side of the ledger: the molecular bands ----
      ! The FUV rows above are the photon INPUT. The three molecular bands are
      ! the new energy flow out of (and into) the same layer, so they belong on
      ! the same ledger. Both columns are per unit area of the star-ward
      ! column, with the same opa_pf weighting.
      !
      ! The check that matters is on the ABSORBED column: an optically thin
      ! layer cannot take more out of the incident field than the field
      ! carries, so sum_X absorbed_X must stay below the incident infrared
      ! flux W sigma T0^4 evaluated at the base, where W is largest. Exceeding
      ! it means the layer is no longer thin at these wavelengths and the
      ! closure of molecular_infrared_cooling.f90 has left its validity range.
      ir_emit = 0.0d0
      ir_abs  = 0.0d0
      if (mol_ir_bands) then
         do j = 1-Ng,N+Ng
            dr_cm = dr_j(j)*R0*opa_pf(j)
            if (base_ir_field) then
               w_ir = 0.5d0*base_sky_fraction(r(j),1.0d0)
            else
               w_ir = 0.0d0
            endif
            ir_e = nmol_eq(j,1)*h2_line_emission_lte(T(j)*T0)
            ir_n = h2_line_net_cooling_rate(T(j)*T0, nmol_eq(j,1), w_ir)
            ir_emit(1) = ir_emit(1) + ir_e*dr_cm
            ir_abs(1)  = ir_abs(1)  + (ir_e - ir_n)*dr_cm
            ir_e = nox_eq(j,2)*h2o_band_emission_lte(T(j)*T0)
            ir_n = h2o_band_net_cooling_rate(T(j)*T0, nox_eq(j,2), w_ir)
            ir_emit(2) = ir_emit(2) + ir_e*dr_cm
            ir_abs(2)  = ir_abs(2)  + (ir_e - ir_n)*dr_cm
            ir_e = nox_eq(j,3)*co_band_emission_lte(T(j)*T0)
            ir_n = co_band_net_cooling_rate(T(j)*T0, nox_eq(j,3), w_ir)
            ir_emit(3) = ir_emit(3) + ir_e*dr_cm
            ir_abs(3)  = ir_abs(3)  + (ir_e - ir_n)*dr_cm
         enddo
      endif

      open(unit = 75, file = './output/FUV_bands.txt')
      write(75,'(A)') '# EXHALE schema 2'
      write(75,'(A)') '# FUV band penetration for the oxygen chemistry.'// &
                     ' H2O and OH absorb these bands as continua,'
      write(75,'(A)') '# so the attenuation is exp(-tau) on the star-'//   &
                     'ward columns and there is no self-shielding'
      write(75,'(A)') '# function (unlike the H2 Lyman-Werner band).'
      do ib = 1,n_fuv_band
         write(75,'(A,A,A,ES12.5,A,ES12.5,A,ES12.5)')                      &
            '# band ', fuv_band_name(ib), ': F = ', fuv_band_flux(ib),     &
            ' erg/cm2/s, sigma(H2O) = ', sigma_H2O_band(ib),               &
            ', sigma(OH) = ', sigma_OH_band(ib)
      enddo
      ! The column list is BUILT from the band table rather than written
      ! out, so a change to n_fuv_band cannot leave the header naming a
      ! different number of bands from the one the loop below writes -- as
      ! it did when the bands went from four to five.
      colhdr = '# columns r[Rp] N_H2O[cm^-2] N_OH[cm^-2]'
      do ib = 1,n_fuv_band
         colhdr = trim(colhdr)//' tau_'//trim(fuv_band_name(ib))
      enddo
      do ib = 1,n_fuv_band
         colhdr = trim(colhdr)//' j_H2O_'//trim(fuv_band_name(ib))
      enddo
      do ib = 1,n_fuv_band
         colhdr = trim(colhdr)//' j_OH_'//trim(fuv_band_name(ib))
      enddo
      colhdr = trim(colhdr)//' heat_FUV[erg/cm3/s]'
      write(75,'(A)') trim(colhdr)
      do j = 1-Ng,N+Ng
         write(75,*) r(j), NH2O_col(j), NOH_col(j),                        &
                    (tau_fuv(j,ib), ib = 1,n_fuv_band),                    &
                    (j_h2o_fuv(j,ib), ib = 1,n_fuv_band),                  &
                    (j_oh_fuv(j,ib), ib = 1,n_fuv_band),                   &
                    heat_fuv(j)
      enddo
      write(75,'(A)') '#'
      write(75,'(A)') '# Band energy ledger (gate G4). Per unit area of'// &
                     ' the column:'
      write(75,'(A)') '#   absorbed_ph   photons absorbed, summed over'//  &
                     ' the grid [cm^-2 s^-1]'
      write(75,'(A)') '#   closed_form   N_band (1 - exp(-tau at the'//    &
                     ' innermost cell)), the exact beam result'
      write(75,'(A)') '#   rel_diff      their relative difference'//     &
                     ' (round-off: the cell rate is the cell MEAN)'
      write(75,'(A)') '#   absorbed_en   absorbed photon energy'//         &
                     ' [erg cm^-2 s^-1], and its split into deposited'
      write(75,'(A)') '#                 heat and bond energy; the two'//  &
                     ' must add back to it exactly'
      write(75,'(A)') '#   incident_en   the band flux itself. absorbed'// &
                     '_en must not exceed it.'
      write(75,'(A)') '# band absorbed_ph closed_form rel_diff'//          &
                     ' absorbed_en heat_en bond_en incident_en'
      worst = 0.0d0
      do ib = 1,n_fuv_band
         nph_b    = fuv_band_photon_flux(fuv_band_flux(ib), ib)
         ! Beam transmission to the innermost cell: the continuum for every
         ! band, times the H2 line transmission for the shared LW band.
         expected = nph_b*(1.0d0 - tr_lines_lw(1-Ng)                       &
                                   *exp(-tau_fuv(1-Ng,ib)))
         if (ib .ne. ib_LW)                                                &
            expected = nph_b*(1.0d0 - exp(-tau_fuv(1-Ng,ib)))
         relerr   = abs(absph(ib) - expected)/max(expected, 1.0d-99)
         if (expected .gt. 0.0d0) worst = max(worst, relerr)
         write(75,'(A,A,7ES14.6)') '# ', fuv_band_name(ib), absph(ib),     &
                    expected, relerr, absen(ib), heat_col(ib),             &
                    bond_col(ib), fuv_band_flux(ib)
      enddo
      e_lw_in = fuv_band_flux(ib_LW)
      write(75,'(A)') '#'
      write(75,'(A)') '# The 912-1110 A band is shared: H2 absorbs it in'//&
                     ' the Lyman-Werner lines while H2O and OH'
      write(75,'(A)') '# absorb it as a continuum, out of ONE beam. The'// &
                     ' split of the LW row above is therefore:'
      write(75,'(A,ES14.6)') '# LW_H2_absorbed_en    ', e_lw_abs
      write(75,'(A,ES14.6)') '# LW_cont_absorbed_en  ', absen(ib_LW)       &
                                                        - e_lw_abs
      write(75,'(A,ES14.6)') '# LW_incident_en       ', e_lw_in
      write(75,'(A,ES14.6)') '# LW_line_removed_frac ', a_lines_lw
      write(75,'(A)') '#'
      write(75,'(A)') '# Infrared side of the ledger: the molecular bands'//&
                     ' ("Molecular IR bands"), column-integrated'
      write(75,'(A)') '# per unit area [erg cm^-2 s^-1]. net = emitted -'// &
                     ' absorbed; a negative net is a heated layer.'
      write(75,'(A)') '# The absorbed total is bounded by the incident'//   &
                     ' infrared flux W sigma T0^4 at the base.'
      write(75,'(A,3ES14.6)') '# IR_emitted  H2/H2O/CO ', ir_emit(1),      &
                                              ir_emit(2), ir_emit(3)
      write(75,'(A,3ES14.6)') '# IR_absorbed H2/H2O/CO ', ir_abs(1),       &
                                              ir_abs(2), ir_abs(3)
      write(75,'(A,3ES14.6)') '# IR_net      H2/H2O/CO ',                  &
                 ir_emit(1) - ir_abs(1), ir_emit(2) - ir_abs(2),           &
                 ir_emit(3) - ir_abs(3)
      ir_bound = 0.0d0
      if (mol_ir_bands .and. base_ir_field)                                &
         ir_bound = 0.5d0*base_sky_fraction(r(1),1.0d0)                    &
                    *5.670374419d-5*T0**4
      write(75,'(A,ES14.6)') '# IR_incident_bound     ', ir_bound
      close(75)

      write(*,'(a,es9.2)') ' (write_output/eq) FUV band ledger: max'//     &
           ' |absorbed/closed-form - 1| = ', worst
      if (e_lw_in .gt. 0.0d0 .and. absen(ib_LW) .gt. e_lw_in*(1.0d0+1.0d-6))&
         write(*,'(a)') ' (write_output/eq) WARNING: the absorbers of'//   &
           ' the shared 912-1110 A band take more energy out of it'//      &
           ' than the band carries -- the shared-beam closure has'//       &
           ' broken (see output/FUV_bands.txt).'
      if (mol_ir_bands) then
         write(*,'(a,3es10.3)') ' (write_output/eq) molecular IR bands,'// &
            ' column net H2/H2O/CO [erg cm^-2 s^-1]: ',                    &
            ir_emit(1) - ir_abs(1), ir_emit(2) - ir_abs(2),                &
            ir_emit(3) - ir_abs(3)
         if (ir_bound .gt. 0.0d0 .and.                                     &
             ir_abs(1) + ir_abs(2) + ir_abs(3) .gt. ir_bound)              &
            write(*,'(a,2es10.3)') ' (write_output/eq) WARNING: the'//     &
              ' molecular bands absorb more than the incident infrared'//  &
              ' flux carries -- the optically thin closure has left its'// &
              ' validity range (absorbed, bound): ',                       &
              ir_abs(1) + ir_abs(2) + ir_abs(3), ir_bound
      endif
      if (maxval(NH2_col_lw) .gt. NH2_richings_max)                        &
         write(*,'(a,es9.2,a,es9.2,a,es9.2,a)') ' (write_output/eq)'//     &
           ' WARNING: the star-ward H2 column reaches ',                   &
           maxval(NH2_col_lw), ' cm^-2. The self-shielding fits are'//     &
           ' checked only to ', NH2_richings_max,                          &
           ' cm^-2 (Richings, Schaye & Oppenheimer 2014, which sets the'// &
           ' photodissociation rate) and ', NH2_db96_max,                  &
           ' cm^-2 (Draine & Bertoldi 1996, which sets the band fraction'//&
           ' the H2 lines remove), so the deepest cells are'//             &
           ' extrapolating at least the first of the two.'

      end subroutine write_oxygen_chemistry

      ! End of module
      end module output_write
