      module output_write
      ! Write the output to the standard output files

      use global_parameters
      use utils, only: write_row_layout_header,                        &
                       write_coupling_state_header,                      &
                       coupling_state_fields,                            &
                       write_provenance_header,                          &
                       state_is_certified, state_certification_reason, &
                       set_state_certified
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
      ! THE CONFIGURATION THE STATE IN THIS FILE IS A STATE OF. The block is
      ! built and parsed in one place, the module that reads a restart file,
      ! so the writer and the loader cannot disagree about its fields.
      use IC_load, only: write_restart_metadata_header
      ! The lower boundary model the state was produced under, and the
      ! prescribed reservoir it stands on (write_base_boundary_header).
      use base_boundary, only: base_boundary_model_id,                     &
                               base_reservoir_prescription_version,        &
                               base_reservoir_p, base_reservoir_T,         &
                               base_reservoir_nhat, r_base_level
      ! A molecular seed run writes the pair a restart READS, not the pair a
      ! solution leaves behind, and records the partition it was built with
      ! in both halves (md/input_schema.md appendix D.3).
      use molecular_seed, only: molecular_seed_on,                        &
                                write_molecular_seed_header
      use species_table, only: n_mion, mion_name, im_OI, melem_i0,       &
                               iel_O, iel_C
      use ionization_equilibrium, only: nmol_eq,   &  ! molecular columns
                                       molecular_carrier_densities_from_state, &
                                       NH2_col_lw, f_shield_lw, k_lw_diss, &
                                       p_lw_single,                        &
                                       tr_lines_lw, a_lines_lw, P_H2_eq,   &
                                       lw_col_over_overlap,                &
                                       nox_eq, n_o1d_eq,                   &
                                       NH2O_col, NOH_col,                  &
                                       NCO_col, k_co_diss,                 &
                                       theta_co_shield,                    &
                                       j_h2o_fuv, j_oh_fuv, tau_fuv,       &
                                       heat_fuv
      use utils_ion_eq, only: fuv_band_flux,                             &
                              fuv_band_absorption_ledger,                 &
                              heat_channel_state, ih_H2_LW_dissoc
      ! The CO share of the Lyman-Werner beam's ledger: the mean photon
      ! energy of a CO dissociation event and the fragment kinetic energy
      ! it leaves.  Both are the assembly's own constants, read here for a
      ! COLUMN-INTEGRATED aggregate of one band; the deposit of one cell is
      ! formed only in heating_of_composition.
      use Cooling_Coefficients, only: base_sky_fraction
      use molecular_infrared_cooling, only: h2_line_emission_lte,          &
                            h2o_band_emission_lte, co_band_emission_lte,   &
                            h2_line_net_cooling_rate,                     &
                            h2o_band_net_cooling_rate,                    &
                            co_band_net_cooling_rate
      use water_photolysis, only: n_fuv_band, fuv_band_name,               &
                                       sigma_H2O_band, sigma_OH_band,      &
                                       e_photon_H2O_band, e_photon_OH_band,&
                                       e_photon_flat_band,                 &
                                       fuv_band_photon_flux,               &
                                       ib_LW, ib_B2
      ! Band B2 is the H I Ly-alpha line, so the transmission of its beam is
      ! the fraction of the stellar line that penetrates the atomic hydrogen
      ! column and not exp(-tau).  ONE definition of it, the one the band
      ! rate is built from (util_ion_eq.f90).
      use lya_rt, only: lya_line_center_optical_depth,                     &
                        lya_stellar_beam_transmission
      use oxygen_rates, only: rk_O1_OH_H2_water, rk_O2_O_H2_hydroxyl,     &
                              rate_from_detailed_balance,                 &
                              ith_H, ith_H2, ith_O, ith_OH, ith_H2O
      use mol_rates, only: rk_R12_H2_thdis, rk_R10_Hp_H2v4,               &
                           rk_R13_Hp_H2_M, rk_R14_H2_edis, rk_R8_H2p_H2,  &
                           rk_R17_Hep_H2_diss,                           &
                           rk_R23_H2_Hep_cx, rk_R18_HeHp_H2,              &
                           rk_R15_3body_H2, rk_R6_H3p_dr_H2,              &
                           rk_R9_H2p_H, rk_R11_H3p_H,                     &
                           h2_association_collider_density
      ! The control key of the collider-resolved third body, from the one
      ! accessor that owns it, so that the R15/R12 rate this file reports is
      ! the rate the run was solved at whichever way the key is set.
      use molecular_reaction_heat, only: reaction_heat_recipients_corrected
      ! n(He 1^1S), the neutral helium of the H2 vibrational cascade and of
      ! the R12/R15 third body, from the one routine that owns the floor.
      use composition, only: he_ground_singlet_density
      use Cooling_Coefficients, only: ioniz_HeI23S_H2
      use diffusive_photochemistry, only: carrier_diffusion_coefficient,   &
                                          carrier_transport_diagnostics,   &
                                          ic_H2,                           &
                                          carrier_co_domain_record,        &
                                          carrier_co_domain_f_dom, ic_CO
      use utils, only: calc_ne, calc_ntot
      ! O I ground-term statistical equilibrium: the same solution the
      ! [O I] fine-structure cooling is built on (Cool_coeff.f90).
      use Cooling_Coefficients, only: n_fsline,                        &
                                      fine_structure_line_transfer,    &
                                      oxygen_ground_term_levels

      ! THE VALIDITY OF EACH ROW OF THE ADVECTION-CORRECTED FILES, in two
      ! fields: adv_T_status for the temperature of the row and
      ! adv_comp_status for its composition.  The advection post-process
      ! (post_process_adv) solves the STEADY energy and ionization equations
      ! along the recorded flow, and those are two separate statements about
      ! a cell: a row whose temperature was corrected can carry the run's own
      ! composition and the reverse, so one integer naming the first refusal
      ! cannot describe the row.
      !
      ! ONE definition, used by the module that assigns the values and by the
      ! writer that labels the columns.  The schema is VERSIONED: a file with
      ! no '# adv_schema' line carries the single field of an earlier writer,
      ! and the validity of its rows is UNKNOWN, not corrected.
      integer, parameter :: adv_schema_version = 2
      ! The steady correction was solved and adopted for this row.  It is
      ! a CONDITIONAL correction: accurate to adv_conditional_tol of itself
      ! in the mass flux, the fraction the row's own measure was compared
      ! with (adv_mass_row, below).
      integer, parameter :: adv_corrected = 0
      ! The run's own value was kept.  Three conditions do that, and all
      ! three are statements about the flow rather than about the solve: the
      ! mass row of the cell, measured by the same face-flux operator the
      ! stationary certification uses, stands above adv_conditional_tol (so
      ! the steady equations of the correction drop terms as large as that
      ! fraction of the ones they keep), the local radiative balance rather
      ! than the flow sets the temperature (thermal Damkohler number above
      ! one), or the gas enters the cell and the step of the recursion has no
      ! upstream state.
      integer, parameter :: adv_retained = 1
      ! The correction was attempted and its solve did not converge, or its
      ! root left the range a state of this gas can occupy.
      integer, parameter :: adv_failed = 2
      ! The closure does not cover this row.  The post-process reconstructs
      ! an H/He + trace-metal gas and omits the molecular (H2 H2+ H3+ HeH+)
      ! and oxygen (OH H2O CO) carriers from its particle and electron
      ! counts, so a cell in which those omitted species carry more of the
      ! particle count than the species it does carry is a cell class it does
      ! not model.  For the TEMPERATURE the same holds of the energy balance:
      ! a row whose heating and cooling in the run are carried more by the
      ! omitted molecular and oxygen channels than by the channels the
      ! post-process assembles keeps the run's own temperature with this
      ! value (the energy class of post_process_adv).
      integer, parameter :: adv_unsupported = 3
      ! The post-process did not reach this row: a ghost row below the first
      ! cell at which the step of its equations can be taken.
      integer, parameter :: adv_not_evaluated = 4

      ! THE CONDITION UNDER WHICH A ROW IS A CORRECTION, and the accuracy
      ! that condition buys.  The advective correction of a cell is a steady
      ! integral along the recorded flow, and the flow it integrates along
      ! carries a fractional change of the face mass flux across the cell,
      ! m_j = |R_1(j)|/s_1(j) (the adv_mass_row column, and the mass row of
      ! the stationary certification).  The correction is first order in
      ! m_j: the terms it drops are the ones the mass divergence puts into
      ! the steady equations, each of them m_j times a term the equations
      ! keep.  So a row whose m_j is at or below this fraction is a
      ! CONDITIONAL correction accurate to that fraction of itself, and a
      ! row above it is not corrected at all and keeps the run's own value.
      !
      ! 1e-2 is the fraction, so a corrected row is accurate to about one
      ! percent of itself in the mass flux.  It is not the tolerance of a
      ! certified stationary state: that is cert_tol_mass, the number the
      ! stationary certification judges the whole state by, reported
      ! separately in the '# adv_input_certified' line.  The two are
      ! different statements about the same measure and the file carries
      ! both.
      real*8, parameter :: adv_conditional_tol = 1.0d-2

      ! The C library rename(3), which replaces the target name in one step
      ! within a file system (POSIX), so a reader of the target sees the old
      ! file or the new one and never a partial one. Standard Fortran has no
      ! rename; this interface is the same for every compiler the Makefile
      ! names.
      interface
         function posix_rename(old_path, new_path) bind(C, name='rename')  &
                  result(rc)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: old_path
            character(kind=c_char), dimension(*), intent(in) :: new_path
            integer(c_int) :: rc
         end function posix_rename
         ! mkdir(2) and remove(3) for the directories of the pass state
         ! (publish_pass_state_generation). mode_t is a 32-bit unsigned
         ! integer on every system the Makefile builds for, passed by value.
         function posix_mkdir(dir_path, dir_mode) bind(C, name='mkdir')    &
                  result(rc)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: dir_path
            integer(c_int), value :: dir_mode
            integer(c_int) :: rc
         end function posix_mkdir
         function posix_remove(entry_path) bind(C, name='remove')          &
                  result(rc)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: entry_path
            integer(c_int) :: rc
         end function posix_remove
      end interface
      private :: posix_mkdir, posix_remove
      private :: pass_state_safe_name, pass_state_generation_name_valid,  &
                 pass_state_current_generation,                           &
                 pass_state_previous_generation,                          &
                 remove_pass_state_directory, stop_if_pass_state_step

      ! THE IDENTITY OF THIS RUN, unique across runs and hosts: the wall
      ! clock of the run start to the millisecond and the count of the
      ! processor clock at that moment, 'yyyymmddThhmmss.sss_c<count>'.
      ! Set once (establish_run_identity) and read by the pass state, whose
      ! generations are named by it, so that two runs writing the same pass
      ! number into one directory cannot publish two states under one name.
      character(len=48), save :: run_identity = ''

      ! THE PASS STATE of the stationary outer iteration is published in
      ! generations under this directory (publish_pass_state_generation):
      ! one subdirectory for each pass, and the one-line file 'current'
      ! naming the generation a restart takes.
      character(len=*), parameter, private :: pass_state_root =           &
                                              './output/pass_state'

      contains

      ! ------------------------------------------------------------------ !

      subroutine establish_run_identity()
      ! Form run_identity from date_and_time (to the millisecond) and the
      ! 64-bit system_clock count; both are standard Fortran. The clock
      ! count differs between hosts started at the same millisecond, the
      ! date between runs of one host. Called once at the start of the run;
      ! a second call keeps the first identity.
      integer :: dt_values(8)
      integer(kind=8) :: clock_count
      if (len_trim(run_identity) .gt. 0) return
      call date_and_time(values = dt_values)
      call system_clock(count = clock_count)
      write(run_identity,'(I4.4,2I2.2,A,3I2.2,A,I3.3,A,I0)')              &
           dt_values(1), dt_values(2), dt_values(3), 'T', dt_values(5),   &
           dt_values(6), dt_values(7), '.', dt_values(8), '_c', clock_count
      run_identity = pass_state_safe_name(run_identity)
      end subroutine establish_run_identity

      ! ------------------------------------------------------------------ !

      function pass_state_safe_name(raw) result(safe)
      ! The name with every character outside [A-Za-z0-9._-] replaced by
      ! '_', so that it can stand as one path component and as one token of
      ! a header line.
      character(len=*), intent(in) :: raw
      character(len=len(raw)) :: safe
      integer :: ic
      character(len=1) :: ch
      safe = raw
      do ic = 1, len_trim(raw)
         ch = raw(ic:ic)
         if (.not. ((ch .ge. 'A' .and. ch .le. 'Z') .or.                  &
                    (ch .ge. 'a' .and. ch .le. 'z') .or.                  &
                    (ch .ge. '0' .and. ch .le. '9') .or.                  &
                    ch .eq. '.' .or. ch .eq. '_' .or. ch .eq. '-'))       &
            safe(ic:ic) = '_'
      enddo
      end function pass_state_safe_name

      subroutine write_derived_provenance_header(unit)
      ! WHOSE CERTIFICATE THE LINE BESIDE A DERIVED PRODUCT IS.
      !
      ! The rows of the advection-corrected files are a composition computed
      ! FROM the state written beside them, and the certification inventory
      ! was evaluated on that state and on no other.  A 'certified=' field
      ! here would attach the verdict of one composition to another one, so
      ! this file carries no coupling header at all; what it carries is the
      ! fields of the state it was derived from -- the same fields as the
      ! '# coupling:' line of the solved Hydro_ioniz.txt, with the
      ! certification pair first (coupling_state_fields: certified,
      ! cert_reason, sec_ion, sec_ion_step, recon, iontrans, mode,
      ! t_phys) -- under a key that says it is
      ! provenance, with the sentence that fixes what the key means.  Both
      ! _adv files carry it (Hydro_ioniz_adv.txt and Ion_species_adv.txt).
      !
      ! The certification pair is stated in words in the adv_input_certified
      ! line of write_adv_validity_header, which the transit tool reads; this
      ! line is the machine-readable form of it, and the two are one
      ! statement.
      integer, intent(in) :: unit
      write(unit,'(A)') '# derived_from: '//                            &
                     trim(coupling_state_fields(certification_first=.true.))
      write(unit,'(A)') '# derived_from is PROVENANCE: the physics and'//   &
                     ' certification fields of the state these rows'
      write(unit,'(A)') '#   were derived from, which is the state written'//&
                     ' beside this file. It is a statement about that state'
      write(unit,'(A)') '#   and not about these rows: the'//               &
                     ' advection-corrected composition is another'
      write(unit,'(A)') '#   composition, and no entry of the'//            &
                     ' certification inventory was evaluated on it.'
      write(unit,'(A)') '#   This file carries no coupling header and'//    &
                     ' makes no certification claim of its own.'
      end subroutine write_derived_provenance_header

      subroutine write_base_boundary_header(unit)
      ! THE LOWER BOUNDARY THIS STATE WAS PRODUCED UNDER, and the reservoir
      ! that boundary was prescribed with.
      !
      ! WHY THE RESERVOIR TRAVELS WITH THE STATE AND THE GHOST ROWS DO NOT.
      ! The boundary is a function of the physical column, the prescribed
      ! reservoir, the radiation context and the model options; of those the
      ! reservoir is the only one the column cannot reconstruct, so it is
      ! what a file has to carry. The ghost rows this file writes are the
      ! boundary's own output and a reader rebuilds them from the column and
      ! the reservoir instead of reading them (load_IC states the same
      ! contract at the read).
      !
      ! THE MEANING AND THE VERSION ARE ON THE LINE: p and T are the pressure
      ! and temperature at the level, in the code's own units; nhat is the
      ! particle count per unit mass of the reservoir gas; r_level is the
      ! radius the three are stated at; and the entropy is not a number but
      ! the isentrope through (p, T) at the base composition. The version
      ! rises when any of those meanings changes.
      !
      ! Both are '#' comments, so no numeric parse and no golden changes.
      integer, intent(in) :: unit
      write(unit,'(A)') '# boundary_model '//base_boundary_model_id()
      ! ES24.16E3 round trips a double exactly, so the four numbers on the
      ! line are the numbers the run held and not a rendering of them (24,
      ! not 23: the field must hold the sign, W >= D + E + 5).
      write(unit,'(A,I0,A,4(1X,ES24.16E3))')                              &
           '# boundary_reservoir version ',                               &
           base_reservoir_prescription_version,                           &
           ' p[p0] T[T0] nhat[n0/rho0] r_level[Rp]',                      &
           base_reservoir_p, base_reservoir_T, base_reservoir_nhat,       &
           r_base_level
      end subroutine write_base_boundary_header

      subroutine write_output(rho,v,p,T,heat,cool,eta,                &
                              nhi,nhii,nhei,nheii,nheiii,nheiTR,      &
                              nm,flag,adv_T_status,adv_comp_status,        &
                              adv_mass_row)
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
      ! Present only for the advection-corrected write, where the pair
      ! carries the validity of each row's temperature and of its
      ! composition (the parameters above). Absent for the equilibrium
      ! write, whose rows are the run's own solution by construction, and
      ! then no such column is written.
      integer, dimension(1-Ng:N+Ng), intent(in), optional :: adv_T_status
      integer, dimension(1-Ng:N+Ng), intent(in), optional :: adv_comp_status
      ! The measure the two fields were decided by, row by row: the
      ! fractional change of the face mass flux across the cell, so that a
      ! reader has the CONDITION of every row and not only its verdict.
      ! One group with the two fields above: the advection-corrected write
      ! passes all three.
      real*8, dimension(1-Ng:N+Ng), intent(in), optional :: adv_mass_row
      ! That measure as written, zero where the post-process did not form
      ! it (a row it never reached).
      real*8, dimension(1-Ng:N+Ng) :: mrow


      character(len=40) :: hyd_path, ion_path

      mrow = 0.0d0
      if (present(adv_mass_row)) mrow = adv_mass_row

      !---- The state pair: thermodynamic and ionization profiles ----!
      if (flag.eq.'eq') then
         if (molecular_seed_on()) then
            ! The product of a seed run IS the restart input: it is written
            ! under the name a restart reads, so that nothing has to be
            ! renamed between the conversion and the run it seeds.
            hyd_path = './output/Hydro_ioniz_IC.txt'
            ion_path = './output/Ion_species_IC.txt'
         else
            hyd_path = './output/Hydro_ioniz.txt'
            ion_path = './output/Ion_species.txt'
         endif
      else	! Change output file after postprocessing
         hyd_path = './output/Hydro_ioniz_adv.txt'
         ion_path = './output/Ion_species_adv.txt'
      endif
      call write_hydro_state_file(trim(hyd_path), flag, rho, v, p, T,     &
                                  heat, cool, mrow, adv_T_status,         &
                                  adv_comp_status)
      call write_species_state_file(trim(ion_path), flag, nhi, nhii,      &
                                    nhei, nheii, nheiii, nheiTR, nm,      &
                                    adv_T_status, adv_comp_status)

      !---- Lyman-Werner photodissociation diagnostic ----!
      ! Written only for a molecular run that carries a Lyman-Werner band
      ! flux, and only for the equilibrium state (the advection
      ! post-process is atomic and does not re-solve the molecules).
      ! It is the record of how deep the band penetrates: f_shield -> 1 in
      ! the thin wind above the H2 -> H front and collapses in the
      ! self-shielded molecular base.  The f_shield column is a DIAGNOSTIC:
      ! h2_self_shielding_level_resolved, the ratio
      ! sigma_diss(N)/sigma_diss(N_bottom) of the level-resolved
      ! overlapping-line table of h2_self_shielding_table at the cell's own
      ! star-ward H2 column (lyman_werner.f90 sec. 2).  The k_LW column next
      ! to it is not this factor times anything: it is the incident band
      ! photon fluence contracted with the same tabulated sigma_diss and the
      ! 912-1201 A continuum, averaged over the cell between its two faces
      ! (lyman_werner_dissociation_rate_cell_mean).  The SHARE of the band
      ! the H2 lines take out of the shared FUV beam is the column integral
      ! of the pump cross section of the SAME table
      ! (h2_lw_band_photon_fraction_absorbed, lyman_werner.f90 sec. 2g); it
      ! matters only when the oxygen chemistry shares the band, and is
      ! reported in output/FUV_bands.txt when it does.
      if (thereis_mol .and. F_LW_star .gt. 0.0d0 .and. flag .eq. 'eq') then
         open(unit = 4, file = './output/Lyman_Werner.txt')
         write(4,'(A)') '# EXHALE schema 2'
         write(4,'(A,ES12.5,A)') '# Lyman-Werner band flux at the planet: ', &
                        F_LW_star, ' erg cm^-2 s^-1 (912-1201 A)'
         write(4,'(A,F6.3,A,ES12.5)') '# dayside dilution applied: ',      &
                        dayside_dilution(), '  -> beam driving k_LW: ',    &
                        fuv_band_flux(ib_LW)
         write(4,'(A)') '# columns r[Rp] T[K] x_H2[2nH2/nH] nH2[cm^-3] '//  &
                        'NH2_star[cm^-2] f_shield k_LW[1/s] '//             &
                        'heat_LW[erg/cm3/s]'
         write(4,'(A)') '# heat_LW is the Lyman-Werner fragment channel'//  &
                        ' of the heating assembly as the sweep deposited'
         write(4,'(A)') '# it, not a second evaluation of it: one energy,'//&
                        ' one place it is formed.'
         call write_row_layout_header(4)
         do j = 1-Ng,N+Ng
            nh_lw = (nhi(j) + nhii(j))*n0                                  &
                  + 2.0d0*(nmol_eq(j,1) + nmol_eq(j,2))                    &
                  + 3.0d0*nmol_eq(j,3) + nmol_eq(j,4)
            write(4,*) r(j), T(j)*T0,                                      &
                       2.0d0*nmol_eq(j,1)/max(nh_lw,1.0d-99),              &
                       nmol_eq(j,1), NH2_col_lw(j), f_shield_lw(j),        &
                       k_lw_diss(j),                                       &
                       heat_channel_state(j,ih_H2_LW_dissoc)
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
         call write_row_layout_header(73)
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

      ! ------------------------------------------------------------------ !

      subroutine write_hydro_state_file(path, flag, rho, v, p, T, heat,   &
                                        cool, mrow, adv_T_status,         &
                                        adv_comp_status, pass_statement,  &
                                        ios_open, state_id)
      ! THE THERMODYNAMIC HALF OF A STATE, to the file named: r, rho, v, p,
      ! T, heat and cool of every row, ghosts included, under the header a
      ! restart reads (the coupling line, provenance, boundary reservoir,
      ! molecular seed statement and the restart metadata block).  One
      ! writer for the solved state, the advection-corrected profiles and
      ! the pass snapshot of the stationary outer iteration, so the three
      ! cannot state the header in different ways.
      !
      ! pass_statement, present only for the pass snapshot, is written as
      ! its own '#' line after the coupling line.  ios_open, when present,
      ! receives the status of the open instead of the run stopping on a
      ! file that cannot be opened: a snapshot that cannot be written must
      ! not end the run that is producing it.  state_id, present only for
      ! the pass snapshot, is the generation the file belongs to, stated as
      ! the last field of the coupling line of both halves
      ! (publish_pass_state_generation); the solved state carries none.
      character(len=*), intent(in) :: path
      character(len=2), intent(in) :: flag
      real*8, dimension(1-Ng:N+Ng), intent(in) :: rho, v, p, T, heat, cool
      real*8, dimension(1-Ng:N+Ng), intent(in) :: mrow
      integer, dimension(1-Ng:N+Ng), intent(in), optional :: adv_T_status
      integer, dimension(1-Ng:N+Ng), intent(in), optional :: adv_comp_status
      character(len=*), intent(in), optional :: pass_statement
      integer, intent(out), optional :: ios_open
      character(len=*), intent(in), optional :: state_id
      integer :: unit_h, j

      if (present(ios_open)) then
         open(newunit = unit_h, file = path, status = 'replace',          &
              action = 'write', iostat = ios_open)
         if (ios_open .ne. 0) return
      else
         open(newunit = unit_h, file = path)
      endif

      ! Schema header ('#' comment lines; readers that predate the header
      ! can skip them, numeric content is unchanged)
      ! Column 2 is rho*n0: the MASS density in units of m_H per cm^3 (metals
      ! included under the eos_metals policy), not a number density. Multiply by
      ! m_H to get g/cm^3.
      write(unit_h,'(A)') '# EXHALE schema 2'
      if (present(adv_T_status)) then
         write(unit_h,'(A)') '# columns r[Rp] rho[mH/cm3] v[cm/s] p[cgs] '//     &
                        'T[K] heat[erg/cm3/s] cool[erg/cm3/s] '//           &
                        'adv_T_status adv_comp_status adv_mass_row'
      else
         write(unit_h,'(A)') '# columns r[Rp] rho[mH/cm3] v[cm/s] p[cgs] '//     &
                        'T[K] heat[erg/cm3/s] cool[erg/cm3/s]'
      endif
      ! Which rows are physical. The loop below writes 1-Ng..N+Ng, so the
      ! first Ng and the last Ng rows are GHOST cells, filled by Apply_BC
      ! from the interior (lower: fixed base state; upper: zero-gradient /
      ! WENO3 extrapolation). Their rho*v*r^2 is a boundary extrapolation,
      ! not part of the solution, and including them in a flux-spread or
      ! residual measure doubles it: the accepted flux spread of the
      ! HD 189733 b solve is 4.64e-3 over the physical cells and 1.05e-2 if
      ! the two upper ghost rows are counted (section 133.6). The line is a
      ! '#' comment, so no reader's numeric parse and no golden changes.

         call write_row_layout_header(unit_h)
         if (flag .eq. 'eq') then
            ! What the state in this file was produced under (see the
            ! routine). Hydro_ioniz.txt is what a restart is fed as
            ! Hydro_ioniz_IC.txt, so this is the file the line has to
            ! travel in.
            call write_coupling_state_header(unit_h, state_id)
            ! The pass snapshot says what it is on a line of its own,
            ! beside the certification pair it carries.
            if (present(pass_statement))                               &
               write(unit_h,'(A)') '# '//trim(pass_statement)
         else
            call write_derived_provenance_header(unit_h)
         endif
         call write_provenance_header(unit_h)
         call write_base_boundary_header(unit_h)
         call write_molecular_seed_header(unit_h)
         ! The reservoir, species schema, physical grid, constants, options
         ! and clock the state was produced under: what a restart of this
         ! file is compared against (see the routine).
         call write_restart_metadata_header(unit_h)
         ! What the two status columns mean, what the stationarity of the
         ! input was judged by, what the product is, and how many rows carry
         ! each value (see the routine).
         if (present(adv_T_status))                                         &
            call write_adv_validity_header(unit_h, adv_T_status, adv_comp_status)
         ! Real columns in ES25.17E3, which reads back to the same binary64
         ! value (see write_species_state_file); the advection rows carry two
         ! integer status columns before the measure.
         do j = 1-Ng,N+Ng
            if (present(adv_T_status)) then
               write(unit_h,'(7(1X,ES25.17E3),2(1X,I0),1X,ES25.17E3)')     &
                        r(j),        &     ! Rad. dist.
                        rho(j)*n0,     &     ! Density
                        v(j)*v0,       &     ! Velocity
                        p(j)*p0,       &     ! Pressure
                        T(j)*T0,       &     ! Temperature
                        heat(j)*q0,    &     ! Rad. heat.
                        cool(j)*q0,    &     ! Rad. cool.
                        adv_T_status(j),  &  ! Validity of its temperature
                        adv_comp_status(j),& ! Validity of its composition
                        mrow(j)              ! The measure both were decided by
            else
               write(unit_h,'(7(1X,ES25.17E3))')                          &
                        r(j),        &     ! Rad. dist.
                        rho(j)*n0,     &     ! Density
                        v(j)*v0,       &     ! Velocity
                        p(j)*p0,       &     ! Pressure
                        T(j)*T0,       &     ! Temperature
                        heat(j)*q0,    &     ! Rad. heat.
                        cool(j)*q0           ! Rad. cool.
            endif
         enddo
      close(unit_h)
      end subroutine write_hydro_state_file

      ! ------------------------------------------------------------------ !

      subroutine write_species_state_file(path, flag, nhi, nhii, nhei,     &
                                          nheii, nheiii, nheiTR, nm,       &
                                          adv_T_status, adv_comp_status,   &
                                          pass_statement, ios_open,        &
                                          state_id)
      ! THE COMPOSITION HALF OF A STATE, to the file named: every species
      ! density of every row, ghosts included. The molecular and oxygen
      ! columns are read from nmol_eq and nox_eq, so a caller refreshes them
      ! from the state it writes first (molecular_carrier_densities_from_state).
      ! pass_statement, ios_open and state_id as in write_hydro_state_file.
      character(len=*), intent(in) :: path
      character(len=2), intent(in) :: flag
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi, nhii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei, nheii, nheiii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nheiTR
      real*8, dimension(1-Ng:N+Ng,n_mion), intent(in) :: nm
      integer, dimension(1-Ng:N+Ng), intent(in), optional :: adv_T_status
      integer, dimension(1-Ng:N+Ng), intent(in), optional :: adv_comp_status
      character(len=*), intent(in), optional :: pass_statement
      integer, intent(out), optional :: ios_open
      character(len=*), intent(in), optional :: state_id
      integer :: unit_i, j, i, n_real
      ! Every real column in ES25.17E3, which reads back to the same binary64
      ! value (17 significant digits); list-directed output leaves the digit
      ! count to the runtime, and the ifx one writes 15, which moves a
      ! reloaded value by up to 12 ulp (measured 2026-09-29) -- more than the
      ! one-ulp scale of the certification rounding floors. The advection
      ! rows end in the two integer status columns.
      character(len=64) :: row_fmt, row_fmt_adv

      ! The column count of the branches below: the oxygen rows carry the
      ! four molecular columns and the three oxygen ones.
      n_real = 7 + n_mion
      if (thereis_oxychem) then
         n_real = n_real + 7
      else if (thereis_mol) then
         n_real = n_real + 4
      endif
      write(row_fmt,'(A,I0,A)')     '(', n_real, '(1X,ES25.17E3))'
      write(row_fmt_adv,'(A,I0,A)') '(', n_real, '(1X,ES25.17E3),2(1X,I0))'

      if (present(ios_open)) then
         open(newunit = unit_i, file = path, status = 'replace',          &
              action = 'write', iostat = ios_open)
         if (ios_open .ne. 0) return
      else
         open(newunit = unit_i, file = path)
      endif

      ! Schema header: species labels in column order, generated from the
      ! species table so they stay correct when species are added.
      write(unit_i,'(A)') '# EXHALE schema 2'
      ! The advection post-process reconstructs an H/He + trace-metal gas and
      ! does not carry the molecular or oxygen species: their columns below
      ! are the EQUILIBRIUM values, written beside advection-corrected H, He
      ! and metal columns. Say so in the file rather than only in the module
      ! that produces it, so that a reader of Ion_species_adv.txt cannot take
      ! them for corrected profiles. Lifting the limitation is separate work.
      if (flag .ne. 'eq' .and. (thereis_mol .or. thereis_oxychem)) then
         write(unit_i,'(A)') '# NOTE the molecular columns (H2 H2p H3p HeHp)'// &
                        ' and, when present, the oxygen columns'
         write(unit_i,'(A)') '#   (OH H2O CO) are NOT advection-corrected:'//   &
                        ' the post-process reconstructs an'
         write(unit_i,'(A)') '#   H/He + trace-metal gas, so those columns'//   &
                        ' are the equilibrium solution and the'
         write(unit_i,'(A)') '#   metal columns inside the molecular layer'//   &
                        ' inherit that approximation.'
      endif
      write(unit_i,'(A)', advance='no') '# columns r[Rp] HI HII HeI HeII '// &
                                   'HeIII HeITR'
      do i = 1,n_mion
         write(unit_i,'(A)', advance='no') ' '//trim(mion_name(i))
      enddo
      ! molecular columns (present only when thereis_mol)
      if (thereis_mol) write(unit_i,'(A)', advance='no') ' H2 H2p H3p HeHp'
      ! oxygen-carrier columns (present only when thereis_oxychem)
      if (thereis_oxychem) write(unit_i,'(A)', advance='no') ' OH H2O CO'
      ! The validity of each row travels with the species columns as well,
      ! so a reader of this file alone can tell a corrected composition from
      ! the run's own (see write_adv_validity_header).
      if (present(adv_T_status)) write(unit_i,'(A)', advance='no')            &
                                   ' adv_T_status adv_comp_status'
      write(unit_i,'(A)') ''

      call write_row_layout_header(unit_i)
      ! THE PASS SNAPSHOT STATES ITS CLAIM IN BOTH HALVES. The solved state
      ! states it on Hydro_ioniz alone; the snapshot writes the coupling
      ! line here as well, with the state_id of its generation, so that a
      ! pair assembled from two generations (or from a generation and a
      ! solved state) states two identities and is refused by load_IC
      ! (adopt_certification_claim) instead of being loaded as one state.
      if (present(pass_statement)) then
         call write_coupling_state_header(unit_i, state_id)
         write(unit_i,'(A)') '# '//trim(pass_statement)
      endif
      ! The derived product states the fields of the state it was derived
      ! from, as Hydro_ioniz_adv.txt does (see the routine).
      if (flag .ne. 'eq') call write_derived_provenance_header(unit_i)
      call write_base_boundary_header(unit_i)
      call write_molecular_seed_header(unit_i)
      ! Both state files carry the block: they are two halves of one state,
      ! and a restart reads both, so a pair whose halves state different
      ! configurations is refused rather than half-loaded.
      call write_restart_metadata_header(unit_i)
      if (present(adv_T_status))                                         &
         call write_adv_validity_header(unit_i, adv_T_status, adv_comp_status)
      do j = 1-Ng,N+Ng

         if (present(adv_T_status)) then
            ! The advection-corrected write. Its own statements, kept apart
            ! from the equilibrium ones below so that the equilibrium file
            ! is written by the statements it always was.
            if (thereis_oxychem) then
               write(unit_i,row_fmt_adv) r(j), nhi(j)*n0, nhii(j)*n0, nhei(j)*n0,       &
                        nheii(j)*n0, nheiii(j)*n0, nheiTR(j)*n0,         &
                        (nm(j,i)*n0, i = 1,n_mion),                      &
                        (nmol_eq(j,i), i = 1,4),                         &
                        (nox_eq(j,i), i = 1,3),                          &
                        adv_T_status(j), adv_comp_status(j)
            else if (thereis_mol) then
               write(unit_i,row_fmt_adv) r(j), nhi(j)*n0, nhii(j)*n0, nhei(j)*n0,       &
                        nheii(j)*n0, nheiii(j)*n0, nheiTR(j)*n0,         &
                        (nm(j,i)*n0, i = 1,n_mion),                      &
                        (nmol_eq(j,i), i = 1,4),                         &
                        adv_T_status(j), adv_comp_status(j)
            else
               write(unit_i,row_fmt_adv) r(j), nhi(j)*n0, nhii(j)*n0, nhei(j)*n0,       &
                        nheii(j)*n0, nheiii(j)*n0, nheiTR(j)*n0,         &
                        (nm(j,i)*n0, i = 1,n_mion),                      &
                        adv_T_status(j), adv_comp_status(j)
            endif
         else if (thereis_oxychem) then
            write(unit_i,row_fmt) r(j), nhi(j)*n0, nhii(j)*n0, nhei(j)*n0,        &
                     nheii(j)*n0, nheiii(j)*n0, nheiTR(j)*n0,          &
                     (nm(j,i)*n0, i = 1,n_mion),                       &
                     (nmol_eq(j,i), i = 1,4),                          &
                     (nox_eq(j,i), i = 1,3)   ! OH H2O CO (already cm^-3)
         else if (thereis_mol) then
            write(unit_i,row_fmt) r(j), nhi(j)*n0, nhii(j)*n0, nhei(j)*n0,        &
                     nheii(j)*n0, nheiii(j)*n0, nheiTR(j)*n0,          &
                     (nm(j,i)*n0, i = 1,n_mion),                       &
                     (nmol_eq(j,i), i = 1,4)  ! H2 H2+ H3+ HeH+ (already cm^-3)
         else
         write(unit_i,row_fmt) r(j),       & ! Rad. dist.
                  nhi(j)*n0,    & ! HI
                  nhii(j)*n0,   & ! HII
                  nhei(j)*n0,   & ! HeI
                  nheii(j)*n0,  & ! HeII
                  nheiii(j)*n0, & ! HeIII
                  nheiTR(j)*n0, & ! HeITR
                  (nm(j,i)*n0, i = 1,n_mion)  ! metal ions (canonical order)
         endif
      enddo
      close(unit_i)
      end subroutine write_species_state_file

      ! ------------------------------------------------------------------ !

      subroutine publish_pass_state_generation(pass_now, pass_statement, &
                                       rho, v, p, T, heat, cool, nhi,     &
                                       nhii, nhei, nheii, nheiii, nheiTR, &
                                       nm, f_sp, written, published_id)
      ! THE STATE ONE OUTER PASS OF THE STATIONARY ITERATION HANDS TO THE
      ! NEXT, published as ONE GENERATION:
      !
      !    output/pass_state/<gen>/Hydro_ioniz.txt
      !    output/pass_state/<gen>/Ion_species.txt
      !    output/pass_state/<gen>/manifest.txt
      !    output/pass_state/current          (one line: <gen>)
      !
      ! <gen> = r<run_identity>_p<pass>, so no two runs and no two passes
      ! name one generation. A run stopped from outside (a wall-clock
      ! timeout, a kill) leaves the state of its last published pass behind
      ! instead of nothing. Copied to Hydro_ioniz_IC.txt / Ion_species_IC.txt
      ! the pair of the generation 'current' names is a "Load IC? True"
      ! restart like any other state file: the writers, header lines and
      ! restart tokens are the ones of the solved state
      ! (write_hydro_state_file, write_species_state_file).
      !
      ! IT IS NEVER A CERTIFIED STATE. The pair carries certified=F
      ! cert_reason=pass_snapshot_p<pass>: the certification of a pass is
      ! taken before its composition update, and the state written here is
      ! the state after it, on which no certification was made. The pair
      ! the run holds for its own final state is put back afterwards, so the
      ! files the run writes at its end are not touched. Both halves carry
      ! the coupling line (the solved state writes it on Hydro_ioniz alone),
      ! ending in state_id=<gen>, so a pair assembled from two generations,
      ! or from a generation and a solved state, is refused by load_IC
      ! (adopt_certification_claim).
      !
      ! THE ORDER OF PUBLICATION, so that an interruption at any point
      ! leaves 'current' naming a complete generation (or, before the first
      ! publication of a directory, no 'current' at all):
      !
      !   1. the two state files and the manifest are written and closed in
      !      output/pass_state/.<gen>.tmp/;
      !   2. that directory is renamed onto output/pass_state/<gen> (one
      !      rename(2));
      !   3. output/pass_state/current.part, holding <gen>, is written,
      !      closed and renamed onto output/pass_state/current (one
      !      rename(2): the single decision that publishes the generation);
      !   4. the generation 'current' named before step 3 is kept, and the
      !      one IT replaced (the previous_state_id of its manifest) is
      !      removed, file by file under its explicit name.
      !
      ! An interruption before step 2 leaves a .<gen>.tmp directory, one
      ! between 2 and 3 a complete generation no pointer names; neither is
      ! read by anything, and neither is removed by a later run (only names
      ! read from a manifest are removed). EXHALE_PASS_STATE_STOP_AT=
      ! <step>[:<pass>] stops the run (stop 3) right after step 1, 2 or 3
      ! of the pass named (of every pass if none is), so what an
      ! interruption at each point leaves can be measured.
      !
      ! nmol_eq and nox_eq, which the molecular and oxygen columns are read
      ! from, are refreshed from (rho, f_sp) for the write and put back
      ! after it, so nothing this routine does is read by the run later.
      integer, intent(in)          :: pass_now
      character(len=*), intent(in) :: pass_statement
      real*8, dimension(1-Ng:N+Ng), intent(in) :: rho, v, p, T, heat, cool
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi, nhii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei, nheii, nheiii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nheiTR
      real*8, dimension(1-Ng:N+Ng,n_mion), intent(in) :: nm
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      logical, intent(out) :: written
      character(len=*), intent(out) :: published_id
      character(len=*), parameter :: hyd_leaf = 'Hydro_ioniz.txt'
      character(len=*), parameter :: ion_leaf = 'Ion_species.txt'
      character(len=*), parameter :: manifest_leaf = 'manifest.txt'
      ! 0755, the mode of a directory the owner writes and others read.
      integer(c_int), parameter :: dir_mode_0755 = 493_c_int
      character(len=96)  :: gen_name, gen_before, gen_retired
      character(len=192) :: tmp_dir, gen_dir, ptr_file
      character(len=len(state_certification_reason)) :: hold_reason
      character(len=32) :: snap_reason
      logical :: hold_cert
      real*8, allocatable :: nmol_hold(:,:), nox_hold(:,:)
      real*8, dimension(1-Ng:N+Ng) :: mrow_none
      integer :: ios_h, ios_i, ios_m, unit_m
      integer(c_int) :: rc_dir, rc_gen, rc_ptr
      integer(kind=8) :: bytes_h, bytes_i

      written = .false.
      published_id = ''
      if (len_trim(run_identity) .eq. 0) call establish_run_identity()
      write(gen_name,'(A,A,I0)') 'r'//trim(run_identity), '_p', pass_now
      gen_name = pass_state_safe_name(gen_name)
      tmp_dir  = pass_state_root//'/.'//trim(gen_name)//'.tmp'
      gen_dir  = pass_state_root//'/'//trim(gen_name)
      ptr_file = pass_state_root//'/current'
      gen_before = pass_state_current_generation()

      ! The directory of the generations exists after the first pass; its
      ! mkdir failing on an existing directory is the ordinary case, and a
      ! directory that cannot be created shows in the mkdir that follows.
      rc_dir = posix_mkdir(pass_state_root//c_null_char, dir_mode_0755)
      rc_dir = posix_mkdir(trim(tmp_dir)//c_null_char, dir_mode_0755)
      if (rc_dir .ne. 0) then
         write(*,'(A,I0,A)') ' (publish_pass_state_generation) the state'// &
              ' of pass ', pass_now, ' is not written: the directory '//  &
              trim(tmp_dir)//' could not be created; the generation'//    &
              ' published before is left as it was'
         return
      endif

      ! ---- 1. the generation, written in the temporary directory ----
      hold_cert   = state_is_certified
      hold_reason = state_certification_reason
      write(snap_reason,'(A,I0)') 'pass_snapshot_p', pass_now
      call set_state_certified(.false., trim(snap_reason))
      if (thereis_mol)     nmol_hold = nmol_eq
      if (thereis_oxychem) nox_hold  = nox_eq
      call molecular_carrier_densities_from_state(rho, f_sp)
      mrow_none = 0.0d0

      call write_hydro_state_file(trim(tmp_dir)//'/'//hyd_leaf, 'eq',     &
                                  rho, v, p, T, heat, cool, mrow_none,    &
                                  pass_statement = pass_statement,        &
                                  ios_open = ios_h,                       &
                                  state_id = trim(gen_name))
      ios_i = -1
      if (ios_h .eq. 0)                                                    &
         call write_species_state_file(trim(tmp_dir)//'/'//ion_leaf,      &
                                       'eq', nhi, nhii, nhei, nheii,      &
                                       nheiii, nheiTR, nm,                &
                                       pass_statement = pass_statement,   &
                                       ios_open = ios_i,                  &
                                       state_id = trim(gen_name))

      if (thereis_mol)     nmol_eq = nmol_hold
      if (thereis_oxychem) nox_eq  = nox_hold
      call set_state_certified(hold_cert, trim(hold_reason))

      ios_m = -1
      if (ios_h .eq. 0 .and. ios_i .eq. 0) then
         bytes_h = -1
         bytes_i = -1
         inquire(file = trim(tmp_dir)//'/'//hyd_leaf, size = bytes_h)
         inquire(file = trim(tmp_dir)//'/'//ion_leaf, size = bytes_i)
         open(newunit = unit_m, file = trim(tmp_dir)//'/'//manifest_leaf, &
              status = 'replace', action = 'write', iostat = ios_m)
         if (ios_m .eq. 0) then
            write(unit_m,'(A)') '# EXHALE pass state generation'//         &
                 ' (publish_pass_state_generation, write_output.f90);'//  &
                 ' the directory of the generation is named by state_id'
            write(unit_m,'(A)') 'state_id='//trim(gen_name)
            write(unit_m,'(A)') 'run_identity='//trim(run_identity)
            write(unit_m,'(A,I0)') 'pass=', pass_now
            if (len_trim(gen_before) .gt. 0) then
               write(unit_m,'(A)') 'previous_state_id='//trim(gen_before)
            else
               write(unit_m,'(A)') 'previous_state_id=none'
            endif
            write(unit_m,'(A)') 'hydro_file='//hyd_leaf
            write(unit_m,'(A,I0)') 'hydro_bytes=', bytes_h
            write(unit_m,'(A)') 'species_file='//ion_leaf
            write(unit_m,'(A,I0)') 'species_bytes=', bytes_i
            write(unit_m,'(A)') 'certified=F'
            write(unit_m,'(A)') 'cert_reason='//trim(snap_reason)
            close(unit_m, iostat = ios_m)
         endif
      endif
      if (ios_h .ne. 0 .or. ios_i .ne. 0 .or. ios_m .ne. 0) then
         write(*,'(A,I0,A)') ' (publish_pass_state_generation) the state'// &
              ' of pass ', pass_now, ' could not be written in '//        &
              trim(tmp_dir)//'; the directory is removed and the'//       &
              ' generation published before is left as it was'
         call remove_pass_state_directory(tmp_dir)
         return
      endif
      call stop_if_pass_state_step(1, pass_now)

      ! ---- 2. the complete generation under its own name ----
      rc_gen = posix_rename(trim(tmp_dir)//c_null_char,                   &
                            trim(gen_dir)//c_null_char)
      if (rc_gen .ne. 0) then
         write(*,'(A,I0,A,I0,A)') ' (publish_pass_state_generation) the'// &
              ' state of pass ', pass_now, ' was written but its'//       &
              ' directory not renamed (rename status ', rc_gen, '); '//   &
              trim(tmp_dir)//' is left, and not published'
         return
      endif
      call stop_if_pass_state_step(2, pass_now)

      ! ---- 3. the publication: one rename of the pointer ----
      open(newunit = unit_m, file = trim(ptr_file)//'.part',              &
           status = 'replace', action = 'write', iostat = ios_m)
      if (ios_m .eq. 0) then
         write(unit_m,'(A)') trim(gen_name)
         close(unit_m, iostat = ios_m)
      endif
      rc_ptr = -1
      if (ios_m .eq. 0) rc_ptr = posix_rename(trim(ptr_file)//'.part'//    &
                                  c_null_char, trim(ptr_file)//c_null_char)
      if (rc_ptr .ne. 0) then
         write(*,'(A,I0,A)') ' (publish_pass_state_generation) the'//     &
              ' state of pass ', pass_now, ' is complete in '//           &
              trim(gen_dir)//' but '//trim(ptr_file)//' could not be'//  &
              ' made to name it; the pointer still names the'//           &
              ' generation published before'
         return
      endif
      call stop_if_pass_state_step(3, pass_now)

      ! ---- 4. keep the generation replaced now, remove the one before ----
      if (len_trim(gen_before) .gt. 0) then
         gen_retired = pass_state_previous_generation(gen_before)
         if (len_trim(gen_retired) .gt. 0 .and.                           &
             trim(gen_retired) .ne. trim(gen_before) .and.                &
             trim(gen_retired) .ne. trim(gen_name))                       &
            call remove_pass_state_directory(pass_state_root//'/'//       &
                                             trim(gen_retired))
      endif
      written = .true.
      published_id = trim(gen_name)
      end subroutine publish_pass_state_generation

      ! ------------------------------------------------------------------ !

      function pass_state_generation_name_valid(gen) result(valid)
      ! A generation name this writer can have made: 'r' followed by
      ! characters of [A-Za-z0-9._-] only. It is therefore one path
      ! component inside output/pass_state and never '.' or '..', which is
      ! what makes a name read from a file safe to remove.
      character(len=*), intent(in) :: gen
      logical :: valid
      valid = .false.
      if (len_trim(gen) .eq. 0) return
      if (gen(1:1) .ne. 'r') return
      valid = (pass_state_safe_name(trim(gen)) .eq. trim(gen))
      end function pass_state_generation_name_valid

      ! ------------------------------------------------------------------ !

      function pass_state_current_generation() result(gen)
      ! The generation output/pass_state/current names, or '' when there is
      ! none (or the line is not a generation name).
      character(len=96) :: gen
      integer :: unit_c, ios_c
      gen = ''
      open(newunit = unit_c, file = pass_state_root//'/current',          &
           status = 'old', action = 'read', iostat = ios_c)
      if (ios_c .ne. 0) return
      read(unit_c,'(A)', iostat = ios_c) gen
      close(unit_c)
      if (ios_c .ne. 0) gen = ''
      gen = adjustl(gen)
      if (.not. pass_state_generation_name_valid(gen)) gen = ''
      end function pass_state_current_generation

      ! ------------------------------------------------------------------ !

      function pass_state_previous_generation(gen) result(gen_prev)
      ! The previous_state_id the manifest of generation gen states, or ''
      ! when it states none or the manifest cannot be read.
      character(len=*), intent(in) :: gen
      character(len=96) :: gen_prev
      character(len=256) :: mline
      integer :: unit_c, ios_c
      gen_prev = ''
      open(newunit = unit_c, file = pass_state_root//'/'//trim(gen)//     &
           '/manifest.txt', status = 'old', action = 'read',              &
           iostat = ios_c)
      if (ios_c .ne. 0) return
      do
         read(unit_c,'(A)', iostat = ios_c) mline
         if (ios_c .ne. 0) exit
         if (index(mline, 'previous_state_id=') .eq. 1) then
            gen_prev = adjustl(mline(len('previous_state_id=')+1:))
            exit
         endif
      enddo
      close(unit_c)
      if (.not. pass_state_generation_name_valid(gen_prev)) gen_prev = ''
      end function pass_state_previous_generation

      ! ------------------------------------------------------------------ !

      subroutine remove_pass_state_directory(dir_path)
      ! Remove a directory of the pass state: the three files a generation
      ! holds, by their names, and then the directory, which remove(3)
      ! removes only when it is empty. Nothing is removed by a pattern, so a
      ! file this writer did not put there stops the removal of the
      ! directory, and that is reported.
      character(len=*), intent(in) :: dir_path
      integer(c_int) :: rc_rm
      rc_rm = posix_remove(trim(dir_path)//'/Hydro_ioniz.txt'//c_null_char)
      rc_rm = posix_remove(trim(dir_path)//'/Ion_species.txt'//c_null_char)
      rc_rm = posix_remove(trim(dir_path)//'/manifest.txt'//c_null_char)
      rc_rm = posix_remove(trim(dir_path)//c_null_char)
      if (rc_rm .ne. 0)                                                    &
         write(*,'(A)') ' (publish_pass_state_generation) '//             &
              trim(dir_path)//' could not be removed (absent, or it'//    &
              ' holds a file this writer did not write); left as it is'
      end subroutine remove_pass_state_directory

      ! ------------------------------------------------------------------ !

      subroutine stop_if_pass_state_step(step_done, pass_now)
      ! EXHALE_PASS_STATE_STOP_AT=<step>[:<pass>] (default: unset, no
      ! effect): stop the run with status 3 right after step <step> of the
      ! publication of pass <pass> (of every pass when no pass is named),
      ! the interruption the order of publication is built against. A
      ! measurement key; a value that cannot be read is ignored.
      integer, intent(in) :: step_done, pass_now
      character(len=32) :: stop_text
      integer :: stop_step, stop_pass, icolon, ios_s
      call get_environment_variable('EXHALE_PASS_STATE_STOP_AT', stop_text)
      if (len_trim(stop_text) .eq. 0) return
      stop_pass = -1
      icolon = index(stop_text, ':')
      if (icolon .gt. 0) then
         read(stop_text(1:icolon-1),*,iostat=ios_s) stop_step
         if (ios_s .ne. 0) return
         read(stop_text(icolon+1:),*,iostat=ios_s) stop_pass
         if (ios_s .ne. 0) return
      else
         read(stop_text,*,iostat=ios_s) stop_step
         if (ios_s .ne. 0) return
      endif
      if (stop_step .ne. step_done) return
      if (stop_pass .ge. 0 .and. stop_pass .ne. pass_now) return
      write(*,'(A,I0,A,I0,A)') ' (publish_pass_state_generation) stopped'//&
           ' after step ', step_done, ' of the publication of pass ',     &
           pass_now, ' (EXHALE_PASS_STATE_STOP_AT)'
      stop 3
      end subroutine stop_if_pass_state_step

      ! ------------------------------------------------------------------ !

      subroutine write_adv_validity_header(unit, T_status, comp_status)
      ! WHAT THE ADVECTION-CORRECTED FILES ARE, AND WHICH ROWS OF THEM CAN
      ! BE READ AS A STEADY SOLUTION.  One author for the block, so that
      ! Hydro_ioniz_adv.txt and Ion_species_adv.txt cannot disagree about
      ! the validity of the same row.
      !
      ! The block carries six things a reader of the file alone cannot
      ! otherwise know:
      !
      !   the schema version, so that a file written before the two fields
      !   existed is read as UNKNOWN validity and not as a corrected one;
      !   the meaning of every value of both fields;
      !   the counts, taken from the columns written next to them, so the
      !   header and the columns cannot drift apart;
      !   the operator each row was judged by, whose value for every row is
      !   a column of the file, and the fraction a corrected row is
      !   accurate to (adv_conditional_tol);
      !   whether the INPUT state passed the stationary certification,
      !   which is a statement about the whole state and not about a row;
      !   what the product is, and the closure's own restrictions.
      integer, intent(in) :: unit
      integer, dimension(1-Ng:N+Ng), intent(in) :: T_status, comp_status
      integer, dimension(0:4) :: nT, nc
      integer :: j, ist

      nT = 0;  nc = 0
      do j = 1-Ng,N+Ng
         ist = T_status(j)
         if (ist .ge. 0 .and. ist .le. 4) nT(ist) = nT(ist) + 1
         ist = comp_status(j)
         if (ist .ge. 0 .and. ist .le. 4) nc(ist) = nc(ist) + 1
      enddo

      write(unit,'(A,I0)') '# adv_schema ', adv_schema_version
      write(unit,'(A)') '# adv_T_status validity of the temperature of'//  &
                     ' the row: 0 corrected (the steady energy'
      write(unit,'(A)') '#   equation was solved and adopted), 1'//        &
                     ' retained (the run''s own T kept: the mass'
      write(unit,'(A)') '#   row of the cell is above the fraction'//      &
                     ' below, or the local radiative balance rather'
      write(unit,'(A)') '#   than the flow sets it, or the gas enters'//   &
                     ' the cell), 2 failed (the cell solve did not'
      write(unit,'(A)') '#   converge or returned an out-of-band'//        &
                     ' root), 3 unsupported (the closure does not'
      write(unit,'(A)') '#   cover this cell class), 4 not_evaluated'//    &
                     ' (the post-process did not reach this row).'
      write(unit,'(A)') '# adv_comp_status the same five values for the'// &
                     ' COMPOSITION of the row: 0 the steady'
      write(unit,'(A)') '#   advection-ionization solution, 1 the'//       &
                     ' equilibrium composition of the run, 2 a cell'
      write(unit,'(A)') '#   solve that did not converge, 3 a cell'//      &
                     ' class the reconstruction does not model, 4 not'
      write(unit,'(A)') '#   reached. 3 and 4 are structural and stand'//  &
                     ' above 0, 1 and 2, which name the flow and the'
      write(unit,'(A)') '#   solve; the two fields are independent, so'//  &
                     ' a corrected T can carry a retained'
      write(unit,'(A)') '#   composition and the reverse.'
      write(unit,'(A,5(1x,i0),A,5(1x,i0))') '# adv_status_counts T',       &
                     (nT(ist), ist = 0,4), ' comp', (nc(ist), ist = 0,4)
      write(unit,'(A)') '# adv_mass_row the measure both fields were'//    &
                     ' decided by, for every row of'
      write(unit,'(A)') '#   Hydro_ioniz_adv.txt: the operator below,'//   &
                     ' against the fraction below it. It is the'
      write(unit,'(A)') '#   CONDITION of the row, so a reader can'//      &
                     ' weigh a corrected row rather than only count it.'
      if (state_is_certified) then
         ! The reason survives a certified state too: 'certified_in_wind'
         ! says the species rows below the wind radius were reported and
         ! did not gate.
         write(unit,'(A)') '# adv_input_certified T the input state'//      &
                     ' passed the stationary certification '//             &
                     trim(state_certification_reason)
      else if (len_trim(state_certification_reason) .gt. 0) then
         write(unit,'(A)') '# adv_input_certified F '//                     &
                     trim(state_certification_reason)//                     &
                     ' (see the certification report of the run)'
      else
         write(unit,'(A)') '# adv_input_certified F no certification was'// &
                     ' made on the input state'
      endif
      write(unit,'(A)') '# adv_stationarity_operator face_flux_mass the'// &
                     ' mass row of the state, |R_1|/s_1 with s_1 the'
      write(unit,'(A)') '#   row''s own largest term, which is the'//      &
                     ' fractional change of the FACE mass flux across'
      write(unit,'(A)') '#   the cell; the same operator the stationary'// &
                     ' certification measures that row by, and its'
      write(unit,'(A)') '#   value for every row is the adv_mass_row'//    &
                     ' column of Hydro_ioniz_adv.txt.'
      write(unit,'(A,ES8.1,A)') '# adv_conditional_tol ',                  &
                     adv_conditional_tol, ' the correction of a'//          &
                     ' corrected row is accurate to this'
      write(unit,'(A)') '#   fraction of itself in the mass flux: the'//   &
                     ' terms the correction drops are the ones the'
      write(unit,'(A)') '#   mass divergence puts into the steady'//       &
                     ' equations, each of them adv_mass_row times a'
      write(unit,'(A)') '#   term they keep, so the product is FIRST'//    &
                     ' ORDER in that measure. A row at or below this'
      write(unit,'(A)') '#   fraction is corrected and is a'//             &
                     ' CONDITIONAL correction to it; a row above it is'
      write(unit,'(A)') '#   retained. This is not the tolerance of a'//   &
                     ' certified stationary state, which is a'
      write(unit,'(A)') '#   statement about the WHOLE state and is'//     &
                     ' reported in the adv_input_certified line.'
      write(unit,'(A)') '# adv_product a ONE-WAY correction on the'//      &
                     ' density and velocity field of the input state:'
      write(unit,'(A)') '#   rho and v are written back unchanged and'//   &
                     ' nothing solved here feeds back into them. A'
      write(unit,'(A)') '#   converged scalar temperature solve'//         &
                     ' therefore does NOT make the row a'
      write(unit,'(A)') '#   self-consistent solution: the momentum and'// &
                     ' mass balances of the corrected pressure were'
      write(unit,'(A)') '#   not re-solved, and a profile assembled'//     &
                     ' from corrected and retained cells need not'
      write(unit,'(A)') '#   satisfy any interface balance.'
      write(unit,'(A)') '# adv_model_restrictions the post-process'//      &
                     ' reconstructs an H/He + trace-metal gas: the'
      write(unit,'(A)') '#   molecular (H2 H2p H3p HeHp) and oxygen'//     &
                     ' (OH H2O CO) columns are NOT advection'
      write(unit,'(A)') '#   corrected and are the equilibrium'//          &
                     ' solution, those species are absent from the'
      write(unit,'(A)') '#   particle and electron counts every'//         &
                     ' equation here is solved with, and the H3+'
      write(unit,'(A)') '#   infrared cooling is absent from its'//        &
                     ' energy balance, so the metal and H/He columns'
      write(unit,'(A)') '#   inside a molecular layer inherit that'//      &
                     ' approximation.'

      end subroutine write_adv_validity_header

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
      ! With "Molecular carrier transport: True" -- the default of the option -- the
      ! carriers are transported, so a cell in which the chemistry is slow
      ! is SOLVED rather than assumed away. The time scales stay a required
      ! output for the other direction: they say where the answer is
      ! chemistry, where it is flow, and where the two are comparable, and
      ! that is what tells a reader which part of the profile the reaction
      ! set is responsible for. With "Molecular carrier transport: False" they say
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
      real*8  :: nH2b, Tb, nOb, nHIb, h2loss(11), h2tot, n3b
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
      real*8 :: nph_b, relerr, worst, e_lw_abs, e_lw_in, ph_lw
      real*8 :: NH2_out_lw, tau_out_lw
      ! The G4 ledger.  cont_ph and cont_beam carry the H2O and OH continuum
      ! alone -- the part of a band that must close cell by cell -- while
      ! absph and absen carry every rated absorber, so in the shared
      ! Lyman-Werner band they also hold the H2 pumps.  beam_loss is what the
      ! beam itself loses over the column, and tr_line_in the transmission of
      ! the band's LINE absorber at the innermost face.
      real*8 :: cont_ph(n_fuv_band), cont_beam(n_fuv_band)
      real*8 :: beam_loss(n_fuv_band), tr_line_in(n_fuv_band)
      real*8 :: unrated(n_fuv_band), unrated_f(n_fuv_band)
      real*8 :: a_cell, dtau_col, dn_h2o, dn_oh, dn_h2, dn_co
      ! How far the H2O density written here has moved from the one the
      ! photon field was built on, as a fraction of the cell's own column
      ! increment, and where the worst cell of the column is.
      real*8 :: drift_worst, drift_r, dstate
      ! Line-centre H I Ly-alpha depth of the state this file writes, for the
      ! B2 beam transmission.
      real*8, dimension(1-Ng:N+Ng) :: tau_lya_col
      ! How much of the oxygen and the carbon the option has moved out of
      ! the atomic coolants, and where.
      ! Column-integrated infrared exchange of the molecular bands, per unit
      ! area: spontaneous emission, absorption of the field from below, and
      ! their difference. Index 1 H2, 2 H2O, 3 CO.
      real*8 :: ir_emit(3), ir_abs(3), ir_e, ir_n, w_ir, ir_bound
      real*8 :: nO_el, nC_el, fO_mol, fC_mol, fO_worst, fC_worst
      ! The domain record of the one-sided CO destruction model.
      integer :: dom_out, dom_hot, dom_hep
      real*8  :: dom_ratio, dom_r, dom_form
      ! The CO row of the Lyman-Werner ledger: photons removed, the energy
      ! they carry and the share of it that reaches the gas, all per unit
      ! area of the star-ward column.
      real*8  :: co_ph, co_en, co_heat, co_share
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
      if (carrier_transport) then
         write(74,'(A)') '# The carriers ARE transported (implicit'//      &
                     ' diffusion-advection with the chemistry).'
      else
         write(74,'(A)') '# The carriers are a LOCAL steady state'//       &
                     ' ("Molecular carrier transport: False"): chemistry alone.'
      endif
      write(74,'(A)') '# columns r[Rp] T[K] n_O n_OII n_OIII n_OH n_H2O'//  &
                     ' n_CO n_O1D x_H2[2nH2/nH] tau_chem[s] tau_adv[s]'//   &
                     ' Da j_H2O[1/s] j_OH[1/s] D_H2[cm2/s] Kzz[cm2/s]'//    &
                     ' tau_diff[s]'
      call write_row_layout_header(74)
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
      if (carrier_transport) then
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
         write(74,'(A)') '# THE DOMAIN OF THE ONE-SIDED CO'//            &
                     ' DESTRUCTION MODEL. The CO row destroys CO --'
         write(74,'(A)') '# He+ + CO -> C+ + O + He (UMIST RATE22'//     &
                     ' 4068) and CO + hv -> C + O on the 912-1201 A'
         write(74,'(A)') '# beam with the Visser et al. (2009)'//        &
                     ' shielding function -- and never forms it. What'
         write(74,'(A)') '# makes that legitimate in a cell is'//        &
                     ' tau_dest << tau_res << tau_form, with'
         write(74,'(A)') '# tau_res = min(r/|v|, dr^2/(D_CO + K_zz)).'//&
                     ' The record is INFORMATIONAL: the rates are on'
         write(74,'(A)') '# everywhere, and where the ordering fails'// &
                     ' it is the omitted FORMATION that fails with'
         write(74,'(A)') '# it, so the transported value stands, which'//&
                     ' is what a transport operator should do in a'
         write(74,'(A)') '# quenched layer.'
         write(74,'(A)') '#   The three counts are CELL VISITS'//     &
                     ' summed over every carrier interval of the run,'
         write(74,'(A)') '#   not distinct cells: what they say is'//    &
                     ' whether the model was ever out of domain and'
         write(74,'(A)') '#   how much of the integration was, and'//    &
                     ' they are not rolled back by a refused step,'
         write(74,'(A)') '#   because whether a state was in domain'//   &
                     ' has an answer whether or not the step that'
         write(74,'(A)') '#   read it was accepted.'
         write(74,'(A)') '#   dom_cells_out   cell visits with'//        &
                     ' tau_dest > f_dom tau_res'
         write(74,'(A)') '#   dom_worst_ratio the largest'//             &
                     ' tau_dest/tau_res reached, and where'
         write(74,'(A)') '#   dom_form_ratio  the largest'//             &
                     ' tau_res/tau_form reached (RATE22 8597,'
         write(74,'(A)') '#                   C + O -> CO + photon);'//  &
                     ' the second inequality asks it to be small'
         write(74,'(A)') '#   dom_cells_hot   cell visits above the'// &
                     ' 512 K excitation-temperature limit of the'
         write(74,'(A)') '#                   shielding table, where'// &
                     ' Visser et al. state a factor of two'
         write(74,'(A)') '#   dom_cells_HeP   cell visits in which'//  &
                     ' the He+ channel removes He+ faster than'
         write(74,'(A)') '#                   recombination does, so'//  &
                     ' this reaction is the leading He+ loss of'
         write(74,'(A)') '#                   the cell; the sweep'//     &
                     ' He+ row carries it, and the count measures'
         write(74,'(A)') '#                   the lag of the operator'// &
                     ' split between the two'
         call carrier_co_domain_record(dom_out, dom_hot, dom_hep,        &
                                       dom_ratio, dom_r, dom_form)
         write(74,'(A,ES12.4)') '# dom_f_dom       ',                    &
                     carrier_co_domain_f_dom()
         write(74,'(A,I0)')     '# dom_cells_out   ', dom_out
         write(74,'(A,ES12.4)') '# dom_worst_ratio ', dom_ratio
         write(74,'(A,ES12.4)') '# dom_worst_r     ', dom_r
         write(74,'(A,ES12.4)') '# dom_form_ratio  ', dom_form
         write(74,'(A,I0)')     '# dom_cells_hot   ', dom_hot
         write(74,'(A,I0)')     '# dom_cells_HeP   ', dom_hep
      else
         write(74,'(A)') '#'
         write(74,'(A)') '# CO IS AN EQUILIBRIUM CLOSURE IN THIS RUN.'// &
                     ' Without the carrier transport there is no CO'
         write(74,'(A)') '# balance row to carry the two destruction'// &
                     ' rates, so CO is the CO <-> C + O chemical'
         write(74,'(A)') '# equilibrium of each cell''s own (n, T).'//   &
                     ' That is a closure with its own domain and not'
         write(74,'(A)') '# a kinetic result: it cannot quench, and'//  &
                     ' it re-forms CO in a cool outer wind that a'
         write(74,'(A)') '# transported CO would have been carried'//   &
                     ' out of the molecular layer into.'
      endif

      ! ---- what actually runs the H2 partition at the base ----
      ! The option exists to make the OXYGEN cycle set the base H2/H
      ! partition, and whether it does is a property of the run rather than
      ! of the code: on a 2331 K base the thermal channel H2 + M -> H + H + M
      ! carries 70-97% of the net and the oxygen cycle a few percent, while
      ! on a 864 K base the oxygen family carries 96-99.6%. A run whose base is
      ! hotter than it should be -- and every molecular run is, until the
      ! H2O and CO infrared bands exist -- is answering a
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
      ! (1) thermal dissociation net of the three-body association R15.
      ! THE THIRD BODY IS THE COLLIDER SUM, k1(H2) n(H2) + k1(H) n(H) +
      ! k1(Ar) n(He) (Cohen & Westberg 1983, p. 559), formed by the one
      ! routine that owns it and carried by BOTH directions, so the pair
      ! stays an exact detailed balance and this diagnostic reports the
      ! rate the composition solver and the energy ledger both ran at.  The
      ! helium collider is the ground-singlet neutral.
      if (reaction_heat_recipients_corrected()) then
         n3b = h2_association_collider_density(Tb, nH2b, nHIb,            &
                          he_ground_singlet_density(nhei(jb), nheiTR(jb))*n0)
      else
         n3b = n_tot_ox(jb)
      endif
      h2loss(1) = rk_R12_H2_thdis(Tb)*n3b*nH2b                            &
                - rk_R15_3body_H2(Tb, n3b)*nHIb*nHIb
      h2loss(2) = k_lw_diss(jb)*nH2b
      h2loss(3) = P_H2_eq(jb)*nH2b
      h2loss(4) = (rk_R10_Hp_H2v4(Tb)                                      &
                 + rk_R13_Hp_H2_M(n_tot_ox(jb)))*nhii(jb)*n0*nH2b
      h2loss(5) = rk_R14_H2_edis(Tb)*ne_ox(jb)*nH2b
      h2loss(6) = rk_R8_H2p_H2()*nmol_eq(jb,2)*nH2b
      ! R17 and R23 only: the HeH+ channel of He+ + H2 (Koskinen R20) was
      ! retired, its cited measurement having bounded it 42 times below the
      ! value Table 1 gave it (mol_rates).
      h2loss(7) = (rk_R17_Hep_H2_diss(Tb)                                  &
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
      ! That is the missing molecular infrared cooling arriving through the
      ! composition instead of through the temperature: A2 is the
      ! composition of the molecular layer and that cooling is its energy.
      ! Measured
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
            ' would really be.'

      ! ---- FUV band penetration and the G4 band ledger ----
      ! THE QUANTITY THAT IS CLOSED. For every band, the photons the beam
      ! loses between the two ends of the column equal the absorptions the
      ! model accounts for -- but only when every absorber of the band has a
      ! rate. Two of the four bands have an absorber that takes photons out
      ! of the beam without one, so the ledger below is written as three
      ! separate statements rather than as one residual that would read as a
      ! defect wherever the physics says otherwise.
      !
      ! (a) CLOSURE OF THE CONTINUUM ABSORBERS, band by band. In cell j the
      !     H2O and OH rates are
      !         j_s = sigma_s N_b tr exp(-tau_out) (1 - exp(-dtau))/dtau
      !     (water_photolysis_rate), with tr the band's line transmission and
      !     tau_out, dtau the CONTINUUM depth at the cell's star-ward face
      !     and across the cell. Applied to the absorbers the columns record
      !     for that cell, dN_H2O and dN_OH, the two rates take
      !         j_H2O dN_H2O + j_OH dN_OH
      !     photons out of the band there, while the beam's own loss between
      !     the cell's two faces is
      !         N_b tr exp(-tau_out) (1 - exp(-dtau))
      !               = (j_H2O/sigma_H2O) dtau .
      !     The two are the same number as long as the band's optical depth
      !     is exactly sigma_H2O dN_H2O + sigma_OH dN_OH, which is how
      !     fuv_band_optical_depth builds it apart from a clamp at negative
      !     column. rel_diff is therefore round-off, and a residual says
      !     either that a column has gone negative or that the depths and
      !     the columns on this file are not one state. It holds whatever
      !     the rate form is -- the two rates differ only by their cross
      !     sections -- so it is NOT the test of the cell mean; (b) is.
      !
      ! (b) THE BEAM BUDGET. beam_loss = N_b (1 - T_line T_cont) at the
      !     innermost face is what the beam loses over the whole column;
      !     rated_ph is every absorption the model's rates account for. Their
      !     difference is the beam that leaves without a rate. Band by band:
      !       B1, B3, B4  no line absorber (T_line = 1) and both continuum
      !                   absorbers rated, so the difference is a closure
      !                   residual -- and it is the test of the cell-mean
      !                   rate form, because
      !                   sum_j N_b exp(-tau_out) (1 - exp(-dtau))
      !                   telescopes to N_b (1 - exp(-tau)) exactly, at any
      !                   grid spacing, only when the rate is the mean over
      !                   the cell. A rate read at one point of the cell
      !                   breaks it one-signed, by 30 percent in the
      !                   Ly-alpha band of an HD 189733 b run
      !                   (water_photolysis.f90).
      !       B2          the band IS the 1215.67 A H I resonance line. The
      !                   H I that scatters the stellar line out of the beam
      !                   dissociates nothing and therefore carries no rate,
      !                   so the difference is the share H I takes; near 1 it
      !                   says the line is thick, which is not a defect.
      !                   T_line is the erfc penetration of the stellar line
      !                   (lya_rt.f90), rebuilt here from the T and n_HI on
      !                   this file, through the same line-centre depth the
      !                   band rate is built from.
      !       LW          over 912-1201 A the same photons are taken by H2 in
      !                   the Lyman-Werner lines and by H2O and OH in a
      !                   continuum, out of ONE beam (water_photolysis.f90
      !                   sec. 3). The H2 share IS rated and is counted in
      !                   rated_ph -- its photon count is the cell mean of
      !                   the PUMP cross section sigma_diss/p_eff, because
      !                   the pumps that do not dissociate also take a photon
      !                   out of the band -- and since 2026-09-06 T_line
      !                   carries the column integral of that same cross
      !                   section, so both sides are one normalization of one
      !                   absorption and this is a CLOSURE residual, as it is
      !                   for B3 and B4. Two things are left in it. The
      !                   discretization: the rate is the cell mean of
      !                   sigma_pump exp(-tau_cont) while the beam identity
      !                   puts the line loss of a cell behind that cell's own
      !                   continuum depth. And, where the star-ward column
      !                   passes the top of the table's column axis, the
      !                   clamped edge cross section drives the line-removed
      !                   fraction past 1, where the beam is capped and the
      !                   rate is not; the run warns about that column
      !                   separately.
      !
      ! (c) THE ENERGY, which is exact by construction: every absorbed photon
      !     carries <hv>_b, of which the threshold goes into the bond and the
      !     rest into the gas, so absorbed = heat + bond cell by cell. It is
      !     computed anyway, because a mismatch would mean the heating and
      !     the rates had drifted apart. In the LW band only the fragment
      !     kinetic energy of the dissociating fraction reaches the gas, so
      !     there the 'bond' column is the fluorescence as well as the bond.
      ! The band ledger is formed in ONE place, beside the heating assembly
      ! (utils_ion_eq::fuv_band_absorption_ledger), because the energy of one
      ! absorption event is the same energy the energy equation is charged.
      ! This file writes what that routine returns and forms none of it.
      call fuv_band_absorption_ledger(T*T0, nhi*n0, nmol_eq(:,1),          &
               nox_eq(:,2), he_ground_singlet_density(nhei, nheiTR)*n0,   &
               (nhi + nhii)*n0 + 2.0d0*(nmol_eq(:,1) + nmol_eq(:,2))      &
                               + 3.0d0*nmol_eq(:,3) + nmol_eq(:,4),       &
               NH2_col_lw, NH2O_col, NOH_col, NCO_col, tau_fuv,           &
               j_h2o_fuv, j_oh_fuv, k_lw_diss, p_lw_single, k_co_diss,    &
               absph, absen, heat_col, bond_col, cont_ph, cont_beam,      &
               co_ph, co_en, co_heat, e_lw_abs, drift_worst, drift_r)

      ! ---- What the beam itself loses over the column, band by band ----
      ! N_b (1 - T_line T_cont) at the innermost face, with T_cont the
      ! continuum transmission exp(-tau) of the band and T_line that of its
      ! LINE absorber: the H2 Lyman-Werner lines in the shared band (1 - A,
      ! the equivalent width of Draine & Bertoldi 1996), the stellar H I
      ! Ly-alpha line in B2, and 1 in the three bands that have no line
      ! absorber. The B2 factor is the transmission of the state on this
      ! file: the same line-centre depth rule and the same erfc penetration
      ! the band rate is built from (lya_rt.f90), evaluated on the T and
      ! n_HI written here.
      tr_line_in = 1.0d0
      tr_line_in(ib_LW) = max(tr_lines_lw(1-Ng), 0.0d0)
      if (F_Lya_star .gt. 0.0d0) then
         call lya_line_center_optical_depth(T*T0, nhi*n0, tau_lya_col)
         tr_line_in(ib_B2) =                                               &
            lya_stellar_beam_transmission(T(1-Ng)*T0, tau_lya_col(1-Ng))
      endif
      do ib = 1,n_fuv_band
         nph_b = fuv_band_photon_flux(fuv_band_flux(ib), ib)
         beam_loss(ib) = nph_b*(1.0d0 - tr_line_in(ib)                     &
                              *exp(-max(tau_fuv(1-Ng,ib), 0.0d0)))
         unrated(ib)   = beam_loss(ib) - absph(ib)
         unrated_f(ib) = 0.0d0
         if (beam_loss(ib) .gt. 0.0d0)                                     &
            unrated_f(ib) = unrated(ib)/beam_loss(ib)
      enddo

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
      colhdr = '# columns r[Rp] N_H2O[cm^-2] N_OH[cm^-2] N_CO[cm^-2]'
      do ib = 1,n_fuv_band
         colhdr = trim(colhdr)//' tau_'//trim(fuv_band_name(ib))
      enddo
      do ib = 1,n_fuv_band
         colhdr = trim(colhdr)//' j_H2O_'//trim(fuv_band_name(ib))
      enddo
      do ib = 1,n_fuv_band
         colhdr = trim(colhdr)//' j_OH_'//trim(fuv_band_name(ib))
      enddo
      ! CO belongs to the Lyman-Werner band alone: its 37 predissociating
      ! lines lie between 912.7 and 1076.1 A, and it takes photons out of
      ! that beam only through its own shielding function, which is why
      ! there is no tau_CO column beside tau_LW (co_photodissociation.f90).
      colhdr = trim(colhdr)//' k_CO[1/s] Theta_CO'
      colhdr = trim(colhdr)//' heat_FUV[erg/cm3/s]'
      write(75,'(A)') trim(colhdr)
      call write_row_layout_header(75)
      do j = 1-Ng,N+Ng
         write(75,*) r(j), NH2O_col(j), NOH_col(j), NCO_col(j),            &
                    (tau_fuv(j,ib), ib = 1,n_fuv_band),                    &
                    (j_h2o_fuv(j,ib), ib = 1,n_fuv_band),                  &
                    (j_oh_fuv(j,ib), ib = 1,n_fuv_band),                   &
                    k_co_diss(j), theta_co_shield(j),                      &
                    heat_fuv(j)
      enddo
      write(75,'(A)') '#'
      write(75,'(A)') '# Band ledger (gate G4). Every quantity is per'//   &
                     ' unit area of the star-ward column.'
      write(75,'(A)') '#'
      write(75,'(A)') '# (a) CLOSURE OF THE CONTINUUM ABSORBERS. In cell'//&
                     ' j the H2O and OH rates take'
      write(75,'(A)') '#     N_b tr exp(-tau_out) (1 - exp(-dtau))'//      &
                     ' photons out of the band, with tr the'
      write(75,'(A)') '#     line transmission of the band and tau_out,'// &
                     ' dtau the CONTINUUM depth at the'
      write(75,'(A)') '#     cell''s star-ward face and across the cell.'
      write(75,'(A)') '#       cont_absorbed_ph  sum over cells of'//      &
                     ' (j_H2O dN_H2O + j_OH dN_OH), the two rates'
      write(75,'(A)') '#                         applied to the'//         &
                     ' absorbers the columns record [cm^-2 s^-1]'
      write(75,'(A)') '#       cont_beam_loss    the same cells as the'//  &
                     ' beam''s loss between their two faces,'
      write(75,'(A)') '#                         (j_H2O/sigma_H2O) dtau,'//&
                     ' which names no density at all'
      write(75,'(A)') '#       rel_diff          |cont_absorbed_ph/'//     &
                     'cont_beam_loss - 1|'
      write(75,'(A)') '#     The two are the same number as long as'//    &
                     ' dtau is exactly sigma_H2O dN_H2O + sigma_OH dN_OH,'
      write(75,'(A)') '#     so rel_diff is round-off: a residual says a'//&
                     ' column has gone negative, or that the'
      write(75,'(A)') '#     depths and the columns on this file are'//    &
                     ' not one state. It holds whatever the rate'
      write(75,'(A)') '#     form is, so it is NOT the test of the cell'// &
                     ' mean; the B1/B3/B4 rows of (b) are.'
      write(75,'(A)') '# band cont_absorbed_ph cont_beam_loss rel_diff'
      worst = 0.0d0
      do ib = 1,n_fuv_band
         relerr = 0.0d0
         if (cont_beam(ib) .gt. 0.0d0)                                     &
            relerr = abs(cont_ph(ib) - cont_beam(ib))/cont_beam(ib)
         worst = max(worst, relerr)
         write(75,'(A,A,3ES14.6)') '# ', fuv_band_name(ib), cont_ph(ib),   &
                    cont_beam(ib), relerr
      enddo
      write(75,'(A)') '#'
      write(75,'(A)') '# (b) BEAM BUDGET.'
      write(75,'(A)') '#       beam_loss_ph  photons the beam loses'//     &
                     ' between the two ends of the column,'
      write(75,'(A)') '#                     N_b (1 - T_line T_cont) at'// &
                     ' the innermost face [cm^-2 s^-1]'
      write(75,'(A)') '#       rated_ph      every absorption the'//       &
                     ' model''s RATES account for over that column'
      write(75,'(A)') '#                     (the H2O and OH continua,'//  &
                     ' and in the LW band the H2 pumps)'
      write(75,'(A)') '#       unrated_ph    their difference, and'//      &
                     ' unrated_frac its share of beam_loss_ph'
      write(75,'(A)') '#     WHAT unrated_frac MEANS, BAND BY BAND:'
      write(75,'(A)') '#       B1 B3 B4  no line absorber (T_line = 1)'//  &
                     ' and every absorber rated, so this is a'
      write(75,'(A)') '#                 CLOSURE residual and is'//        &
                     ' round-off. It is the test of the cell-mean rate'
      write(75,'(A)') '#                 form: the cell sums telescope'// &
                     ' to N_b (1 - exp(-tau)) exactly, at any grid'
      write(75,'(A)') '#                 spacing, only because the rate'//&
                     ' is the mean over the cell.'
      write(75,'(A)') '#       B2        the band IS the 1215.67 A H I'//  &
                     ' resonance line. The H I that scatters'
      write(75,'(A)') '#                 the stellar line out of the'//    &
                     ' beam dissociates nothing and carries no'
      write(75,'(A)') '#                 rate, so this is the SHARE H I'// &
                     ' takes. A value near 1 says the line is'
      write(75,'(A)') '#                 thick; it is not a defect.'//     &
                     ' T_line is the erfc penetration of the'
      write(75,'(A)') '#                 stellar line, on the T and'//     &
                     ' n_HI this file writes.'
      write(75,'(A)') '#       LW        the H2 Lyman-Werner lines ARE'//  &
                     ' rated and are counted in rated_ph, and'
      write(75,'(A)') '#                 T_line is 1 - A with A the'//     &
                     ' column integral of the pump cross section'
      write(75,'(A)') '#                 of the same level-resolved'//     &
                     ' table the rate comes from, so this is a'
      write(75,'(A)') '#                 CLOSURE residual like B3 and'//   &
                     ' B4: one normalization of one absorption.'
      write(75,'(A)') '#                 What is left in it is the'//      &
                     ' discretization of the shared beam inside'
      write(75,'(A)') '#                 a cell, and it vanishes with'//   &
                     ' the cell continuum depth.'
      write(75,'(A)') '# band rated_ph beam_loss_ph unrated_ph'//          &
                     ' unrated_frac T_line'
      do ib = 1,n_fuv_band
         write(75,'(A,A,5ES14.6)') '# ', fuv_band_name(ib), absph(ib),     &
                    beam_loss(ib), unrated(ib), unrated_f(ib),             &
                    tr_line_in(ib)
      enddo
      write(75,'(A)') '#'
      write(75,'(A)') '# (b2) CO ON THE SAME BEAM. CO predissociates in'//&
                     ' 37 lines between 912.7 and 1076.1 A, all inside'
      write(75,'(A)') '#      the LW interval, and it shields itself in'//&
                     ' them (Visser, van Dishoeck & Black 2009).'
      write(75,'(A)') '#      Its absorptions are RATED -- k_CO is a'//   &
                     ' column of this file -- but they are printed'
      write(75,'(A)') '#      here and not in rated_ph above, because'//  &
                     ' beam_loss_ph carries no CO term: CO adds'
      write(75,'(A)') '#      nothing to tau_cont and its equivalent'//   &
                     ' width lives inside its own shielding'
      write(75,'(A)') '#      function, so the other three absorbers'//   &
                     ' see a beam undepleted by CO. Counting CO'
      write(75,'(A)') '#      against that beam would compare two'//      &
                     ' different beams. co_frac IS the size of that'
      write(75,'(A)') '#      approximation: the share of the LW beam''s'//&
                     ' own loss that CO takes a second time.'
      write(75,'(A)') '#      co_absorbed_en is charged at the CO'//      &
                     ' events'' own mean photon energy, 12.87 eV,'
      write(75,'(A)') '#      not the band mean 11.74 eV, because CO'//   &
                     ' absorbs at the blue end of the beam.'
      write(75,'(A)') '#      Only the dissociating share is counted;'//  &
                     ' the oscillator-strength weighted'
      write(75,'(A)') '#      dissociation efficiency of Visser Table 1'//&
                     ' is 0.96, so about 4 percent of the CO'
      write(75,'(A)') '#      absorptions are not here either.'
      co_share = 0.0d0
      if (beam_loss(ib_LW) .gt. 0.0d0)                                     &
         co_share = co_ph/beam_loss(ib_LW)
      write(75,'(A,4ES14.6)') '# co_absorbed_ph co_frac co_absorbed_en'// &
                     ' co_heat_en ', co_ph, co_share, co_en, co_heat
      write(75,'(A)') '#'
      write(75,'(A)') '# (c) ENERGY. Every absorbed photon carries'//      &
                     ' <hv>_b, of which the threshold goes into the'
      write(75,'(A)') '#     bond and the rest into the gas, so'//         &
                     ' absorbed_en = heat_en + bond_en exactly.'
      write(75,'(A)') '#     incident_en is the band flux;'//              &
                     ' absorbed_en must not exceed it. [erg cm^-2 s^-1]'
      write(75,'(A)') '# band absorbed_en heat_en bond_en incident_en'
      do ib = 1,n_fuv_band
         write(75,'(A,A,4ES14.6)') '# ', fuv_band_name(ib), absen(ib),     &
                    heat_col(ib), bond_col(ib), fuv_band_flux(ib)
      enddo
      write(75,'(A)') '#'
      write(75,'(A)') '# (d) THE STATE. The ledger above sums the'//       &
                     ' absorbers the optical-depth columns record,'
      write(75,'(A)') '#     which is the state the photon field was'//    &
                     ' built on. The H2O density written in'
      write(75,'(A)') '#     Oxygen_chemistry.txt is the state at the'//   &
                     ' end of the sweep, and the two need not be'
      write(75,'(A)') '#     the same cell by cell. state_drift is the'//  &
                     ' largest |n_H2O dr / dN_H2O - 1| over the'
      write(75,'(A)') '#     column and the radius at which it stands;'//  &
                     ' it goes to zero as the sweep converges.'
      write(75,'(A,ES14.6,A,F10.6)') '# state_drift ', drift_worst,        &
                     '   at r[Rp] ', drift_r
      e_lw_in = fuv_band_flux(ib_LW)
      write(75,'(A)') '#'
      write(75,'(A)') '# The 912-1201 A band is shared: H2 absorbs it in'//&
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

      ! The closure of paragraph (a): the continuum rates against the same
      ! cells re-formed from the optical-depth columns. This is the number
      ! that must be round-off. The beam budget of paragraph (b) is a
      ! physical statement about each band and not a residual, so the two
      ! bands that carry an absorber without a rate are reported separately.
      write(*,'(a,es9.2)') ' (write_output/eq) FUV band ledger: max'//     &
           ' |continuum rate sum/beam loss - 1| = ', worst
      if (F_Lya_star .gt. 0.0d0)                                           &
         write(*,'(a,f7.4,a)') ' (write_output/eq) band B2: the H I'//     &
           ' Ly-alpha line scatters ', unrated_f(ib_B2), ' of that'//      &
           ' band out of the beam, which no photolysis rate carries'
      if (beam_loss(ib_LW) .gt. 0.0d0)                                     &
         write(*,'(a,es10.3,a)') ' (write_output/eq) band LW: the H2'//    &
           ' rates and the beam loss of the same lines close to ',         &
           unrated_f(ib_LW), ' of the beam loss'
      if (e_lw_in .gt. 0.0d0 .and. absen(ib_LW) .gt. e_lw_in*(1.0d0+1.0d-6))&
         write(*,'(a)') ' (write_output/eq) WARNING: the absorbers of'//   &
           ' the shared 912-1201 A band take more energy out of it'//      &
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
      if (lw_col_over_overlap .gt. 1.0d0)                                 &
         write(*,'(a,es9.2,a,f6.2,a)') ' (write_output/eq)'//              &
           ' WARNING: the star-ward H2 column reaches ',                   &
           maxval(NH2_col_lw), ' cm^-2, which is ',                        &
           lw_col_over_overlap, ' times the top of the column axis of'//  &
           ' the overlapping-line self-shielding table, above which its'//&
           ' edge value is returned rather than a calculated one, for'//   &
           ' the dissociation rate and the band share the H2 lines'//      &
           ' remove alike.'

      end subroutine write_oxygen_chemistry

      ! End of module
      end module output_write
