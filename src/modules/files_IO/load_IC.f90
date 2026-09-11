      module IC_load
      ! Module to load previous ICs, stored in the following files:
      ! - Hydro_ioniz_IC.txt
      ! - Ion_species_IC.txt
      !
      ! Two file generations are supported:
      !  * schema-2 files (written by write_output with '# ...' header
      !    lines): species columns are identified by the labels on the
      !    "# columns" line, so ALL species present in the file --
      !    including the metal ions -- are restored. Restarts therefore
      !    preserve the metal ionization state.
      !  * legacy headerless files: the original fixed 7-column read
      !    (r + H/He/HeITR); metals are initialized from the abundance.
      !
      ! An element whose stages are not carried by the file -- either no
      ! column at all, or columns that are identically zero because the run
      ! that wrote them had metals off -- is initialized from the abundance
      ! exactly as a cold start does, and the substitution is reported.
      !
      ! An element the file DOES carry, and whose reservoir a lower-atmosphere
      ! handoff states, has its loaded column renormalized onto that reservoir
      ! by one factor common to every ionization stage, so a closure iteration
      ! restarts on the El/H it just moved to. Everything else about the
      ! restart -- including an element the handoff is silent about -- is
      ! unchanged.

      use global_parameters
      use species_table, only: isp_HI, isp_HII, isp_HeI, isp_HeII,    &
                               isp_H2, isp_H2p, isp_H3p, isp_HeHp,     &
                               isp_HeIII, isp_HeTR,                    &
                               n_mion, n_melem, mion_fsp, mion_name,   &
                               mion_elem, melem_i0, melem_top,         &
                               melem_name, isp_OH, isp_H2O, isp_CO,     &
                               iel_O, iel_C
      use utils, only: calc_rho
      ! Which build produced a state file: the source revision stamped into
      ! the executable (the run never calls git).
      use build_stamp, only: build_git, build_dirty
      ! Oxygen-chemistry seed for a restart file written before the option
      ! existed (see the block near the end of load_IC).
      use oxygen_rates, only: co_equilibrium_density,                   &
                              oxygen_chemical_equilibrium_fractions

      implicit none

      ! HOW CLOSE A RESTART'S OWN He/H MUST BE TO THE INPUT'S before the
      ! loader carries the loaded hydrogen and helium onto the input
      ! composition. The input file is the authority on composition, so a
      ! state written at another He/H is rescaled (or, when it carries HeH+
      ! and helium is not diffusing, refused): the threshold separates a
      ! restart at the composition it was written with, whose reconstruction
      ! round trip is a few units in the last place, from a state whose
      ! composition really is another one. It is a round-off allowance and
      ! not a physical tolerance; a state that exceeds it differs in its
      ! element ratio, and the answer is the producing path or an explicit
      ! conservative mapping, not a larger number here.
      real*8, parameter :: heh_dev_tol = 1.0d-6

      ! Elements whose density had to be built from the abundance because the
      ! restart file did not carry it (see the metals block below). Allocated
      ! by load_IC, so a run that starts cold leaves it unallocated; the setup
      ! report tests for that.
      logical, allocatable :: melem_from_abundance(:)

   ! What the restart file says it was produced under, parsed from its
   ! '# coupling:' line (write_coupling_state_header, utilities.f90).
   ! ic_coupling_present is .false. for a file written before that line
   ! existed, and then every field below is meaningless and no caller may act
   ! on it -- such a file restarts exactly as it did before.
   logical :: ic_coupling_present   = .false.
   logical :: ic_sec_ion_active     = .false.
   integer :: ic_sec_ion_armed_step = -1
   real*8  :: ic_valve_eps          = -1.0d0
   character(len=8) :: ic_rec_method = ''
   ! Was the state in the restart file produced with the hydrogen ionization
   ! state carried? Absent from the header means no, which is what every
   ! file written before the option existed means as well.
   logical :: ic_ionization_transport = .false.
   real*8  :: ic_base_flux_const    = -1.0d0
   ! WHICH RUN STATE PRODUCED THE STATE IN THE FILE, and its clock
   ! (docs/a0_run_mode_contract_20260906.md sections 3 and 5). A file written
   ! before the field existed carries no statement, and the only thing that
   ! can be said about it is that no elapsed time was recorded with it, which
   ! is what run_mode_init means: ic_run_mode therefore starts at
   ! run_mode_init and ic_run_mode_present says whether the file stated it.
   integer :: ic_run_mode         = run_mode_init
   logical :: ic_run_mode_present = .false.
   real*8  :: ic_t_phys           = 0.0d0
   ! Did the header actually carry a readable, finite t_phys? A file that
   ! states mode=phys is claiming to stand on a trajectory, and a trajectory
   ! without the time it reached is not one a continuation can be started
   ! from; this separates "the field was read" from "the field defaulted".
   logical :: ic_t_phys_present   = .false.
   ! Was the state in the restart file certified as a stationary solution by
   ! the run that wrote it ('certified=' of the same line)?  An evaluation of
   ! that state re-answers exactly that question, so the claim the file makes
   ! is what its re-evaluation is held to.
   logical :: ic_certified          = .false.

   ! ------------------------------------------------------------------ !
   ! THE RESTART METADATA BLOCK (docs/restart_contract_design_20260909.md
   ! section 3).
   !
   ! WHAT IT IS FOR. The state files carry the state; what they carried until
   ! now about the RUN that produced it was the coupling line and the
   ! provenance line, neither of which states the composition reservoir, the
   ! physical grid, the constant set or the option set. A state is a solution
   ! OF a configuration: reloaded under another one it is not the same
   ! problem's state, and nothing in the files said so. The block states the
   ! configuration the state solves, one field per line, and the loader
   ! compares every field with the run's own.
   !
   ! WHAT IS COMPARED, AND HOW STRICTLY.
   !   reservoir         the elemental composition the state was solved at.
   !                     Compared within heh_dev_tol, the same round-off
   !                     allowance the He/H check uses: the file's value is
   !                     the printed decimal of the input's.
   !   grid, constants   compared EXACTLY, as the text of the field. A cell
   !                     center from one construction and a face from another
   !                     do not belong to one discretization, and a changed
   !                     Jupiter radius is a changed physical grid; a
   !                     conservative remap is a separate workflow.
   !   options           compared EXACTLY. The tokens are the switches that
   !                     decide WHICH EQUATIONS the state solves, so a
   !                     difference means the file's state is a solution of
   !                     another system.
   !   species_columns   reported, not refused: the loader maps species by
   !                     the label on the '# columns' line, and an element
   !                     the file does not carry is built from the abundance
   !                     and reported (the block above this one).
   !   source            reported. It identifies the executable, and the same
   !                     equations solved by two builds are the same
   !                     equations.
   !
   ! A file with no 'restart_schema' line is a file written before this block
   ! existed: it is loaded exactly as it was before, and the fact that
   ! nothing is known about the configuration it was produced under travels
   ! with the run -- printed on load, and written into the files this run
   ! produces (ic_provenance_unknown).
   integer, parameter :: restart_schema_version = 1
   ! One field per line, in this order, each line's first token being its
   ! name. The writer emits them, the loader stores them by name, so the two
   ! cannot disagree about which line is which field.
   integer, parameter :: n_meta = 8
   integer, parameter :: imeta_schema   = 1
   integer, parameter :: imeta_reservoir = 2
   integer, parameter :: imeta_species  = 3
   integer, parameter :: imeta_grid     = 4
   integer, parameter :: imeta_const    = 5
   integer, parameter :: imeta_options  = 6
   integer, parameter :: imeta_t_phys   = 7
   integer, parameter :: imeta_source   = 8
   character(len=16), parameter :: meta_tag(n_meta) = [ character(len=16) :: &
        'restart_schema', 'reservoir', 'species_columns', 'grid',            &
        'constants', 'options', 't_phys[s]', 'source' ]
   integer, parameter :: meta_len = 1024
   ! The block the restart files carried, one entry per field ('' = the file
   ! did not carry that field).
   character(len=meta_len) :: ic_meta(n_meta) = ''
   ! Did the restart pair carry the block at all? .false. is a file written
   ! before it existed, whose configuration is unknown.
   logical :: ic_restart_schema_present = .false.
   logical :: ic_provenance_unknown     = .false.

   ! ------------------------------------------------------------------ !
   ! THE TOKENS OF THE 'options' FIELD, AND WHICH OF THEM A RESTART MAY BE
   ! ALLOWED TO CHANGE (decision 21 of
   ! docs/To_be_determined_by_user_20260906.md, option a).
   !
   ! Comparing the whole field exactly is right for a state that claims to
   ! be stationary and wrong for the way this project reaches its
   ! solutions: converge without an option, restart with it on, converge
   ! again (the option ladder). The input key "Restart option change:" names
   ! the tokens that are ALLOWED to differ; every other difference is still
   ! refused by name, and the change that was allowed is written into the
   ! new state's block, so a rung of a ladder states what it was reached
   ! from instead of the change being silent.
   integer, parameter :: n_opt = 21
   character(len=16), parameter :: opt_name(n_opt) = [ character(len=16) :: &
        'He23S', 'metals', 'eos_metals', 'mol', 'molbase', 'oxychem',       &
        'carrier', 'carrier_newton', 'iontrans', 'he_diff',                 &
        'he_metal_diff', 'sec_ion', 'caloric_mono', 'excH', 'base_ir',      &
        'mol_ir', 'mol_heat', 'visc', 'cond', 'jlya', 'wellbal' ]
   ! WHICH TOKENS MAY NEVER BE NAMED, and why: these six decide HOW MANY
   ! UNKNOWNS the state has and which rows the state files carry. metals
   ! adds the metal ionization stages, mol the four molecular carriers,
   ! oxychem the three oxygen carriers; carrier makes the carriers
   ! transported and carrier_newton puts their rows in the Newton unknown
   ! vector; iontrans carries the hydrogen ionization state as a
   ! transported row. A state whose rows are not this run's rows is not
   ! this run's state, so such a change is a cold start and not a restart,
   ! whatever the input names. The grid is not a token of this field at all
   ! (it is the 'grid' field, N included) and is refused the same way.
   logical, parameter :: opt_changes_layout(n_opt) = [                      &
        .false., .true.,  .false., .true.,  .false., .true.,                &
        .true.,  .true.,  .true.,  .false.,                                 &
        .false., .false., .false., .false., .false.,                        &
        .false., .false., .false., .false., .false., .false. ]
   ! The tokens the input named as allowed to differ, set by input_read
   ! from "Restart option change:" (which is where an unknown token and a
   ! layout token are refused, the input file being what states them).
   logical :: restart_option_change_named(n_opt) = .false.
   logical :: restart_option_change_given = .false.
   ! THE STATE'S OPTION-CHANGE HISTORY, one line per restart at which a
   ! named option was allowed to differ: the line this run adds and the
   ! lines the loaded file already carried, in that order, so the rungs of
   ! a ladder can be read back from the last state it produced. The cap is
   ! a history length and not a physical limit; when it is reached the
   ! oldest line is replaced by one line stating how many were dropped, so
   ! a truncated history says that it is truncated.
   integer, parameter :: max_option_change = 32
   character(len=meta_len) :: ic_option_change(max_option_change) = ''
   integer :: n_ic_option_change      = 0
   integer :: n_option_change_dropped = 0
   ! Did the load actually change a named option, and did it find a named
   ! token that does not differ? Reported by the setup report.
   logical :: ic_option_change_applied = .false.
   logical :: ic_option_change_inert   = .false.

      contains

      subroutine load_IC(rho,v,p,T,f_sp,W)

      ! Integer variables
      integer :: j, k, ios, nlab, c, e, i0, im, nrec

      ! Loaded number densities for every f_sp column (zero = not in file)
      real*8, dimension(1-Ng:N+Ng,n_species) :: nsp_l
      ! Metal and molecular densities assembled for the calc_rho mass policy
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm_l
      real*8, dimension(1-Ng:N+Ng,4)      :: nmol_l
      ! Oxygen carriers of the loaded state (OH, H2O, CO); zero for a file
      ! written before the oxygen chemistry existed.
      real*8, dimension(1-Ng:N+Ng,3)      :: nox_l
      ! Oxygen-chemistry seeding of a restart file that predates the option.
      logical :: ox_seeded
      integer :: ic0
      real*8  :: T_K_l, nOtot, nCtot, nCO_l, nOfam, nOH_l, nH2O_l
      real*8  :: f_oh_l, f_h2o_l, sO_l, sC_l, dH_l
      real*8, dimension(1-Ng:N+Ng)        :: rho_dim
      ! Cell centers the restart file carries. The grid is a property of the
      ! run: define_grid builds r, r_edg, dr_j, j_min and j_flux from
      ! "Grid type", "Base grid", "Grid cells" and "Outer radius", and only
      ! r appears in the state file. Assigning the file's column into r would
      ! march centers from one construction against faces, widths and window
      ! indices from another, so the column is read here and compared with r
      ! instead (docs/development_plan_20260905_rev3.md section 10.2 item 3).
      real*8, dimension(1-Ng:N+Ng)        :: r_file
      ! Hydrogen nuclei density of the loaded state (free + bound in molecules)
      real*8, dimension(1-Ng:N+Ng)        :: nH_l
      ! Helium nuclei density of the loaded state, and the two factors that
      ! carry the loaded composition onto the input one
      real*8, dimension(1-Ng:N+Ng)        :: nHe_l, sH_l, sHe_l
      real*8 :: heh_loaded, heh_dev, gotH_l, gotHe_l
      ! The cell whose He/H is farthest from the input's, its ratio
      ! there and that cell's own departure: the check is a maximum
      ! over the physical cells and the messages name the cell it
      ! attained (report_deviating_heh_cell).
      integer :: j_heh
      real*8  :: heh_at_j, heh_dev_j
      logical :: col_present(n_species), elem_ok
      ! Auxiliary temporary variable
      real*8 :: tmp
      ! Loaded El/H at the base cell, and the single factor that carries the
      ! loaded metal column onto a handoff reservoir (metals block below)
      real*8 :: elh_l, r_el
      real*8 :: vals(80)

      character(len=8192) :: line
      ! The restart metadata block of each of the two state files, stored by
      ! field name as it is read (parse_restart_metadata_line) and compared
      ! with this run's own configuration afterwards.
      character(len=meta_len) :: meta_h(n_meta), meta_i(n_meta)
      character(len=16)   :: labels(80)
      integer :: col2fsp(80)
      logical :: has_header

      ! Output variables
      real*8, dimension(1-Ng:N+Ng),   intent(out) :: rho,v,p,T
      real*8, dimension(1-Ng:N+Ng,n_species), intent(out) :: f_sp
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: W


	   !-------------------------------------!

      ! The IC files carry one record per cell, N + 2*Ng rows. Since N is a
      ! runtime value ("Grid cells:"), a restart file written at a different
      ! N would otherwise die below with a bare end-of-file error; count the
      ! data records first and state the actual mismatch.
      nrec = 0
      ic_coupling_present = .false.
      meta_h = ''
      meta_i = ''
      ! The loaded state's own option-change history starts empty and is
      ! filled from the file, so this run's line is appended after it.
      ic_option_change        = ''
      n_ic_option_change      = 0
      n_option_change_dropped = 0
      open(unit = 1, file = 'output/Hydro_ioniz_IC.txt')
      do
         read(1,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (is_comment(line)) then
            if (index(line,'coupling:') .gt. 0) call parse_coupling_header(line)
            call parse_restart_metadata_line(line, meta_h)
            call collect_option_change_line(line)
         else
            nrec = nrec + 1
         endif
      enddo
      close(1)
      ! WHAT THE STATE WAS PRODUCED UNDER, AGAINST WHAT THIS RUN WILL DO TO
      ! IT. The H+ column of a file written with the ionization state
      ! carried is a transported quantity; restarted with the option off,
      ! the first sweep replaces it by the local root of its own cell, which
      ! for the wind this option exists for is a factor of several. That is
      ! a legitimate thing to ask for and it is not refused -- but it is not
      ! something to discover afterwards from the profile.
      if (ic_ionization_transport .and. .not. ionization_transport) then
         write(*,*) ' (load_IC) NOTE: this state was written with'
         write(*,*) '   "Ionization transport: True" and is being restarted'
         write(*,*) '   without it. Its H+ column is a transported'
         write(*,*) '   ionization state and the first sweep will replace'
         write(*,*) '   it by the local equilibrium of each cell.'
      endif
      if (ionization_transport .and. ic_coupling_present .and.               &
          .not. ic_ionization_transport) then
         write(*,*) ' (load_IC) NOTE: "Ionization transport: True", but the'
         write(*,*) '   state being restarted was produced without it, so'
         write(*,*) '   its H+ column is a local equilibrium. The transport'
         write(*,*) '   starts from it and relaxes over the ionization time.'
      endif
      if (nrec .ne. N + 2*Ng) then
         write(*,'(A,I0,A,I0,A)')                                            &
            ' (load_IC) ERROR: output/Hydro_ioniz_IC.txt has ', nrec,        &
            ' data rows, but the grid needs N + 2*Ng = ', N + 2*Ng,          &
            ' (the IC was written at a different "Grid cells:" N,'//         &
            ' or the file is truncated).'
         error stop 1
      endif

      ! Load thermodynamic variables (skip any '#' header lines)
      open(unit = 1, file = 'output/Hydro_ioniz_IC.txt')
         read(1,'(A)') line
         do while (is_comment(line))
            read(1,'(A)') line
         enddo
         read(line,*) r_file(1-Ng), tmp, v(1-Ng), p(1-Ng), T(1-Ng),      &
                      tmp, tmp
         do j = 2-Ng,N+Ng
            read(1,*) r_file(j), tmp, v(j), p(j), T(j), tmp, tmp
         enddo
      close(1)
      call verify_restart_radii(r_file, 'output/Hydro_ioniz_IC.txt')

      ! Adimensionalize
      v = v/v0
      p = p/p0
      T = T/T0

      ! Load ionization profiles
      nsp_l       = 0.0d0
      col_present = .false.

      open(unit = 2, file = 'output/Ion_species_IC.txt')
      read(2,'(A)') line
      has_header = is_comment(line)

      if (has_header) then
         ! ---- schema-2 file: map columns by label ----
         nlab = 0
         do while (is_comment(line))
            ! THE FIELD OF A COMMENT LINE IS ITS FIRST TOKEN, not a
            ! substring of it: 'species_columns' and the NOTE prose of an
            ! advection-corrected file both contain 'columns', and either
            ! read as a label line would replace the species labels with
            ! whatever words followed. (Found while adding the metadata
            ! block, which carries such a field name.)
            if (comment_field_is(line, 'columns'))                        &
               call parse_labels(line, labels, nlab)
            call parse_restart_metadata_line(line, meta_i)
            read(2,'(A)',iostat=ios) line
            if (ios .ne. 0) exit
         enddo
         if (nlab .lt. 2) error stop '(load_IC) header found but no "# columns" line'

         ! Build the label -> f_sp column map (0 = ignore). labels(1) is r.
         col2fsp = 0
         do k = 2, nlab
            col2fsp(k) = species_column(labels(k))
            if (col2fsp(k) .gt. 0) col_present(col2fsp(k)) = .true.
         enddo

         ! First data record is already in 'line'
         read(line,*) (vals(k), k = 1,nlab)
         call scatter_row(1-Ng, vals, col2fsp, nlab, nsp_l, r_file)
         do j = 2-Ng,N+Ng
            read(2,*) (vals(k), k = 1,nlab)
            call scatter_row(j, vals, col2fsp, nlab, nsp_l, r_file)
         enddo

      else
         ! ---- legacy headerless file: fixed 7-column layout ----
         rewind(2)
         do j = 1-Ng,N+Ng
            read(2,*) r_file(j),             &
                      nsp_l(j,isp_HI),       &
                      nsp_l(j,isp_HII),      &
                      nsp_l(j,isp_HeI),      &
                      nsp_l(j,isp_HeII),     &
                      nsp_l(j,isp_HeIII),    &
                      nsp_l(j,isp_HeTR)
         enddo
         col_present(isp_HI:isp_HeTR) = .true.
      endif
      close(2)
      call verify_restart_radii(r_file, 'output/Ion_species_IC.txt')

      ! THE CONFIGURATION THE STATE IS A STATE OF, before anything is done
      ! with the state itself.
      call verify_restart_metadata(meta_h, meta_i)

      ! Zero the whole f_sp before any assignment so columns not restored
      ! below -- in particular the molecular columns, absent from most IC
      ! files -- are defined.
      f_sp = 0.0d0

      ! ---- carry the loaded H/He onto the input composition ----
      ! The file holds absolute species densities, so a state written by a run
      ! at a different "He/H number ratio" would otherwise be used as it
      ! stands: the restart would run the donor's composition while the setup
      ! report echoes the input one, and nothing downstream would notice --
      ! measured, a He/H = 10 restart seeded from a He/H = 1 solution converged
      ! back to He/H = 1 everywhere. The input file is the authority on
      ! composition, so the loaded hydrogen and helium are rescaled to it,
      ! holding each cell's H + He nuclei count and each element's
      ! ionization-stage split; the trace metals are defined per hydrogen
      ! nucleus and so follow hydrogen. A restart at the composition it was
      ! written with leaves every density untouched.
      !
      ! HeH+ carries one nucleus of each, so a single factor cannot set both
      ! counts; a molecular state whose composition disagrees with the input
      ! is refused rather than approximated.
      if (thereis_He) then
         nH_l  = nsp_l(:,isp_HI)  + nsp_l(:,isp_HII)                      &
               + 2.0d0*(nsp_l(:,isp_H2) + nsp_l(:,isp_H2p))               &
               + 3.0d0*nsp_l(:,isp_H3p) + nsp_l(:,isp_HeHp)               &
               + nsp_l(:,isp_OH) + 2.0d0*nsp_l(:,isp_H2O)
         ! The HeI column of an output file is the TOTAL He I density, He
         ! 2^3S included, so the triplet column is not added again here.
         nHe_l = nsp_l(:,isp_HeI) + nsp_l(:,isp_HeII)                     &
               + nsp_l(:,isp_HeIII) + nsp_l(:,isp_HeHp)
         ! OVER THE PHYSICAL CELLS ONLY.  The ghosts are boundary data: the
         ! lower pair is the inflow reservoir, whose molecular partition the
         ! handoff imposes and the sweep re-pins every step, and the upper
         ! pair mirrors the top cell.  They are rebuilt on the first step of
         ! the restart, and the input's He/H is the authority on the
         ! atmosphere, not on them.  How this was found (2026-09-05): a
         ! carrier-transport snapshot carried its lower ghost with He/H 2.1e-5
         ! off the input while every physical cell agreed to 1e-6, and that
         ! ghost alone sent the restart into the rescale branch, which the
         ! HeH+ clause below then refuses.  The ghost itself was a writer
         ! defect -- the molecular columns came from the sweep's arrays while
         ! the atomic ones were f_sp*rho at the write, and the ghost's rho
         ! had moved in between (fixed the same day, docs/Update_EXHALE_stage1.md
         ! section 169; files written before it still carry such ghosts).
         ! The range stays physical either way: a check on the input's He/H
         ! has nothing to say about boundary data.
         ! THE CELL THE CHECK IS ABOUT IS THE ONE FARTHEST FROM THE INPUT,
         ! and it is the cell the messages below name. The quantity that
         ! decides is heh_dev, the largest relative departure over the
         ! physical cells; a molecular restart whose He/H varies with radius
         ! can agree with the input at the top cell to four digits and
         ! disagree by percent at a front, so a message naming any other cell
         ! states that the file was written at the composition the input asks
         ! for and then refuses it (MEASURED on a coupled-carrier state: top
         ! cell 7.9307E-02 against the input 7.9307E-02, largest departure
         ! 3.96E-02 relative at cell 199, r = 1.0906).
         heh_dev  = 0.0d0
         j_heh    = 0
         heh_at_j = 0.0d0
         do j = 1, N
            if (nH_l(j) .gt. 0.0d0) then
               heh_loaded = nHe_l(j)/nH_l(j)
               heh_dev_j  = abs(heh_loaded - HeH)/max(HeH,1.0d-30)
               if (heh_dev_j .gt. heh_dev) then
                  heh_dev  = heh_dev_j
                  j_heh    = j
                  heh_at_j = heh_loaded
               endif
            endif
         enddo
         if (heh_dev .gt. heh_dev_tol) then
            if (maxval(nHe_l) .le. 0.0d0) then
               write(*,'(A,ES11.4,A)')                                    &
                  ' (load_IC) ERROR: the input asks for He/H =', HeH,     &
                  ' but the restart file carries no helium at all;'//     &
                  ' there is nothing to rescale. Start this composition'//&
                  ' cold, or restart from a helium-bearing state.'
               error stop 1
            endif
            if (maxval(abs(nsp_l(:,isp_HeHp))) .gt. 0.0d0 .and.            &
                .not. he_diffusion) then
               write(*,'(A)')                                             &
                  ' (load_IC) ERROR: the restart file does not carry'//   &
                  ' the He/H the input asks for, and the state carries'
               write(*,'(A)') '   HeH+, whose nucleus of each element'//  &
                  ' cannot be rescaled by one factor. Restart a'//        &
                  ' molecular state at'
               write(*,'(A)') '   its own composition.'
               call report_deviating_heh_cell(j_heh, heh_at_j, heh_dev)
               error stop 1
            endif
            sH_l  = (nH_l + nHe_l)/(1.0d0 + HeH)/max(nH_l, 1.0d-30)
            sHe_l = HeH*(nH_l + nHe_l)/(1.0d0 + HeH)/max(nHe_l, 1.0d-30)
            ! With He_diffusion the cell-by-cell element split IS the state
            ! being restarted: the diffusion operator produced it, and a
            ! separated profile is the physics, not a defect of the file.
            ! Only the base cells are set to the reservoir composition --
            ! that is where the operator itself holds a Dirichlet HeH -- and
            ! every cell above keeps the helium fraction it was written with.
            ! Without the flag the input file remains the authority on the
            ! composition and the whole column is rescaled, as before.
            if (he_diffusion) then
               do j = 2, N+Ng
                  sH_l(j)  = 1.0d0
                  sHe_l(j) = 1.0d0
               enddo
            endif
            nsp_l(:,isp_HI)    = nsp_l(:,isp_HI)   *sH_l
            nsp_l(:,isp_HII)   = nsp_l(:,isp_HII)  *sH_l
            nsp_l(:,isp_H2)    = nsp_l(:,isp_H2)   *sH_l
            nsp_l(:,isp_H2p)   = nsp_l(:,isp_H2p)  *sH_l
            nsp_l(:,isp_H3p)   = nsp_l(:,isp_H3p)  *sH_l
            ! OH and H2O carry only H among the two elements being rescaled
            ! (their oxygen follows hydrogen through melem_ab, like every
            ! trace metal), so they scale with the hydrogen factor.
            nsp_l(:,isp_OH)    = nsp_l(:,isp_OH)   *sH_l
            nsp_l(:,isp_H2O)   = nsp_l(:,isp_H2O)  *sH_l
            nsp_l(:,isp_CO)    = nsp_l(:,isp_CO)   *sH_l
            nsp_l(:,isp_HeI)   = nsp_l(:,isp_HeI)  *sHe_l
            nsp_l(:,isp_HeII)  = nsp_l(:,isp_HeII) *sHe_l
            nsp_l(:,isp_HeIII) = nsp_l(:,isp_HeIII)*sHe_l
            nsp_l(:,isp_HeTR)  = nsp_l(:,isp_HeTR) *sHe_l
            do im = 1, n_mion
               nsp_l(:,mion_fsp(im)) = nsp_l(:,mion_fsp(im))*sH_l
            enddo
            ! HeH+ carries a nucleus of each element, so no single factor can
            ! rescale it. With He_diffusion the only cells rescaled at all are
            ! the base and its inner ghosts -- the column above keeps the
            ! diffused split -- so those cells are projected onto their two
            ! element totals exactly as binary_element_diffusion writes back:
            ! HeH+ takes the smaller of the two factors and the nuclei that
            ! leaves short are deposited into the neutral ground species.
            if (he_diffusion) then
               do j = 1-Ng, 1
                  if (nsp_l(j,isp_HeHp) .le. 0.0d0) cycle
                  nsp_l(j,isp_HeHp) = nsp_l(j,isp_HeHp)                    &
                                      *min(sH_l(j), sHe_l(j))
                  gotH_l  = nsp_l(j,isp_HI) + nsp_l(j,isp_HII)             &
                          + 2.0d0*(nsp_l(j,isp_H2) + nsp_l(j,isp_H2p))     &
                          + 3.0d0*nsp_l(j,isp_H3p) + nsp_l(j,isp_HeHp)     &
                          + nsp_l(j,isp_OH) + 2.0d0*nsp_l(j,isp_H2O)
                  gotHe_l = nsp_l(j,isp_HeI) + nsp_l(j,isp_HeII)           &
                          + nsp_l(j,isp_HeIII) + nsp_l(j,isp_HeHp)
                  if (nH_l(j)*sH_l(j) .gt. gotH_l)                         &
                     nsp_l(j,isp_HI)  = nsp_l(j,isp_HI)                    &
                                        + (nH_l(j)*sH_l(j) - gotH_l)
                  if (nHe_l(j)*sHe_l(j) .gt. gotHe_l)                      &
                     nsp_l(j,isp_HeI) = nsp_l(j,isp_HeI)                   &
                                        + (nHe_l(j)*sHe_l(j) - gotHe_l)
               enddo
            endif
            if (he_diffusion) then
               write(*,'(A,ES11.4,A,ES11.4,A)')                           &
                  ' (load_IC) He_diffusion: the diffused He/H profile of'//&
                  ' the restart file is kept (top cell He/H =',           &
                  nHe_l(N)/nH_l(N), '); base cells set to the input He/H =',&
                  HeH, '.'
            else
               write(*,'(A,ES13.6,A)')                                    &
                  ' (load_IC) restart file rescaled to the input He/H =', &
                  HeH, '; its own He/H at the cell farthest from that:'
               call report_deviating_heh_cell(j_heh, heh_at_j, heh_dev)
            endif
         endif
      endif

      ! ---- metals the file does not carry: build them from the abundance ----
      ! Two ways a restart file can fail to carry an element: it has no column
      ! for some ionization stage (legacy headerless files, or a file written
      ! by a run with fewer elements), or it has the columns but they are
      ! identically zero because the run that wrote them had metals OFF. The
      ! schema-2 writer emits the metal columns unconditionally, so the column
      ! test alone accepts a metals-off file and the restart then runs with
      ! zero metal density everywhere -- while the base boundary condition
      ! still counts metals in the mass and particle budget (eos_metals),
      ! which leaves an O(1) residual in the first cells.
      !
      ! Such an element is initialized exactly as a cold start does (set_IC):
      ! n_X = melem_ab(e) * n_H with all of it in the neutral stage. It is put
      ! into the loaded density array BEFORE the mass density is reconstructed,
      ! so calc_rho sees it and the restart keeps the mass closure
      ! sum_i f_i A_i = 1 that the cold start has by construction.
      nH_l = nsp_l(:,isp_HI)  + nsp_l(:,isp_HII)                         &
           + 2.0d0*(nsp_l(:,isp_H2) + nsp_l(:,isp_H2p))                  &
           + 3.0d0*nsp_l(:,isp_H3p) + nsp_l(:,isp_HeHp)                  &
           + nsp_l(:,isp_OH) + 2.0d0*nsp_l(:,isp_H2O)
      if (.not. allocated(melem_from_abundance))                         &
         allocate(melem_from_abundance(n_melem))
      melem_from_abundance = .false.
      do e = 1, n_melem
         i0 = melem_i0(e)
         elem_ok = .true.
         do k = 0, melem_top(e)
            if (.not. col_present(mion_fsp(i0+k))) elem_ok = .false.
         enddo
         if (elem_ok) then
            tmp = 0.0d0
            do k = 0, melem_top(e)
               tmp = tmp + sum(abs(nsp_l(1:N,mion_fsp(i0+k))))
            enddo
            ! An oxygen-chemistry state can hold nearly all of its oxygen and
            ! all of its carbon in OH, H2O and CO, with the ion stages at
            ! zero. Counting only the stages would call the element absent
            ! and rebuild it from the abundance, doubling it.
            if (e .eq. iel_O) tmp = tmp                                    &
               + sum(abs(nsp_l(1:N,isp_OH))) + sum(abs(nsp_l(1:N,isp_H2O))) &
               + sum(abs(nsp_l(1:N,isp_CO)))
            if (e .eq. iel_C) tmp = tmp + sum(abs(nsp_l(1:N,isp_CO)))
            if (tmp .le. 0.0d0) elem_ok = .false.
         endif
         if (elem_ok) cycle
         do k = 0, melem_top(e)
            c = mion_fsp(i0+k)
            if (k .eq. 0) then
               nsp_l(:,c) = melem_ab(e)*nH_l
            else
               nsp_l(:,c) = 0.0d0
            endif
         enddo
         if (melem_ab(e) .gt. 0.0d0) then
            melem_from_abundance(e) = .true.
            write(*,'(A)') ' (load_IC) WARNING: the restart file carries no'// &
                 ' density for element '//trim(melem_name(e))//               &
                 '; initializing it as neutral at the input abundance.'
         endif
      enddo

      ! ---- elements the handoff states: carry the loaded column onto that
      ! reservoir --------------------------------------------------------
      ! A restart file carries the metal densities of the state it was written
      ! from. When the reservoir itself has moved since -- which is what a
      ! flux-closure iteration does to the elemental ratios of the lower
      ! atmosphere -- the base boundary condition uses the new El/H while the
      ! loaded column above the base still holds the old one, and the elemental
      ! budget n_El/n_H = (El/H)_resolved cannot close above the first cell.
      !
      ! Only an element the handoff itself states is touched: melem_from_handoff
      ! is set by set_element_abundance, the one door the "<El>_H_base" keys of
      ! base.inp and the elemental ratios of a "Lower atmosphere profile:" both
      ! go through. An abundance that came from metals.inp alone, and every
      ! restart with no handoff at all, leaves this loop doing nothing, so such
      ! a restart is bit-for-bit the one the previous code produced.
      !
      ! The whole column of the element is multiplied by the single factor
      !    r_El = (handoff El/H) / (loaded El/H at the base cell),
      ! the same factor for every ionization stage, so the loaded ionization
      ! split and the settling shape of the profile survive untouched and only
      ! the reservoir normalization changes. El/H at the base cell is counted in
      ! nuclei, summing the element over its stages against the hydrogen nuclei
      ! of nH_l (free plus the H bound in H2, H2+, H3+ and HeH+), which is the
      ! convention of element_ratio_HeH and of src/utils/element_budget.py.
      ! Helium keeps its own convention -- the base cells are set to the input
      ! He/H, the column above may carry a diffused split -- and is not touched
      ! here.
      ! (the flag is allocated by input_read; the test keeps the unit tests of
      !  src/tests, which build a state without it, on the untouched path)
      do e = 1, n_melem
         if (.not. allocated(melem_from_handoff)) exit
         if (.not. melem_from_handoff(e))  cycle
         ! Rebuilt from the abundance just above: it already IS the reservoir.
         if (melem_from_abundance(e))      cycle
         if (melem_ab(e) .le. 0.0d0)       cycle
         if (nH_l(1)     .le. 0.0d0)       cycle
         i0 = melem_i0(e)
         tmp = 0.0d0
         do k = 0, melem_top(e)
            tmp = tmp + nsp_l(1,mion_fsp(i0+k))
         enddo
         if (tmp .le. 0.0d0) cycle
         elh_l = tmp/nH_l(1)
         r_el  = melem_ab(e)/elh_l
         do k = 0, melem_top(e)
            c = mion_fsp(i0+k)
            nsp_l(:,c) = nsp_l(:,c)*r_el
         enddo
         ! Enough digits that the factor can be compared with the reservoir
         ! change it is supposed to equal: a closure iteration moves El/H by
         ! parts in 1e4, which ES10.3 would print as 1.000E+00.
         write(*,'(A,ES13.6,A,ES13.6,A,ES13.6)')                           &
            ' (load_IC) '//trim(melem_name(e))//'/H at the base cell:'//   &
            ' restart file ', elh_l, ', handoff ', melem_ab(e),            &
            '; column rescaled by ', r_el
      enddo

      ! ---- an oxygen-chemistry run restarted from a state without it ----
      ! A pre-oxygen-chemistry file has no OH / H2O / CO columns, so the
      ! loader leaves them at zero. Starting the option from zero is not
      ! neutral: zero is itself a root of the water cycle, and hybrd1 is
      ! already known to be bistable from a zero molecular seed
      ! (ionization_equilibrium), so the run would sit in the empty basin.
      ! The state is therefore seeded with the CHEMICAL EQUILIBRIUM of the
      ! loaded (T, H2/H) and the loaded oxygen and carbon totals -- the same
      ! partition the molecular-basin retry of the cell solve uses -- and the
      ! seeding is printed, never silent.
      !
      ! Both elements stay conserved: what goes into the carriers is taken
      ! out of the ion stages of the SAME element, in proportion, and the H
      ! nuclei the carriers hold are taken out of atomic H (and out of H2 if
      ! atomic H runs short). A file that already carries the columns is left
      ! alone.
      if (thereis_oxychem) then
         ox_seeded = (sum(abs(nsp_l(1:N,isp_OH)))                          &
                    + sum(abs(nsp_l(1:N,isp_H2O)))                         &
                    + sum(abs(nsp_l(1:N,isp_CO)))) .le. 0.0d0
         if (ox_seeded) then
            i0 = melem_i0(iel_O)
            ic0 = melem_i0(iel_C)
            do j = 1-Ng, N+Ng
               T_K_l  = T(j)*T0
               nOtot  = nsp_l(j,mion_fsp(i0))   + nsp_l(j,mion_fsp(i0+1))  &
                      + nsp_l(j,mion_fsp(i0+2))
               nCtot  = nsp_l(j,mion_fsp(ic0))  + nsp_l(j,mion_fsp(ic0+1)) &
                      + nsp_l(j,mion_fsp(ic0+2))
               if (nOtot .le. 0.0d0) cycle
               nCO_l  = co_equilibrium_density(nCtot, nOtot, T_K_l)
               nOfam  = nOtot - nCO_l
               call oxygen_chemical_equilibrium_fractions(T_K_l,           &
                       nsp_l(j,isp_H2), nsp_l(j,isp_HI), f_oh_l, f_h2o_l)
               nOH_l  = f_oh_l *nOfam
               nH2O_l = f_h2o_l*nOfam
               ! Take the carriers out of the ion stages of their element.
               sO_l = max(nOtot - nOH_l - nH2O_l - nCO_l, 0.0d0)/nOtot
               do k = 0, melem_top(iel_O)
                  nsp_l(j,mion_fsp(i0+k)) = nsp_l(j,mion_fsp(i0+k))*sO_l
               enddo
               if (nCtot .gt. 0.0d0) then
                  sC_l = max(nCtot - nCO_l, 0.0d0)/nCtot
                  do k = 0, melem_top(iel_C)
                     nsp_l(j,mion_fsp(ic0+k)) = nsp_l(j,mion_fsp(ic0+k))   &
                                                *sC_l
                  enddo
               endif
               ! Take the H nuclei of the carriers out of atomic H, then out
               ! of H2 if atomic H cannot supply them.
               dH_l = nOH_l + 2.0d0*nH2O_l
               if (nsp_l(j,isp_HI) .ge. dH_l) then
                  nsp_l(j,isp_HI) = nsp_l(j,isp_HI) - dH_l
               else
                  dH_l = dH_l - nsp_l(j,isp_HI)
                  nsp_l(j,isp_HI) = 0.0d0
                  nsp_l(j,isp_H2) = max(nsp_l(j,isp_H2) - 0.5d0*dH_l,      &
                                        0.0d0)
               endif
               nsp_l(j,isp_OH)  = nOH_l
               nsp_l(j,isp_H2O) = nH2O_l
               nsp_l(j,isp_CO)  = nCO_l
            enddo
            write(*,'(A)') ' (load_IC) the restart file carries no OH /'// &
               ' H2O / CO columns; seeding them from the chemical'
            write(*,'(A)') '   equilibrium of the loaded (T, H2/H) and'//  &
               ' the loaded oxygen and carbon totals.'
         endif
      endif

      ! Reconstruct the mass density (adimensional) from the LOADED densities
      ! with the SAME mass policy as the run (calc_rho): the trace-metal mass
      ! under the eos_metals policy and the molecular mass are included,
      ! exactly as calc_rho does in the main loop. He 2^3S is NOT a separate
      ! mass term -- it is an excited level of He I, whose density (nheiTR)
      ! is already inside the He I column calc_rho sums (bsp_is_excited_level,
      ! section 75). A restart therefore preserves the conserved mass exactly.
      ! The old H/He-only formula (rho = (nHI+nHII+4*(nHeI+nHeII+nHeIII))/n0)
      ! is gone deliberately: it dropped the metal/molecular mass and left a
      ! mass discontinuity on reload. Columns absent from the file are zero in
      ! nsp_l -- except the metals the block above rebuilt from the abundance
      ! -- so they add nothing (calc_rho honors thereis_He and
      ! eos_include_metals .and. thereis_metals).
      do im = 1, n_mion
         nm_l(:,im) = nsp_l(:,mion_fsp(im))
      enddo
      nmol_l(:,1) = nsp_l(:,isp_H2)
      nmol_l(:,2) = nsp_l(:,isp_H2p)
      nmol_l(:,3) = nsp_l(:,isp_H3p)
      nmol_l(:,4) = nsp_l(:,isp_HeHp)
      nox_l(:,1)  = nsp_l(:,isp_OH)
      nox_l(:,2)  = nsp_l(:,isp_H2O)
      nox_l(:,3)  = nsp_l(:,isp_CO)
      call calc_rho(nsp_l(:,isp_HI),  nsp_l(:,isp_HII),   nsp_l(:,isp_HeI),   &
                    nsp_l(:,isp_HeII), nsp_l(:,isp_HeIII),                    &
                    rho_dim, nm_l, nmol_l, nox_l)
      rho = rho_dim/n0

      ! H/He(+HeITR) and molecular fractions (f = n/(rho*n0)).
      f_sp(:,isp_HI)    = nsp_l(:,isp_HI)/(rho*n0)
      f_sp(:,isp_HII)   = nsp_l(:,isp_HII)/(rho*n0)
      f_sp(:,isp_HeI)   = nsp_l(:,isp_HeI)/(rho*n0)
      f_sp(:,isp_HeII)  = nsp_l(:,isp_HeII)/(rho*n0)
      f_sp(:,isp_HeIII) = nsp_l(:,isp_HeIII)/(rho*n0)
      f_sp(:,isp_HeTR)  = nsp_l(:,isp_HeTR)/(rho*n0)
      f_sp(:,isp_H2)    = nsp_l(:,isp_H2)/(rho*n0)
      f_sp(:,isp_H2p)   = nsp_l(:,isp_H2p)/(rho*n0)
      f_sp(:,isp_H3p)   = nsp_l(:,isp_H3p)/(rho*n0)
      f_sp(:,isp_HeHp)  = nsp_l(:,isp_HeHp)/(rho*n0)
      f_sp(:,isp_OH)    = nsp_l(:,isp_OH)/(rho*n0)
      f_sp(:,isp_H2O)   = nsp_l(:,isp_H2O)/(rho*n0)
      f_sp(:,isp_CO)    = nsp_l(:,isp_CO)/(rho*n0)

      ! Metal fractions. One rule for every element: the density array holds
      ! either the loaded state or the abundance-built one from the block
      ! above, and the fraction is n/(rho*n0) in both cases. For an element
      ! built from the abundance this is f_X = melem_ab/mass_per_H whenever the
      ! loaded H/He respects the input He/H ratio, i.e. the cold-start value.
      do e = 1, n_melem
         i0 = melem_i0(e)
         do k = 0, melem_top(e)
            c = mion_fsp(i0+k)
            f_sp(:,c) = nsp_l(:,c)/(rho*n0)
         enddo
      enddo

      ! Construct matrix of primitive profiles
      W(1,:) = rho
      W(2,:) = v
      W(3,:) = p

      ! End of subroutine
      end subroutine load_IC

      !-------------------------------------!

      !-------------------------------------!

      subroutine build_restart_metadata(blk)
      ! THE CONFIGURATION THIS RUN'S STATE IS A STATE OF, one field per
      ! entry, in the order of meta_tag. Every number is written with
      ! ES23.16, which round trips a double exactly, so a field compared as
      ! text is compared as the number it stands for.
      character(len=meta_len), intent(out) :: blk(n_meta)
      character(len=meta_len) :: s
      integer :: e, i, ncol
      blk = ''

      write(blk(imeta_schema),'(A,I0)') trim(meta_tag(imeta_schema))//' ', &
           restart_schema_version

      ! The elemental reservoir the state was solved at: He/H always, and
      ! every trace element the run carries a reservoir for. An element with
      ! zero abundance is absent from the line, so a metals-off state and a
      ! metals-on one do not carry the same reservoir field.
      s = trim(meta_tag(imeta_reservoir))//' He/H '//trim(meta_num(HeH))
      if (thereis_metals .and. allocated(melem_ab)) then
         do e = 1, n_melem
            if (melem_ab(e) .le. 0.0d0) cycle
            s = trim(s)//' '//trim(melem_name(e))//'/H '//                 &
                trim(meta_num(melem_ab(e)))
         enddo
      endif
      blk(imeta_reservoir) = s

      ! The species the state carries, BUILT FROM THE SPECIES TABLE in the
      ! order write_output puts them on the '# columns' line: one table, so a
      ! species added to it appears in both lines or in neither. Reported by
      ! the loader and never a refusal -- the loader maps species by label,
      ! and an element the file does not carry is built from the abundance.
      s = 'r'
      ncol = 1
      s = trim(s)//' HI HII HeI HeII HeIII HeITR'
      ncol = ncol + 6
      do i = 1, n_mion
         s = trim(s)//' '//trim(mion_name(i))
         ncol = ncol + 1
      enddo
      if (thereis_mol) then
         s = trim(s)//' H2 H2p H3p HeHp'
         ncol = ncol + 4
      endif
      if (thereis_oxychem) then
         s = trim(s)//' OH H2O CO'
         ncol = ncol + 3
      endif
      write(blk(imeta_species),'(A,I0,A)')                                 &
           trim(meta_tag(imeta_species))//' ', ncol, ' '//trim(s)

      ! The physical grid: the cell count, the length scale the radii are in
      ! units of, the first and last physical cell center, and the
      ! construction. Compared exactly, because centers from one
      ! construction and faces, widths and window indices from another are
      ! not one discretization.
      write(blk(imeta_grid),'(A,I0,A)') trim(meta_tag(imeta_grid))//' N ', &
           N, ' R0[cm] '//trim(meta_num(R0))//                                   &
           ' r_min[Rp] '//trim(meta_num(r(1)))//                                 &
           ' r_max[Rp] '//trim(meta_num(r(N)))//                                 &
           ' mode '//trim(grid_type)

      ! The constants the state was solved with. The set is named so that
      ! one token identifies a revision of parameters.f90's constants, and
      ! the two that scale the problem (the radius unit through R0 and the
      ! Boltzmann constant through every temperature) are written out.
      blk(imeta_const) = trim(meta_tag(imeta_const))//                     &
           ' set IAU2015+CODATA2018 RJ[cm] '//trim(meta_num(RJ))//               &
           ' kB[erg/K] '//trim(meta_num(kb_erg))

      ! WHICH EQUATIONS THE STATE SOLVES. Every token here is a switch that
      ! adds or removes an equation, an unknown or a source term of the
      ! system, so a difference means the file's state solves another
      ! system and is refused unless the input named that token (decision
      ! 21). The numeric inputs (fluxes, tolerances, abundance-independent
      ! rates) are NOT here: they change the coefficients of the same
      ! system, and a state carried across such a change is a legitimate
      ! starting point whose residual simply is not zero.
      ! The names come from opt_name and the values from opt_value, so the
      ! vocabulary the comparison and the input key use is the vocabulary
      ! the field is written with.
      blk(imeta_options) = trim(meta_tag(imeta_options))
      do i = 1, n_opt
         blk(imeta_options) = trim(blk(imeta_options))//' '//              &
              trim(opt_name(i))//'='//trim(opt_value(opt_name(i)))
      enddo

      ! The clock of the state. Zero, and meaningless, for a state written
      ! by an initialization or continuation run: the '# coupling:' line's
      ! mode field says which of the two this is.
      blk(imeta_t_phys) = trim(meta_tag(imeta_t_phys))//' '//              &
           trim(meta_num(merge(t_phys, 0.0d0, run_mode .eq. run_mode_phys)))

      ! What produced the file. Informational: two builds of one source
      ! solve the same equations. The last token is written only when this
      ! run's own initial state came from a restart file that carried no
      ! metadata block at all, so that a state descended from an unknown
      ! configuration says so in every file it appears in.
      blk(imeta_source) = trim(meta_tag(imeta_source))//' git='//          &
           trim(build_git)//' tree='//trim(build_dirty)
      if (ic_provenance_unknown) blk(imeta_source) =                       &
           trim(blk(imeta_source))//' restart_input=provenance_unknown'
      end subroutine build_restart_metadata

      !-------------------------------------!

      subroutine write_restart_metadata_header(unit)
      ! The block, as '#' comment lines. Every reader of these files skips
      ! '#' lines (the regression comparator drops them with grep -v), so no
      ! numeric parse and no golden verdict changes.
      integer, intent(in) :: unit
      character(len=meta_len) :: blk(n_meta)
      integer :: i
      call build_restart_metadata(blk)
      do i = 1, n_meta
         if (len_trim(blk(i)) .eq. 0) cycle
         write(unit,'(A)') '# '//trim(blk(i))
      enddo
      ! The state's option-change history (decision 21): one line per
      ! restart at which a named physics option was allowed to differ,
      ! oldest first. Informational and not compared, like the source
      ! field: it says what the state was reached from.
      do i = 1, n_ic_option_change
         write(unit,'(A)') '# '//trim(ic_option_change(i))
      enddo
      end subroutine write_restart_metadata_header

      !-------------------------------------!

      character(len=24) function meta_num(x)
      ! ES23.16 carries seventeen significant digits, so the decimal written
      ! here reads back as the same double.
      real*8, intent(in) :: x
      write(meta_num,'(ES23.16)') x
      meta_num = adjustl(meta_num)
      end function meta_num

      !-------------------------------------!

      character(len=1) function tf(flag)
      logical, intent(in) :: flag
      tf = 'F'
      if (flag) tf = 'T'
      end function tf

      !-------------------------------------!

      character(len=8) function opt_value(name)
      ! THE VALUE OF ONE OPTION TOKEN, one statement per token. The names
      ! live in opt_name and the values here, so a token that has a name
      ! and no value stops the run the first time a state file is written,
      ! which is every run: the two cannot drift apart unnoticed.
      character(len=*), intent(in) :: name
      opt_value = ''
      select case (trim(name))
      case ('He23S');          opt_value = tf(thereis_HeITR)
      case ('metals');         opt_value = tf(thereis_metals)
      case ('eos_metals');     opt_value = tf(eos_include_metals)
      case ('mol');            opt_value = tf(thereis_mol)
      case ('molbase');        opt_value = tf(molecular_base)
      case ('oxychem');        opt_value = tf(thereis_oxychem)
      case ('carrier');        opt_value = tf(carrier_transport)
      case ('carrier_newton'); opt_value = tf(carrier_in_newton)
      case ('iontrans');       opt_value = tf(ionization_transport)
      case ('he_diff');        opt_value = tf(he_diffusion)
      case ('he_metal_diff');  opt_value = tf(he_metal_diffusion)
      case ('sec_ion');        opt_value = tf(use_sec_ion)
      case ('caloric_mono');   opt_value = tf(caloric_eos_monatomic)
      case ('excH');           opt_value = tf(use_excited_H)
      case ('base_ir');        opt_value = tf(base_ir_field)
      case ('mol_ir');         opt_value = tf(mol_ir_bands)
      case ('mol_heat');       opt_value = tf(mol_reaction_heat)
      case ('visc');           opt_value = tf(visc_on)
      case ('cond');           opt_value = tf(cond_on)
      case ('jlya');           write(opt_value,'(I0)') jlya_mode
      case ('wellbal');        opt_value = tf(well_balanced)
      case default
         write(*,'(A)') ' (load_IC) ERROR: the option token "'//           &
              trim(name)//'" is named in opt_name and has no value in'//   &
              ' opt_value.'
         error stop 1
      end select
      end function opt_value

      !-------------------------------------!

      character(len=16) function opt_field_value(s, name)
      ! The value of '<name>=' in an options field, or '<absent>' when the
      ! field does not carry that token at all. Every token of the field is
      ! preceded by a blank (the field begins with its own tag), so the
      ! search is for ' <name>=' and 'mol=' cannot be found inside
      ! 'mol_ir=' or 'carrier=' inside 'carrier_newton='.
      character(len=*), intent(in) :: s, name
      integer :: p, q
      opt_field_value = '<absent>'
      p = index(s, ' '//trim(name)//'=')
      if (p .le. 0) return
      p = p + len_trim(name) + 2
      if (p .gt. len_trim(s)) return
      q = index(trim(s(p:)), ' ')
      if (q .le. 0) then
         opt_field_value = adjustl(s(p:len_trim(s)))
      else
         opt_field_value = adjustl(s(p:p+q-2))
      endif
      end function opt_field_value

      !-------------------------------------!

      subroutine append_option_change(text)
      ! One line onto the state's option-change history. A full history
      ! loses its oldest line and says how many it has lost, so a history
      ! read back from a file is either complete or states that it is not.
      character(len=*), intent(in) :: text
      integer :: k
      if (n_ic_option_change .lt. max_option_change) then
         n_ic_option_change = n_ic_option_change + 1
         ic_option_change(n_ic_option_change) = text
         return
      endif
      n_option_change_dropped = n_option_change_dropped + 1
      do k = 2, max_option_change - 1
         ic_option_change(k) = ic_option_change(k+1)
      enddo
      ic_option_change(max_option_change) = text
      write(ic_option_change(1),'(A,I0,A)') 'option_change ',              &
           n_option_change_dropped,                                        &
           ' earlier changes are not recorded (history length reached)'
      end subroutine append_option_change

      !-------------------------------------!

      subroutine collect_option_change_line(line)
      ! An 'option_change' line the loaded state already carried: the
      ! changes ITS restarts were allowed, kept in the order the file
      ! carries them and ahead of the line this restart adds, so the
      ! history travels with the state (decision 21).
      character(len=*), intent(in) :: line
      character(len=len(line)+1) :: t
      t = adjustl(line)
      if (t(1:1) .eq. '#') t = adjustl(t(2:))
      if (index(t, 'option_change ') .ne. 1) return
      call append_option_change(trim(t))
      end subroutine collect_option_change_line

      !-------------------------------------!

      subroutine parse_restart_metadata_line(line, blk)
      ! Store one '# <field> ...' line under its field name. A line whose
      ! first token names no field is not part of the block and is ignored,
      ! which is what keeps the coupling, provenance and row-layout lines
      ! out of it; a field name this version does not know is ignored the
      ! same way, so a file written by a later version is still readable.
      character(len=*), intent(in)    :: line
      character(len=meta_len), intent(inout) :: blk(n_meta)
      character(len=len(line)) :: t
      character(len=32) :: tok
      integer :: i, pos
      t = adjustl(line)
      if (len_trim(t) .eq. 0) return
      if (t(1:1) .eq. '#') t = adjustl(t(2:))
      pos = index(trim(t), ' ')
      if (pos .le. 1) then
         tok = trim(t)
      else
         tok = t(1:pos-1)
      endif
      do i = 1, n_meta
         if (trim(tok) .eq. trim(meta_tag(i))) then
            blk(i) = trim(t)
            return
         endif
      enddo
      end subroutine parse_restart_metadata_line

      !-------------------------------------!

      subroutine verify_restart_metadata(meta_h, meta_i)
      ! THE CONFIGURATION THE RESTART FILES STATE, AGAINST THE ONE THIS RUN
      ! RESOLVED (docs/restart_contract_design_20260909.md section 3).
      !
      ! Three outcomes, and no fourth: the pair carries no block and is
      ! loaded as a state of an unknown configuration; the pair carries a
      ! block that agrees with this run and the load proceeds; a field
      ! disagrees and the run stops naming that field. Nothing is repaired
      ! and nothing is inferred.
      character(len=meta_len), intent(in) :: meta_h(n_meta), meta_i(n_meta)
      character(len=meta_len) :: blk(n_meta)
      integer :: i, iv, ios
      real*8  :: t_file

      ic_restart_schema_present = (len_trim(meta_h(imeta_schema)) .gt. 0)
      ic_provenance_unknown     = .not. ic_restart_schema_present

      ! The two files are two halves of ONE state, so a pair whose halves
      ! state different configurations is not a state at all.
      if (ic_restart_schema_present .neqv.                                 &
          (len_trim(meta_i(imeta_schema)) .gt. 0)) then
         write(*,'(A)') ' (load_IC) ERROR: one of the two restart files'// &
              ' carries a restart metadata block and the other does not,'
         write(*,'(A)') '   so the pair was not written by one run.'//     &
              ' Restart from a matching pair of state files.'
         error stop 1
      endif

      if (.not. ic_restart_schema_present) then
         write(*,'(A)') ' (load_IC) provenance unknown: the restart'//     &
              ' files carry no restart metadata block (written before'
         write(*,'(A)') '   the block existed), so the composition'//      &
              ' reservoir, the grid, the constants and the physics'
         write(*,'(A)') '   options the state was produced under cannot'// &
              ' be compared with this run''s. The state is loaded as it'
         write(*,'(A)') '   was before, and every file this run writes'//  &
              ' carries restart_input=provenance_unknown.'
         return
      endif

      do i = 1, n_meta
         if (trim(meta_h(i)) .ne. trim(meta_i(i))) then
            write(*,'(A)') ' (load_IC) ERROR: the two restart files'//     &
                 ' disagree in the metadata field "'//trim(meta_tag(i))//  &
                 '":'
            write(*,'(A)') '   Hydro_ioniz_IC.txt: '//trim(meta_h(i))
            write(*,'(A)') '   Ion_species_IC.txt: '//trim(meta_i(i))
            error stop 1
         endif
      enddo

      call build_restart_metadata(blk)

      ! The schema version. A file of a later version may carry fields this
      ! version cannot compare, and a silent partial comparison is the one
      ! thing the block exists to prevent.
      read(meta_h(imeta_schema)(len_trim(meta_tag(imeta_schema))+1:),      &
           *, iostat=ios) iv
      if (ios .ne. 0 .or. iv .ne. restart_schema_version) then
         write(*,'(A)') ' (load_IC) ERROR: metadata field'//               &
              ' "restart_schema": the restart files state'
         write(*,'(A)') '   '//trim(meta_h(imeta_schema))//', this'//      &
              ' version writes and compares version'
         write(*,'(A,I0,A)') '   ', restart_schema_version, '.'
         error stop 1
      endif

      call compare_reservoir_field(meta_h(imeta_reservoir),                &
                                   blk(imeta_reservoir))
      call refuse_unequal_field(imeta_grid,    meta_h, blk)
      call refuse_unequal_field(imeta_const,   meta_h, blk)
      ! The options field is compared TOKEN BY TOKEN, so that the tokens
      ! the input named as allowed to differ can differ and every other
      ! difference is still refused by its own name (decision 21).
      call compare_options_field(meta_h(imeta_options), blk(imeta_options), &
                                 meta_h(imeta_source))

      ! Reported, not refused (see the block's header): the label map is
      ! what restores the species, and a missing element is built from the
      ! abundance and reported by the metals block of load_IC.
      if (trim(meta_h(imeta_species)) .ne. trim(blk(imeta_species))) then
         write(*,'(A)') ' (load_IC) NOTE: the restart files and this run'//&
              ' carry different species columns.'
         write(*,'(A)') '   file: '//trim(meta_h(imeta_species))
         write(*,'(A)') '   run:  '//trim(blk(imeta_species))
         write(*,'(A)') '   Species are restored by the label on the'//    &
              ' "# columns" line; an element the file does not carry is'
         write(*,'(A)') '   built from the abundance and reported above.'
      endif
      ! A STATE DESCENDED FROM AN UNKNOWN CONFIGURATION STAYS MARKED. The
      ! mark is a property of the state, not of one reload: if the file this
      ! run starts from was itself produced by a run that read a
      ! block-less state, nothing more is known about the configuration
      ! that state was reached under here than there, and dropping the mark
      ! at the second reload would lose it for good.
      if (index(meta_h(imeta_source), 'restart_input=provenance_unknown')  &
          .gt. 0) ic_provenance_unknown = .true.
      if (index(meta_h(imeta_source), 'git='//trim(build_git)) .le. 0)     &
         write(*,'(A)') ' (load_IC) the state was produced by another'//   &
              ' build: '//trim(meta_h(imeta_source))

      ! The clock, if the coupling line did not already carry it. ONE
      ! authority: write_coupling_state_header states t_phys for a physical
      ! state and is parsed first, so this line is read only for a file
      ! whose coupling line predates that field.
      if (.not. ic_t_phys_present) then
         read(meta_h(imeta_t_phys)(len_trim(meta_tag(imeta_t_phys))+1:),   &
              *, iostat=ios) t_file
         if (ios .eq. 0 .and. t_file .eq. t_file .and.                     &
             t_file .ge. 0.0d0 .and. abs(t_file) .le. huge(t_file)) then
            ic_t_phys         = t_file
            ic_t_phys_present = .true.
         endif
      endif
      end subroutine verify_restart_metadata

      !-------------------------------------!

      subroutine refuse_unequal_field(ifield, meta_file, meta_run)
      ! A field whose comparison is EXACT: the text of the field is the
      ! statement, and two different statements are two different problems.
      ! The grid and the constant set are compared this way; the options
      ! field has its own token-by-token comparison, because a named token
      ! may differ there.
      integer, intent(in) :: ifield
      character(len=meta_len), intent(in) :: meta_file(n_meta)
      character(len=meta_len), intent(in) :: meta_run(n_meta)
      if (trim(meta_file(ifield)) .eq. trim(meta_run(ifield))) return
      write(*,'(A)') ' (load_IC) ERROR: metadata field "'//                &
           trim(meta_tag(ifield))//'" differs between the restart files'// &
           ' and this run.'
      write(*,'(A)') '   file: '//trim(meta_file(ifield))
      write(*,'(A)') '   run:  '//trim(meta_run(ifield))
      write(*,'(A)') '   The state in the file was built on another'//     &
           ' discretization or another constant set. A state is not'
      write(*,'(A)') '   carried across such a change without a'//         &
           ' conservative remap, which is a separate workflow.'
      error stop 1
      end subroutine refuse_unequal_field

      !-------------------------------------!

      subroutine compare_options_field(s_file, s_run, s_source)
      ! THE OPTIONS FIELD, TOKEN BY TOKEN (decision 21 of
      ! docs/To_be_determined_by_user_20260906.md, option a).
      !
      ! Three things happen to a token, and nothing else can:
      !   it differs and the input did not name it  -> the load is refused
      !                                                and the token named;
      !   it differs and the input named it         -> the change is allowed
      !                                                and written into the
      !                                                new state's block;
      !   it was named and does not differ          -> said so, and nothing.
      ! A file carrying a token this version has no name for is refused as a
      ! whole: a comparison that silently skips part of the field is the one
      ! thing the block exists to prevent.
      character(len=*), intent(in) :: s_file, s_run, s_source
      character(len=16) :: vf, vr
      character(len=meta_len) :: chg_from, chg_to, refused
      integer :: i, nch, nref, nsame, nef, ner
      ic_option_change_applied = .false.
      ic_option_change_inert   = .false.
      nch   = 0
      nref  = 0
      nsame = 0
      chg_from = ''
      chg_to   = ''
      refused  = ''
      ! One token carries one '=' sign, so the count of them is the count of
      ! tokens: two fields of the same schema version with different counts
      ! do not carry the same vocabulary.
      nef = count_char(s_file, '=')
      ner = count_char(s_run,  '=')
      if (nef .ne. ner) then
         write(*,'(A)') ' (load_IC) ERROR: metadata field "options":'//    &
              ' the restart files carry a different set of option'
         write(*,'(A,I0,A,I0,A)') '   tokens than this version knows (',   &
              nef, ' tokens in the file, ', ner, ' here), so the'
         write(*,'(A)') '   comparison would silently skip part of the'//  &
              ' field.'
         write(*,'(A)') '   file: '//trim(s_file)
         write(*,'(A)') '   run:  '//trim(s_run)
         error stop 1
      endif
      do i = 1, n_opt
         vf = opt_field_value(s_file, opt_name(i))
         vr = opt_field_value(s_run,  opt_name(i))
         if (trim(vf) .eq. trim(vr)) then
            if (restart_option_change_named(i)) nsame = nsame + 1
            cycle
         endif
         if (restart_option_change_named(i)) then
            nch = nch + 1
            chg_from = trim(chg_from)//' '//trim(opt_name(i))//'='//trim(vf)
            chg_to   = trim(chg_to)//' '//trim(opt_name(i))//'='//trim(vr)
         else
            nref = nref + 1
            refused = trim(refused)//' '//trim(opt_name(i))//': '//        &
                 trim(vf)//' -> '//trim(vr)//';'
         endif
      enddo

      if (nref .gt. 0) then
         write(*,'(A)') ' (load_IC) ERROR: metadata field "options"'//     &
              ' differs between the restart files and this run in'
         write(*,'(A,I0,A)') '   ', nref, ' token(s) the input did not'//  &
              ' name as allowed to change:'
         write(*,'(A)') '  '//trim(refused)
         write(*,'(A)') '   The state in the file is a state of another'// &
              ' equation set. Restart it under the options it was'
         write(*,'(A)') '   produced with, start this configuration'//     &
              ' cold, or name the tokens on a'
         write(*,'(A)') '   "Restart option change:" line, which allows'// &
              ' exactly the tokens it names to differ.'
         error stop 1
      endif

      if (nsame .gt. 0) then
         ic_option_change_inert = .true.
         write(*,'(A,I0,A)') ' (load_IC) NOTE: "Restart option change"'//  &
              ' names ', nsame, ' token(s) that do not differ between'
         write(*,'(A)') '   the state and this run. Nothing was allowed'// &
              ' that was needed; the key is harmless here.'
      endif

      if (nch .gt. 0) then
         ic_option_change_applied = .true.
         write(*,'(A)') ' (load_IC) the restart CHANGES the physics'//     &
              ' options it was allowed to change:'
         write(*,'(A)') '  '//trim(adjustl(chg_from))//' ->'//trim(chg_to)
         write(*,'(A)') '   The loaded state is a solution of the'//       &
              ' former and is not stationary under the latter; it is a'
         write(*,'(A)') '   starting point, and the change is recorded'//  &
              ' in every file this run writes.'
         call append_option_change('option_change '//                      &
              trim(adjustl(chg_from))//' ->'//trim(chg_to)//              &
              ' at restart of '//trim(source_body(s_source)))
      endif
      end subroutine compare_options_field

      !-------------------------------------!

      integer function count_char(s, c)
      character(len=*), intent(in) :: s
      character(len=1), intent(in) :: c
      integer :: k
      count_char = 0
      do k = 1, len_trim(s)
         if (s(k:k) .eq. c) count_char = count_char + 1
      enddo
      end function count_char

      !-------------------------------------!

      character(len=meta_len) function source_body(s_source)
      ! What the loaded file says produced it, without the field's own tag:
      ! the state an option change was made at is identified by the build
      ! that wrote it.
      character(len=*), intent(in) :: s_source
      integer :: p
      p = len_trim(meta_tag(imeta_source)) + 1
      source_body = 'an unstated source'
      if (len_trim(s_source) .le. p) return
      source_body = adjustl(s_source(p+1:))
      end function source_body

      !-------------------------------------!

      subroutine compare_reservoir_field(s_file, s_run)
      ! The elemental reservoir, element by element. The file carries the
      ! decimal the writing run's input stated, so the comparison is the
      ! same round-off allowance the He/H check uses (heh_dev_tol) and not
      ! an exact one; an element present in one reservoir and absent from
      ! the other is a different composition whatever the tolerance.
      character(len=*), intent(in) :: s_file, s_run
      character(len=32) :: nm_f(n_melem+1), nm_r(n_melem+1)
      real*8            :: v_f(n_melem+1), v_r(n_melem+1)
      integer :: nf, nr, i, k
      real*8  :: dev
      call parse_reservoir_field(s_file, nm_f, v_f, nf)
      call parse_reservoir_field(s_run,  nm_r, v_r, nr)
      do i = 1, nf
         do k = 1, nr
            if (trim(nm_f(i)) .eq. trim(nm_r(k))) exit
         enddo
         if (k .gt. nr) call refuse_reservoir_element(nm_f(i), v_f(i),     &
                             'the restart files carry',                    &
                             'this run does not')
         dev = abs(v_f(i) - v_r(k))/max(abs(v_r(k)), 1.0d-30)
         if (dev .gt. heh_dev_tol) then
            write(*,'(A)') ' (load_IC) ERROR: metadata field'//            &
                 ' "reservoir": the elemental ratio '//trim(nm_f(i))//     &
                 ' of the'
            write(*,'(A,ES16.9)') '   restart files is', v_f(i)
            write(*,'(A,ES16.9)') '   this run asks for', v_r(k)
            write(*,'(A,ES11.4,A,ES11.4)') '   relative departure',        &
                 dev, ', threshold', heh_dev_tol
            write(*,'(A)') '   The state is a solution at its own'//       &
                 ' composition. Restart it there, or start this'
            write(*,'(A)') '   composition cold.'
            error stop 1
         endif
      enddo
      do k = 1, nr
         do i = 1, nf
            if (trim(nm_r(k)) .eq. trim(nm_f(i))) exit
         enddo
         if (i .gt. nf) call refuse_reservoir_element(nm_r(k), v_r(k),     &
                             'this run carries',                           &
                             'the restart files do not')
      enddo
      end subroutine compare_reservoir_field

      !-------------------------------------!

      subroutine refuse_reservoir_element(name, value, has, lacks)
      character(len=*), intent(in) :: name, has, lacks
      real*8,           intent(in) :: value
      write(*,'(A)') ' (load_IC) ERROR: metadata field "reservoir": '//    &
           trim(has)//' a reservoir for '//trim(name)
      write(*,'(A,ES16.9,A)') '   (', value, ') and '//trim(lacks)//       &
           '. The two are different compositions, and the'
      write(*,'(A)') '   state in the file is a solution at its own.'//    &
           ' Restart it there, or start this composition cold.'
      error stop 1
      end subroutine refuse_reservoir_element

      !-------------------------------------!

      subroutine parse_reservoir_field(s, names, values, n_out)
      ! 'reservoir He/H <v> [<El>/H <v> ...]' -> the pairs it carries.
      character(len=*),  intent(in)  :: s
      character(len=32), intent(out) :: names(:)
      real*8,            intent(out) :: values(:)
      integer,           intent(out) :: n_out
      character(len=len(s)) :: rest
      character(len=64) :: tok
      integer :: pos, l, ios
      n_out = 0
      rest  = adjustl(s)
      ! Drop the field name.
      pos = index(trim(rest), ' ')
      if (pos .le. 0) return
      rest = adjustl(rest(pos:))
      do
         l = len_trim(rest)
         if (l .eq. 0) exit
         pos = index(trim(rest), ' ')
         if (pos .le. 1) then
            tok = trim(rest);  rest = ''
         else
            tok = rest(1:pos-1);  rest = adjustl(rest(pos:))
         endif
         if (index(tok, '/H') .le. 0) cycle
         ! The value is the next token.
         l = len_trim(rest)
         if (l .eq. 0) exit
         if (n_out .ge. size(names)) exit
         n_out = n_out + 1
         names(n_out) = trim(tok)
         pos = index(trim(rest), ' ')
         if (pos .le. 1) then
            read(rest,*,iostat=ios) values(n_out);  rest = ''
         else
            read(rest(1:pos-1),*,iostat=ios) values(n_out)
            rest = adjustl(rest(pos:))
         endif
         if (ios .ne. 0) values(n_out) = -1.0d0
      enddo
      end subroutine parse_reservoir_field

      !-------------------------------------!

      subroutine report_deviating_heh_cell(j_dev, heh_dev_cell, dev)
      ! The cell that decided the He/H comparison: its index, its radius,
      ! the He/H the file carries there, the He/H the input asks for and
      ! the relative departure against the threshold. Printed with enough
      ! digits to separate a round-off difference from a composition
      ! difference -- ES11.4 renders 0.0793 and 0.079307 identically, and
      ! two numbers that print the same cannot explain a refusal.
      integer, intent(in) :: j_dev
      real*8,  intent(in) :: heh_dev_cell, dev
      if (j_dev .le. 0) then
         write(*,'(A)') '   no physical cell carries hydrogen, so'//    &
            ' no cell states an He/H at all.'
         return
      endif
      write(*,'(A,I0,A,ES13.6,A)') '   farthest cell ', j_dev,          &
         ' at r =', r(j_dev), ' Rp:'
      write(*,'(A,ES16.9)') '     He/H in the file  ', heh_dev_cell
      write(*,'(A,ES16.9)') '     He/H in the input ', HeH
      write(*,'(A,ES11.4,A,ES11.4)') '     relative departure',         &
         dev, ', threshold', heh_dev_tol
      end subroutine report_deviating_heh_cell

      !-------------------------------------!

      subroutine parse_coupling_header(line)
      ! Read the '# coupling: key=value ...' line of a restart file. Unknown
      ! keys are skipped, so a file written by a later version that carries
      ! more of them is still readable; a key this version knows but the file
      ! omits keeps its "not stated" default.
      character(len=*), intent(in) :: line
      character(len=len(line)) :: rest
      character(len=64) :: tok
      character(len=64) :: key, val
      integer :: pos, l, ieq
      rest = adjustl(line)
      do
         l = len_trim(rest)
         if (l .eq. 0) exit
         pos = index(rest, ' ')
         if (pos .le. 1) then
            tok = rest(1:min(l,len(tok)));  rest = ''
         else
            tok = rest(1:min(pos-1,len(tok)));  rest = adjustl(rest(pos:))
         endif
         ieq = index(tok, '=')
         if (ieq .le. 1) cycle
         key = tok(1:ieq-1)
         val = tok(ieq+1:)
         select case (trim(key))
         case ('sec_ion')
            ic_sec_ion_active   = (trim(val) .eq. 'T')
            ic_coupling_present = .true.
         case ('sec_ion_step')
            read(val,*,iostat=pos) ic_sec_ion_armed_step
         case ('valve')
            read(val,*,iostat=pos) ic_valve_eps
         case ('recon')
            ic_rec_method = trim(val)
         case ('iontrans')
            ic_ionization_transport = (trim(val) .eq. 'T')
         case ('certified')
            ! The stationary claim the writing run made ABOUT THIS STATE. A
            ! re-evaluation of the state answers the same question, so this
            ! is what the answer is held to: a state written as certified
            ! and re-evaluated as uncertified is a refused claim.
            ic_certified = (trim(val) .eq. 'T')
         case ('mode')
            ic_run_mode_present = .true.
            if (trim(val) .eq. 'phys') then
               ic_run_mode = run_mode_phys
            else
               ic_run_mode = run_mode_init
            endif
         case ('t_phys')
            read(val,*,iostat=pos) ic_t_phys
            ic_t_phys_present = (pos .eq. 0) .and. (ic_t_phys .eq. ic_t_phys) &
                                .and. (abs(ic_t_phys) .le. huge(ic_t_phys))   &
                                .and. (ic_t_phys .ge. 0.0d0)
            if (.not. ic_t_phys_present) ic_t_phys = 0.0d0
         case ('fluxconst')
            read(val,*,iostat=pos) ic_base_flux_const
         end select
      enddo
      end subroutine parse_coupling_header

      !-------------------------------------!

      logical function comment_field_is(line, name)
      ! Is this comment line's FIRST token the field 'name'? The header
      ! lines of a state file are 'name value ...' records behind a '#', so
      ! the field is the first token and nothing else.
      character(len=*), intent(in) :: line, name
      character(len=len(line)) :: t
      integer :: pos
      comment_field_is = .false.
      t = adjustl(line)
      if (len_trim(t) .eq. 0) return
      if (t(1:1) .eq. '#') t = adjustl(t(2:))
      pos = index(trim(t), ' ')
      if (pos .le. 1) then
         comment_field_is = (trim(t) .eq. trim(name))
      else
         comment_field_is = (t(1:pos-1) .eq. trim(name))
      endif
      end function comment_field_is

      !-------------------------------------!

      logical function is_comment(line)
      ! True if the (left-adjusted) line starts with '#'
      character(len=*), intent(in) :: line
      character(len=len(line)) :: t
      t = adjustl(line)
      is_comment = (len_trim(t) .gt. 0 .and. t(1:1) .eq. '#')
      end function is_comment

      !-------------------------------------!

      subroutine parse_labels(line, labels, nlab)
      ! Split a "# columns r[Rp] HI HII ..." line into its column labels
      ! ('#' and 'columns' tokens are dropped; labels(1) is the r column).
      character(len=*), intent(in)  :: line
      character(len=16), intent(out) :: labels(:)
      integer, intent(out) :: nlab
      integer :: pos, l, start
      character(len=len(line)) :: rest
      character(len=64) :: tok

      rest = adjustl(line)
      nlab = 0
      do
         l = len_trim(rest)
         if (l .eq. 0) exit
         pos = index(rest, ' ')
         if (pos .le. 1) then
            tok = rest(1:l); rest = ''
         else
            tok = rest(1:pos-1); rest = adjustl(rest(pos:))
         endif
         if (trim(tok) .eq. '#')        cycle
         if (trim(tok) .eq. 'columns')  cycle
         ! The caller's labels array and its vals row are the same length, so
         ! a header longer than that has no place to be stored and the read
         ! of the data row below would write past the row buffer.
         if (nlab .ge. size(labels)) then
            write(*,'(A,I0,A)') ' (load_IC) ERROR: the "# columns" header'// &
               ' carries more than ', size(labels), ' labels, which is'//    &
               ' more columns than the reader holds.'
            error stop 1
         endif
         nlab = nlab + 1
         labels(nlab) = trim(tok)
      enddo
      end subroutine parse_labels

      !-------------------------------------!

      integer function species_column(label)
      ! f_sp column for a species label (0 = unknown/ignored, e.g. r[Rp])
      character(len=*), intent(in) :: label
      integer :: i
      select case (trim(label))
         case ('HI');    species_column = isp_HI
         case ('HII');   species_column = isp_HII
         case ('HeI');   species_column = isp_HeI
         case ('HeII');  species_column = isp_HeII
         case ('HeIII'); species_column = isp_HeIII
         case ('HeITR'); species_column = isp_HeTR
         ! molecular columns
         case ('H2');    species_column = isp_H2
         case ('H2p');   species_column = isp_H2p
         case ('H3p');   species_column = isp_H3p
         case ('HeHp');  species_column = isp_HeHp
         ! oxygen-chemistry columns
         case ('OH');    species_column = isp_OH
         case ('H2O');   species_column = isp_H2O
         case ('CO');    species_column = isp_CO
         case default
            species_column = 0
            do i = 1, n_mion
               if (trim(label) .eq. trim(mion_name(i))) then
                  species_column = mion_fsp(i)
                  return
               endif
            enddo
      end select
      end function species_column

      !-------------------------------------!

      subroutine verify_restart_radii(r_file, fname)
      ! THE GRID BELONGS TO THE RUN, NOT TO THE STATE FILE.
      ! define_grid builds the cell centers r together with the faces r_edg,
      ! the widths dr_j and the window indices j_min and j_flux; a state file
      ! carries only the centers. A file written under a different
      ! "Grid type", "Base grid" or "Outer radius" passes the row-count
      ! guard above whenever it has the same "Grid cells", so the centers are
      ! compared here and a mismatch is refused
      ! (docs/development_plan_20260905_rev3.md section 10.2 item 3).
      !
      ! Compared over the physical cells 1..N only. The ghost rows are
      ! boundary data, rebuilt on the first step of the restart, and the
      ! writer's ghost convention changed (docs/Update_EXHALE_stage1.md
      ! section 169), so files written before it carry ghosts of their own.
      !
      ! Tolerance 1e-10 relative: the writer emits list-directed reals with
      ! 17 significant digits, so a restart on the grid the state was written
      ! on reproduces every center to round-off, while any other grid moves a
      ! center by many orders of magnitude more (MEASURED 2026-09-05: the
      ! smallest real grid difference in the tree, 3.3e-4). The coarsest
      ! writer of a state file is the Wind-AE initial condition
      ! (wae_ic_writer.f90, ES17.10 = 11 significant digits), whose worst
      ! relative round-off is 5.0e-11, half of this tolerance. This reader
      ! does not interpolate; a state on another grid is mapped onto this one
      ! before it is handed back.
      real*8,  intent(in) :: r_file(1-Ng:N+Ng)
      character(len=*), intent(in) :: fname
      integer :: j, j_worst
      real*8  :: d, d_max
      d_max   = 0.0d0
      j_worst = 1
      do j = 1, N
         d = abs(r_file(j) - r(j))/max(abs(r(j)), 1.0d-30)
         if (d .gt. d_max) then
            d_max   = d
            j_worst = j
         endif
      enddo
      if (d_max .gt. 1.0d-10) then
         write(*,'(A)') ' (load_IC) ERROR: '//trim(fname)//' holds cell'//   &
            ' centers of a different grid from the one this run built.'
         write(*,'(A,A,A,I0,A)') '   this run: "Grid type: ',                &
            trim(grid_type), '" with ', N, ' cells;'
         write(*,'(A,I0,A,ES23.16)') '   cell ', j_worst,                    &
            ': the run has r = ', r(j_worst)
         write(*,'(A,ES23.16)') '                    the file has r = ',     &
            r_file(j_worst)
         write(*,'(A,ES9.2,A)') '   worst relative difference ', d_max,      &
            ', tolerance 1.0e-10.'
         write(*,'(A)') '   The faces, widths and window indices are this'// &
            ' run''s and are not re-read, so the state cannot be marched'
         write(*,'(A)') '   on the file''s centers. Restart on the grid the'//&
            ' state was written on, or map the state onto this grid first'
         write(*,'(A)') '   (docs/development_plan_20260905_rev3.md'//       &
            ' section 10.2 item 3).'
         error stop 1
      endif
      end subroutine verify_restart_radii

      !-------------------------------------!

      subroutine scatter_row(j, vals, col2fsp, nlab, nsp_l, r_file)
      ! Store one data record: vals(1) is the cell center the file states,
      ! kept apart from the run's own r for the comparison in
      ! verify_restart_radii; vals(k>=2) go to their mapped f_sp columns.
      integer, intent(in) :: j, nlab
      real*8,  intent(in) :: vals(:)
      integer, intent(in) :: col2fsp(:)
      real*8,  intent(inout) :: nsp_l(1-Ng:N+Ng, n_species)
      real*8,  intent(inout) :: r_file(1-Ng:N+Ng)
      integer :: k
      r_file(j) = vals(1)
      do k = 2, nlab
         if (col2fsp(k) .gt. 0) nsp_l(j,col2fsp(k)) = vals(k)
      enddo
      end subroutine scatter_row

      ! End of module
      end module IC_load
