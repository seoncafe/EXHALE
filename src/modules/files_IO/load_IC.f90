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
      ! Oxygen-chemistry seed for a restart file written before the option
      ! existed (see the block near the end of load_IC).
      use oxygen_rates, only: co_equilibrium_density,                   &
                              oxygen_chemical_equilibrium_fractions

      implicit none

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
      ! Hydrogen nuclei density of the loaded state (free + bound in molecules)
      real*8, dimension(1-Ng:N+Ng)        :: nH_l
      ! Helium nuclei density of the loaded state, and the two factors that
      ! carry the loaded composition onto the input one
      real*8, dimension(1-Ng:N+Ng)        :: nHe_l, sH_l, sHe_l
      real*8 :: heh_loaded, heh_dev, gotH_l, gotHe_l
      logical :: col_present(n_species), elem_ok
      ! Auxiliary temporary variable
      real*8 :: tmp
      ! Loaded El/H at the base cell, and the single factor that carries the
      ! loaded metal column onto a handoff reservoir (metals block below)
      real*8 :: elh_l, r_el
      real*8 :: vals(80)

      character(len=8192) :: line
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
      open(unit = 1, file = 'output/Hydro_ioniz_IC.txt')
      do
         read(1,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (is_comment(line)) then
            if (index(line,'coupling:') .gt. 0) call parse_coupling_header(line)
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
         read(line,*) tmp, tmp, v(1-Ng), p(1-Ng), T(1-Ng), tmp, tmp
         do j = 2-Ng,N+Ng
            read(1,*) tmp, tmp, v(j), p(j), T(j), tmp, tmp
         enddo
      close(1)

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
            if (index(line,'columns') .gt. 0) call parse_labels(line, labels, nlab)
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
         call scatter_row(1-Ng, vals, col2fsp, nlab, nsp_l)
         do j = 2-Ng,N+Ng
            read(2,*) (vals(k), k = 1,nlab)
            call scatter_row(j, vals, col2fsp, nlab, nsp_l)
         enddo

      else
         ! ---- legacy headerless file: fixed 7-column layout ----
         rewind(2)
         do j = 1-Ng,N+Ng
            read(2,*) r(j),                  &
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
         ! had moved in between (fixed the same day, docs/Update_EXHALE.md
         ! section 169; files written before it still carry such ghosts).
         ! The range stays physical either way: a check on the input's He/H
         ! has nothing to say about boundary data.
         heh_dev = 0.0d0
         do j = 1, N
            if (nH_l(j) .gt. 0.0d0) then
               heh_loaded = nHe_l(j)/nH_l(j)
               heh_dev = max(heh_dev, abs(heh_loaded - HeH)/max(HeH,1.0d-30))
            endif
         enddo
         if (heh_dev .gt. 1.0d-6) then
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
               write(*,'(A,ES11.4,A,ES11.4,A)')                           &
                  ' (load_IC) ERROR: the restart file was written at'//   &
                  ' He/H =', nHe_l(N)/nH_l(N), ', the input asks for',    &
                  HeH, ', and the state carries HeH+, whose nucleus of'// &
                  ' each element cannot be rescaled by one factor.'//     &
                  ' Restart a molecular state at its own composition.'
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
               write(*,'(A,ES11.4,A,ES11.4)')                             &
                  ' (load_IC) restart file written at He/H =',            &
                  nHe_l(N)/nH_l(N), ' rescaled to the input He/H =', HeH
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
         case ('fluxconst')
            read(val,*,iostat=pos) ic_base_flux_const
         end select
      enddo
      end subroutine parse_coupling_header

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

      subroutine scatter_row(j, vals, col2fsp, nlab, nsp_l)
      ! Store one data record: vals(1) is r, vals(k>=2) go to their
      ! mapped f_sp columns.
      integer, intent(in) :: j, nlab
      real*8,  intent(in) :: vals(:)
      integer, intent(in) :: col2fsp(:)
      real*8,  intent(inout) :: nsp_l(1-Ng:N+Ng, n_species)
      integer :: k
      r(j) = vals(1)
      do k = 2, nlab
         if (col2fsp(k) .gt. 0) nsp_l(j,col2fsp(k)) = vals(k)
      enddo
      end subroutine scatter_row

      ! End of module
      end module IC_load
