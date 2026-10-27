   module utils
	! Collection of auxiliary subroutines

   use global_parameters
   use build_stamp, only: build_git, build_dirty
   use species_table, only: n_mion, n_mphot, mion_isphot, mion_iphot,  &
	                         mion_stage, mion_elem, melem_A,  &
	                         bsp_charge, bsp_mass

   implicit none

   ! WHETHER THE STATE NOW HELD PASSED THE STATIONARY CERTIFICATION
   ! (certification.f90). It is a verdict about the state a file carries,
   ! written into the '# coupling:' header so that a reader of the file --
   ! and not only a reader of the run log -- can tell a certified stationary
   ! solution from a state that was merely written. .false. until a
   ! certification is made, so a run that never reaches a declaration point
   ! reports what it is: not certified.
   logical, save :: state_is_certified = .false.
   ! Why it is .false., where the answer is not "an equation of the inventory
   ! refused it": a run that ends on a step cap or any other bound makes NO
   ! stationarity claim at all, and a state written by such a run is a
   ! relaxation snapshot. The distinction is the run state, and it is
   ! written into the header so that a reader of the file can tell the two
   ! apart. Empty
   ! for a certified state and for a refusal on the entries themselves.
   character(len=32), save :: state_certification_reason = ''

	contains

   subroutine set_state_certified(flag, reason)
   ! Set by the certification contexts (certification.f90) and by the
   ! run's own adoption of a loaded claim; every other reader takes it. The
   ! pass state writer (publish_pass_state_generation) sets the pair of the
   ! pass snapshot for its own write and puts the run's pair back after it.
   logical,          intent(in) :: flag
   character(len=*), intent(in) :: reason
   state_is_certified         = flag
   state_certification_reason = reason
   end subroutine set_state_certified



   subroutine write_row_layout_header(unit)
   ! The one line every profile file carries to say what its rows ARE.
   !
   ! Every output written on the full 1-Ng:N+Ng range begins with Ng ghost
   ! rows and ends with Ng more, and they are not part of the solution:
   ! including them in a flux-spread or a residual measure doubles it (the
   ! accepted flux spread of the HD 189733 b solve is 4.64e-3 over the
   ! physical cells and 1.05e-2 with the two upper ghost rows counted,
   ! section 133.6 of docs/Update_EXHALE_stage1.pdf).  A reader that does not know
   ! this silently averages two rows of boundary data into every profile.
   !
   ! It was written by Hydro_ioniz(_adv).txt alone, which made it a property
   ! of one file rather than of the convention; it is now written from here
   ! by every profile writer, so the sentence has one author.  The line is a
   ! '#' comment, so no numeric parse and no golden changes.
      integer, intent(in) :: unit
   write(unit,'(A,I0,A,I0,A,I0,A,I0)') '# rows ', N+2*Ng, ': ', Ng,       &
        ' ghost cells at each end; physical cells are rows ', Ng+1,       &
        ' to ', Ng+N
   end subroutine write_row_layout_header



   integer(8) function file_rolling_checksum(fname) result(ck)
   ! A CHECKSUM OF A FILE'S BYTES, for provenance only.
   !
   ! It is NOT a cryptographic digest and is deliberately not called one: two
   ! independent modular rolling sums,
   !     h_i <- mod( h_i * a_i + byte, 2^31 - 1 ),   a = (131, 8191),
   ! packed into one 62-bit integer. Each rolling multiply-add is below
   ! 2^44 and the packed result is below 2^62, so arithmetic is exact in integer(8)
   ! on every compiler -- which a 64-bit FNV or CRC written in standard
   ! Fortran is not, since they rely on unsigned wraparound the language does
   ! not define.
   !
   ! Returns -1 for a file that cannot be opened, which is how an absent
   ! metals.inp or base.inp is reported.
   character(len=*), intent(in) :: fname
   integer, parameter :: m = 2147483647            ! 2^31 - 1
   integer(8) :: h1, h2
   integer :: u, ios
   integer(1) :: b
   logical :: there
   ck = -1_8
   inquire(file=fname, exist=there)
   if (.not. there) return
   open(newunit=u, file=fname, access='stream', form='unformatted',       &
        status='old', action='read', iostat=ios)
   if (ios .ne. 0) return
   h1 = 1_8;  h2 = 1_8
   do
      read(u, iostat=ios) b
      if (ios .ne. 0) exit
      h1 = mod(h1*131_8   + int(iand(int(b,4), 255), 8), int(m,8))
      h2 = mod(h2*8191_8  + int(iand(int(b,4), 255), 8), int(m,8))
   enddo
   close(u)
   ck = h1*int(m,8) + h2
   end function file_rolling_checksum


   subroutine write_provenance_header(unit)
   ! WHAT PRODUCED THIS FILE, so that a profile found on disk in a year can
   ! be tied to an executable, an input and a set of physics options without
   ! asking anyone. The external review of 2026-09-03
   ! asks for exactly this: scientific products marked with the executable
   ! and source revision and the physics options they were made under.
   !
   !   git / dirty           the revision the EXECUTABLE was built from and
   !                         whether the working tree carried uncommitted
   !                         changes at that moment. Stamped in by the
   !                         Makefile; the run never calls git. NO BUILD TIME
   !                         is compiled in: a timestamp in a constant makes
   !                         two builds of one source differ, which destroys
   !                         the only cheap test for a stale object.
   !   run                   when THIS run wrote the file, from date_and_time.
   !   ck_input / ck_base / ck_metals
   !                         checksums of the input files AS READ, computed
   !                         now (file_rolling_checksum above; -1 = absent).
   !   recon / base_bc / carrier / restart_schema / resid_def / N
   !                         the discretization, the lower boundary model,
   !                         the molecular-carrier model, the version of the
   !                         restart file's own layout, the version of the
   !                         residual definition the gates used, and the grid.
   !
   ! It is a '#' comment, so no numeric parse and no golden changes: the
   ! regression compares with grep -v '^ *#'.
      integer, intent(in) :: unit
      character(len=16) :: bmod, cmod
      character(len=8)  :: d_ymd
      character(len=10) :: d_hms
      ! THERE IS ONE LOWER BOUNDARY, so this field names it rather than
      ! selecting among options. Section 152 replaced the three component-wise
      ! ghost closures this used to distinguish -- `hydrostatic_base`,
      ! `Base ghost temperature: continuous` and the isothermal default -- with
      ! a characteristic condition imposed at the face, and retired their keys.
      ! The field is kept, and kept in this position, so that a reader of an
      ! older profile can still tell which boundary produced it.
      bmod = 'characteristic'
      if (.not. thereis_mol) then
         cmod = 'none'
      else if (carrier_transport) then
         cmod = 'transported'
      else
         cmod = 'local'
      endif
      call date_and_time(date=d_ymd, time=d_hms)
      write(unit,'(A)') '# provenance: git='//trim(build_git)//             &
           ' tree='//trim(build_dirty)//                                    &
           ' run='//d_ymd(1:4)//'-'//d_ymd(5:6)//'-'//d_ymd(7:8)//'T'//     &
           d_hms(1:2)//':'//d_hms(3:4)//':'//d_hms(5:6)
      write(unit,'(A,I0,A,I0,A,I0)') '# provenance: ck_input=',             &
           file_rolling_checksum(trim(inp_file)),                              &
           ' ck_base=',   file_rolling_checksum(companion_path('base.inp')),                &
           ' ck_metals=', file_rolling_checksum(companion_path('metals.inp'))
      write(unit,'(A,I0)') '# provenance: recon='//                        &
           trim(reconstruction_operator_label())//                          &
           ' base_bc='//trim(bmod)//' carrier='//trim(cmod)//              &
           ' restart_schema=3 resid_def=145 N=', N
   end subroutine write_provenance_header

   subroutine write_coupling_state_header(unit, state_id)
   ! Write the # coupling line that records the staged physics of this
   ! state, independent of the selected input's current defaults.
   ! sec_ion and sec_ion_step record secondary-ionization activation;
   ! recon names the active reconstruction, including its continuation;
   ! certified, cert_reason and mode describe the solve status.
   ! coupling_state_fields also supplies the derived-product provenance.
   ! A pass-state generation appends state_id to both halves of its pair;
   ! load_IC refuses unequal or one-sided identities. Older lines omit it.
   ! Header fields and restart behavior: ../EXHALE_v1.00/md/input_schema.md appendix D.2.
      integer, intent(in) :: unit
      character(len=*), intent(in), optional :: state_id
      if (present(state_id)) then
         write(unit,'(A)') '# coupling: '//trim(coupling_state_fields())// &
                           ' state_id='//trim(state_id)
      else
         write(unit,'(A)') '# coupling: '//trim(coupling_state_fields())
      endif
   end subroutine write_coupling_state_header

   function coupling_state_fields(certification_first) result(fields)
   ! The key=value fields of the '# coupling:' line, without its label, in
   ! the order that line has always carried them (write_coupling_state_header
   ! states what each field means). With certification_first present and
   ! true, the same fields with the certification pair (certified,
   ! cert_reason) moved to the front: the order of the '# derived_from:'
   ! line of the _adv products, whose readers key on 'derived_from:
   ! certified='.
      logical, intent(in), optional :: certification_first
      character(len=512) :: fields
      character(len=1) :: s, c
      character(len=48) :: why
      character(len=32) :: ptr
      character(len=64) :: mtok
      character(len=32) :: tbuf
      character(len=16) :: sbuf
      s = 'F';  if (sec_ion_active) s = 'T'
      ! iontrans: the state in this file was produced with the hydrogen
      ! ionization state CARRIED (Ionization transport), so its H+ column is a
      ! transported quantity and not the local root of its own cell. A
      ! reader that restarts it with the option off will overwrite that
      ! column on the first sweep, which is a real change of physics and
      ! worth saying out loud. Written only when the option is on, so every
      ! file a run without it produces is unchanged, byte for byte.
      ptr = ''
      if (ionization_transport) ptr = ' iontrans=T'
      ! he23strans: the He 2^3S column is a population the flow carried
      ! (He 2^3S transport), not the local root of its cell's balance.
      ! Written only when the option is on, for the reason above.
      if (he23s_transport) ptr = trim(ptr)//' he23strans=T'
      ! certified: the stationary certification of A2 was made on this state
      ! and every active equation it could evaluate was within its tolerance,
      ! none of them was unavailable, no cell of the state was without a
      ! chemical root, and no unbudgeted accepted correction stands in its
      ! history. F is the honest answer for a state no certification was made
      ! on at all, which is what a run that stops on du reports.
      c = 'F';  if (state_is_certified) c = 'T'
      why = ''
      if (len_trim(state_certification_reason) .gt. 0)                    &
         why = ' cert_reason='//trim(state_certification_reason)
      ! mode / t_phys: WHICH OF THE THREE RUN STATES PRODUCED THIS STATE.
      ! A state written by an initialization or continuation run is a
      ! relaxation
      ! snapshot: it has no elapsed time and none is written for it, so a
      ! restart cannot invent one. A state written by a physical integration
      ! carries the time it was reached at, in seconds, as the sum of the
      ! global dt of the accepted steps behind it -- which is what makes a
      ! continuation of that trajectory possible at all.
      if (run_mode .eq. run_mode_phys) then
         write(tbuf,'(ES23.16)') t_phys
         mtok = ' mode=phys t_phys='//trim(adjustl(tbuf))
      else
         mtok = ' mode=init'
      endif
      write(sbuf,'(I0)') sec_ion_armed_step
      fields = 'sec_ion='//s//' sec_ion_step='//trim(sbuf)//              &
               ' recon='//trim(reconstruction_operator_label())//         &
               trim(ptr)//' certified='//c//trim(why)//trim(mtok)
      if (present(certification_first)) then
         if (certification_first)                                         &
            fields = 'certified='//c//trim(why)//' sec_ion='//s//         &
                     ' sec_ion_step='//trim(sbuf)//                       &
                     ' recon='//trim(reconstruction_operator_label())//   &
                     trim(ptr)//trim(mtok)
      endif
   end function coupling_state_fields

	subroutine calc_ne(nhii,nheii,nheiii,ne,nm,nmol)
	! Calculate the free electron density.
	! The optional nm (metal densities for each ion, same units as nhii) adds
	! the metal electrons when the eos_metals policy is on; omitting it
	! (or eos_metals 0) reproduces the legacy H/He-only electron count.
	! The base H/He/molecular electrons are accumulated in the canonical bsp
	! order of species_table, each species weighted by its net charge
	! (bsp_charge = free electrons released); neutral species (charge 0) are
	! skipped as no-ops.  Accumulating in bsp order deliberately fixes the FP
	! add order (a golden re-snapshot decision).

	integer :: im
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nhii
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nheii,nheiii
	! Optional molecular ions: cols 1 H2 (neutral), 2 H2+, 3 H3+,
	! 4 HeH+ -- each molecular ion carries one electron.
	real*8, dimension(1-Ng:N+Ng,4), intent(in), optional :: nmol
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in), optional :: nm
	real*8, dimension(1-Ng:N+Ng), intent(out) :: ne

	ne = 0.0d0

	! Base H/He electrons (bsp 2,4,5; bsp 1 HI and 3 HeI are neutral, and
	! bsp 6 HeTR is both neutral and an excited level of HeI, already inside
	! the HeI column -- see bsp_is_excited_level in species_table).
	call accum(nhii,   dble(bsp_charge(2)))        ! HII   (bsp 2)
	if (thereis_He) then
		call accum(nheii,  dble(bsp_charge(4)))     ! HeII  (bsp 4)
		call accum(nheiii, dble(bsp_charge(5)))     ! HeIII (bsp 5)
	endif

	if (present(nm) .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			if (mion_stage(im) .gt. 0)                                  &
				ne = ne + dble(mion_stage(im))*nm(:,im)
		enddo
	endif

	! Molecular ions (bsp 8,9,10; bsp 7 H2 is neutral).
	if (present(nmol)) then
		call accum(nmol(:,2), dble(bsp_charge(8)))     ! H2+  (bsp 8)
		call accum(nmol(:,3), dble(bsp_charge(9)))     ! H3+  (bsp 9)
		call accum(nmol(:,4), dble(bsp_charge(10)))    ! HeH+ (bsp 10)
	endif

	contains
		subroutine accum(vec, w)
		! Accumulate w*vec into ne in canonical bsp order.
		real*8, dimension(1-Ng:N+Ng), intent(in) :: vec
		real*8, intent(in) :: w
		ne = ne + w*vec
		end subroutine accum

	end subroutine calc_ne
	

	subroutine calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol,nox)
	! Calculate gas-particle number density, excluding free electrons.
	! Every species counts as ONE gas particle, so the base H/He/molecular
	! contribution is accumulated with unit weight in the canonical bsp order
	! of species_table; the optional nm adds the metal nuclei (all stages)
	! when the eos_metals policy is on.  bsp-order accumulation deliberately
	! fixes the FP add order (a golden re-snapshot decision).
	! The He 2^3S column is NOT an argument: it is an excited level of He I
	! (bsp_is_excited_level), so its gas particle is already the He I particle
	! counted through nhei, and adding it counted the triplet twice.
	! The optional nox holds the oxygen-chemistry carriers OH, H2O and CO
	! (bsp 11..13); each is one more gas particle.  Its oxygen and carbon
	! nuclei have already been removed from nm by the ionization solve, so
	! nothing is counted twice.

	integer :: im
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi,nhii
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiii
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in), optional :: nm
	real*8, dimension(1-Ng:N+Ng,4), intent(in), optional :: nmol  ! molecular
	real*8, dimension(1-Ng:N+Ng,3), intent(in), optional :: nox   ! OH H2O CO
	real*8, dimension(1-Ng:N+Ng), intent(out) :: n_tot

	n_tot = 0.0d0

	! Base H/He particles (bsp 1..5), one particle each; bsp 6 (HeTR) is an
	! excited level of bsp 3 (HeI) and is already inside it.
	call accum(nhi)                     ! HI    (bsp 1)
	call accum(nhii)                    ! HII   (bsp 2)
	if (thereis_He) then
		call accum(nhei)                ! HeI   (bsp 3, He 2^3S included)
		call accum(nheii)               ! HeII  (bsp 4)
		call accum(nheiii)              ! HeIII (bsp 5)
	endif

	if (present(nm) .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			n_tot = n_tot + nm(:,im)
		enddo
	endif

	! Each molecule is ONE gas particle (bsp 7..10).
	if (present(nmol)) then
		call accum(nmol(:,1))           ! H2   (bsp 7)
		call accum(nmol(:,2))           ! H2+  (bsp 8)
		call accum(nmol(:,3))           ! H3+  (bsp 9)
		call accum(nmol(:,4))           ! HeH+ (bsp 10)
	endif

	! Oxygen-chemistry carriers, one gas particle each (bsp 11..13).
	if (present(nox)) then
		call accum(nox(:,1))            ! OH   (bsp 11)
		call accum(nox(:,2))            ! H2O  (bsp 12)
		call accum(nox(:,3))            ! CO   (bsp 13)
	endif

	contains
		subroutine accum(vec)
		! Accumulate vec into n_tot (unit weight; one particle per species).
		real*8, dimension(1-Ng:N+Ng), intent(in) :: vec
		n_tot = n_tot + vec
		end subroutine accum

	end subroutine calc_ntot


	subroutine hydrogen_helium_nuclei_density(nhi,nhii,nhei,nheii,nheiii,   &
	                                          nh,nhe,nmol,nox)
	! Total hydrogen and helium NUCLEI densities [same units as the inputs].
	! Every carrier of an H or He nucleus is counted with its nucleus
	! multiplicity, so nh and nhe are conserved by the chemistry: H2 and H2+
	! carry two H nuclei, H3+ three, HeH+ one H and one He, OH one H and H2O
	! two.  CO carries none and is absent from the sum.  The He 2^3S column is
	! not an argument -- it is an excited level of He I (bsp_is_excited_level),
	! already inside nhei.
	! This is the ONE definition of the H/He nucleus totals: the ionization
	! equilibrium solve, its heating dump and the advection-corrected
	! post-process all call it, so the ionized fraction
	!   xion = n_e/(nh + nhe)
	! that Dalgarno, Yan & Liu (1999) section 7 define ("the number density
	! ratio of the electrons to the hydrogen and helium nuclei") and the H
	! nucleus density that scales the FUV/Lyman-Werner beam cannot drift apart
	! between them.
	! The optional nmol (cols 1 H2, 2 H2+, 3 H3+, 4 HeH+) and nox (cols 1 OH,
	! 2 H2O, 3 CO) are ignored unless the run carries that chemistry; a caller
	! that does not track those carriers simply omits them.  The statement
	! order is fixed: it is the FP add order the goldens were snapshotted with.

	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi,nhii
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiii
	real*8, dimension(1-Ng:N+Ng,4), intent(in), optional :: nmol  ! molecular
	real*8, dimension(1-Ng:N+Ng,3), intent(in), optional :: nox   ! OH H2O CO
	real*8, dimension(1-Ng:N+Ng), intent(out) :: nh,nhe

	nh  = nhi  + nhii
	nhe = nhei + nheii + nheiii

	if (present(nmol) .and. thereis_mol) then
		nh  = nh  + 2.0d0*(nmol(:,1) + nmol(:,2))                        &
		          + 3.0d0*nmol(:,3) + nmol(:,4)
		nhe = nhe + nmol(:,4)
	endif

	! OH carries one H nucleus and H2O two (bsp_nH of the species table).
	! Written as its own statement so the sum above stays bit-for-bit the one
	! a run without the oxygen chemistry evaluates.
	if (present(nox) .and. thereis_oxychem)                               &
		nh = nh + nox(:,1) + 2.0d0*nox(:,2)

	end subroutine hydrogen_helium_nuclei_density

	
	subroutine calc_rho(nhi,nhii,nhei,nheii,nheiii,n_out,nm,nmol,nox)
	! Calculate the total mass density (adimensional).
	! The base H/He/molecular mass is accumulated in the canonical bsp order
	! of species_table, each species weighted by bsp_mass [m_H units]; this
	! bsp-order weighted accumulation deliberately reorders the FP adds versus
	! the old factored 4.0*(nhei+...) form (a golden re-snapshot
	! decision).  The optional nm adds the metal mass (melem_A per
	! nucleus, all stages) when the eos_metals policy is on.  (HeH+ carries
	! 4.9715 m_H: its He nucleus is NOT in the nhei..nheiii free-He arrays.)
	! The He 2^3S column is NOT an argument: it is an excited level of He I
	! (bsp_is_excited_level), so its 3.9715 m_H are already the He I atom's mass
	! counted through nhei, and adding it put the triplet mass in rho twice.
	! The optional nox adds the oxygen-chemistry carriers OH (16.875 m_H),
	! H2O (17.875) and CO (27.793).  Those weights are the H mass plus the
	! metal block's own melem_A, and the oxygen and carbon they carry have
	! been removed from nm by the ionization solve, so the mass per nucleus
	! is the same either way and none of it is counted twice.

	integer :: im
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi,nhii
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiii
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in), optional :: nm
	real*8, dimension(1-Ng:N+Ng,4), intent(in), optional :: nmol  ! molecular
	real*8, dimension(1-Ng:N+Ng,3), intent(in), optional :: nox   ! OH H2O CO
	real*8, dimension(1-Ng:N+Ng), intent(out) :: n_out

	n_out = 0.0d0

	! Base H/He mass (bsp 1..5), weighted by bsp_mass; bsp 6 (HeTR) is an
	! excited level of bsp 3 (HeI) and is already inside it.
	call accum(nhi,    bsp_mass(1))     ! HI    (bsp 1)
	call accum(nhii,   bsp_mass(2))     ! HII   (bsp 2)
	if (thereis_He) then
		call accum(nhei,   bsp_mass(3))     ! HeI   (bsp 3, He 2^3S included)
		call accum(nheii,  bsp_mass(4))     ! HeII  (bsp 4)
		call accum(nheiii, bsp_mass(5))     ! HeIII (bsp 5)
	endif

	if (present(nm) .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			n_out = n_out + melem_A(mion_elem(im))*nm(:,im)
		enddo
	endif

	! Molecular mass (bsp 7..10): H2/H2+ = 2, H3+ = 3, HeH+ = 4.9715 m_H.
	if (present(nmol)) then
		call accum(nmol(:,1), bsp_mass(7))     ! H2   (bsp 7)
		call accum(nmol(:,2), bsp_mass(8))     ! H2+  (bsp 8)
		call accum(nmol(:,3), bsp_mass(9))     ! H3+  (bsp 9)
		call accum(nmol(:,4), bsp_mass(10))    ! HeH+ (bsp 10)
	endif

	! Oxygen-chemistry carriers (bsp 11..13).
	if (present(nox)) then
		call accum(nox(:,1), bsp_mass(11))     ! OH   (bsp 11)
		call accum(nox(:,2), bsp_mass(12))     ! H2O  (bsp 12)
		call accum(nox(:,3), bsp_mass(13))     ! CO   (bsp 13)
	endif

	contains
		subroutine accum(vec, w)
		! Accumulate w*vec into n_out in canonical bsp order.
		real*8, dimension(1-Ng:N+Ng), intent(in) :: vec
		real*8, intent(in) :: w
		n_out = n_out + w*vec
		end subroutine accum

	end subroutine calc_rho


	subroutine calc_column_dens(nhi,nhei,nheii,nheiTR,N1,N15,N2,NTR)
	! Calculates the column densities for given ionization profiles
	! 	by method of rectangles
	!
	! THE COLUMN STARTS AT THE OUTER FACE OF THE DOMAIN, r_edg(N). Nothing
	! is assumed to absorb beyond it: the stellar beam reaches the outer
	! face unattenuated, the statement every escape model with a finite
	! domain makes. The outer ghost cells N+1 .. N+Ng are a boundary
	! closure for the reconstruction stencil, not gas on the line of
	! sight; their column is zero (they see the unattenuated beam), and
	! they add no depth to the cells below. Counting them made the depth of
	! the domain depend on the ghost rule: a hydrostatic ghost past the
	! Roche L1 ridge of 1e41 cm^-3 (MEASURED, WASP-52 b at 10 R_p, frozen
	! f8a38866) shadowed every cell from the XUV. Cell j's column includes
	! its weighted width: N1(j)-N1(j+1) = n_HI(j) dr_j(j) R0 opa_pf(j).
    
	integer :: j
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiTR
   real*8 :: dr 
	real*8, dimension(1-Ng:N+Ng), intent(out) :: N1
	real*8, dimension(1-Ng:N+Ng), intent(out) :: N15,N2,NTR
	
	! Initialize outputs
	N1  = 0.0
	N15 = 0.0
	N2  = 0.0
	NTR = 0.0
	
	! Outer ghosts: zero column (initialized above). From the last physical
	! cell inward (opa_pf weights the opacity for the 'P' model; =1
	! otherwise); N1(N+1) = 0 starts the sum at the outer face.
	do j = N,1-Ng,-1

	    ! Spacing
	    dr = dr_j(j)*R0*opa_pf(j)

	      ! Evaluate new column densities by integration
	    N1(j)  = N1(j+1)  + nhi(j)*dr         ! HI

		if (thereis_He) then
			N15(j) = N15(j+1) + nhei(j)*dr	  					  ! HeI
			N2(j)  = N2(j+1)  + nheii(j)*dr       				  ! HeII
			if(thereis_HeITR) NTR(j) = NTR(j+1) + nheiTR(j)*dr	  ! HeI triplet
		endif

	enddo
	
	end subroutine calc_column_dens


	subroutine calc_column_dens_one(nsp, Ncol)
	! Column density of a single species (same rectangle rule, opa_pf
	! weighting and start at the outer face r_edg(N) as calc_column_dens,
	! whose header states why the outer ghosts carry no column). Used for
	! N_H2 (molecular).
	integer :: j
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nsp
	real*8, dimension(1-Ng:N+Ng), intent(out) :: Ncol
	Ncol = 0.0
	do j = N,1-Ng,-1
		Ncol(j) = Ncol(j+1) + nsp(j)*dr_j(j)*R0*opa_pf(j)
	enddo
	end subroutine calc_column_dens_one


	subroutine calc_column_dens_metals(nm, Nm_col)
	! Column densities for metals (same scheme as calc_column_dens,
	! starting at the outer face r_edg(N); the outer ghosts carry none).
	! nm(:,1:n_mion) holds the metal ion densities in canonical
	! species_table order; Nm_col(:,1:n_mphot) returns the column
	! density of each photo-ionizable metal ion in iphot order.

	integer :: j,i,k
	real*8, dimension(1-Ng:N+Ng,n_mion),  intent(in)  :: nm
	real*8 :: dr
	real*8, dimension(1-Ng:N+Ng,n_mphot), intent(out) :: Nm_col

	Nm_col = 0.0

	! From the last physical cell inward (opa_pf weights the opacity for
	! the 'P' model)
	do j = N,1-Ng,-1
	   dr = dr_j(j)*R0*opa_pf(j)
	   do i = 1,n_mion
	      if (.not. mion_isphot(i)) cycle
	      k = mion_iphot(i)
	      Nm_col(j,k) = Nm_col(j+1,k) + nm(j,i)*dr
	   enddo
	enddo

	end subroutine calc_column_dens_metals


	! Evaluate A(x)=(1-exp(-x))/x for x>=0; return its limit 1 at
	! x<=0. This one numerical form serves cell-mean absorption, locally
	! absorbed recombination radiation, one-face slab escape, and Sobolev
	! Ly-alpha escape; callers retain their distinct physical meanings.
	! For x<1e-8 use 1-x/2+x^2/6; for 1e-8<=x<1e-2 use the
	! Taylor series through x^8 in Horner form; otherwise use the closed
	! form. The series avoid cancellation of 1-exp(-x) and its staircase
	! in density-dependent optical depths. At x=1e-2 the closed-form
	! relative rounding error is <=2.2e-14 and the series truncation is
	! <3e-25. The measured steady-residual failure is in
	! ../EXHALE_v1.00/md/Update_EXHALE_stage3.md section 58.
	elemental double precision function one_minus_exp_over_x(x) result(g)
	real*8, intent(in) :: x
	if (x .le. 0.0d0) then
		g = 1.0d0
	else if (x .lt. 1.0d-8) then
		g = 1.0d0 - 0.5d0*x + x*x/6.0d0
	else if (x .lt. 1.0d-2) then
		g = 1.0d0 - x*(1.0d0/2.0d0 - x*(1.0d0/6.0d0 - x*(1.0d0/24.0d0    &
		  - x*(1.0d0/120.0d0 - x*(1.0d0/720.0d0 - x*(1.0d0/5040.0d0      &
		  - x*(1.0d0/40320.0d0 - x/362880.0d0)))))))
	else
		g = (1.0d0 - exp(-x))/x
	endif
	end function one_minus_exp_over_x


	! (1 - exp(-d))/d, the fraction of the photons entering a cell that the
	! cell absorbs, per unit of its own optical depth d.  This is the ONLY
	! definition of that quantity: the XUV beam takes it through
	! cell_mean_attenuation below, and the FUV bands take it from here for
	! the H2O and OH rates (water_photolysis.f90). The quotient is evaluated
	! by one_minus_exp_over_x above (branches, switch points and errors
	! there); the limit at a cell of no optical depth, d <= 0, is 1.
	elemental double precision function absorbed_fraction_per_unit_depth(d) &
	                                   result(fr)
	real*8, intent(in) :: d
	fr = one_minus_exp_over_x(d)
	end function absorbed_fraction_per_unit_depth


	! MEAN over one cell of the attenuation of the stellar beam, for a beam
	! whose optical depth is tau_out at the cell's star-ward face and
	! tau_out + dtau at its inner face.  This is the field the RATE of the
	! cell sees.
	!
	! WHY A MEAN AND NOT A FACE VALUE.  calc_column_dens* accumulate from the
	! top down, so N(j) already holds the WHOLE of cell j and the depth built
	! from it is the depth at that cell's INNER face.  The photoionization and
	! photoheating rates of the cell are what a particle experiences anywhere
	! inside it, averaged over it, and the beam falls across the cell by
	! exactly that cell's own optical depth.  A one-point rule at the inner
	! face is low by dtau/2 to first order, in one direction at every cell,
	! and by 42 per cent at dtau = 1, which is what a cell carries at the
	! ionization front.  The H2O and OH rates (water_photolysis.f90 sec. 3,
	! where the one-point rule was measured 30 per cent low in the Ly-alpha
	! band) and the Lyman-Werner rate (lyman_werner.f90) of the FUV beam
	! already take this mean; this is the same discretization for the XUV
	! beam and its absorbers.
	!
	! THE IDENTITY.  Inside a cell the absorber densities are uniform, which
	! is the rectangle rule the column integration itself uses, so the depth
	! runs linearly across the cell, tau(s) = tau_out + s dtau, and for the
	! pure exponential field
	!
	!   <exp(-tau)> = int_0^1 exp(-tau_out - s dtau) ds
	!               = exp(-tau_out) (1 - exp(-dtau))/dtau ,
	!
	! exact for a piecewise-constant absorber density at any grid spacing.
	! Multiplied by the cell's own dtau and summed over the column it
	! telescopes to 1 - exp(-tau_total): the photons the cells absorb are the
	! photons the beam loses, cell by cell and over the whole column.  The
	! inner-face rule satisfies neither statement.
	!
	! With a_tau>0, integrate exp(-tau)/(1+a_tau*tau) across the cell.
	! Only the first min(dtau,40) in optical depth is needed: monotonicity
	! of 1/(1+a_tau*tau) bounds the discarded tail relative to the retained
	! integral by exp(-40)/(1-exp(-40)) < 4.3e-18.
	! Eight-point Gauss-Legendre uses width h<=0.25/(1+b), where
	! b=a_tau/(1+a_tau*tau_left). For f=exp(-tau)/(1+a_tau*tau),
	! |f^(16)|/f <= 16!*(1+b)^16. The Gauss error coefficient is
	! (8!)^4/[17*(16!)^3]; hence each segment's relative quadrature error
	! is <= C*16!*0.25^16*exp(0.25)*1.25 < 1.4e-19, before roundoff.
	! These bounds require nonnegative optical depths and a_tau>0. The
	! existing clipping supplies nonnegative depths. a_tau<=0 retains the
	! exact exponential cell mean. Underflow and summation roundoff remain
	! floating-point limits, not errors covered by the quadrature bound.
	elemental double precision function cell_mean_attenuation(tau_out, dtau) &
	                                   result(f)
	real*8, intent(in) :: tau_out, dtau
	integer, parameter :: n_gl = 8
	! Gauss-Legendre nodes and weights on [0,1], rounded to binary64.
	real*8, parameter :: gl_s(n_gl) = (/                              &
	  0.01985507175123188416d0, 0.1016667612931866302d0,              &
	  0.2372337950418355071d0,  0.4082826787521750975d0,              &
	  0.5917173212478249025d0,  0.7627662049581644929d0,              &
	  0.8983332387068133698d0,  0.9801449282487681158d0 /)
	real*8, parameter :: gl_w(n_gl) = (/                              &
	  0.05061426814518812958d0, 0.1111905172266872353d0,              &
	  0.1568533229389436437d0,  0.1813418916891809915d0,              &
	  0.1813418916891809915d0,  0.1568533229389436437d0,              &
	  0.1111905172266872353d0,  0.05061426814518812958d0 /)
	real*8 :: tau_face, dt, extent, offset, width, tau_s, slope, geometry
	integer :: g

	tau_face = max(tau_out,0.0d0)
	dt = max(dtau,0.0d0)
	if (a_tau .le. 0.0d0) then
	   f = exp(-tau_face)*absorbed_fraction_per_unit_depth(dt)
	   return
	endif
	if (a_tau .le. 1.0d0) then
	   geometry = 1.0d0/(1.0d0+a_tau*tau_face)
	else
	   geometry = (1.0d0/a_tau)/(1.0d0/a_tau+tau_face)
	endif
	if (dt .eq. 0.0d0) then
	   f = exp(-tau_face)*geometry
	   return
	endif
	f = 0.0d0
	if (exp(-tau_face) .eq. 0.0d0) return
	extent = min(dt,40.0d0)
	offset = 0.0d0
	do while (offset .lt. extent)
	   if (a_tau .le. 1.0d0) then
	      slope = a_tau/(1.0d0+a_tau*(tau_face+offset))
	   else
	      slope = 1.0d0/(1.0d0/a_tau+tau_face+offset)
	   endif
	   width = min(0.25d0/(1.0d0+slope),extent-offset)
	   do g = 1,n_gl
	      tau_s = tau_face+offset+gl_s(g)*width
	      if (a_tau .le. 1.0d0) then
	         geometry = 1.0d0/(1.0d0+a_tau*tau_s)
	      else
	         geometry = (1.0d0/a_tau)/(1.0d0/a_tau+tau_s)
	      endif
	      f = f + gl_w(g)*(width/dt)*exp(-tau_s)*geometry
	   enddo
	   offset = offset+width
	enddo
	end function cell_mean_attenuation

	subroutine calc_mmw(nh,nhe,ne,mmw,nm)
	! Calculate the mean molecular weight for a certain ionization profile.
	! The H/He nucleus masses come from the species_table metadata
	! (bsp_mass(1) = HI = 1, bsp_mass(3) = HeI = 3.9715259 in m_H units; for
	! the atomic species the bsp position equals the f_sp column), so this
	! diagnostic weighs a helium atom exactly as calc_rho does.
	! The optional nm adds the metal mass and
	! metal nuclei under the same eos_metals policy as calc_rho/calc_ntot,
	! so the post-process temperature solve uses the same composition as
	! the main loop (the metal mass raises mmw by ~1%; omitting nm keeps
	! the legacy H/He-only diagnostic).

	integer :: im
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nh,nhe,ne
	real*8, dimension(1-Ng:N+Ng,n_mion), intent(in), optional :: nm
	real*8, dimension(1-Ng:N+Ng), intent(out) :: mmw
	real*8, dimension(1-Ng:N+Ng) :: mass_l, npart_l

	if (thereis_He) then
		mass_l  = bsp_mass(1)*nh + bsp_mass(3)*nhe
		npart_l = nh + nhe + ne
	else
		mass_l  = nh
		npart_l = nh + ne
	endif

	if (present(nm) .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			mass_l  = mass_l + melem_A(mion_elem(im))*nm(:,im)
			npart_l = npart_l + nm(:,im)
		enddo
	endif

	mmw = mass_l/npart_l

	end subroutine calc_mmw

	end module utils 
