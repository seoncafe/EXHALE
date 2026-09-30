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


   ! ------------------------------------------------------------------ !

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

   ! ------------------------------------------------------------------ !

   ! ------------------------------------------------------------------ !

   integer(8) function file_rolling_checksum(fname) result(ck)
   ! A CHECKSUM OF A FILE'S BYTES, for provenance only.
   !
   ! It is NOT a cryptographic digest and is deliberately not called one: two
   ! independent modular rolling sums,
   !     h_i <- mod( h_i * a_i + byte, 2^31 - 1 ),   a = (131, 8191),
   ! packed into one 62-bit integer. Both moduli are prime and every
   ! intermediate stays below 2^38, so the arithmetic is exact in integer(8)
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

   ! ------------------------------------------------------------------ !

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
           file_rolling_checksum('input.inp'),                              &
           ' ck_base=',   file_rolling_checksum('base.inp'),                &
           ' ck_metals=', file_rolling_checksum('metals.inp')
      write(unit,'(A,I0)') '# provenance: recon='//                        &
           trim(reconstruction_operator_label())//                          &
           ' base_bc='//trim(bmod)//' carrier='//trim(cmod)//              &
           ' restart_schema=3 resid_def=145 N=', N
   end subroutine write_provenance_header

   subroutine write_coupling_state_header(unit, state_id)
   ! THE LINE THAT SAYS WHAT PHYSICS THE STATE IN THIS FILE WAS PRODUCED
   ! UNDER, for every switch of the run that is NOT fixed by input.inp but
   ! changes while the run converges.
   !
   ! WHY. A restart re-reads the state but re-derives the switches from
   ! input.inp, and for a STAGED switch that is a different setting from the
   ! one the state was converged under. Measured: the WASP-121b steady
   ! solution was reached with the SvS85
   ! secondary ionization armed at step 2242, the restart re-staged it, the
   ! first ionization sweep moved the particle count by 8.1e-3 rather than
   ! 2.6e-4, and the layer left the root by a peak-to-peak 4.3 times the wind
   ! mass flux instead of 0.017.
   !
   ! WHAT IS IN IT. The census of every switch EXHALE_main changes inside the
   ! marching loop, so a reader of this line knows the whole run state and not
   ! just the part today's load_IC acts on:
   !
   !   sec_ion    the SvS85/Dalgarno photoelectron secondary-ionization
   !              coupling, as applied when the file was written (T/F), and
   !              the step it was armed at (-1 = never staged in)
   !   recon      the discrete reconstruction operator in force, which is not
   !              rec_method alone while the PLM -> WENO3 continuation of
   !              section 148 is armed: the label comes from
   !              reconstruction_operator_label() so that a file written
   !              inside a ramp says so
   !
   ! WHAT LEFT IT, section 152. `valve` (valve_eps) and `fluxconst`
   ! (base_flux_const) were the two switches the OLD base boundary chose for
   ! itself, and the characteristic boundary has neither: there is no one-way
   ! valve to smooth and no moving average standing in for a mass flux, so a
   ! header that reported them would be reporting settings the run cannot
   ! have. `load_IC`'s parser reads the fields it recognizes one at a time and
   ! ignores the rest, so a file written before section 152 still restores its
   ! sec_ion state and an older reader still finds `sec_ion` in a file written
   ! after it.
   !
   ! It is a '#' comment, so no numeric parse and no golden changes: the
   ! regression compares with grep -v '^ *#'.
   !
   ! The fields are formed by coupling_state_fields, which the provenance
   ! line of the derived (_adv) products writes as well, so the two lines
   ! cannot state different sets of fields.
   !
   ! state_id, present only for a generation of the pass state
   ! (publish_pass_state_generation, write_output.f90), is appended as the
   ! last field, state_id=<generation>, on both halves of the pair; load_IC
   ! refuses a pair whose halves state different identities or only one of
   ! them states one. Absent, the line is the one every other file carries.
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
      character(len=16) :: ptr
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
	! ------------------------------------------------------!

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
	
	! ------------------------------------------------------!

	subroutine calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol,nox)
	! Calculate the total atomic number density.
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

	! ------------------------------------------------------!

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

	! ------------------------------------------------------!
	
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

	! ------------------------------------------------------!

	subroutine calc_column_dens(nhi,nhei,nheii,nheiTR,N1,N15,N2,NTR)
	! Calculates the column densities for given ionization profiles
	! 	by method of rectangles
    
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
	
	! Outer point (opa_pf weights the opacity for the 'P' model; =1 otherwise)
	N1(N+Ng)  = dr_j(N+Ng)*R0*nhi(N+Ng)*opa_pf(N+Ng)
	if (thereis_He) then

		N15(N+Ng) = dr_j(N+Ng)*R0*nhei(N+Ng)*opa_pf(N+Ng)
		N2(N+Ng)  = dr_j(N+Ng)*R0*nheii(N+Ng)*opa_pf(N+Ng)
		if(thereis_HeITR) NTR(N+Ng) = dr_j(N+Ng)*R0*nheiTR(N+Ng)*opa_pf(N+Ng)

	endif

	do j = N+Ng-1,1-Ng,-1

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
	
	! End of subroutine
	end subroutine calc_column_dens

	! ------------------------------------------------------!

	subroutine calc_column_dens_one(nsp, Ncol)
	! Column density of a single species (same rectangle rule and
	! opa_pf weighting as calc_column_dens).  Used for N_H2 (molecular).
	integer :: j
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nsp
	real*8, dimension(1-Ng:N+Ng), intent(out) :: Ncol
	Ncol = 0.0
	Ncol(N+Ng) = dr_j(N+Ng)*R0*nsp(N+Ng)*opa_pf(N+Ng)
	do j = N+Ng-1,1-Ng,-1
		Ncol(j) = Ncol(j+1) + nsp(j)*dr_j(j)*R0*opa_pf(j)
	enddo
	end subroutine calc_column_dens_one

	! ------------------------------------------------------!

	subroutine calc_column_dens_metals(nm, Nm_col)
	! Column densities for metals (same scheme as calc_column_dens).
	! nm(:,1:n_mion) holds the metal ion densities in canonical
	! species_table order; Nm_col(:,1:n_mphot) returns the column
	! density of each photo-ionizable metal ion in iphot order.

	integer :: j,i,k
	real*8, dimension(1-Ng:N+Ng,n_mion),  intent(in)  :: nm
	real*8 :: dr
	real*8, dimension(1-Ng:N+Ng,n_mphot), intent(out) :: Nm_col

	Nm_col = 0.0

	! Outer point (opa_pf weights the opacity for the 'P' model)
	do i = 1,n_mion
	   if (.not. mion_isphot(i)) cycle
	   k = mion_iphot(i)
	   Nm_col(N+Ng,k) = dr_j(N+Ng)*R0*nm(N+Ng,i)*opa_pf(N+Ng)
	enddo

	do j = N+Ng-1,1-Ng,-1
	   dr = dr_j(j)*R0*opa_pf(j)
	   do i = 1,n_mion
	      if (.not. mion_isphot(i)) cycle
	      k = mion_iphot(i)
	      Nm_col(j,k) = Nm_col(j+1,k) + nm(j,i)*dr
	   enddo
	enddo

	end subroutine calc_column_dens_metals

	! ------------------------------------------------------!

	! (1 - exp(-x))/x for x >= 0, the ONE evaluation of that quotient in the
	! code. It is the same mathematical function behind several physical
	! quantities, each of which keeps its own name, limits and meaning and
	! is evaluated through this: the fraction of the photons entering a
	! cell that the cell absorbs per unit of its own optical depth
	! (absorbed_fraction_per_unit_depth below, and through it the XUV cell
	! mean and the FUV H2O/OH band rates), the same fraction for the locally
	! absorbed recombination radiation (util_ion_eq.f90,
	! absorbed_fraction_over_tau and absorbed_photon_fraction), the one-face
	! escape probability of a uniformly emitting slab (Cool_coeff.f90,
	! line_escape_probability_one_face) and the Sobolev escape probability
	! of Ly-alpha (lya_rt.f90, beta_sob). x <= 0 returns the limit 1.
	!
	! THREE BRANCHES, AND WHY THE CLOSED FORM STARTS ONLY AT x = 1e-2.
	! The closed form subtracts two numbers that approach each other as
	! x -> 0: exp(-x) is rounded to a double near 1, whose spacing is
	! eps/2 = 1.1e-16, so 1 - exp(-x) carries an absolute error of that
	! size and the quotient a relative error of about eps/x. That error is
	! not only an inaccuracy. As x moves continuously the rounded exp(-x)
	! jumps from one double to the next, so the closed form is a STAIRCASE
	! in x with relative steps of about eps/x, i.e. a discontinuous function
	! of the densities x is formed from. MEASURED at x ~ 1.5e-7 (the H I
	! threshold depth of an LHS 1140 b wind cell): steps of 4-5e-11 in the
	! photoionization rate every 1.7e-11 in n(H I), enough for the
	! composition iteration of the steady residual to find no fixed point
	! and alternate between two states on either side of a step.
	!   x < 1e-8          the three-term Taylor series 1 - x/2 + x^2/6,
	!                     truncated at x^3/24, i.e. 4.2e-26 at x = 1e-8;
	!   1e-8 <= x < 1e-2  the Taylor series through x^8, in Horner form,
	!                     truncated at x^9/10!, i.e. below 3e-25 at
	!                     x = 1e-2; its rounding is a few ulps and it is
	!                     smooth to the last bit;
	!   x >= 1e-2         the closed form, whose error eps/x is then at most
	!                     2.2e-14 and falls as x grows; exp(-x) underflows
	!                     harmlessly to 0 for large x, leaving 1/x.
	! The switch at 1e-2 is where the closed form's error has fallen to
	! 2.2e-14 while the series needs only nine terms to stay below 3e-25;
	! the three branches agree to 1e-14 at both switch points.
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

	! ------------------------------------------------------!

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

	! ------------------------------------------------------!

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
	! THE 2D RATE CORRECTION.  With a_tau > 0 the field carries the further
	! factor 1/(1 + a_tau tau), the geometry approximation selected by
	! appx_mth (parameters.f90), which is not exponential in the column, so
	! the product has no closed form.  The mean is then taken by composite
	! three-point Gauss-Legendre in the depth across the cell, on segments at
	! most seg_dtau wide; the rule is exact for polynomials of degree five, and
	! its error on a segment of width h is h^7 f^(6)/2016000, i.e. a relative
	! truncation of h^6/2016000, which is 1.2e-10 at seg_dtau = 0.25.
	! a_tau = 0, the default and the only value for which the field is the
	! pure exponential, takes the closed form above.
	elemental double precision function cell_mean_attenuation(tau_out, dtau) &
	                                   result(f)
	real*8, intent(in) :: tau_out, dtau
	! Three-point Gauss-Legendre on [0,1]: nodes (1 -+ sqrt(3/5))/2 and 1/2,
	! weights 5/18, 8/18, 5/18.
	integer, parameter :: n_gl = 3
	real*8, parameter :: gl_s(n_gl) = (/ 0.1127016653792583d0,             &
	                                     0.5d0,                            &
	                                     0.8872983346207417d0 /)
	real*8, parameter :: gl_w(n_gl) = (/ 5.0d0/18.0d0, 8.0d0/18.0d0,       &
	                                     5.0d0/18.0d0 /)
	real*8, parameter  :: seg_dtau  = 0.25d0
	integer, parameter :: nseg_max  = 256
	real*8  :: tau_face, dt, ds, s_lo, tau_s, acc
	integer :: nseg, i, g

	tau_face = max(tau_out, 0.0d0)
	dt       = max(dtau,    0.0d0)

	if (a_tau .le. 0.0d0) then
		f = exp(-tau_face)*absorbed_fraction_per_unit_depth(dt)
		return
	endif

	nseg = min(max(int(dt/seg_dtau) + 1, 1), nseg_max)
	ds   = 1.0d0/dble(nseg)
	acc  = 0.0d0
	do i = 0,nseg-1
		s_lo = dble(i)*ds
		do g = 1,n_gl
			tau_s = tau_face + (s_lo + gl_s(g)*ds)*dt
			acc   = acc + gl_w(g)*ds*exp(-tau_s)/(1.0d0 + a_tau*tau_s)
		enddo
	enddo
	f = acc

	end function cell_mean_attenuation

	! ------------------------------------------------------!

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

	! End of subroutine
	end subroutine calc_mmw

	! End of module
	end module utils 
