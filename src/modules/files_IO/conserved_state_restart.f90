      module conserved_state_restart
      ! THE EXACT CODE-UNIT STATE BESIDE THE DIMENSIONAL PAIR.
      !
      ! output/conserved_state.txt holds the conserved variables
      ! u = (rho, rho v, E) and every species fraction f_sp of every cell,
      ! ghosts included, in code units and binary64, written as ES25.17E3
      ! text (17 significant digits, which reproduces a binary64 value
      ! exactly on reading). Hydro_ioniz.txt and Ion_species.txt remain the
      ! state files: they are what analysis reads and what a restart reads
      ! when no exact file is present. The exact file only removes the
      ! rounding of the dimensional route (rho n0, v v0, p p0 written to 16
      ! digits, then the energy rebuilt from p through the caloric equation
      ! of state), so a restart from it starts at the bits the run ended on.
      !
      ! READ AS output/conserved_state_IC.txt, beside Hydro_ioniz_IC.txt and
      ! Ion_species_IC.txt. It is used only when every recorded condition
      ! matches the present run: format, shape (N, Ng, n_species), the
      ! checksums of the two dimensional files it was written with, the
      ! physical configuration fields of the restart metadata (reservoir and
      ! options), the base boundary model, the species order, the
      ! normalization constants, the physical input files, the gravity and
      ! irradiation fields, the cell and face radii, and a checksum of its own
      ! data. Any mismatch, a damaged file, or a load that changed an option
      ! or migrated a legacy token leaves the file unused with a printed
      ! notice, and the restart proceeds from the dimensional pair, which
      ! load_IC has already accepted. The checksums are rolling sums for
      ! detecting accidental mismatches (file_rolling_checksum, utilities.f90),
      ! not cryptographic digests.
      !
      ! WHAT IS RESTORED. Physical cells 1..N and the upper ghosts are taken
      ! from the file. The lower ghosts are set to cell 1 as placeholders and
      ! the composition of the lower ghosts is the one load_IC built from the
      ! run's reservoir: the lower boundary is derived again by Apply_BC from
      ! the restored interior.
      !
      ! WRITING never stops a run. A state that cannot be written exactly
      ! (a nonfinite value, a nonpositive density or internal energy in a
      ! physical cell, a negative species fraction, an unreadable input file,
      ! an I/O failure) is not written; any file already at that path is
      ! removed, because it would describe another state, and the reason is
      ! printed. The dimensional pair is then the restart.
      use, intrinsic :: iso_fortran_env, only: int64
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use global_parameters
      use species_table, only: n_mion, mion_name, mion_fsp,            &
           isp_HI, isp_HII, isp_HeI, isp_HeII, isp_HeIII, isp_HeTR,   &
           isp_H2, isp_H2p, isp_H3p, isp_HeHp, isp_OH, isp_H2O,       &
           isp_CO
      use utils, only: file_rolling_checksum
      use Conversion, only: admissible_conserved_cell
      use IC_load, only: build_restart_metadata, n_meta, meta_len,    &
           imeta_reservoir, imeta_options, ic_restart_schema_present,  &
           ic_option_change_applied, ic_legacy_factor_migrated
      use base_boundary, only: base_boundary_model_id,                 &
           base_ghost_composition_seed_id
      use lower_atmosphere_profile, only: lap_in_use, lap_file
      implicit none
      private
      public :: write_conserved_state, read_conserved_state
      public :: remove_conserved_state, conserved_state_not_read
      integer, parameter :: format_version = 2
      integer(int64), parameter :: hash_mod = 2147483647_int64
      integer, parameter :: line_capacity = 8192
      integer, parameter :: why_len = 256

      contains

      ! ------------------------------------------------------------------ !

      subroutine write_conserved_state(path, hydro_path, ion_path, u,    &
                                       f_sp, written)
      ! Write the exact state of (u, f_sp) paired with the dimensional files
      ! hydro_path and ion_path, which must already be written and closed.
      character(len=*), intent(in) :: path, hydro_path, ion_path
      real*8, intent(in) :: u(3,1-Ng:N+Ng)
      real*8, intent(in) :: f_sp(1-Ng:N+Ng,n_species)
      logical, intent(out), optional :: written
      character(len=line_capacity) :: line
      character(len=why_len) :: why
      character(len=80) :: fmt
      integer :: iu, ios, j, k
      integer(int64) :: h1, h2, pair_h, pair_i
      real*8 :: norm(13)
      character(len=meta_len) :: meta(n_meta)
      character(len=8) :: names(n_species)
      integer(int64) :: files(9), fields(5)

      if (present(written)) written = .false.
      why = ''
      call require_binary64(why)
      call require_writable_state(u, f_sp, why)
      call pair_checksums(hydro_path, ion_path, pair_h, pair_i, why)
      call current_normalization(norm, why)
      call current_file_checksums(files, why)
      call current_field_checksums(fields, why)
      if (len_trim(why) .gt. 0) then
         call skip_exact_write(path, why)
         return
      endif
      call build_restart_metadata(meta)
      call current_species_names(names)

      open(newunit=iu, file=path, status='replace', action='write',     &
           iostat=ios)
      if (ios .ne. 0) then
         call skip_exact_write(path, 'the file cannot be opened')
         return
      endif
      write(iu,'(A,I0)') 'EXHALE_CONSERVED_STATE ', format_version
      write(iu,'(A,3(1X,I0))') 'SHAPE', N, Ng, n_species
      write(iu,'(A,2(1X,I0))') 'PAIR_CHECKSUM', pair_h, pair_i
      do k = imeta_reservoir, imeta_options
         write(iu,'(A,I0,1X,A)') 'CONFIG ', k, trim(meta(k))
      enddo
      write(iu,'(A)') 'BOUNDARY '//base_boundary_model_id()//' '//     &
           base_ghost_composition_seed_id
      do k = 1, n_species
         write(iu,'(A,I0,1X,A)') 'SPECIES ', k, trim(names(k))
      enddo
      write(iu,'(A,13(1X,ES25.17E3))') 'NORMALIZATION', norm
      write(iu,'(A,9(1X,I0))') 'SOURCE_CHECKSUM', files
      write(iu,'(A,5(1X,I0))') 'FIELD_CHECKSUM', fields
      do j = 1-Ng, N+Ng
         write(iu,'(A2,1X,I8,1X,ES25.17E3)') 'R ', j, r(j)
      enddo
      do j = 0, N
         write(iu,'(A2,1X,I8,1X,ES25.17E3)') 'E ', j, r_edg(j)
      enddo
      h1 = 1_int64; h2 = 1_int64
      write(fmt,'("(A2,1X,I8,1X,",I0,"ES25.17E3)")') n_species
      do j = 1-Ng, N+Ng
         write(iu,'(A2,1X,I8,1X,3ES25.17E3)') 'U ', j, u(:,j)
         write(line,fmt) 'F ', j, f_sp(j,:)
         write(iu,'(A)') trim(line)
         do k = 1, 3
            call hash_real(h1, h2, u(k,j))
         enddo
         do k = 1, n_species
            call hash_real(h1, h2, f_sp(j,k))
         enddo
      enddo
      write(iu,'(A,2(1X,I0))') 'DATA_CHECKSUM', h1, h2
      write(iu,'(A)',iostat=ios) 'END_EXHALE_CONSERVED_STATE'
      if (ios .eq. 0) close(iu, iostat=ios)
      if (ios .ne. 0) then
         close(iu, iostat=ios)
         call skip_exact_write(path, 'the file could not be completed')
         return
      endif
      if (present(written)) written = .true.
      end subroutine write_conserved_state

      ! ------------------------------------------------------------------ !

      subroutine read_conserved_state(path, hydro_path, ion_path, u,     &
                                      f_sp, found)
      ! Restore the exact state from path when it belongs to the pair
      ! hydro_path / ion_path (already read by load_IC) and to this run's
      ! configuration. found = .false. leaves u and f_sp as they were on
      ! entry; the caller then converts the dimensional pair.
      character(len=*), intent(in) :: path, hydro_path, ion_path
      real*8, intent(inout) :: u(3,1-Ng:N+Ng)
      real*8, intent(inout) :: f_sp(1-Ng:N+Ng,n_species)
      logical, intent(out) :: found
      character(len=line_capacity) :: line, expected
      character(len=why_len) :: why
      character(len=80) :: fmt
      character(len=2) :: tag
      integer :: iu, ios, j, k, got_j
      integer(int64) :: h1, h2, pair_h, pair_i
      real*8 :: norm(13), ur(3), fr(n_species)
      real*8 :: u_read(3,1-Ng:N+Ng)
      real*8 :: f_read(1-Ng:N+Ng,n_species)
      character(len=meta_len) :: meta(n_meta)
      character(len=8) :: names(n_species)
      integer(int64) :: files(9), fields(5)
      logical :: exists

      found = .false.
      inquire(file=path, exist=exists)
      if (.not. exists) return

      ! A load that changed a named option, migrated a legacy token or read
      ! a pair without restart metadata produced a state of another
      ! configuration; the dimensional route is the one that carries it.
      if (ic_option_change_applied) then
         call conserved_state_not_read(path, 'the load applied an'//     &
              ' allowed option change ("Restart option change")')
         return
      endif
      if (ic_legacy_factor_migrated) then
         call conserved_state_not_read(path, 'the load migrated a'//     &
              ' legacy factor token')
         return
      endif
      if (.not. ic_restart_schema_present) then
         call conserved_state_not_read(path, 'the dimensional pair'//    &
              ' carries no restart metadata')
         return
      endif

      why = ''
      call require_binary64(why)
      call pair_checksums(hydro_path, ion_path, pair_h, pair_i, why)
      call current_normalization(norm, why)
      call current_file_checksums(files, why)
      call current_field_checksums(fields, why)
      if (len_trim(why) .gt. 0) then
         call conserved_state_not_read(path, why)
         return
      endif
      call build_restart_metadata(meta)
      call current_species_names(names)

      open(newunit=iu, file=path, status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
         call conserved_state_not_read(path, 'the file cannot be opened')
         return
      endif

      compare: block
      write(expected,'(A,I0)') 'EXHALE_CONSERVED_STATE ', format_version
      if (.not. next_line_is(iu, expected, 'format', why)) exit compare
      write(expected,'(A,3(1X,I0))') 'SHAPE', N, Ng, n_species
      if (.not. next_line_is(iu, expected, 'shape', why)) exit compare
      write(expected,'(A,2(1X,I0))') 'PAIR_CHECKSUM', pair_h, pair_i
      if (.not. next_line_is(iu, expected, 'dimensional pair', why))    &
         exit compare
      do k = imeta_reservoir, imeta_options
         write(expected,'(A,I0,1X,A)') 'CONFIG ', k, trim(meta(k))
         if (.not. next_line_is(iu, expected, 'physical configuration', &
                                why)) exit compare
      enddo
      expected = 'BOUNDARY '//base_boundary_model_id()//' '//         &
           base_ghost_composition_seed_id
      if (.not. next_line_is(iu, expected, 'boundary model', why))      &
         exit compare
      do k = 1, n_species
         write(expected,'(A,I0,1X,A)') 'SPECIES ', k, trim(names(k))
         if (.not. next_line_is(iu, expected, 'species order', why))    &
            exit compare
      enddo
      write(expected,'(A,13(1X,ES25.17E3))') 'NORMALIZATION', norm
      if (.not. next_line_is(iu, expected, 'normalization', why))       &
         exit compare
      write(expected,'(A,9(1X,I0))') 'SOURCE_CHECKSUM', files
      if (.not. next_line_is(iu, expected, 'physical input files', why)) &
         exit compare
      write(expected,'(A,5(1X,I0))') 'FIELD_CHECKSUM', fields
      if (.not. next_line_is(iu, expected,                              &
                             'gravity or irradiation fields', why))     &
         exit compare
      do j = 1-Ng, N+Ng
         write(expected,'(A2,1X,I8,1X,ES25.17E3)') 'R ', j, r(j)
         if (.not. next_line_is(iu, expected, 'cell grid', why))        &
            exit compare
      enddo
      do j = 0, N
         write(expected,'(A2,1X,I8,1X,ES25.17E3)') 'E ', j, r_edg(j)
         if (.not. next_line_is(iu, expected, 'face grid', why))        &
            exit compare
      enddo
      h1 = 1_int64; h2 = 1_int64
      write(fmt,'("(A2,1X,I8,1X,",I0,"ES25.17E3)")') n_species
      do j = 1-Ng, N+Ng
         read(iu,'(A)',iostat=ios) line
         if (ios .ne. 0 .or. len_trim(line) .ne. 12 + 3*25) then
            why = 'conserved row truncated or of wrong width'
            exit compare
         endif
         read(line,'(A2,1X,I8,1X,3ES25.17E3)',iostat=ios) tag, got_j, ur
         if (ios .ne. 0 .or. tag .ne. 'U ' .or. got_j .ne. j) then
            why = 'invalid conserved row'
            exit compare
         endif
         read(iu,'(A)',iostat=ios) line
         if (ios .ne. 0 .or. len_trim(line) .ne. 12 + n_species*25) then
            why = 'species row truncated or of wrong width'
            exit compare
         endif
         read(line,fmt,iostat=ios) tag, got_j, fr
         if (ios .ne. 0 .or. tag .ne. 'F ' .or. got_j .ne. j) then
            why = 'invalid species row'
            exit compare
         endif
         if (.not. all(ieee_is_finite(ur)) .or.                         &
             .not. all(ieee_is_finite(fr))) then
            why = 'nonfinite value'
            exit compare
         endif
         if (j .ge. 1 .and. j .le. N) then
            if (.not. admissible_conserved_cell(ur)) then
               why = 'nonpositive density or internal energy in a'//    &
                     ' physical cell'
               exit compare
            endif
         endif
         if (any(fr .lt. 0.0d0)) then
            why = 'negative species fraction'
            exit compare
         endif
         u_read(:,j) = ur
         f_read(j,:) = fr
         do k = 1, 3
            call hash_real(h1, h2, ur(k))
         enddo
         do k = 1, n_species
            call hash_real(h1, h2, fr(k))
         enddo
      enddo
      write(expected,'(A,2(1X,I0))') 'DATA_CHECKSUM', h1, h2
      if (.not. next_line_is(iu, expected, 'data checksum', why))       &
         exit compare
      if (.not. next_line_is(iu, 'END_EXHALE_CONSERVED_STATE',          &
                             'end marker', why)) exit compare
      read(iu,'(A)',iostat=ios) line
      if (ios .ge. 0) why = 'data after the end marker'
      end block compare
      close(iu)
      if (len_trim(why) .gt. 0) then
         call conserved_state_not_read(path, why)
         return
      endif

      u(:,1:N+Ng)    = u_read(:,1:N+Ng)
      f_sp(1:N+Ng,:) = f_read(1:N+Ng,:)
      do j = 1-Ng, 0
         u(:,j) = u(:,1)
      enddo
      found = .true.
      write(*,'(A)') ' (conserved_state_restart) exact code-unit state'//  &
           ' read from '//trim(path)//': physical cells and upper'//       &
           ' ghosts restored, lower ghosts rebuilt by the boundary.'
      end subroutine read_conserved_state

      ! ------------------------------------------------------------------ !

      subroutine conserved_state_not_read(path, why)
      ! The notice of an exact state file that is present and not used.
      character(len=*), intent(in) :: path, why
      logical :: exists
      inquire(file=path, exist=exists)
      if (.not. exists) return
      write(*,'(A)') ' (conserved_state_restart) ********************'//   &
           '****************************************'
      write(*,'(A)') ' (conserved_state_restart) NOTICE: '//trim(path)//   &
           ' is IGNORED: '//trim(why)//'.'
      write(*,'(A)') ' (conserved_state_restart) The state is restarted'// &
           ' from the dimensional pair (Hydro_ioniz_IC.txt,'//            &
           ' Ion_species_IC.txt), not bitwise.'
      write(*,'(A)') ' (conserved_state_restart) ********************'//   &
           '****************************************'
      end subroutine conserved_state_not_read

      ! ------------------------------------------------------------------ !

      subroutine remove_conserved_state(path)
      ! Remove an exact state file that no longer describes the dimensional
      ! pair beside it. A file that cannot be removed is reported; a later
      ! restart still refuses it, since its pair checksum names other files.
      character(len=*), intent(in) :: path
      logical :: exists
      integer :: iu, ios
      inquire(file=path, exist=exists)
      if (.not. exists) return
      open(newunit=iu, file=path, status='old', iostat=ios)
      if (ios .eq. 0) close(iu, status='delete', iostat=ios)
      if (ios .ne. 0)                                                    &
         write(*,'(A)') ' (conserved_state_restart) WARNING: '//         &
              trim(path)//' describes another state and could not be'//  &
              ' removed; a restart will not use it (pair checksum).'
      end subroutine remove_conserved_state

      ! ------------------------------------------------------------------ !

      subroutine skip_exact_write(path, why)
      character(len=*), intent(in) :: path, why
      call remove_conserved_state(path)
      write(*,'(A)') ' (conserved_state_restart) WARNING: exact state'//   &
           ' not written to '//trim(path)//': '//trim(why)//'. The'//     &
           ' dimensional pair is the restart of this state.'
      end subroutine skip_exact_write

      ! ------------------------------------------------------------------ !

      subroutine pair_checksums(hydro_path, ion_path, ck_h, ck_i, why)
      character(len=*), intent(in) :: hydro_path, ion_path
      integer(int64), intent(out) :: ck_h, ck_i
      character(len=*), intent(inout) :: why
      ck_h = file_rolling_checksum(hydro_path)
      ck_i = file_rolling_checksum(ion_path)
      if (len_trim(why) .gt. 0) return
      if (ck_h .lt. 0_int64 .or. ck_i .lt. 0_int64)                      &
         why = 'the dimensional pair is missing or unreadable'
      end subroutine pair_checksums

      ! ------------------------------------------------------------------ !

      subroutine current_normalization(x, why)
      ! The code-unit scales and the planetary parameters a code-unit state
      ! is expressed in.
      real*8, intent(out) :: x(13)
      character(len=*), intent(inout) :: why
      x = [R0,n0,v0,p0,T0,mu,RJ,kb_erg,Mp,Mstar,a_orb,b0,J_XUV]
      if (len_trim(why) .gt. 0) return
      if (.not. all(ieee_is_finite(x)))                                  &
         why = 'nonfinite normalization or planetary parameter'
      end subroutine current_normalization

      ! ------------------------------------------------------------------ !

      subroutine current_file_checksums(ck, why)
      ! The input files that define the physical problem: the stellar
      ! spectrum, the lower-atmosphere profile or base.inp, metals.inp, the
      ! opacity configuration (its parameters, not the path of its tables)
      ! and the opacity tables themselves.
      integer(int64), intent(out) :: ck(9)
      character(len=*), intent(inout) :: why
      ck = -1_int64
      if (do_read_sed) then
         if (allocated(sed_file)) then
            if (len_trim(sed_file) .gt. 0)                               &
               ck(1) = file_rolling_checksum(sed_file)
         endif
         if (ck(1) .lt. 0_int64 .and. len_trim(why) .eq. 0)              &
            why = 'the spectrum file is unreadable'
      endif
      if (lap_in_use) then
         ck(2) = file_rolling_checksum(trim(lap_file))
         if (ck(2) .lt. 0_int64 .and. len_trim(why) .eq. 0)              &
            why = 'the lower-atmosphere profile is unreadable'
      else
         ck(3) = file_rolling_checksum('base.inp')
      endif
      ck(4) = file_rolling_checksum('metals.inp')
      ck(5) = opacity_configuration_checksum(why)
      if (allocated(opa_file_HI)) ck(6) = file_rolling_checksum(opa_file_HI)
      if (allocated(opa_file_HeI)) ck(7) = file_rolling_checksum(opa_file_HeI)
      if (allocated(opa_file_HeII)) ck(8) = file_rolling_checksum(opa_file_HeII)
      if (allocated(opa_file_HeITR))                                     &
         ck(9) = file_rolling_checksum(opa_file_HeITR)
      end subroutine current_file_checksums

      ! ------------------------------------------------------------------ !

      integer(int64) function opacity_configuration_checksum(why) result(ck)
      character(len=*), intent(inout) :: why
      integer(int64) :: h1, h2, model_byte
      real*8 :: coefficients(7)
      integer :: k
      h1 = 1_int64; h2 = 1_int64
      model_byte = int(iachar(opacity_model),int64)
      h1 = mod(h1*131_int64 + model_byte, hash_mod)
      h2 = mod(h2*8191_int64 + model_byte, hash_mod)
      coefficients = [opa_const_HI,opa_const_HeI,opa_const_HeII,          &
           opa_const_HeITR,opa_pb_factor,opa_pb_exponent,opa_pb_pivot]
      if (.not. all(ieee_is_finite(coefficients)) .and.                  &
          len_trim(why) .eq. 0) why = 'nonfinite opacity coefficient'
      do k = 1, size(coefficients)
         call hash_real(h1,h2,coefficients(k))
      enddo
      ck = h1*hash_mod + h2
      end function opacity_configuration_checksum

      ! ------------------------------------------------------------------ !

      subroutine current_field_checksums(ck, why)
      ! The gravitational potential at centers and interfaces and the
      ! irradiation field with its energy grid.
      integer(int64), intent(out) :: ck(5)
      character(len=*), intent(inout) :: why
      ck = -1_int64
      if (allocated(Gphi_c)) ck(1) = array_checksum(Gphi_c, why)
      if (allocated(Gphi_i)) ck(2) = array_checksum(Gphi_i, why)
      if (allocated(F_XUV))  ck(3) = array_checksum(F_XUV, why)
      if (allocated(e_v))    ck(4) = array_checksum(e_v, why)
      if (allocated(de_v))   ck(5) = array_checksum(de_v, why)
      end subroutine current_field_checksums

      ! ------------------------------------------------------------------ !

      integer(int64) function array_checksum(x, why) result(ck)
      real*8, intent(in) :: x(:)
      character(len=*), intent(inout) :: why
      integer(int64) :: h1, h2
      integer :: k
      h1 = 1_int64; h2 = 1_int64
      if (.not. all(ieee_is_finite(x)) .and. len_trim(why) .eq. 0)       &
         why = 'nonfinite gravity or irradiation field'
      do k = 1, size(x)
         call hash_real(h1,h2,x(k))
      enddo
      ck = h1*hash_mod + h2
      end function array_checksum

      ! ------------------------------------------------------------------ !

      subroutine current_species_names(names)
      ! The species order of f_sp, by name, as the run's species table
      ! defines it.
      character(len=8), intent(out) :: names(n_species)
      integer :: k
      names = ''
      names(isp_HI)='HI'; names(isp_HII)='HII'
      names(isp_HeI)='HeI'; names(isp_HeII)='HeII'
      names(isp_HeIII)='HeIII'; names(isp_HeTR)='HeITR'
      do k = 1, n_mion
         names(mion_fsp(k)) = trim(mion_name(k))
      enddo
      names(isp_H2)='H2'; names(isp_H2p)='H2p'
      names(isp_H3p)='H3p'; names(isp_HeHp)='HeHp'
      names(isp_OH)='OH'; names(isp_H2O)='H2O'; names(isp_CO)='CO'
      ! Every slot of f_sp is named above; an unnamed one is a species
      ! table this routine does not describe, which is a code error.
      if (any(names .eq. '')) error stop                                 &
         'conserved_state_restart: a species of f_sp has no name'
      end subroutine current_species_names

      ! ------------------------------------------------------------------ !

      subroutine require_writable_state(u, f_sp, why)
      ! The state a reader would accept: finite everywhere, admissible
      ! (positive density and internal energy) in the physical cells, and
      ! no negative species fraction in any row.
      real*8, intent(in) :: u(3,1-Ng:N+Ng)
      real*8, intent(in) :: f_sp(1-Ng:N+Ng,n_species)
      character(len=*), intent(inout) :: why
      integer :: j
      if (len_trim(why) .gt. 0) return
      if (.not. all(ieee_is_finite(u)) .or.                             &
          .not. all(ieee_is_finite(f_sp))) then
         why = 'nonfinite value in the state'
         return
      endif
      do j = 1, N
         if (.not. admissible_conserved_cell(u(:,j))) then
            write(why,'(A,I0)') 'nonpositive density or internal'//     &
                 ' energy in physical cell ', j
            return
         endif
      enddo
      do j = 1-Ng, N+Ng
         if (any(f_sp(j,:) .lt. 0.0d0)) then
            write(why,'(A,I0)') 'negative species fraction in row ', j
            return
         endif
      enddo
      end subroutine require_writable_state

      ! ------------------------------------------------------------------ !

      logical function next_line_is(unit, expected, description, why)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: expected, description
      character(len=*), intent(inout) :: why
      character(len=line_capacity) :: got
      integer :: ios
      read(unit,'(A)',iostat=ios) got
      next_line_is = (ios .eq. 0 .and. trim(got) .eq. trim(expected))
      if (next_line_is) return
      if (ios .ne. 0) then
         why = 'file truncated at '//description
      else
         why = 'differs from this run in '//description
      endif
      end function next_line_is

      ! ------------------------------------------------------------------ !

      subroutine hash_real(h1, h2, x)
      ! The two rolling sums of file_rolling_checksum (utilities.f90) over
      ! the eight bytes of one binary64 value.
      integer(int64), intent(inout) :: h1, h2
      real*8, intent(in) :: x
      integer(int64) :: bits, byte
      integer :: k
      bits = transfer(x, bits)
      do k = 0, 7
         byte = int(ibits(bits, 8*k, 8),int64)
         h1 = mod(h1*131_int64 + byte, hash_mod)
         h2 = mod(h2*8191_int64 + byte, hash_mod)
      enddo
      end subroutine hash_real

      ! ------------------------------------------------------------------ !

      subroutine require_binary64(why)
      character(len=*), intent(inout) :: why
      if (len_trim(why) .gt. 0) return
      if (storage_size(0.0d0) .ne. 64 .or. radix(0.0d0) .ne. 2 .or.     &
          digits(0.0d0) .ne. 53) why = 'real*8 is not IEEE binary64'
      end subroutine require_binary64

      end module conserved_state_restart
