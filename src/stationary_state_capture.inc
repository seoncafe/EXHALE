      subroutine write_stationary_stage_restart(pass, stage)
      ! The exact (u, f_sp) held at one named alternation stage. This is a
      ! diagnostic restart state, never a stationarity claim. The three
      ! files are complete only when manifest.txt exists.
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
      use conserved_state_restart, only: write_conserved_state
      integer, intent(in) :: pass
      character(len=*), intent(in) :: stage
      real*8 :: W_stage(3,1-Ng:N+Ng), T_stage(1-Ng:N+Ng)
      real*8 :: nhi_stage(1-Ng:N+Ng), nhii_stage(1-Ng:N+Ng)
      real*8 :: nhei_stage(1-Ng:N+Ng), nheii_stage(1-Ng:N+Ng)
      real*8 :: nheiii_stage(1-Ng:N+Ng), nheiTR_stage(1-Ng:N+Ng)
      real*8 :: ne_stage(1-Ng:N+Ng), ntot_stage(1-Ng:N+Ng)
      real*8 :: nm_stage(1-Ng:N+Ng,n_mion)
      real*8 :: mrow_none(1-Ng:N+Ng), npart1_hold
      real*8, allocatable :: nmol_hold(:,:), nox_hold(:,:)
      type(caloric_mixture_record) :: mixture_hold
      character(len=96) :: stage_id
      character(len=256) :: dir, hydro_path, ion_path, exact_path
      character(len=200) :: statement
      character(len=len(state_certification_reason)) :: reason_hold
      logical :: certified_hold, exact_written
      integer :: ios_h, ios_i, ios_m, unit_m, j
      integer(c_int) :: rc
      interface
         function stage_mkdir(path, mode) bind(C, name='mkdir') result(status)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: path
            integer(c_int), value :: mode
            integer(c_int) :: status
         end function stage_mkdir
         function stage_rename(old_path, new_path) bind(C, name='rename') &
              result(status)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: old_path
            character(kind=c_char), dimension(*), intent(in) :: new_path
            integer(c_int) :: status
         end function stage_rename
      end interface

      ! This also prevents a future caller from placing arbitrary text in
      ! a path. All current calls use these four stage names.
      select case (trim(stage))
      case ('hydro_entry', 'hydro_return', 'certification',              &
            'composition_return')
      case default
         write(*,'(A)') ' (EXHALE_main) unknown stage; restart omitted'
         return
      end select
      do j = 1-Ng, N+Ng
         if (.not. admissible_conserved_cell(u(:,j))) then
            write(*,'(A,I0,A,A)') ' (EXHALE_main) state at pass ', pass, &
                 ' is not a gas; restart omitted at ', trim(stage)
            return
         endif
      enddo

      if (len_trim(run_identity) .eq. 0) call establish_run_identity()
      write(stage_id,'(A,A,A,I0,A,I0,A,A)') 'r', trim(run_identity),     &
           '_s', n_stationary_solves_run, '_p', pass, '_', trim(stage)
      dir = trim(output_dir)//'stage_restart/'//trim(stage_id)
      rc = stage_mkdir(trim(output_dir)//'stage_restart'//c_null_char, 493_c_int)
      rc = stage_mkdir(trim(dir)//c_null_char, 493_c_int)
      if (rc .ne. 0) then
         write(*,'(A,A)') ' (EXHALE_main) cannot create stage restart ', &
              trim(dir)
         return
      endif
      hydro_path = trim(dir)//'/Hydro_ioniz.txt'
      ion_path   = trim(dir)//'/Ion_species.txt'
      exact_path = trim(dir)//'/conserved_state.txt'

      ! get_species_densities updates the caloric mixture and the base
      ! particle count. The species writer updates nmol_eq and nox_eq.
      ! Restore each value after forming the local dimensional view.
      call hold_caloric_mixture(mixture_hold)
      npart1_hold = n_part_cell1
      certified_hold = state_is_certified
      reason_hold = state_certification_reason
      if (thereis_mol)     nmol_hold = nmol_eq
      if (thereis_oxychem) nox_hold  = nox_eq
      nhei_stage = 0.0d0; nheii_stage = 0.0d0
      nheiii_stage = 0.0d0; nheiTR_stage = 0.0d0
      call get_species_densities(u(1,:), f_sp, nhi_stage, nhii_stage,   &
           nhei_stage, nheii_stage, nheiii_stage, nheiTR_stage,         &
           nm_stage, ne_stage, ntot_stage)
      call U_to_W(u, W_stage)
      call comp_T_from_p(W_stage(3,:), ntot_stage, ne_stage, T_stage)
      call molecular_carrier_densities_from_state(W_stage(1,:), f_sp)
      mrow_none = 0.0d0
      write(statement,'(A,I0,A,I0,A,A,A)') 'stage_snapshot solve=',      &
           n_stationary_solves_run, ' pass=', pass, ' stage=',           &
           trim(stage), ' (not certified)'
      call set_state_certified(.false., 'stage_snapshot')
      call write_hydro_state_file(trim(hydro_path), 'eq', W_stage(1,:), &
           W_stage(2,:), W_stage(3,:), T_stage, heat, cool, mrow_none,  &
           pass_statement=trim(statement), ios_open=ios_h,             &
           state_id=trim(stage_id))
      ios_i = -1
      if (ios_h .eq. 0) call write_species_state_file(trim(ion_path),   &
           'eq', nhi_stage, nhii_stage, nhei_stage, nheii_stage,        &
           nheiii_stage, nheiTR_stage, nm_stage,                         &
           pass_statement=trim(statement), ios_open=ios_i,             &
           state_id=trim(stage_id))
      exact_written = .false.
      if (ios_h .eq. 0 .and. ios_i .eq. 0) call write_conserved_state(   &
           trim(exact_path), trim(hydro_path), trim(ion_path), u, f_sp, &
           exact_written)
      if (thereis_mol)     nmol_eq = nmol_hold
      if (thereis_oxychem) nox_eq  = nox_hold
      call set_state_certified(certified_hold, trim(reason_hold))
      call restore_caloric_mixture(mixture_hold)
      n_part_cell1 = npart1_hold

      if (.not. exact_written) then
         write(*,'(A,A)') ' (EXHALE_main) incomplete stage restart ',  &
              trim(dir)
         return
      endif
      open(newunit=unit_m, file=trim(dir)//'/manifest.part',           &
           status='new', action='write', iostat=ios_m)
      if (ios_m .ne. 0) then
         write(*,'(A,A)') ' (EXHALE_main) cannot complete stage restart ', &
              trim(dir)
         return
      endif
      write(unit_m,'(A)') '# EXHALE stage restart; not certified'
      write(unit_m,'(A)') 'state_id='//trim(stage_id)
      write(unit_m,'(A)') 'run_identity='//trim(run_identity)
      write(unit_m,'(A,I0)') 'solve=', n_stationary_solves_run
      write(unit_m,'(A,I0)') 'pass=', pass
      write(unit_m,'(A)') 'stage='//trim(stage)
      write(unit_m,'(A)') 'hydro_file=Hydro_ioniz.txt'
      write(unit_m,'(A)') 'species_file=Ion_species.txt'
      write(unit_m,'(A)') 'conserved_state=conserved_state.txt'
      write(unit_m,'(A)') 'certified=F'
      close(unit_m, iostat=ios_m)
      if (ios_m .ne. 0) then
         write(*,'(A,A)') ' (EXHALE_main) stage manifest close failed ', &
              trim(dir)
      else
         rc = stage_rename(trim(dir)//'/manifest.part'//c_null_char,    &
                           trim(dir)//'/manifest.txt'//c_null_char)
         if (rc .eq. 0) then
            write(*,'(A,A)') ' (EXHALE_main) stage restart ', trim(dir)
         else
            write(*,'(A,A)') ' (EXHALE_main) manifest rename failed ',  &
                 trim(dir)
         endif
      endif
      end subroutine write_stationary_stage_restart
