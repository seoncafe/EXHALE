      subroutine write_carrier_balance_record(label)
      ! Preserve the carrier terms written by the immediately preceding
      ! verdict-recording certification. The one-sweep comparison uses
      ! record_verdict=.false. and does not write the generic file.
      ! This copies bytes only; it does not evaluate a carrier row or
      ! change the conserved, composition, or carrier module state.
      use, intrinsic :: iso_fortran_env, only: int64
      character(len=*), intent(in) :: label
      character(len=8) :: enabled
      character(len=64) :: destination
      character(len=8192) :: payload
      integer :: input_unit, output_unit, ios, ios_close, count_now
      integer(int64) :: bytes_left
      logical :: source_exists

      ! EXHALE_CARRIER_ROW_TERMS: Set to 1 to write carrier-row terms for this recorded state. Unset
      ! omits the record; diagnostic only.
      call get_environment_variable('EXHALE_CARRIER_ROW_TERMS', enabled)
      if (trim(enabled) .ne. '1') return
      select case (trim(label))
      case ('stored_pair')
         destination = trim(output_dir)//'carrier_row_terms_stored_pair.txt'
      case ('final_certification')
         destination = trim(output_dir)//'carrier_row_terms_final_certification.txt'
      case default
         write(*,'(A,A)') ' (EXHALE_main) unknown carrier balance record: ', &
              trim(label)
         stop 1
      end select

      bytes_left = 0_int64
      inquire(file=trim(output_dir)//'carrier_row_terms.txt', exist=source_exists,    &
              size=bytes_left, iostat=ios)
      if (ios .ne. 0 .or. .not. source_exists .or. bytes_left .le. 0) then
         write(*,'(A,A)') ' (EXHALE_main) missing carrier terms for ',     &
              trim(label)
         stop 1
      endif
      open(newunit=input_unit, file=trim(output_dir)//'carrier_row_terms.txt',        &
           access='stream', form='unformatted', status='old',              &
           action='read', iostat=ios)
      if (ios .ne. 0) then
         write(*,'(A)') ' (EXHALE_main) cannot read carrier row terms'
         stop 1
      endif
      open(newunit=output_unit, file=trim(destination),                    &
           access='stream', form='unformatted', status='replace',          &
           action='write', iostat=ios)
      if (ios .ne. 0) then
         close(input_unit)
         write(*,'(A,A)') ' (EXHALE_main) cannot write ', trim(destination)
         stop 1
      endif

      do while (bytes_left .gt. 0_int64)
         count_now = int(min(bytes_left, int(len(payload),int64)))
         read(input_unit, iostat=ios) payload(1:count_now)
         if (ios .ne. 0) exit
         write(output_unit, iostat=ios) payload(1:count_now)
         if (ios .ne. 0) exit
         bytes_left = bytes_left - int(count_now,int64)
      enddo
      close(input_unit)
      close(output_unit, iostat=ios_close)
      if (ios .ne. 0 .or. ios_close .ne. 0 .or.                      &
          bytes_left .ne. 0_int64) then
         write(*,'(A,A)') ' (EXHALE_main) incomplete carrier balance ',   &
              trim(destination)
         stop 1
      endif
      end subroutine write_carrier_balance_record
