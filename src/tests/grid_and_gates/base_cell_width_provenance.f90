      program base_cell_width_provenance
      ! WHAT IS ASSERTED
      !   1. round_trip_decimal writes a binary64 value as a decimal string
      !      that a list-directed read returns to the SAME BITS, for both
      !      signs, for zero and over the representable exponent range. The
      !      comparison is on the bit pattern (transfer to integer*8), not on
      !      a relative difference: a serializer that loses one ulp is not a
      !      round-trip serializer.
      !   2. The string never overflows its field. A field one character too
      !      narrow (ES24.17E3, which holds a positive value but not the sign
      !      of a negative one) fills with asterisks, and asterisks parse as
      !      nothing at all; that failure is shown as a diagnostic below so
      !      the width of the field is a measured requirement.
      !   3. dr_base_default, the width an input.inp without the
      !      "Base grid [dr,cells]:" key builds its grid on, carries exactly
      !      the bits a list-directed read of "2.0e-4" returns, so the key
      !      spelled the way a reader writes the default states the default.
      !   4. The pinned width "1.9999999494757503e-4", which every input run
      !      before 2026-09-19 without the key now states, reads back to the
      !      bits of the default-real literal 2.0e-4 widened to double: the
      !      width those runs were made on. Both are checked on the bits and
      !      not assumed, because each holds only if the read rounds
      !      correctly.
      !
      ! The serializer under test is the PRODUCTION one (module setup_report),
      ! the one that writes base_cell_width_Rp into EXHALE_resolved.out.

      use global_parameters, only: dr_base_default
      use setup_report,      only: round_trip_decimal

      implicit none

      integer, parameter :: nval = 13
      real*8    :: vals(nval)
      real*8    :: v, w
      integer*8 :: iv, iw
      integer   :: k, nfail, ios
      character(len=26) :: s
      character(len=24) :: narrow
      character(len=64) :: tag
      character(len=32) :: spelled

      nfail = 0

      vals = (/  0.0d0,                                                      &
                 1.0d0,              -1.0d0,                                 &
                 1.0d-4,             -1.0d-4,                                &
                 2.0d-4,             -2.0d-4,                                &
                 1.0d-300,           -1.0d-300,                              &
                 1.0d300,            -1.0d300,                               &
                 dr_base_default,    -dr_base_default /)

      do k = 1, nval
         v = vals(k)
         s = round_trip_decimal(v)
         if (index(s, '*') .gt. 0) then
            write(tag,'(A,I0,A)') 'round_trip_field_holds_value[', k, ']'
            call verdict_i(trim(tag), 1, 0, nfail)
            write(*,'(A,A,A)') '     printed [', trim(s), ']'
            cycle
         endif
         w   = 0.0d0
         ios = 0
         read(s, *, iostat=ios) w
         iv = transfer(v, iv)
         iw = transfer(w, iw)
         write(tag,'(A,I0,A)') 'round_trip_same_bits[', k, ']'
         if (ios .ne. 0) then
            call verdict_i(trim(tag), ios, 0, nfail)
         else
            call verdict_i8(trim(tag), iw, iv, nfail)
         endif
         write(*,'(A,A,A,Z16.16,A,Z16.16)') '     [', trim(s),               &
              ']  read=', iw, '  written=', iv
      enddo

      ! Why the field is 26 wide and not 24: the same value, same digits, in
      ! a field that cannot hold the sign.
      write(narrow,'(ES24.17E3)') -dr_base_default
      write(*,'(A,A,A)') 'DIAGNOSTIC ES24.17E3 of a negative value: [',      &
           narrow, ']'

      ! The default width is the double a list-directed read of 2.0e-4
      ! returns, the read input_read applies to the key.
      spelled = '2.0e-4'
      read(spelled, *) w
      iv = transfer(w, iv)
      iw = transfer(dr_base_default, iw)
      call verdict_i8('dr_base_default_is_the_read_of_2.0e-4',               &
                      iw, iv, nfail)
      write(*,'(A,A)') '     dr_base_default = ',                            &
           trim(round_trip_decimal(dr_base_default))

      ! The pinned width reads back to the width of the runs made before
      ! 2026-09-19: the default-real literal 2.0e-4, widened on assignment.
      v  = 2.0e-4
      iv = transfer(v, iv)
      spelled = '1.9999999494757503e-4'
      read(spelled, *) w
      iw = transfer(w, iw)
      call verdict_i8('pinned_width_is_the_pre_20260919_default',            &
                      iw, iv, nfail)
      write(*,'(A,ES12.5)') '     (default - pinned)/pinned = ',             &
           (dr_base_default - w)/w

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'base_cell_width_provenance: ', nfail,          &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'base_cell_width_provenance: all assertions PASSED'

      contains

      subroutine verdict_i8(name, measured, reference, nfail)
      character(len=*), intent(in)    :: name
      integer*8,        intent(in)    :: measured, reference
      integer,          intent(inout) :: nfail
      if (measured .eq. reference) then
         write(*,'(A,A,A,Z16.16,A,Z16.16,A)') 'PASS ', name,                 &
              ' measured=', measured, ' reference=', reference, ' tol=0'
      else
         write(*,'(A,A,A,Z16.16,A,Z16.16,A)') 'FAIL ', name,                 &
              ' measured=', measured, ' reference=', reference, ' tol=0'
         nfail = nfail + 1
      endif
      end subroutine verdict_i8

      subroutine verdict_i(name, measured, reference, nfail)
      character(len=*), intent(in)    :: name
      integer,          intent(in)    :: measured, reference
      integer,          intent(inout) :: nfail
      if (measured .eq. reference) then
         write(*,'(A,A,A,I0,A,I0,A)') 'PASS ', name,                         &
              ' measured=', measured, ' reference=', reference, ' tol=0'
      else
         write(*,'(A,A,A,I0,A,I0,A)') 'FAIL ', name,                         &
              ' measured=', measured, ' reference=', reference, ' tol=0'
         nfail = nfail + 1
      endif
      end subroutine verdict_i

      end program base_cell_width_provenance
