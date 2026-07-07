      module opacity_input
      ! Optional opacity.inp reader. Format is KEY = VALUE, one entry
      ! per line; blank lines and lines starting with '#' are ignored;
      ! every key is optional. Missing file or missing keys leave the
      ! defaults from parameters.f90 untouched (i.e., the legacy ATES
      ! analytic cross sections).
      !
      ! Recognized keys:
      !   OPACITY_MODEL     = A | C | P | T
      !   OPA_CONST_HI      = <real>      (factor; default 1.0)
      !   OPA_CONST_HEI     = <real>
      !   OPA_CONST_HEII    = <real>
      !   OPA_CONST_HEITR   = <real>
      !   OPA_PB_FACTOR     = <real>      (Robinson&Catling a)
      !   OPA_PB_EXPONENT   = <real>      (n)
      !   OPA_PB_PIVOT      = <real>      (dyne/cm^2; default 1e5)
      !   OPA_FILE_HI       = <path>      (.opa table for HI)
      !   OPA_FILE_HEI      = <path>
      !   OPA_FILE_HEII     = <path>
      !   OPA_FILE_HEITR    = <path>

      use global_parameters

      implicit none

      private
      public :: read_opacity_input

      character(len = *), parameter :: opa_inp_file = 'opacity.inp'

      contains

      !----------------------------------------------------------!

      subroutine read_opacity_input
      integer :: io, eq_pos, n_set
      character(len = 256) :: line
      character(len = :), allocatable :: trimmed, key, val
      logical :: file_exists

      inquire(file = opa_inp_file, exist = file_exists)
      if (.not. file_exists) then
         write(*,*) '(opacity_input) No opacity.inp found; using defaults (model=', &
                    opacity_model, ').'
         return
      endif

      write(*,*) '(opacity_input) Reading opacity.inp..'
      open(unit = 38, file = opa_inp_file, status = 'old', action = 'read')

      n_set = 0
      do
         read(38, '(A)', iostat = io) line
         if (io /= 0) exit

         trimmed = trim(adjustl(line))
         if (len(trimmed) == 0) cycle
         if (trimmed(1:1) == '#') cycle

         eq_pos = index(trimmed, '=')
         if (eq_pos <= 1) then
            write(*,*) '  WARN: skipping malformed line: ', trim(trimmed)
            cycle
         endif

         key = trim(adjustl(trimmed(1:eq_pos-1)))
         val = trim(adjustl(trimmed(eq_pos+1:)))
         call to_upper(key)

         select case (key)
            case ('OPACITY_MODEL')
               opacity_model = val(1:1)
               call to_upper_char(opacity_model)
            case ('OPA_CONST_HI');     read(val,*) opa_const_HI
            case ('OPA_CONST_HEI');    read(val,*) opa_const_HeI
            case ('OPA_CONST_HEII');   read(val,*) opa_const_HeII
            case ('OPA_CONST_HEITR');  read(val,*) opa_const_HeITR
            case ('OPA_PB_FACTOR');    read(val,*) opa_pb_factor
            case ('OPA_PB_EXPONENT');  read(val,*) opa_pb_exponent
            case ('OPA_PB_PIVOT');     read(val,*) opa_pb_pivot
            case ('OPA_FILE_HI');      opa_file_HI    = val
            case ('OPA_FILE_HEI');     opa_file_HeI   = val
            case ('OPA_FILE_HEII');    opa_file_HeII  = val
            case ('OPA_FILE_HEITR');   opa_file_HeITR = val
            case default
               write(*,*) '  WARN: unknown key: ', trim(key)
               cycle
         end select

         n_set = n_set + 1
      enddo

      close(38)

      if (.not. (opacity_model == 'A' .or. opacity_model == 'C' .or. &
                 opacity_model == 'P' .or. opacity_model == 'T')) then
         write(*,*) '(opacity_input) ERROR: OPACITY_MODEL must be one of A/C/P/T, got "', &
                    opacity_model, '".'
         stop 1
      endif

      write(*,'(a,i0,a,a)') ' (opacity_input) Done. ', n_set, &
                            ' key(s) set; model = ', opacity_model

      end subroutine read_opacity_input

      !----------------------------------------------------------!

      subroutine to_upper(s)
      character(len = :), allocatable, intent(inout) :: s
      integer :: i, ic
      do i = 1, len(s)
         ic = iachar(s(i:i))
         if (ic >= iachar('a') .and. ic <= iachar('z')) then
            s(i:i) = achar(ic - 32)
         endif
      enddo
      end subroutine to_upper

      subroutine to_upper_char(c)
      character(len = 1), intent(inout) :: c
      integer :: ic
      ic = iachar(c)
      if (ic >= iachar('a') .and. ic <= iachar('z')) then
         c = achar(ic - 32)
      endif
      end subroutine to_upper_char

      ! End of module
      end module opacity_input
