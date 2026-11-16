      subroutine gas_admissible_perturbation(y_base, direction, amplitude, &
                                             y_trial, label)
      ! Shorten a diagnostic direction until every physical cell has
      ! positive density and internal energy. Never clip the state itself.
      real*8, intent(in) :: y_base(:), direction(:)
      real*8, intent(inout) :: amplitude
      real*8, intent(out) :: y_trial(:)
      character(len=*), intent(in) :: label
      real*8 :: initial_amplitude
      integer :: reduction, j_gas, offset
      logical :: gas_ok
      initial_amplitude = amplitude
      do reduction=0,60
         y_trial = y_base + amplitude*direction
         gas_ok = .true.
         do j_gas=1,N
            offset = nvar_jac*(j_gas-1)
            if (.not. admissible_conserved_cell(                         &
                         y_trial(offset+1:offset+3))) gas_ok=.false.
         enddo
         if (gas_ok .and. any(y_trial .ne. y_base)) then
            write(*,'(A,A,A,ES12.4,A,I0)') ' (diagnostic perturbation) ', &
                 trim(label), ' amplitude factor ',                     &
                 amplitude/initial_amplitude, ' halvings ', reduction
            return
         endif
         amplitude = 0.5d0*amplitude
      enddo
      write(*,'(A,A)') ' ERROR: no nonzero gas-admissible probe: ',trim(label)
      stop 1
      end subroutine gas_admissible_perturbation
