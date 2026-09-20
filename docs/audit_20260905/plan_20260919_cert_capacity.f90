program certification_capacity_probe
   use certification, only: cert_report, cert_max_entries, &
                            ionization_stage_sum_entry
   implicit none
   type(cert_report) :: rep

   ! Exercise the public entry writer on a report with no free slot.
   ! The linked production routine must not silently rewrite the last row.
   rep%n = cert_max_entries
   rep%e(rep%n)%name = 'existing equation sentinel'
   rep%e(rep%n)%row_max = 42.0d0
   call ionization_stage_sum_entry(rep, 1, 0.0d0, 7, .true.)
   write(*,'(A,I0)') 'entry_count=', rep%n
   write(*,'(A,A)') 'last_name=', trim(rep%e(rep%n)%name)
   write(*,'(A,ES24.16)') 'last_measure=', rep%e(rep%n)%row_max
   write(*,'(A,L1)') 'silent_overwrite=', &
        rep%n == cert_max_entries .and. rep%e(rep%n)%row_max /= 42.0d0
end program certification_capacity_probe
