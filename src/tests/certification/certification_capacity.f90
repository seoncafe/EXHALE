      program certification_capacity
      ! ONE EQUATION PAST A FULL REPORT, through the public entry writer
      ! (PLAN_20260919_rev1, item P5a).
      !
      ! THE CONTRACT: add_entry refuses an entry that does not fit, names
      ! the equation, the capacity and the entry that stands last, and stops
      ! the run BEFORE writing any field. Nothing the report already holds
      ! moves. The refusal cannot be observed from inside the process it
      ! stops, so this driver is a program of its own and run.sh reads its
      ! exit status and its text.
      !
      ! The entry text returned silently instead, and the entry writer then
      ! wrote the new equation's measure, tolerance and verdict into the
      ! entry that already stood last, under that entry's name. Against that
      ! text this program runs to its end and prints the overwritten
      ! sentinel.
      use certification,     only: cert_report, cert_max_entries,          &
                                   ionization_stage_sum_entry
      use element_inventory, only: ien_H, ien_He
      implicit none
      type(cert_report) :: rep
      integer :: i
      rep%n = 0
      do i = 1, cert_max_entries - 1
         call ionization_stage_sum_entry(rep, ien_He, 0.0d0, i, .true.)
      enddo
      ! THE SENTINEL is the last equation that fits: another element, a
      ! measure no other entry carries and a cell of its own.
      call ionization_stage_sum_entry(rep, ien_H, 42.0d0,                  &
                                      cert_max_entries, .true.)
      write(*,'(A,I0)')     'filled entry_count=', rep%n
      write(*,'(A,A)')      'filled last_name=', trim(rep%e(rep%n)%name)
      write(*,'(A,ES13.6)') 'filled last_measure=', rep%e(rep%n)%row_max
      flush(6)
      ! ONE MORE EQUATION. Everything below is reached only by a text that
      ! writes into an entry it did not create.
      call ionization_stage_sum_entry(rep, ien_He, 7.0d0, 1, .true.)
      write(*,'(A,I0)')     'after entry_count=', rep%n
      write(*,'(A,A)')      'after last_name=',                           &
                            trim(rep%e(cert_max_entries)%name)
      write(*,'(A,ES13.6)') 'after last_measure=',                         &
                            rep%e(cert_max_entries)%row_max
      write(*,'(A,L1)')     'silent_overwrite=',                           &
                            rep%e(cert_max_entries)%row_max .ne. 42.0d0
      end program certification_capacity
