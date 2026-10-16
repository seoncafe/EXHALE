      subroutine write_stationary_residual_record(label)
      ! Full precision rows and signed area-weighted mass fluxes of the
      ! immediately preceding assembly. No source or boundary is recomputed.
      character(len=*), intent(in) :: label
      real*8 :: Frho(1-Ng:N+Ng)
      integer :: unit_r, j, k
      if (.not. boundary_trace_armed()) return
      call face_mass_flux_of_state(u(1,:), Frho)
      open(newunit=unit_r, file='output/residual_'//trim(label)//'.txt', &
           status='replace', action='write')
      write(unit_r,'(A)') '# code units; physical cells 1..N'
      write(unit_r,'(A)') '# j r R_mass R_momentum R_energy'//          &
           ' scale_mass scale_momentum scale_energy heat cool'
      do j=1,N
         write(unit_r,'(I6,9ES26.16E3)') j,r(j),Rres(:,j),             &
              (residual_row_scale(k,j,u),k=1,3),heat(j),cool(j)
      enddo
      close(unit_r)
      open(newunit=unit_r, file='output/face_mass_'//trim(label)//'.txt', &
           status='replace', action='write')
      write(unit_r,'(A)') '# f r_face raw_mass_flux signed_r2_mass_flux'
      write(unit_r,'(A)') '# all physical faces 0..N; code units'
      do j=0,N
         write(unit_r,'(I6,3ES26.16E3)') j,r_edg(j),Frho(j),            &
              r_edg(j)**2*Frho(j)
      enddo
      close(unit_r)
      end subroutine write_stationary_residual_record
