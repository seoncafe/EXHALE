      module output_write
      ! Write the output to the standard output files

      use global_parameters
      use species_table, only: n_mion, mion_name

      contains

      subroutine write_output(rho,v,p,T,heat,cool,eta,                &
                              nhi,nhii,nhei,nheii,nheiii,nheiTR,      &
                              nm,flag)
      ! Metal ion densities are passed as the 2D array nm(:, 1:n_mion),
      ! one column per metal ion stage in the canonical species_table
      ! order (CI, CII, CIII, OI, ..., MgIII). This keeps the argument
      ! list fixed as metals are added.

      character(len = 2) :: flag
      integer :: j,i
      real*8, dimension(1-Ng:N+Ng), intent(in) :: rho,v,p,T
      real*8, dimension(1-Ng:N+Ng), intent(in) :: heat,cool
      real*8, dimension(1-Ng:N+Ng), intent(in) :: eta
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi,nhii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei,nheii,nheiii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nheiTR
      real*8, dimension(1-Ng:N+Ng,n_mion), intent(in) :: nm

      
      !---- Write thermodynamic profiles ----! 
          
      if (flag.eq.'eq') then
      	open(unit = 2, file = './output/Hydro_ioniz.txt')
	   else	! Change output file after postprocessing
		   open(unit = 2, file = './output/Hydro_ioniz_adv.txt')
	   endif

      ! Schema header ('#' comment lines; readers that predate the header
      ! can skip them, numeric content is unchanged)
      write(2,'(A)') '# ATES-metal schema 2'
      write(2,'(A)') '# columns r[Rp] n[cm-3] v[cm/s] p[cgs] T[K] '//   &
                     'heat[erg/cm3/s] cool[erg/cm3/s]'

         do j = 1-Ng,N+Ng
            write(2,*) r(j),        &     ! Rad. dist.
                     rho(j)*n0,     &     ! Density
                     v(j)*v0,       &     ! Velocity
                     p(j)*p0,       &     ! Pressure
                     T(j)*T0,       &     ! Temperature
                     heat(j)*q0,    &     ! Rad. heat.
                     cool(j)*q0           ! Rad. cool.
         enddo
      close(2)
      
      !---- Write ionization profiles ----!  
      if (flag.eq.'eq') then
      	open(unit = 3, file = './output/Ion_species.txt')
	   else	! Change output file after postprocessing
		   open(unit = 3, file = './output/Ion_species_adv.txt')
	   endif

      ! Schema header: species labels in column order, generated from the
      ! species table so they stay correct when species are added.
      write(3,'(A)') '# ATES-metal schema 2'
      write(3,'(A)', advance='no') '# columns r[Rp] HI HII HeI HeII '// &
                                   'HeIII HeITR'
      do i = 1,n_mion
         write(3,'(A)', advance='no') ' '//trim(mion_name(i))
      enddo
      write(3,'(A)') ''

      do j = 1-Ng,N+Ng
         
         write(3,*) r(j),       & ! Rad. dist.
                  nhi(j)*n0,    & ! HI
                  nhii(j)*n0,   & ! HII
                  nhei(j)*n0,   & ! HeI
                  nheii(j)*n0,  & ! HeII
                  nheiii(j)*n0, & ! HeIII
                  nheiTR(j)*n0, & ! HeITR
                  (nm(j,i)*n0, i = 1,n_mion)  ! metal ions (canonical order)
      enddo
      close(3)
                  
      ! End of subroutine
      end subroutine write_output
      
      ! End of module
      end module output_write
