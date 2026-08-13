      module output_write
      ! Write the output to the standard output files

      use global_parameters
      use species_table, only: n_mion, mion_name
      use ionization_equilibrium, only: nmol_eq,   &  ! molecular columns
                                       NH2_col_lw, f_shield_lw, k_lw_diss
      use lyman_werner_photodissociation, only: e_lw_fragment_erg

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
      real*8  :: nh_lw     ! H nuclei density [cm^-3], LW diagnostic
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
      ! Column 2 is rho*n0: the MASS density in units of m_H per cm^3 (metals
      ! included under the eos_metals policy), not a number density. Multiply by
      ! m_H to get g/cm^3.
      write(2,'(A)') '# EXHALE schema 2'
      write(2,'(A)') '# columns r[Rp] rho[mH/cm3] v[cm/s] p[cgs] T[K] '//   &
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
      write(3,'(A)') '# EXHALE schema 2'
      write(3,'(A)', advance='no') '# columns r[Rp] HI HII HeI HeII '// &
                                   'HeIII HeITR'
      do i = 1,n_mion
         write(3,'(A)', advance='no') ' '//trim(mion_name(i))
      enddo
      ! molecular columns (present only when thereis_mol)
      if (thereis_mol) write(3,'(A)', advance='no') ' H2 H2p H3p HeHp'
      write(3,'(A)') ''

      do j = 1-Ng,N+Ng

         if (thereis_mol) then
            write(3,*) r(j), nhi(j)*n0, nhii(j)*n0, nhei(j)*n0,        &
                     nheii(j)*n0, nheiii(j)*n0, nheiTR(j)*n0,          &
                     (nm(j,i)*n0, i = 1,n_mion),                       &
                     (nmol_eq(j,i), i = 1,4)  ! H2 H2+ H3+ HeH+ (already cm^-3)
         else
         write(3,*) r(j),       & ! Rad. dist.
                  nhi(j)*n0,    & ! HI
                  nhii(j)*n0,   & ! HII
                  nhei(j)*n0,   & ! HeI
                  nheii(j)*n0,  & ! HeII
                  nheiii(j)*n0, & ! HeIII
                  nheiTR(j)*n0, & ! HeITR
                  (nm(j,i)*n0, i = 1,n_mion)  ! metal ions (canonical order)
         endif
      enddo
      close(3)

      !---- Lyman-Werner photodissociation diagnostic ----!
      ! Written only for a molecular run that carries a Lyman-Werner band
      ! flux, and only for the equilibrium state (the advection
      ! post-process is atomic and does not re-solve the molecules).
      ! It is the record of how deep the band penetrates: f_shield -> 1 in
      ! the thin wind above the H2 -> H front and collapses in the
      ! self-shielded molecular base.
      if (thereis_mol .and. F_LW_star .gt. 0.0d0 .and. flag .eq. 'eq') then
         open(unit = 4, file = './output/Lyman_Werner.txt')
         write(4,'(A)') '# EXHALE schema 2'
         write(4,'(A,ES12.5,A)') '# Lyman-Werner band flux at the planet: ', &
                        F_LW_star, ' erg cm^-2 s^-1 (912-1110 A)'
         write(4,'(A)') '# columns r[Rp] T[K] x_H2[2nH2/nH] nH2[cm^-3] '//  &
                        'NH2_star[cm^-2] f_shield k_LW[1/s] '//             &
                        'heat_LW[erg/cm3/s]'
         do j = 1-Ng,N+Ng
            nh_lw = (nhi(j) + nhii(j))*n0                                  &
                  + 2.0d0*(nmol_eq(j,1) + nmol_eq(j,2))                    &
                  + 3.0d0*nmol_eq(j,3) + nmol_eq(j,4)
            write(4,*) r(j), T(j)*T0,                                      &
                       2.0d0*nmol_eq(j,1)/max(nh_lw,1.0d-99),              &
                       nmol_eq(j,1), NH2_col_lw(j), f_shield_lw(j),        &
                       k_lw_diss(j),                                       &
                       k_lw_diss(j)*nmol_eq(j,1)*e_lw_fragment_erg
         enddo
         close(4)
      endif

      ! End of subroutine
      end subroutine write_output
      
      ! End of module
      end module output_write
