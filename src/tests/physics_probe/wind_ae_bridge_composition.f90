      program wind_ae_bridge_composition
      ! The composition an EXHALE input.inp hands the Wind-AE solver.
      !
      ! Production routine exercised: wae_read_exhale_input
      ! (src/modules/wind_ae/wae_exhale_input.f90), the front end that maps
      ! an EXHALE input.inp onto the Wind-AE parameter list, reached by
      ! `IC mode: windae` and by the standalone wind_ae_ic.x.
      !
      ! REFERENCE.  Wind-AE never stores a number ratio: it forms every
      ! number density from the mass fractions and the atomic masses it was
      ! given, n_j = rho HX(j)/atomic_mass(j) (wae_soe.f90 138, 166, 216;
      ! wae_glq_rates.f90 67).  The helium-to-hydrogen number ratio the
      ! solver therefore works at is
      !     (HX(2)/m_He) / (HX(1)/m_H),
      ! and it must equal the `He/H number ratio` of the input file exactly,
      ! for every value of it: a run asks for a helium abundance and the
      ! solver must relax at that abundance and no other.  The mass
      ! fractions must also close, HX(1) + HX(2) = 1, since hydrogen and
      ! helium are the only two species Wind-AE carries.
      !
      ! Both tolerances are round-off (1e-12 relative): these are exact
      ! identities of a unit conversion, not fits.
      use wae_params,       only: wae_par
      use wae_exhale_input, only: wae_read_exhale_input
      use assertion_report
      implicit none

      real*8, parameter :: heh_test(4) =                                  &
           (/ 0.0d0, 0.0793d0, 0.1d0, 1.0d0 /)
      character(len=*), parameter :: fname =                              &
           'wind_ae_bridge_composition_input.inp'
      real*8 :: heh_solver, mass_sum
      character(len=64) :: label
      character(len=16) :: tag
      integer :: i, u

      do i = 1, size(heh_test)
         call write_bridge_input(fname, heh_test(i))
         call wae_read_exhale_input(fname)
         open(newunit=u, file=fname, status='old')
         close(u, status='delete')

         heh_solver = (wae_par%HX(2)/wae_par%atomic_mass(2))              &
                    / (wae_par%HX(1)/wae_par%atomic_mass(1))
         mass_sum   = wae_par%HX(1) + wae_par%HX(2)

         write(tag,'(f7.4)') heh_test(i)
         tag = adjustl(tag)
         label = 'wind_ae_helium_hydrogen_number_ratio_at_'//trim(tag)
         call check_relative(trim(label), heh_solver, heh_test(i), 1.0d-12)
         label = 'wind_ae_mass_fraction_sum_at_'//trim(tag)
         call check_relative(trim(label), mass_sum, 1.0d0, 1.0d-12)
      enddo

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'wind_ae_bridge_composition: ',              &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'wind_ae_bridge_composition: all assertions passed'

      contains

      subroutine write_bridge_input(path, heh)
      ! A minimal input.inp carrying the keys the bridge reads, so that no
      ! field of the parameter list is built from a missing value.  Only the
      ! helium abundance is varied.
      character(len=*), intent(in) :: path
      real*8,           intent(in) :: heh
      integer :: iu
      open(newunit=iu, file=path, status='replace', action='write')
      write(iu,'(a)')          'Planet radius [R_J]:               1.0'
      write(iu,'(a)')          'Planet mass [M_J]:                 1.0'
      write(iu,'(a)')          'Equilibrium temperature [K]:       1000.0'
      write(iu,'(a)')          'Orbital distance [AU]:             0.05'
      write(iu,'(a)')          'Parent star mass [M_sun]:          1.0'
      write(iu,'(a)')          'X-ray luminosity [log10 erg/s]:    27.0'
      write(iu,'(a)')          'EUV luminosity [log10 erg/s]:      28.0'
      write(iu,'(a)')          'Stellar Teff [K]:                  5000.0'
      write(iu,'(a)')          'Stellar radius [R_sun]:            1.0'
      write(iu,'(a,es23.15)')  'He/H number ratio:                 ', heh
      close(iu)
      end subroutine write_bridge_input

      end program wind_ae_bridge_composition
