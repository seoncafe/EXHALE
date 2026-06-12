      module composition
      ! Single point that turns (rho, f_sp) into all species number
      ! densities plus the free-electron (ne) and total-particle (n_tot)
      ! densities, and converts between pressure and temperature.
      !
      ! This collapses the three near-identical extraction blocks that were
      ! repeated in ATES_main (and the inline copies in energy_semi_implicit
      ! / post_process), so that the electron-density and total-density
      ! POLICY lives in exactly one place.
      !
      ! The eos_include_metals policy (global_parameters, default .true.;
      ! runtime key 'eos_metals 0|1' in metals.inp) is applied here by
      ! passing nm through to calc_ne / calc_ntot, which add the metal
      ! electrons and metal nuclei. With eos_metals 0 (or metals off) the
      ! legacy H/He-only behavior is reproduced exactly.

      use global_parameters
      use species_table, only: n_mion, mion_fsp,                       &
                               isp_HI, isp_HII, isp_HeI, isp_HeII,      &
                               isp_HeIII, isp_HeTR
      use utils, only: calc_ne, calc_ntot

      implicit none
      private
      public :: get_species_densities, comp_T_from_p, comp_p_from_T

      contains

      ! ------------------------------------------------------!

      subroutine get_species_densities(rho, f_sp, nhi, nhii, nhei,     &
                                       nheii, nheiii, nheiTR, nm,       &
                                       ne, n_tot)
      ! (rho, f_sp) -> all number densities + ne + n_tot, reproducing the
      ! legacy ATES_main extraction blocks exactly. The He arrays are
      ! intent(inout): like the original module-scope locals they keep
      ! their previous values when helium is absent (they are then never
      ! read by calc_ne/calc_ntot, which guard on thereis_He).

      real*8, dimension(1-Ng:N+Ng),          intent(in)    :: rho
      real*8, dimension(1-Ng:N+Ng,n_species),intent(in)    :: f_sp
      real*8, dimension(1-Ng:N+Ng),          intent(out)   :: nhi, nhii
      real*8, dimension(1-Ng:N+Ng),          intent(inout) :: nhei, nheii
      real*8, dimension(1-Ng:N+Ng),          intent(inout) :: nheiii, nheiTR
      real*8, dimension(1-Ng:N+Ng,n_mion),   intent(out)   :: nm
      real*8, dimension(1-Ng:N+Ng),          intent(out)   :: ne, n_tot
      integer :: im

      nhi  = rho*f_sp(:,isp_HI)
      nhii = rho*f_sp(:,isp_HII)
      if (thereis_He) then
         nhei   = rho*f_sp(:,isp_HeI)
         nheii  = rho*f_sp(:,isp_HeII)
         nheiii = rho*f_sp(:,isp_HeIII)
         if (thereis_HeITR) nheiTR = rho*f_sp(:,isp_HeTR)
      endif
      do im = 1, n_mion
         nm(:,im) = rho*f_sp(:,mion_fsp(im))
      enddo

      call calc_ne(nhii, nheii, nheiii, ne, nm)
      call calc_ntot(nhi, nhii, nhei, nheii, nheiii, nheiTR, n_tot, nm)

      end subroutine get_species_densities

      ! ------------------------------------------------------!

      subroutine comp_T_from_p(p, n_tot, ne, T)
      ! Adimensional ideal-gas inverse: T = p / (n_tot + ne).
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: p, n_tot, ne
      real*8, dimension(1-Ng:N+Ng), intent(out) :: T
      T = p/(n_tot + ne)
      end subroutine comp_T_from_p

      ! ------------------------------------------------------!

      subroutine comp_p_from_T(T, n_tot, ne, p)
      ! Adimensional ideal-gas law: p = (n_tot + ne) * T.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: T, n_tot, ne
      real*8, dimension(1-Ng:N+Ng), intent(out) :: p
      p = (n_tot + ne)*T
      end subroutine comp_p_from_T

      ! End of module
      end module composition
