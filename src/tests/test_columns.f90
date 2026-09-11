      module test_columns
      ! SYNTHETIC COLUMNS WHOSE SPECIES CARRY THE MASS DENSITY THEY ARE
      ! MEASURED AGAINST, and the measurement of that closure.  Compiled
      ! only by the test drivers that use it (it is not in the production
      ! SRC of the Makefile).
      !
      ! THE INVARIANT: sum_i m_i n_i = rho.  An operator that is handed a
      ! composition and a density is handed one gas, and what its species
      ! weigh under the code's own mass policy (calc_rho) has to be what
      ! the hydrodynamics evolves.  A column built by scaling helium or the
      ! metals by a radial factor while leaving hydrogen alone breaks that:
      ! its element ratios and its density are two different atmospheres,
      ! and every closure measured on it afterwards is the fixture's
      ! departure and not the operator's.  The construction below takes
      ! every mass fraction it gives helium and the trace elements OUT OF
      ! THE HYDROGEN, so the invariant holds to arithmetic round-off:
      ! MEASURED 1.0e-15 at the atomic entry state and 1.0e-14 on these
      ! synthetic columns (P1 report, item P1, 2026-09-11), which is the
      ! 1.0d-14 bound the drivers hold their rows to.
      !
      ! ONE CONSTRUCTION, TWO SUITES.  src/tests/element_operator measures
      ! the element transport operator on a column with a helium gradient;
      ! src/tests/steady_species_rows measures the drift the fixed-wind
      ! relaxation reports on a column whose trace metals swing about their
      ! own reservoir.  Both are the same column with different arguments,
      ! so the invariant is stated and held in one place.

      use global_parameters
      use species_table
      use composition, only: mass_per_H_nucleus_without_He
      use utils, only: calc_rho
      implicit none

      contains

      ! ================================================================= !

      subroutine column_carrying_its_own_density(f_c, q_h2, he_taper,     &
                                                 metal_swing)
      ! A column written as mass fractions that sum to one: helium carries
      ! X_He, each trace element its own Z_e, hydrogen the remainder, and
      ! the species vector holds a mass fraction divided by the species
      ! mass.  Nothing is added on top of a fixed hydrogen fraction, so the
      ! species reconstruct the density they are handed exactly.
      !
      !   q_h2        the fraction of the hydrogen NUCLEI bound into H2; at
      !               zero the column is atomic.
      !   he_taper    true makes the helium mass fraction fall by a factor
      !               two from the base to the outermost cell, so that the
      !               element operator has a gradient to work on; false
      !               puts helium at its reservoir value everywhere, so it
      !               has almost no distance to travel.
      !   metal_swing the amplitude of the swing every trace element makes
      !               about its own reservoir mixing ratio, in its neutral
      !               stage.  Cell 1 is the Dirichlet reservoir and is left
      !               at the reservoir value, which is what the relaxation
      !               measures each element against.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(out) :: f_c
      real*8,                                 intent(in)  :: q_h2
      logical,                                intent(in)  :: he_taper
      real*8,                                 intent(in)  :: metal_swing
      real*8  :: mpH, xhe0, xhe, zsum, ze(n_melem), fH, swing
      integer :: j, ie
      mpH  = mass_per_H_nucleus_without_He() + m_He_over_m_H*HeH
      xhe0 = m_He_over_m_H*HeH/mpH
      f_c  = 0.0d0
      do j = 1-Ng, N+Ng
         if (he_taper) then
            xhe = xhe0/(1.0d0 + (r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30))
         else
            xhe = xhe0
         endif
         swing = 1.0d0
         if (j .ne. 1) swing = 1.0d0                                      &
                             + metal_swing*sin(6.0d0*(r(j) - 1.0d0))
         zsum = 0.0d0
         do ie = 1, n_melem
            ze(ie) = melem_A(ie)*melem_ab(ie)*swing/mpH
            zsum   = zsum + ze(ie)
         enddo
         ! hydrogen NUCLEI per unit mass, then split between H I and H2
         fH = (1.0d0 - xhe - zsum)/bsp_mass(isp_HI)
         f_c(j,isp_HI)  = (1.0d0 - q_h2)*fH
         f_c(j,isp_H2)  = 0.5d0*q_h2*fH
         f_c(j,isp_HeI) = xhe/m_He_over_m_H
         do ie = 1, n_melem
            f_c(j,mion_fsp(melem_i0(ie))) = ze(ie)/melem_A(ie)
         enddo
      enddo
      end subroutine column_carrying_its_own_density

      ! ================================================================= !

      function column_mass_closure(rho_c, f_c) result(closure)
      ! |sum_i m_i n_i - rho|/rho over the physical cells, with the code's
      ! own mass policy (calc_rho): whether a state's species still carry
      ! the density the hydrodynamics evolves.  The molecular columns are
      ! passed whatever the configuration: with no molecules they are zero
      ! and add nothing to the sum.
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho_c
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_c
      real*8  :: closure
      real*8, dimension(1-Ng:N+Ng)        :: mrho
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm_c
      real*8, dimension(1-Ng:N+Ng,4)      :: nmol_c
      integer :: j, k
      do k = 1, n_mion
         nm_c(:,k) = f_c(:,mion_fsp(k))*rho_c*n0
      enddo
      nmol_c(:,1) = f_c(:,isp_H2) *rho_c*n0
      nmol_c(:,2) = f_c(:,isp_H2p)*rho_c*n0
      nmol_c(:,3) = f_c(:,isp_H3p)*rho_c*n0
      nmol_c(:,4) = f_c(:,isp_HeHp)*rho_c*n0
      call calc_rho(f_c(:,isp_HI)*rho_c*n0, f_c(:,isp_HII)*rho_c*n0,      &
                    f_c(:,isp_HeI)*rho_c*n0, f_c(:,isp_HeII)*rho_c*n0,    &
                    f_c(:,isp_HeIII)*rho_c*n0, mrho, nm = nm_c,           &
                    nmol = nmol_c)
      closure = 0.0d0
      do j = 1, N
         closure = max(closure, abs(mrho(j) - rho_c(j)*n0)                &
                                /max(rho_c(j)*n0, 1.0d-300))
      enddo
      end function column_mass_closure

      end module test_columns
