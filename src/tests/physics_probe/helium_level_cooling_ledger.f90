      program helium_level_cooling_ledger
      ! Which helium level each cooling channel of eval_cool is charged to,
      ! and at which recombination coefficient the He II channel is charged.
      !
      ! Origin: the two helium findings of
      ! docs/audit_20260905/source_audit/group_B_report.md.  Production
      ! routine exercised: eval_cool (src/modules/radiation/util_ion_eq.f90)
      ! and the coefficient routines of Cool_coeff.f90 it assembles.
      !
      ! REFERENCES AND WHAT IS ASSERTED.
      !
      !  1. Collisional ionization and collisional excitation of neutral
      !     helium.  The state vector carries ONE neutral-helium column,
      !     n(He I), which CONTAINS the 2^3S metastable
      !     (composition.f90: he_ground_singlet_density is where the
      !     difference is taken).  The two levels are 19.8 eV apart and have
      !     nothing in common: the 24.6 eV ionization potential and the Cen
      !     (1992) collisional-excitation coefficient belong to the 1^1S
      !     ground singlet, the 4.8 eV potential (Black 1981), the 10830 A
      !     excitation and the 2^3S -> 2^1S / 2^1P conversions belong to the
      !     metastable.  So each of the two densities must appear exactly
      !     once, multiplied by its own level's coefficient:
      !
      !         coio_He = 24.6 eV * a_ion_HeI * n(1^1S)
      !                 +  4.8 eV * ci_HeI23S * n(2^3S)
      !         coex_He = lambda_coex_HeI      * n(1^1S)
      !                 + (10830 + 0.80 q31a + 1.40 q31b) * n(2^3S)
      !
      !     These are identities of the assembly, not fits, so the tolerance
      !     is round-off: 1e-12 relative.  The configuration removes H I and
      !     He II so that the two breakdown columns contain the helium terms
      !     alone and the comparison carries no cancellation.
      !
      !     The split is asserted at three temperatures because the size of
      !     the defect it targets depends on temperature: the ground-state
      !     coefficients are exponentially suppressed below about 2e4 K
      !     (24.6 eV against kT), so charging them to the metastable as well
      !     is invisible at 1e4 K and is of order n(2^3S)/n(1^1S) at 5e4 K.
      !
      !  2. He II recombination cooling.  The electron gas loses kT of
      !     thermal energy per recombination the balance actually performs,
      !     so the cooling coefficient must be kT times the coefficient with
      !     which the balance removes He+.  With the metastable tracked that
      !     coefficient is the sum of the two capture channels the balance
      !     writes into its rows, rcheiiB + rcheiTR
      !     (ionization_equilibrium.f90), i.e. rec_HeII_11S + rec_HeII_23S;
      !     without it, the case-B coefficient rec_HeII_B.  Asserted against
      !     the same routines the balance calls, to round-off, and together
      !     with the identity that the recombination column of the breakdown
      !     is n_e times the three (coefficient x ion density) terms.
      !
      ! At 35d9dd5 the He I coefficients were charged to the summed neutral
      ! helium and the 2^3S terms added on top, and the He II recombination
      ! cooling was kT times the case-B coefficient in every configuration.

      use global_parameters
      use species_table,  only: n_mion
      use Cooling_Coefficients, only: ion_coeff_HeI, ci_HeI23S,            &
                                      coex_rate_HeI, coex_rate_HeI23S_10830,&
                                      coex_HeI_23S_21S, coex_HeI_23S_21P,  &
                                      rec_cool_HII, rec_cool_HeII,         &
                                      rec_cool_HeIII, rec_HeII_B
      use utils_ion_eq,   only: eval_cool, HeITR_coeffs
      use assertion_report

      implicit none

      ! Level populations of the slab [cm^-3].  n(2^3S)/n(He I) = 1e-3, the
      ! ratio the finding is stated at (the largest measured on wasp_full is
      ! 2.2e-4).
      real*8, parameter :: n_singlet = 1.0d5
      real*8, parameter :: n_meta    = 1.0d2
      real*8, parameter :: n_HII     = 1.0d6
      real*8, parameter :: n_HeII    = 1.0d5
      real*8, parameter :: n_HeIII   = 1.0d3

      real*8, dimension(:), allocatable :: T_K, nhi,nhii,nhei,nheii,nheiii
      real*8, dimension(:), allocatable :: nheiTR
      real*8, dimension(:), allocatable :: rchiiB,rcheiiB,rcheiiiB
      real*8, dimension(:), allocatable :: aHI,aHeI,aHeII,aHeTR
      real*8, dimension(:), allocatable :: cool
      real*8, dimension(:,:), allocatable :: nm, rec_m, aion_m, cool_chan
      ! Reference coefficients of the two helium levels.
      real*8, dimension(:), allocatable :: a_ion_HeI_ref, ci_HeTR_ref
      real*8, dimension(:), allocatable :: coex_HeI_ref, coex_10830_ref
      real*8, dimension(:), allocatable :: q31a_ref, q31b_ref
      real*8, dimension(:), allocatable :: lam_HII,lam_HeII,lam_HeIII
      real*8, dimension(:), allocatable :: rcheiTR_bal,rcheii_bal,alphaB
      real*8, dimension(:), allocatable :: q13_d,q31a_d,q31b_d,Q31_d
      real*8 :: A31_d, ne_slab, reference
      character(len=56) :: label
      integer :: j

      ! --- run-wide configuration ----------------------------------------
      ! Helium with the metastable tracked, no metals, no molecules: the
      ! cooling reduces to the H/He channels this test reads.
      thereis_He      = .true.
      thereis_HeITR   = .true.
      thereis_metals  = .false.
      thereis_mol     = .false.
      thereis_oxychem = .false.
      use_2lev_cool   = .false.
      base_ir_field   = .false.
      mol_ir_bands    = .false.

      N  = 3
      R0 = 1.0d10
      call allocate_grid_arrays

      allocate(T_K(1-Ng:N+Ng), nhi(1-Ng:N+Ng), nhii(1-Ng:N+Ng),           &
               nhei(1-Ng:N+Ng), nheii(1-Ng:N+Ng), nheiii(1-Ng:N+Ng),      &
               nheiTR(1-Ng:N+Ng))
      allocate(rchiiB(1-Ng:N+Ng), rcheiiB(1-Ng:N+Ng), rcheiiiB(1-Ng:N+Ng),&
               aHI(1-Ng:N+Ng), aHeI(1-Ng:N+Ng), aHeII(1-Ng:N+Ng),         &
               aHeTR(1-Ng:N+Ng), cool(1-Ng:N+Ng))
      allocate(nm(1-Ng:N+Ng,n_mion), rec_m(1-Ng:N+Ng,n_mion),             &
               aion_m(1-Ng:N+Ng,n_mion), cool_chan(1-Ng:N+Ng,10+n_mion))
      allocate(a_ion_HeI_ref(1-Ng:N+Ng), ci_HeTR_ref(1-Ng:N+Ng),          &
               coex_HeI_ref(1-Ng:N+Ng), coex_10830_ref(1-Ng:N+Ng),        &
               q31a_ref(1-Ng:N+Ng), q31b_ref(1-Ng:N+Ng))
      allocate(lam_HII(1-Ng:N+Ng), lam_HeII(1-Ng:N+Ng),                   &
               lam_HeIII(1-Ng:N+Ng))
      allocate(rcheiTR_bal(1-Ng:N+Ng), rcheii_bal(1-Ng:N+Ng),             &
               alphaB(1-Ng:N+Ng), q13_d(1-Ng:N+Ng), q31a_d(1-Ng:N+Ng),    &
               q31b_d(1-Ng:N+Ng), Q31_d(1-Ng:N+Ng))

      ! Three temperatures across the onset of the He I ground-state
      ! channels; the ghost cells carry the coolest of them.
      T_K    = 1.0d4
      T_K(1) = 1.0d4
      T_K(2) = 2.0d4
      T_K(3) = 5.0d4

      nm     = 0.0d0
      nhei   = n_singlet + n_meta
      nheiTR = n_meta
      nhii   = n_HII

      ! --- 1. the two neutral-helium levels -------------------------------
      ! No H I and no He II, so the collisional-ionization column and the
      ! He I collisional-excitation column of the breakdown hold the two
      ! helium-level terms and nothing else.  n_e is the proton density.
      nhi    = 0.0d0
      nheii  = 0.0d0
      nheiii = 0.0d0
      ne_slab = n_HII

      call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,                  &
                     rchiiB,rcheiiB,rcheiiiB, rec_m,                      &
                     aHI,aHeI,aHeII, aion_m, cool, cool_chan,             &
                     nheiTR=nheiTR, a_ion_HeITR=aHeTR)

      call ion_coeff_HeI(T_K, a_ion_HeI_ref)
      call ci_HeI23S(T_K, ci_HeTR_ref)
      call coex_rate_HeI(T_K, coex_HeI_ref)
      call coex_rate_HeI23S_10830(T_K, coex_10830_ref)
      call coex_HeI_23S_21S(T_K, q31a_ref)
      call coex_HeI_23S_21P(T_K, q31b_ref)

      do j = 1, N
         ! Collisional ionization: the He I ground-state potential out of the
         ! singlet, the 2^3S potential out of the metastable, each density
         ! once. Both energies are the global thresholds (e_th_HeI_erg,
         ! e_th_HeTR_erg) since batch 2c unified the duplicated literals;
         ! the single-precision 3.940e-11 erg eval_cool once carried is gone.
         reference = ne_slab*(e_th_HeI_erg*a_ion_HeI_ref(j)*n_singlet      &
                     + e_th_HeTR_erg*ci_HeTR_ref(j)*n_meta)
         write(label,'(a,i0,a)') 'helium_collisional_ionization_by_level_T',&
              int(T_K(j)), 'K'
         call check_relative(label, cool_chan(j,2), reference, 1.0d-12)

         ! Collisional excitation: the Cen (1992) He I coefficient out of the
         ! singlet, the 10830 A and singlet-conversion channels out of the
         ! metastable.
         reference = ne_slab*(coex_HeI_ref(j)*n_singlet                    &
                     + (coex_10830_ref(j)                                  &
                        + (0.80d0*q31a_ref(j) + 1.40d0*q31b_ref(j))        &
                          /erg2eV)*n_meta)
         write(label,'(a,i0,a)') 'helium_collisional_excitation_by_level_T',&
              int(T_K(j)), 'K'
         call check_relative(label, cool_chan(j,4), reference, 1.0d-12)
      enddo

      ! --- 2. He II recombination cooling ---------------------------------
      ! kT per recombination, at the coefficient the balance removes He+
      ! with.  First the coefficient itself, then the recombination column
      ! of the breakdown built from it.
      call HeITR_coeffs(T_K,rcheiTR_bal,rcheii_bal,A31_d,q13_d,q31a_d,     &
                        q31b_d,Q31_d)
      call rec_cool_HeII(T_K, lam_HeII)
      do j = 1, N
         write(label,'(a,i0,a)')                                           &
              'he_ii_recombination_cooling_per_event_metastable_T',        &
              int(T_K(j)), 'K'
         call check_relative(label, lam_HeII(j)/(kb_erg*T_K(j)),           &
                             rcheii_bal(j) + rcheiTR_bal(j), 1.0d-12)
      enddo

      nheii  = n_HeII
      nheiii = n_HeIII
      ne_slab = n_HII + n_HeII + 2.0d0*n_HeIII
      call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,                  &
                     rchiiB,rcheiiB,rcheiiiB, rec_m,                      &
                     aHI,aHeI,aHeII, aion_m, cool, cool_chan,             &
                     nheiTR=nheiTR, a_ion_HeITR=aHeTR)
      call rec_cool_HII(T_K, lam_HII)
      call rec_cool_HeIII(T_K, lam_HeIII)
      do j = 1, N
         reference = ne_slab*(lam_HII(j)*n_HII + lam_HeII(j)*n_HeII        &
                              + lam_HeIII(j)*n_HeIII)
         write(label,'(a,i0,a)') 'recombination_cooling_column_T',         &
              int(T_K(j)), 'K'
         call check_relative(label, cool_chan(j,1), reference, 1.0d-12)
      enddo

      ! Without the metastable the balance keeps the case-B coefficient and
      ! the cooling must follow it back.
      thereis_HeITR = .false.
      call rec_HeII_B(T_K, alphaB)
      call rec_cool_HeII(T_K, lam_HeII)
      do j = 1, N
         write(label,'(a,i0,a)')                                           &
              'he_ii_recombination_cooling_per_event_case_b_T',            &
              int(T_K(j)), 'K'
         call check_relative(label, lam_HeII(j)/(kb_erg*T_K(j)),           &
                             alphaB(j), 1.0d-12)
      enddo
      thereis_HeITR = .true.

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'helium_level_cooling_ledger: ',              &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'helium_level_cooling_ledger: all assertions passed'

      end program helium_level_cooling_ledger
