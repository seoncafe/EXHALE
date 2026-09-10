      module species_table
      ! Canonical metadata table for the trace-metal ion stages, plus
      ! composition metadata (mass/charge/nuclei) for the base H/He, the
      ! molecular and the oxygen-chemistry species (bsp_* below).
      !
      ! This module exists so that adding a metal becomes "add rows to a
      ! table" rather than "thread N new arguments through M subroutines."
      ! The ordering of the metal ions here MUST match the f_sp species
      ! columns 7..30 set up elsewhere (ionization_equilibrium, write_output,
      ! set_IC), i.e.
      !
      !     index  ion    f_sp col   element  stage
      !       1    C I        7         C       0
      !       2    C II       8         C       1
      !       3    C III      9         C       2
      !       4    O I       10         O       0
      !       5    O II      11         O       1
      !       6    O III     12         O       2
      !       7    N I       13         N       0
      !       8    N II      14         N       1
      !       9    N III     15         N       2
      !      10    Mg I      16        Mg       0
      !      11    Mg II     17        Mg       1
      !      12    Mg III    18        Mg       2
      !      13    Si I      19        Si       0
      !      14    Si II     20        Si       1
      !      15    Si III    21        Si       2
      !      16    Ca I      22        Ca       0
      !      17    Ca II     23        Ca       1
      !      18    Ca III    24        Ca       2
      !      19    Na I      25        Na       0
      !      20    Na II     26        Na       1
      !      21    K I       27         K       0
      !      22    K II      28         K       1
      !      23    S I       29         S       0
      !      24    S II      30         S       1
      !      25    Fe I      31        Fe       0
      !      26    Fe II     32        Fe       1
      !      27    Fe III    33        Fe       2
      !
      ! "stage" is also the net ionic charge, used for the bremsstrahlung
      ! Z^2 weight and for the free-electron sum inside the MINPACK system.
      ! Top stages are inert: not photo-ionized further and not line coolants
      ! in this set. For C/O/N/Mg/Si/Ca the top stage solved is the doubly
      ! ionized (++); for Na/K/S the top stage is the singly ionized (+),
      ! i.e. only two stages are carried (melem_top = 1).

      use global_parameters, only: e_th_MgI, e_th_MgII,                 &
                                   m_He_over_m_H, amu_over_m_H

      implicit none

      public

      ! ---- sizes ----
      integer, parameter :: n_mion  = 27   ! tracked metal ion stages
      integer, parameter :: n_melem = 10   ! C, O, N, Mg, Si, Ca, Na, K, S, Fe
      integer, parameter :: n_mphot = 17   ! photo-ionizable metal ions

      ! ---- base (H/He) species columns in f_sp ----
      ! Named constants for the fixed H/He/HeITR layout of f_sp(:,1:6), so the
      ! "which column is which species" question is answered here and the rest
      ! of the code stops using bare literals 1..6. Values are unchanged.
      integer, parameter :: isp_HI    = 1
      integer, parameter :: isp_HII   = 2
      integer, parameter :: isp_HeI   = 3
      integer, parameter :: isp_HeII  = 4
      integer, parameter :: isp_HeIII = 5
      integer, parameter :: isp_HeTR  = 6   ! He 2^3S metastable triplet
      ! ---- molecular species columns (molecular; zero unless thereis_mol) ----
      integer, parameter :: isp_H2    = 34
      integer, parameter :: isp_H2p   = 35
      integer, parameter :: isp_H3p   = 36
      integer, parameter :: isp_HeHp  = 37
      ! ---- oxygen-chemistry columns (zero unless thereis_oxychem) ----
      integer, parameter :: isp_OH    = 38
      integer, parameter :: isp_H2O   = 39
      integer, parameter :: isp_CO    = 40

      ! ---- element indices ----
      integer, parameter :: iel_C  = 1
      integer, parameter :: iel_O  = 2
      integer, parameter :: iel_N  = 3
      integer, parameter :: iel_Mg = 4
      integer, parameter :: iel_Si = 5
      integer, parameter :: iel_Ca = 6
      integer, parameter :: iel_Na = 7
      integer, parameter :: iel_K  = 8
      integer, parameter :: iel_S  = 9
      integer, parameter :: iel_Fe = 10

      ! ---- element masses ----
      ! melem_A_u: the standard atomic weight of each element in the unified
      ! atomic mass unit u (CIAAW conventional values), which is the unit
      ! these numbers are tabulated in and the only unit in which they are a
      ! measured quantity.
      real*8,  parameter :: melem_A_u(n_melem) = &
           [ 12.011d0, 15.999d0, 14.007d0, 24.305d0, 28.085d0,        &
             40.078d0, 22.990d0, 39.098d0, 32.06d0,  55.845d0 ]
      ! melem_A: the same weights in the code's mass unit, the hydrogen ATOM
      ! (global_parameters mu), which is the unit calc_rho, calc_mmw,
      ! comp_mass_per_H, the element census closure and the element
      ! diffusion all add masses in.  One u is 0.99223573 hydrogen atoms, so
      ! a weight left in u overstates every metal mass by 0.78 per cent.
      real*8,  parameter :: melem_A(n_melem) = amu_over_m_H*melem_A_u

      ! ---- metadata for each base (H/He/molecular) species ----
      ! Composition metadata for the non-metal species, so the EOS helpers
      ! (calc_rho / calc_ne / calc_ntot in utilities.f90) can eventually be
      ! driven from one table instead of hardcoded branches (§5.3 of
      ! docs/refactor_plan_system_composition_parser.md). Every species counts
      ! as ONE gas particle in the pressure/EOS particle sum.
      ! UNIT OF bsp_mass: the hydrogen ATOM, mu = 1.67353284e-24 g
      ! (parameters.f90), which is the mass unit of the density
      ! normalization. It is NOT the atomic mass unit u, so a species mass
      ! here is m_species/m_H and not its atomic weight in u.
      !   He I, He II, He III and the He 2^3S metastable all weigh
      !   m_He_over_m_H = 3.9715259, the helium-4 atom over the hydrogen
      !   atom (global_parameters), an ionization stage changing the mass
      !   by 1.4e-4, below the precision at which it is used; and
      !   HeH+ = 1 + m_He_over_m_H.
      ! Decision 14 of docs/development_plan_20260905_rev3.md section 10.5.
      ! H2, H2+ and H3+ are NUCLEUS COUNTS in the same unit, 2, 2 and 3 m_H:
      ! the molecular binding energies are 4.5, 2.7 and 4.4 eV, i.e. a mass
      ! defect of order 1e-9 of the molecule, which is below every digit
      ! this table carries and below the precision of the mass unit itself.
      ! Rows 11-13 (OH, H2O, CO) are the oxygen-chemistry carriers, built
      ! from the same H mass and the same element weights melem_A the metal
      ! ion columns use, so one nucleus weighs the same whether it is
      ! counted through a metal ion or through a molecule and calc_rho
      ! cannot count it twice. The oxygen and carbon bound in these
      ! molecules is removed from the metal ion totals by the ionization
      ! solve (the O I column means FREE ATOMIC oxygen when the oxygen
      ! chemistry is on), which is where that non-double-counting is made.
      integer, parameter :: n_bsp = 13   ! HI HII HeI HeII HeIII HeTR H2 H2+ H3+ HeH+ OH H2O CO
      ! f_sp species column for each base species
      integer, parameter :: bsp_fsp(n_bsp) = &
           [ isp_HI, isp_HII, isp_HeI, isp_HeII, isp_HeIII, isp_HeTR,   &
             isp_H2, isp_H2p, isp_H3p, isp_HeHp,                        &
             isp_OH, isp_H2O, isp_CO ]
      ! mass [m_H units] as used by calc_rho (see the unit note above)
      real*8,  parameter :: bsp_mass(n_bsp) = &
           [ 1.0d0, 1.0d0,                                              &
             m_He_over_m_H, m_He_over_m_H, m_He_over_m_H,               &
             m_He_over_m_H,                                             &
             2.0d0, 2.0d0, 3.0d0, 1.0d0 + m_He_over_m_H,                &
             1.0d0 + melem_A(iel_O), 2.0d0 + melem_A(iel_O),            &
             melem_A(iel_C) + melem_A(iel_O) ]
      ! net charge = free electrons contributed (calc_ne: each molecular
      ! ion carries +1; the He 2^3S triplet is neutral; OH, H2O and CO are
      ! neutral)
      integer, parameter :: bsp_charge(n_bsp) = &
           [ 0, 1, 0, 1, 2, 0,  0, 1, 1, 1,  0, 0, 0 ]
      ! H nuclei carried by one particle of the species
      integer, parameter :: bsp_nH(n_bsp) = &
           [ 1, 1, 0, 0, 0, 0,  2, 2, 3, 1,  1, 2, 0 ]
      ! He nuclei carried (the He nucleus in HeH+ is NOT in the free-He
      ! arrays; see the calc_rho comment)
      integer, parameter :: bsp_nHe(n_bsp) = &
           [ 0, 0, 1, 1, 1, 1,  0, 0, 0, 1,  0, 0, 0 ]
      ! O nuclei carried. The counterpart of bsp_nH for the oxygen element:
      ! without it every O nucleus bound in OH, H2O or CO disappears from
      ! the element total, and the total-oxygen closure
      !   n_O,tot = n(OI) + n(OII) + n(OIII) + n(OH) + n(H2O) + n(CO)
      ! silently loses half the element (docs/a2_oxygen_option_design.md
      ! sec. 4.1). Zero for every species that predates the oxygen option,
      ! so it is inert for a run without it.
      integer, parameter :: bsp_nO(n_bsp) = &
           [ 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  1, 1, 1 ]
      ! C nuclei carried; only CO has any.
      integer, parameter :: bsp_nC(n_bsp) = &
           [ 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 1 ]
      ! .true. where the column is an EXCITED LEVEL of another species, i.e.
      ! a sub-population already contained in that species' density, not an
      ! independent chemical species.  He 2^3S is the metastable triplet
      ! level of He I and the He I density the code carries (nhei, f_sp
      ! column isp_HeI) is the TOTAL He I population, triplet included --
      ! the singlet is formed by he_ground_singlet_density (composition),
      ! which is the one place that difference is taken and floored at zero.  Its mass, its gas particle and its helium nucleus
      ! are therefore already counted through He I, so EVERY BUDGET SUM
      ! (mass density, particle count, element nuclei, free electrons,
      ! collision partners) must skip the flagged column.  The level
      ! populations themselves are still transported and still scale with
      ! their parent, so only the sums skip them, never the updates.
      logical, parameter :: bsp_is_excited_level(n_bsp) = &
           [ .false., .false., .false., .false., .false., .true.,       &
             .false., .false., .false., .false.,                        &
             .false., .false., .false. ]
      ! human-readable species label (diagnostics only)
      character(len=5), parameter :: bsp_name(n_bsp) = &
           [ 'HI   ', 'HII  ', 'HeI  ', 'HeII ', 'HeIII', 'HeTR ',      &
             'H2   ', 'H2+  ', 'H3+  ', 'HeH+ ',                        &
             'OH   ', 'H2O  ', 'CO   ' ]

      ! ---- indices for each ion (canonical mion order; see table above) ----
      ! Used for index-based rate dispatch (rec/ion/cool_coeff_by_ion), so
      ! rate routines are selected without a string comparison on each call.
      integer, parameter :: im_CI    =  1, im_CII   =  2, im_CIII  =  3
      integer, parameter :: im_OI    =  4, im_OII   =  5, im_OIII  =  6
      integer, parameter :: im_NI    =  7, im_NII   =  8, im_NIII  =  9
      integer, parameter :: im_MgI   = 10, im_MgII  = 11, im_MgIII = 12
      integer, parameter :: im_SiI   = 13, im_SiII  = 14, im_SiIII = 15
      integer, parameter :: im_CaI   = 16, im_CaII  = 17, im_CaIII = 18
      integer, parameter :: im_NaI   = 19, im_NaII  = 20
      integer, parameter :: im_KI    = 21, im_KII   = 22
      ! Sulfur uses an underscore: Fortran identifiers are case-insensitive,
      ! so "im_SII" (S+) would collide with "im_SiI" (neutral Si).
      integer, parameter :: im_S_I   = 23, im_S_II  = 24
      integer, parameter :: im_FeI   = 25, im_FeII  = 26, im_FeIII = 27

      ! ---- metadata for each ion (length n_mion) ----
      ! f_sp species column for this ion
      integer, parameter :: mion_fsp(n_mion) = &
           [  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,          &
             19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,          &
             31, 32, 33 ]
      ! parent element index
      integer, parameter :: mion_elem(n_mion) = &
           [ iel_C,  iel_C,  iel_C,   iel_O,  iel_O,  iel_O,          &
             iel_N,  iel_N,  iel_N,   iel_Mg, iel_Mg, iel_Mg,         &
             iel_Si, iel_Si, iel_Si,  iel_Ca, iel_Ca, iel_Ca,        &
             iel_Na, iel_Na,  iel_K,  iel_K,   iel_S,  iel_S,         &
             iel_Fe, iel_Fe, iel_Fe ]
      ! ionization stage = net ionic charge
      integer, parameter :: mion_stage(n_mion) = &
           [ 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,                   &
             0, 1, 2,  0, 1, 2,  0, 1,  0, 1,  0, 1,  0, 1, 2 ]
      ! bremsstrahlung charge^2 weight (= stage^2; neutral -> 0, a no-op term)
      integer, parameter :: mion_z2(n_mion) = &
           [ 0, 1, 4,  0, 1, 4,  0, 1, 4,  0, 1, 4,                   &
             0, 1, 4,  0, 1, 4,  0, 1,  0, 1,  0, 1,  0, 1, 4 ]
      ! photo-ionizable? (.false. for the inert top stage)
      logical, parameter :: mion_isphot(n_mion) = &
           [ .true., .true., .false.,  .true., .true., .false.,      &
             .true., .true., .false.,  .true., .true., .false.,      &
             .true., .true., .false.,  .true., .true., .false.,      &
             .true., .false.,  .true., .false.,  .true., .false.,    &
             .true., .true., .false. ]
      ! column index into the photo cross-section / photo-rate tables
      ! (1..n_mphot for photo-ionizable ions, 0 otherwise)
      integer, parameter :: mion_iphot(n_mion) = &
           [ 1, 2, 0,  3, 4, 0,  5, 6, 0,  7, 8, 0,                   &
             9, 10, 0, 11, 12, 0, 13, 0,  14, 0,  15, 0,              &
             16, 17, 0 ]
      ! photoionization threshold [eV] (0 for the inert top stage).
      ! C/N/O thresholds are the literals used in PH_heat_HHe; Mg uses the
      ! named parameters so the values stay bit-identical. Si/Ca/Na/K/S
      ! thresholds are the VFKY96/VY95 Table 1 thresholds (= adopted NIST IP).
      real*8,  parameter :: mion_ethr(n_mion) = &
           [ 11.26d0, 24.38d0, 0.0d0,  13.62d0, 35.12d0, 0.0d0,      &
             14.53d0, 29.60d0, 0.0d0,  e_th_MgI, e_th_MgII, 0.0d0,   &
             8.152d0, 16.35d0, 0.0d0,  6.113d0, 11.87d0, 0.0d0,      &
             5.139d0, 0.0d0,   4.341d0, 0.0d0,   10.36d0, 0.0d0,     &
             7.902d0, 16.199d0, 0.0d0 ]
      ! contributes to metal line cooling (cool_M)?  (Huang+2023)
      ! adds CHIANTI Ca II H&K, Na I D, and Fe II to the C/N/O/Mg coolants
      ! already present (cool_coeff_metal dispatches each by name), plus Fe I
      ! (NIST f-values + Van Regemorter, Huang Fig 5). Si/K/S still have no
      ! line-cooling fit and stay excluded.
      ! Fe II cooling is density-dependent: eval_cool overrides the coronal
      ! coefficient with the multilevel-SE Lambda_eff(T,ne) (cool_FeII_ne), so
      ! the forbidden a6D fine-structure / metastable lines saturate (LTE) at
      ! the dense base instead of being overcounted ~1e4x.
      logical, parameter :: mion_iscool(n_mion) = &
           [ .true., .true., .false.,  .true., .true., .false.,      &
             .true., .true., .false.,  .true., .true., .false.,      &
             .false., .false., .false., .false., .true.,  .false.,   &
             .true.,  .false., .false., .false., .false., .false.,   &
             .true.,  .true.,  .false. ]
      ! human-readable ion label (diagnostics only)
      character(len=5), parameter :: mion_name(n_mion) = &
           [ 'CI   ', 'CII  ', 'CIII ', 'OI   ', 'OII  ', 'OIII ',    &
             'NI   ', 'NII  ', 'NIII ', 'MgI  ', 'MgII ', 'MgIII',    &
             'SiI  ', 'SiII ', 'SiIII', 'CaI  ', 'CaII ', 'CaIII',    &
             'NaI  ', 'NaII ', 'KI   ', 'KII  ', 'SI   ', 'SII  ',    &
             'FeI  ', 'FeII ', 'FeIII' ]

      ! ---- metadata for each element (length n_melem) ----
      ! nuclear charge (Z), used for the metal Gaunt factor
      integer, parameter :: melem_Z(n_melem)   = &
           [ 6, 8, 7, 12, 14, 20, 11, 19, 16, 26 ]
      ! The atomic weights melem_A_u [u] and melem_A [m_H], which give the
      ! metal mass contribution to the gas mass density and to the mean
      ! molecular weight (eos_include_metals), are declared with the element
      ! indices at the top of this module, because the oxygen-chemistry rows
      ! of bsp_mass are built from them.
      ! index of the neutral ion in the mion list
      integer, parameter :: melem_i0(n_melem)  = &
           [ 1, 4, 7, 10, 13, 16, 19, 21, 23, 25 ]
      ! highest ionization stage solved for this element
      ! (2 = up to ++ for C/O/N/Mg/Si/Ca/Fe; 1 = up to + for Na/K/S)
      integer, parameter :: melem_top(n_melem) = &
           [ 2, 2, 2, 2, 2, 2, 1, 1, 1, 2 ]
      ! element label (diagnostics only)
      character(len=2), parameter :: melem_name(n_melem) = &
           [ 'C ', 'O ', 'N ', 'Mg', 'Si', 'Ca', 'Na', 'K ', 'S ', 'Fe' ]

      end module species_table
