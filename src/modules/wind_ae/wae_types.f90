      module wae_types
      ! Derived types ported from wind-ae wind.h (PARAMLIST, EQNVARS,
      ! I_EQNVARS, VARLIST). Field names and meanings are kept identical
      ! to the C structs so the soe/difeq port reads 1:1 against the C.
      use wae_config, only: nsp => wae_nspecies, totalpts => wae_totalpts
      implicit none

      ! PARAMLIST: planetary system + BC + term-indicator + technical +
      ! physics parameters (set_parameters() in io.c).
      type :: wae_paramlist
         ! planetary system properties
         real*8 :: Mp, Rp, Mstar, semimajor, Ftot, Lstar, H0
         ! boundary conditions
         real*8 :: Rmin, Rmax, T_rmin
         real*8 :: Ys_rmin(nsp)
         real*8 :: rho_rmin
         real*8 :: Ncol_sp(nsp)
         real*8 :: erf_drop(2)
         ! term indicators
         integer :: lyacool
         real*8  :: tidalforce
         real*8  :: bolo_heat_cool
         integer :: integrate_outward
         ! technical parameters
         real*8  :: breezeparam, rapidity, erfn, mach_limit
         ! physics parameters
         real*8       :: HX(nsp)            ! mass fraction per species
         character(len=10) :: species(nsp) ! Roman-numeral name
         real*8       :: atomic_mass(nsp)   ! [g]
         integer      :: Z(nsp)             ! atomic number
         integer      :: N_e(nsp)           ! electrons in atom
         real*8       :: molec_adjust
      end type wae_paramlist

      ! EQNVARS: full-grid solution arrays (units per wind.h comments).
      type :: wae_eqnvars
         real*8 :: r(totalpts)
         real*8 :: rho(totalpts)
         real*8 :: v(totalpts)
         real*8 :: z(totalpts)
         real*8 :: Ys(totalpts, nsp)
         real*8 :: Ncol(totalpts, nsp)
         real*8 :: T(totalpts)
         real*8 :: q(totalpts)
      end type wae_eqnvars

      ! I_EQNVARS: one grid point's variables.
      type :: wae_i_eqnvars
         real*8 :: r, rho, v, z
         real*8 :: Ys(nsp), Ncol(nsp)
         real*8 :: T, q
      end type wae_i_eqnvars

      ! VARLIST: one value per relaxed variable (numerical-derivative steps).
      type :: wae_varlist
         real*8 :: rho, v, T
         real*8 :: Ys(nsp), Ncol(nsp)
         real*8 :: z
      end type wae_varlist

      end module wae_types
