      module wae_config
      ! Compile-time configuration for the Wind-AE Fortran port, mirroring
      ! wind-ae src/defs.h. Stage 1 keeps NSPECIES fixed at 2 (H/He); the
      ! dynamic-species generalization is out of scope (metals are EXHALE's
      ! job -- see md_backup/wind_ae_fortran_port_plan.md).
      implicit none

      integer, parameter :: wae_nspecies = 2          ! defs.h NSPECIES
      integer, parameter :: wae_ne  = 4 + 2*wae_nspecies
      integer, parameter :: wae_nb  = 2 + wae_nspecies
      integer, parameter :: wae_m   = 1501            ! defs.h M
      integer, parameter :: wae_inpts  = 0            ! defs.h INPTS
      integer, parameter :: wae_addpts = 1000         ! defs.h ADDPTS
      integer, parameter :: wae_totalpts = wae_inpts + wae_m + wae_addpts

      real*8,  parameter :: wae_NCOL0 = 1.0d17        ! defs.h NCOL0
      ! Lower limit of the x-ray energy range, in erg (glq_rates.c).
      real*8,  parameter :: wae_xray_lo = 6.408707d-11

      ! Physical constants / unit scales (defs.h). CS0 = sqrt(K*T0/MH).
      real*8,  parameter :: wae_T0   = 1.0d4
      real*8,  parameter :: wae_RHO0 = 1.0d-15
      real*8,  parameter :: wae_gamma_atomic = 5.0d0/3.0d0

      ! Relaxation controls (defs.h)
      integer, parameter :: wae_itmax = 100
      real*8,  parameter :: wae_conv  = 1.0d-10
      real*8,  parameter :: wae_slowc = 1.0d-2
      real*8,  parameter :: wae_derivdiv = 1.0d7

      ! Convergence scales (relax.h). NOTE: RHOSCALE is 100 for most
      ! planets but 10 for strongly-bound ones (the Python wrapper rewrites
      ! defs.h); carried as a runtime variable, defaulted to 100 here.
      real*8 :: wae_rhoscale  = 100.0d0
      real*8,  parameter :: wae_vscale   = 2.0d0
      real*8,  parameter :: wae_tempscale= 1.0d0
      real*8,  parameter :: wae_ncolscale= 100.0d0
      real*8,  parameter :: wae_fpscale  = 1.0d0
      real*8,  parameter :: wae_zscale   = 3.0d0

      end module wae_config
