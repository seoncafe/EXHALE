   module J_incident
   ! Energy-dependent incident spectrum.
   !
   ! ONE SPECTRUM TYPE BUILDS EVERY BAND. "Spectrum type:" in input.inp
   ! states the
   ! spectrum of the WHOLE photon grid, the XUV and the part below 13.6 eV
   ! alike -- the part where the He 2^3S metastable (4.80 eV) and the low-IP
   ! metals absorb:
   !   Power-law      J_inc below: the power law of "Power-law index",
   !                  normalized on [e_low, e_mid] and evaluated wherever the
   !                  grid reaches. Below e_low it is an extrapolation of an
   !                  EUV fit, which is what the type means, not an accident.
   !   Planck         planck_stellar_flux_eV below: the photospheric
   !                  blackbody pi B_nu(T_eff) (R_star/a)^2 of "Stellar Teff"
   !                  and "Stellar radius".
   !   Load           the SED table (sed_read), which must cover the grid.
   ! No band is filled from a type the input did not select.

   use global_parameters

   implicit none

   ! Frequency of a one-electronvolt photon, 1 eV/h [Hz]. The single
   ! definition: excited_hydrogen's Balmer continuum reads it from here.
   real*8, parameter :: eV2Hz = 2.417989242d14

   ! Ionization threshold of hydrogen in n = 2, the head of the Balmer
   ! continuum: e_th_HI/n^2 with n = 2, i.e. 3.3996 eV = 3647 A. Hydrogenic
   ! scaling is exact here to the fine structure: NIST ASD puts the 1s -> 2s
   ! excitation at 82258.9543 cm^-1 = 10.198834 eV, so the n = 2 binding
   ! energy is 13.598435 - 10.198834 = 3.399601 eV against e_th_HI/4 =
   ! 3.399609 eV.
   !
   ! ONE DEFINITION, three consumers, as for every other threshold
   ! (parameters.f90 "ONE DEFINITION EACH"): the band excited_hydrogen
   ! integrates the H(n=2) photoionization over, and the floor of the photon
   ! grid in sed_read (a loaded table must reach it) and in
   ! set_energy_vectors (the analytic types), whenever the excited-hydrogen
   ! coupling is armed. H(n=2) is an absorber like the metastable helium and
   ! the low-IP metals, and its threshold enters the grid floor the same way
   ! theirs do.
   real*8, parameter :: e_th_HI_n2 = e_th_HI/4.0d0

   contains

   logical function spectrum_is_planck()
   ! Whether the run's spectrum type is the photospheric blackbody. sp_type
   ! (input.inp "Spectrum type:") is the single record of the choice, so the
   ! comparison lives in one place. It is deferred-length and is allocated by
   ! input_read: before that -- a test driver that sets the spectrum flags
   ! itself and never reads an input file -- there is no type to compare, and
   ! the answer is no.
   spectrum_is_planck = .false.
   if (.not. allocated(sp_type)) return
   spectrum_is_planck = (trim(sp_type) .eq. 'Planck')
   end function spectrum_is_planck

   double precision function planck_stellar_flux_nu(Tstar, R_over_a, nu)
   ! Flux density [erg cm^-2 s^-1 Hz^-1] at the planet's orbit from a star
   ! whose photosphere radiates as a blackbody of effective temperature
   ! Tstar, seen at R_over_a = R_star/a:
   !   F_nu = pi B_nu(Tstar) (R_star/a)^2 ,
   !   B_nu = (2 h nu^3/c^2) / (exp(h nu / k Tstar) - 1)
   ! (Rybicki & Lightman 1979, eqs. 1.51 and 1.13). Integrated over all
   ! frequencies this is sigma Tstar^4 (R_star/a)^2.
   !
   ! VALIDITY. A real photosphere is not a blackbody: line blanketing and
   ! the Balmer jump depress the near-ultraviolet of an F or G star below
   ! this curve, so the 4.8-13.6 eV field of "Spectrum type: Planck" is an
   ! upper bound there. Above 13.6 eV the chromospheric and coronal EUV of
   ! an active star exceeds it by orders of magnitude; a run that needs both
   ! bands right needs a measured spectrum ("Spectrum type: Load").
   real*8, intent(in) :: Tstar, R_over_a, nu
   real*8 :: x, Bnu

   planck_stellar_flux_nu = 0.0d0
   if (Tstar .le. 0.0d0 .or. nu .le. 0.0d0) return

   x = hp_erg*nu/(kb_erg*Tstar)
   if (x .gt. 7.0d2) then
      ! Wien tail: 1/(exp(x) - 1) = exp(-x)/(1 - exp(-x)) and exp(-x) is
      ! below 1e-304 here, so the two forms agree to double precision while
      ! exp(+x) would overflow. The XUV top of the grid, 1.24 keV, is
      ! x = 2228 for a 6459 K star.
      Bnu = (2.0d0*hp_erg*nu**3.0/c_light**2.0)*exp(-x)
   else
      Bnu = (2.0d0*hp_erg*nu**3.0/c_light**2.0)/(exp(x) - 1.0d0)
   endif
   planck_stellar_flux_nu = pi*Bnu*R_over_a**2.0

   end function planck_stellar_flux_nu

   double precision function planck_stellar_flux_eV(Tstar, R_over_a, E)
   ! The same field per unit photon ENERGY [erg cm^-2 s^-1 eV^-1], which is
   ! the unit of F_XUV on the photon grid: F_E = F_nu dnu/dE with nu = E/h,
   ! so dnu/dE = eV2Hz for E in eV.
   real*8, intent(in) :: Tstar, R_over_a, E

   planck_stellar_flux_eV =                                                &
      planck_stellar_flux_nu(Tstar, R_over_a, E*eV2Hz)*eV2Hz

   end function planck_stellar_flux_eV

   logical function spectrum_covers_eV(E)
   ! Whether the run's spectrum type states a field at the photon energy E
   ! [eV]. The two analytic types are defined at every energy; a loaded
   ! table is defined only over the band it was read on, from
   ! loaded_table_floor_eV to the top of the photon grid; a monochromatic
   ! run states a field at one
   ! energy and nowhere else. A band of the code that needs a field outside
   ! the covered range has to stop the run rather than fill the band from a
   ! type the input did not select.
   real*8, intent(in) :: E

   spectrum_covers_eV = .false.
   if (is_PL_sed .or. spectrum_is_planck()) then
      spectrum_covers_eV = .true.
   else if (do_read_sed) then
      if (allocated(e_sed_node) .and. size(e_sed_node) .ge. 2)             &
         spectrum_covers_eV = (E .ge. loaded_table_floor_eV() .and.        &
                               E .le. e_sed_node(size(e_sed_node)))
   endif

   end function spectrum_covers_eV

   double precision function loaded_table_floor_eV()
   ! Lowest photon energy [eV] a loaded table states a field at.
   !
   ! read_sed selects the rows of the file down to the photon-grid floor
   ! e_low and stops the run when the file ends above it, so a table that
   ! got this far covers [e_low, e_top]. The lowest SELECTED row can still
   ! sit up to one tabulated row above e_low, because the selection stops at
   ! the first row below the floor; the field between the two is the table's
   ! own, over less than one of its rows, and it is read by extrapolating
   ! the first segment in stellar_flux_eV rather than declared absent.
   !
   ! The TABLE ROWS decide this, not the photon grid: e_v holds bin centres,
   ! and the centre of the lowest bin sits above the lowest row.
   real*8 :: e_bot

   loaded_table_floor_eV = 0.0d0
   if (.not. allocated(e_sed_node)) return
   if (size(e_sed_node) .lt. 2) return
   e_bot = e_sed_node(1)
   if (e_low .gt. 0.0d0 .and. e_low .lt. e_bot) e_bot = e_low
   loaded_table_floor_eV = e_bot

   end function loaded_table_floor_eV

   double precision function stellar_flux_eV(E)
   ! The incident stellar flux per unit photon energy [erg cm^-2 s^-1 eV^-1]
   ! at the planet's orbit, from the run's spectrum type, at ANY photon
   ! energy: the single statement of what field the run has at E, which
   ! every band of the code reads. Undiluted, as J_inc and
   ! planck_stellar_flux_eV are: the caller applies dayside_dilution() where
   ! the beam is used.
   !   Power-law   the power law of PLind, normalized on [e_low, e_mid] and
   !               [e_mid, e_top] and evaluated wherever it is asked; below
   !               e_low that is an extrapolation, which is what the type
   !               means.
   !   Planck      pi B_nu(T_eff) (R_star/a)^2 dnu/dE.
   !   Load        the loaded table, log-log interpolated on its own rows,
   !               down to loaded_table_floor_eV. Outside
   !               that band the table states nothing and this returns zero;
   !               ask spectrum_covers_eV first, and stop the run there
   !               rather than integrate a zero the input never stated.
   real*8, intent(in) :: E
   integer :: k, klo, khi, n_node
   real*8  :: t, f_lo, f_hi

   stellar_flux_eV = 0.0d0
   if (E .le. 0.0d0) return

   if (is_PL_sed) then

      stellar_flux_eV = J_inc(E)

   else if (spectrum_is_planck()) then

      stellar_flux_eV = planck_stellar_flux_eV(T_star_eff,                 &
                           R_star/max(a_orb,1.0d-30), E)

   else if (do_read_sed) then

      if (.not. allocated(e_sed_node) .or. .not. allocated(F_sed_node))    &
         return
      n_node = size(e_sed_node)
      if (n_node .lt. 2) return
      if (E .lt. loaded_table_floor_eV() .or. E .gt. e_sed_node(n_node))   &
         return

      ! The TABLE ROWS are interpolated, not the photon grid: the grid holds
      ! bin centres and is cut at the ionization thresholds, so its points
      ! are a property of the absorbers while the field is a property of the
      ! file. e_sed_node is ascending in energy and undiluted.
      ! Below the first row the bracket is the first segment, and the
      ! interpolation weight below is then negative: the field of a table
      ! that reaches below the grid floor, over the less than one row
      ! between the floor and the lowest selected row.
      klo = 1
      khi = n_node
      do while (khi - klo .gt. 1)
         k = (khi + klo)/2
         if (e_sed_node(k) .gt. E) then
            khi = k
         else
            klo = k
         endif
      enddo

      f_lo = F_sed_node(klo)
      f_hi = F_sed_node(khi)
      if (f_lo .gt. 0.0d0 .and. f_hi .gt. 0.0d0 .and.                      &
          e_sed_node(khi) .gt. e_sed_node(klo)) then
         ! A stellar spectrum spans decades over a few eV; interpolate the
         ! logarithm, which is exact for the power-law segments a tabulated
         ! continuum is made of.
         t = log(E/e_sed_node(klo))/log(e_sed_node(khi)/e_sed_node(klo))
         stellar_flux_eV = exp((1.0d0-t)*log(f_lo) + t*log(f_hi))
      else if (e_sed_node(khi) .gt. e_sed_node(klo)) then
         t = (E - e_sed_node(klo))/(e_sed_node(khi) - e_sed_node(klo))
         stellar_flux_eV = (1.0d0-t)*f_lo + t*f_hi
      else
         stellar_flux_eV = f_lo
      endif

   endif

   end function stellar_flux_eV

   double precision function J_inc(E)
   real*8, intent(in) ::  E
   real*8 :: Lrapp,J_X,J_EUV
   real*8 :: JEUVnorm,JXnorm
   real*8 :: P1
   
   ! Substitution
   P1 = PLind + 1.0
   
   ! Ratio of luminosities
   Lrapp = 10.0**(LX-LEUV)
   
   ! X-ray Flux
   J_X = Lrapp*J_XUV/(1.0+Lrapp)
   
   ! EUV flux
   J_EUV = J_XUV/(1.0+Lrapp)
   
   ! Normalizations for different power-law index
   if (PLind.eq.(-1.0)) then
      
      JEUVnorm = J_EUV/log(e_mid/e_low)
      if (thereis_Xray) then
         JXnorm   = J_X/log(e_top/e_mid)
      else      		
         JXnorm   = 0.0
      endif
               
   else
   
      JEUVnorm = J_EUV*P1/(e_mid**P1 - e_low**P1)
      if (thereis_Xray) then
         JXnorm   = J_X*P1/(e_top**P1 - e_mid**P1)
      else      	
         JXnorm   = 0.0
      endif
   
   endif      
   
   ! Parametrization of incident spectrum
   if(E.lt.e_mid) then
      J_inc = JEUVnorm*E**PLind
   else
      J_inc = JXnorm*E**PLind
   endif
   
   end function
   
   ! End of module
   end module J_incident
   
