      module co_photodissociation
      ! Photodissociation of CO on the 912-1201 A Lyman-Werner beam:
      !
      !     CO + hv  ->  CO(predissociating state)  ->  C + O .
      !
      ! This is the destruction channel that takes over from He+ charge
      ! transfer once the star-ward CO and H2 columns thin out, and together
      ! with that channel it is the whole of the one-sided CO destruction
      ! model.
      !
      ! ---------------------------------------------------------------
      ! 1. THE RATE, AND WHY IT IS NOT WRITTEN AS VISSER'S k0 THETA
      !
      ! Visser, van Dishoeck and Black (2009), A&A 503, 323, sec. 3.3
      ! eq. (2) write the rate of an isotopologue as
      !
      !     k_i = chi k0_i Theta_i exp(-gamma A_V) ,
      !
      ! with k0_i the unattenuated rate IN A STATED RADIATION FIELD.  Their
      ! sec. 3.2 says how much that statement matters: "Clearly, the rate
      ! depends on the choice of radiation field.  If we adopt Habing
      ! (1968), Gondhalekar et al. (1980) or Mathis et al. (1983) instead of
      ! Draine (1978), the photodissociation rate becomes 2.0, 2.0 or
      ! 2.3 x 10^-10 s^-1, respectively", against their own 2.6e-10.
      !
      ! EXHALE has a band flux and not an interstellar field, so k0 is not
      ! the quantity to carry.  What is carried is the band-mean CO
      ! dissociation cross section, multiplied by the band PHOTON flux of
      ! the run's own beam:
      !
      !     k_CO = ( F_LW / <hv>_band ) sigma_CO_band Theta(N_CO, N_H2)
      !            exp(-tau_cont) ,
      !
      ! which is the same contraction lyman_werner.f90 makes for H2 and
      ! water_photolysis.f90 for the H2O and OH continua.  Four absorbers of
      ! one beam, one convention.
      !
      ! THE ONE-SIDED REDUCTION IS ALREADY APPLIED, and it is applied in one
      ! place for all four absorbers.  Visser's sec. 3.3: "Equation (2)
      ! assumes the radiation is coming from all directions.  If this is not
      ! the case, such as for a cloud irradiated only from one side, k0,i
      ! should be reduced accordingly."  F_LW here is fuv_band_flux(ib_LW),
      ! which is the stated band flux times dayside_dilution(); no second
      ! geometric factor is applied, because applying one would give the
      ! four absorbers of one beam two different geometries.
      !
      ! ---------------------------------------------------------------
      ! 2. THE BAND-MEAN CROSS SECTION, AND THE WINDOW CO ACTUALLY ABSORBS IN
      !
      ! The beam is 912-1201 A (lyman_werner.f90 sec. 1).  CO does not
      ! absorb across all of it: its 37 predissociating bands lie between
      ! 912.7037 and 1076.0796 A (Visser Table 1) and its dissociation
      ! continuum ends at the C=O bond energy, 11.1157 eV = 1115.4 A, so the
      ! CO window is 912-1118 A and the remaining 83 A of the beam are
      ! transparent to it.  sigma_CO_band below is therefore the mean over
      ! the WHOLE beam of a cross section that is zero above 1118 A:
      !
      !     sigma_CO_band = int_912^1118 sigma(lam) lam dlam
      !                     / int_912^1201 lam dlam ,
      !
      ! the lam weighting being the photon flux of a flat F_lambda band, the
      ! same spectral assumption <hv>_band = 2hc/(912 A + 1201 A) makes.
      ! MEASURED here (2026-09-06) from the CO photodissociation cross
      ! section of VULCAN/thermo/photo_cross/CO/CO_cross.csv (Heays et al.
      ! 2017 as concatenated by Photochem; 0.1 nm grid, peak 2.599e-16 cm^2
      ! at 1076 A, which is the E0 band Visser call the strongest
      ! contributor at the cloud edge):
      !
      !     sigma_CO_band = 1.0160e-17 cm^2
      !
      ! and the same integral over 912-1110 A alone, the band as it stood
      ! before the Lyman-Werner interval was renormalized, is 1.5496e-17
      ! cm^2.  The ratio of the two is the flat-F_lambda share of the beam's
      ! photons that lie in the CO window.
      !
      ! FLAT F_LAMBDA IS ASSUMED FOR BOTH SPECTRUM TYPES, and that is a
      ! stated approximation rather than an oversight.  With
      ! "Spectrum type: Load" the run integrates the file over 912-1201 A to
      ! get F_LW (sed_read.f90) but keeps no shape inside the band, and
      ! every other absorber of this beam -- the H2 self-shielding table,
      ! the H2O and OH band cross sections -- is normalized per photon of a
      ! flat-F_lambda band.  Giving CO alone a shape-aware cross section
      ! would put the four absorbers of one beam on two conventions; the
      ! error it would remove is the tilt of the stellar spectrum across
      ! 912-1118 A, which is the same tilt the other three already carry.
      !
      ! CROSS-CHECK AGAINST VISSER'S OWN NORMALIZATION.  Their k0 belongs to
      ! the Draine (1978) field, whose 912-1110 A photon flux at chi = 1 is
      ! 1.232e7 cm^-2 s^-1 (Draine & Bertoldi 1996 Table 1, the F/chi column,
      ! READ through lyman_werner.f90 sec. 1).  Their Table 6 header gives
      ! k0 = 2.590e-10 s^-1, so the cross section implied by their data over
      ! that band is 2.102e-17 cm^2 against the 1.550e-17 measured above on
      ! the same band: OUR VALUE IS 26 PER CENT LOWER (MEASURED).  Visser
      ! state their absolute rate accuracy as "about 20%" (their sec. 3.5),
      ! so the two do not quite meet inside that figure and the difference
      ! is recorded here rather than absorbed.  Two things it can be: the
      ! Draine field is redder across 912-1118 A than a flat F_lambda and so
      ! weights the strong 1076 A band more heavily, and the 0.1 nm
      ! continuum table is a different reduction of the same molecular data
      ! from the band-resolved line list the shielding function is built on.
      ! The cross-section route is kept, because it is the run's own beam
      ! that multiplies it.  The test co_shielding_table.f90 measures this
      ! ratio, so it cannot drift unrecorded.
      !
      ! ---------------------------------------------------------------
      ! 3. THE PHOTON ENERGY OF ONE DISSOCIATION EVENT
      !
      ! The fragments leave with hv - D0(CO), and D0(CO) comes from the one
      ! formation-energy table (oxygen_rates::photolysis_threshold_erg), so
      ! the two CO destruction channels cannot disagree about how much
      ! energy the C=O bond holds.  hv is NOT the flat-band mean: CO absorbs
      ! in 37 discrete bands, not across the beam, and the flat-band mean
      ! 11.7354 eV is set by the band edges 912 and 1201 A, including the
      ! 83 A longward of the CO window where CO does not absorb; it lies
      ! only 0.62 eV above the 11.1157 eV threshold, so a deposit taken from
      ! it would be an artefact of those edges.
      !
      ! What is used is the OSCILLATOR-STRENGTH WEIGHTED MEAN over Visser's
      ! Table 1, sum_i f_i hv_i / sum_i f_i over the 37 bands.  MEASURED
      ! here from that table: 12.8674 eV, so the deposit is
      ! 12.8674 - 11.1157 = 1.7517 eV per event.  Weighting instead by
      ! f_i eta_i (eta the dissociation efficiency of the band, also in
      ! their Table 1) gives 12.9149 eV and by f_i eta_i lambda_i (the
      ! photon flux of a flat F_lambda) 12.8823 eV, so the choice of weight
      ! is worth 0.05 eV out of 1.75 and the constant below is not sensitive
      ! to it.
      !
      ! ---------------------------------------------------------------
      ! 4. THE CELL MEAN, AND WHY A FACE VALUE WILL NOT DO
      !
      ! Theta falls by more than three decades along the CO column axis of
      ! the table and by more than six along the H2 axis.  A cell at the CO
      ! front carries a large fraction of a decade of column, so a rate
      ! taken at one face and applied across the cell is wrong in one
      ! direction everywhere in that cell.  This is the same argument, and
      ! the same quadrature, that lyman_werner_dissociation_rate_cell_mean
      ! makes for H2 and that water_photolysis makes for the continua of the
      ! same beam: composite three-point Gauss-Legendre on segments cut
      ! geometrically in the column, at most seg_dex decades of column and
      ! seg_dtau of continuum depth wide.
      !
      ! Two columns vary across the cell here where H2 has one.  The
      ! segments are cut geometrically in whichever of the two spans more
      ! decades across the cell, and the segment count satisfies the
      ! seg_dex criterion for BOTH columns and the seg_dtau criterion for
      ! the continuum, so the finest of the three wins.
      !
      ! ---------------------------------------------------------------
      ! References: Visser, van Dishoeck & Black (2009), A&A 503, 323
      ! (publisher PDF, references/Visser_2009A&A_503_323.pdf), Table 1
      ! (the 37 bands, their wavelengths, oscillator strengths and
      ! dissociation efficiencies), secs. 3.2, 3.3 and 3.5; Draine &
      ! Bertoldi (1996) ApJ 468, 269 Table 1 (the Draine-field band photon
      ! flux used for the cross-check); Heays, Bosman & van Dishoeck (2017)
      ! A&A 602, A105 through VULCAN/thermo/photo_cross/CO/CO_cross.csv.

      use co_self_shielding_table, only: co_self_shielding
      use oxygen_rates, only: photolysis_threshold_erg, ich_CO_C_O

      implicit none
      private
      public :: co_photodissociation_rate,                                &
                co_photodissociation_rate_cell_mean,                      &
                sigma_CO_band, e_co_photon_erg,                           &
                heat_per_co_dissociation

      integer, parameter :: dp = kind(1.0d0)

      ! Mean photon energy of a flat-F_lambda 912-1201 A band [erg], the
      ! SAME constant lyman_werner.f90 carries: F_LW/e_lw_photon_erg is the
      ! band photon flux, and sigma_CO_band below is normalized per photon
      ! of that band.  It is repeated here rather than imported because
      ! importing lyman_werner would make the CO rate depend on the H2
      ! self-shielding table; the test co_shielding_table.f90 asserts the
      ! two are the same number.
      real(dp), parameter :: e_lw_photon_erg = 1.88021d-11

      ! Band-mean CO dissociation cross section per photon of the
      ! 912-1201 A beam [cm^2] (header sec. 2).
      real(dp), parameter :: sigma_CO_band = 1.0160d-17

      ! Oscillator-strength weighted mean photon energy of a CO
      ! dissociation event, 12.8674 eV in erg (header sec. 3).
      real(dp), parameter :: e_co_photon_erg = 2.061583d-11

      contains

      ! LOCAL CO photodissociation rate [s^-1] at one point of the column:
      ! the band flux F_LW [erg cm^-2 s^-1] at the planet, the star-ward CO
      ! and H2 columns [cm^-2] and the star-ward continuum depth tau_cont of
      ! the same 912-1201 A interval, all at that point.  What a grid cell
      ! needs is the mean of this over the cell, which is the function
      ! below.
      !
      ! tau_cont is the H2O and OH continuum of the oxygen chemistry
      ! (water_photolysis.f90).  CO itself contributes NOTHING to it: its
      ! lines shield CO, and H2's lines shield CO through N_H2 inside
      ! Theta, which is exactly the construction Visser's eq. (2) makes.
      ! Adding a CO term to tau_cont would shield the H2O and OH continua
      ! with a line absorber whose equivalent width is already inside
      ! Theta, and would count the same photons twice.
      double precision function co_photodissociation_rate(F_LW, N_CO,     &
                                N_H2, tau_cont) result(k)
      real(dp), intent(in) :: F_LW, N_CO, N_H2, tau_cont
      if (F_LW .le. 0.0d0) then
         k = 0.0d0
         return
      endif
      k = F_LW/e_lw_photon_erg*sigma_CO_band                              &
          *co_self_shielding(N_CO, N_H2)
      if (tau_cont .gt. 0.0d0) k = k*exp(-tau_cont)
      end function co_photodissociation_rate

      ! MEAN of that rate over one grid cell [s^-1], for the star-ward face
      ! values (N_CO_out, N_H2_out, tau_out) and the inner face values
      ! (N_CO_in, N_H2_in, tau_in) of the two columns and of the 912-1201 A
      ! continuum depth.  Header sec. 4 states why the mean and how it is
      ! taken; within a cell every absorber density is uniform, which is the
      ! rectangle rule the column integration itself uses, so both columns
      ! and the depth run linearly across the cell.
      double precision function co_photodissociation_rate_cell_mean(F_LW, &
                                N_CO_out, N_CO_in, N_H2_out, N_H2_in,     &
                                tau_out, tau_in) result(k)
      real(dp), intent(in) :: F_LW, N_CO_out, N_CO_in, N_H2_out, N_H2_in
      real(dp), intent(in) :: tau_out, tau_in
      ! Three-point Gauss-Legendre on [0,1]: nodes (1 -+ sqrt(3/5))/2 and
      ! 1/2, weights 5/18, 8/18, 5/18.  Exact for polynomials of degree 5.
      integer,  parameter :: n_gl = 3
      real(dp), parameter :: gl_s(n_gl) = (/ 0.1127016653792583d0,        &
                                             0.5d0,                       &
                                             0.8872983346207417d0 /)
      real(dp), parameter :: gl_w(n_gl) = (/ 5.0d0/18.0d0, 8.0d0/18.0d0,  &
                                             5.0d0/18.0d0 /)
      real(dp), parameter :: seg_dex  = 0.05d0
      real(dp), parameter :: seg_dtau = 0.5d0
      real(dp), parameter :: col_head_ratio = 1.0d-10
      integer,  parameter :: nseg_max = 256
      real(dp) :: se(0:nseg_max+1)
      real(dp) :: Nco_lo, Nco_hi, dNco, Nh2_lo, Nh2_hi, dNh2
      real(dp) :: tau_lo, dtau, dex_co, dex_h2, dstart, ratio, Nx
      real(dp) :: s0, s1, ds, ss, acc
      integer  :: nseg, nedge, i, g
      logical  :: cut_on_co

      if (F_LW .le. 0.0d0) then
         k = 0.0d0
         return
      endif

      Nco_lo = max(N_CO_out, 0.0d0)
      Nco_hi = max(N_CO_in,  Nco_lo)
      dNco   = Nco_hi - Nco_lo
      Nh2_lo = max(N_H2_out, 0.0d0)
      Nh2_hi = max(N_H2_in,  Nh2_lo)
      dNh2   = Nh2_hi - Nh2_lo
      tau_lo = max(tau_out, 0.0d0)
      dtau   = max(tau_in - tau_lo, 0.0d0)

      ! Decades of each column across the cell.  A cell whose star-ward
      ! column is zero has no geometric starting point, so the count starts
      ! at col_head_ratio of the inner-face column and the remainder below
      ! that is one further segment; the table returns its edge value below
      ! its bottom knot, so the rate is constant over that head and one
      ! segment integrates it exactly.
      dex_co = 0.0d0
      dex_h2 = 0.0d0
      if (dNco .gt. 0.0d0)                                                &
         dex_co = log10(Nco_hi/max(Nco_lo, Nco_hi*col_head_ratio))
      if (dNh2 .gt. 0.0d0)                                                &
         dex_h2 = log10(Nh2_hi/max(Nh2_lo, Nh2_hi*col_head_ratio))

      nseg = 1
      nseg = max(nseg, int(dex_co/seg_dex) + 1)
      nseg = max(nseg, int(dex_h2/seg_dex) + 1)
      nseg = max(nseg, int(dtau/seg_dtau) + 1)
      nseg = min(nseg, nseg_max)

      ! Segment edges as positions s across the cell.  They are cut
      ! geometrically in whichever column spans more decades, because that
      ! is the variable the rate is a power law of over the widest range;
      ! the other column and the continuum depth are resolved by the
      ! segment COUNT, which already satisfies their own criteria.
      cut_on_co = (dex_co .ge. dex_h2) .and. (dNco .gt. 0.0d0)
      se(0) = 0.0d0
      nedge = 0
      if (cut_on_co .or. (dNh2 .gt. 0.0d0 .and. dex_h2 .gt. 0.0d0)) then
         if (cut_on_co) then
            dstart = max(Nco_lo, Nco_hi*col_head_ratio)
            if (dstart .gt. Nco_lo) then
               nedge = 1
               se(1) = (dstart - Nco_lo)/dNco
            endif
            ratio = (Nco_hi/dstart)**(1.0d0/dble(nseg))
            Nx    = dstart
            do i = 1,nseg-1
               Nx = Nx*ratio
               se(nedge+i) = (Nx - Nco_lo)/dNco
            enddo
         else
            dstart = max(Nh2_lo, Nh2_hi*col_head_ratio)
            if (dstart .gt. Nh2_lo) then
               nedge = 1
               se(1) = (dstart - Nh2_lo)/dNh2
            endif
            ratio = (Nh2_hi/dstart)**(1.0d0/dble(nseg))
            Nx    = dstart
            do i = 1,nseg-1
               Nx = Nx*ratio
               se(nedge+i) = (Nx - Nh2_lo)/dNh2
            enddo
         endif
         nedge = nedge + nseg
      else
         do i = 1,nseg-1
            se(i) = dble(i)/dble(nseg)
         enddo
         nedge = nseg
      endif
      se(nedge) = 1.0d0

      acc = 0.0d0
      do i = 1,nedge
         s0 = se(i-1)
         s1 = se(i)
         ds = s1 - s0
         if (ds .le. 0.0d0) cycle
         do g = 1,n_gl
            ss  = s0 + ds*gl_s(g)
            acc = acc + ds*gl_w(g)                                        &
                *co_photodissociation_rate(F_LW,                          &
                        Nco_lo + ss*dNco, Nh2_lo + ss*dNh2,               &
                        tau_lo + ss*dtau)
         enddo
      enddo
      k = acc
      end function co_photodissociation_rate_cell_mean

      ! Kinetic energy left to the C and O fragments by one CO
      ! photodissociation [erg]: the mean photon energy of the event less
      ! the 0 K bond energy, which the photon pays and the gas does not.
      ! The threshold is the difference of the SAME formation energies the
      ! He+ charge-transfer channel is charged with, so the two channels of
      ! one molecule cannot state two bond energies.
      double precision function heat_per_co_dissociation() result(e)
      e = max(e_co_photon_erg - photolysis_threshold_erg(ich_CO_C_O),     &
              0.0d0)
      end function heat_per_co_dissociation

      ! End of module
      end module co_photodissociation
