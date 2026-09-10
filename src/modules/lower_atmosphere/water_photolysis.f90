      module water_photolysis
      ! Photodissociation of H2O and OH by the stellar far-ultraviolet
      ! continuum, in the four bands of the A2 oxygen option:
      !
      !     H2O + hv  ->  OH + H          (O3)
      !               ->  H2 + O(1D)      (O4)
      !               ->  O  + H + H      (O5)
      !     OH  + hv  ->  O  + H          (O7)
      !
      ! O3 alone carries 54.5% of the H2O loss at the HD 189733 b base
      ! (docs/a2_oxygen_option_design.md sec. 2.2), which is what makes the
      ! OH/H2O cycle catalytic rather than a one-way sink; without it the
      ! oxygen chemistry cannot set the H2/H partition.  The reaction set,
      ! the quantum yields and the band edges are fixed by
      ! docs/a2_reaction_audit.md sec. 7 and are carried in oxygen_rates,
      ! which this module reads rather than restating.
      !
      ! ---------------------------------------------------------------
      ! 1. WHY BANDS AND NOT THE PHOTON ENERGY GRID
      !
      ! The code's energy grid is built from the ionization thresholds of
      ! the species the run carries: it starts at e_th_HI = 13.6 eV and is
      ! lowered only to the He 2^3S edge or to the lowest active low-IP
      ! metal threshold (set_energy_vectors, sed_read).  Water photolysis
      ! lives at 912-2304 A, i.e. 5.4-13.6 eV, almost entirely BELOW that
      ! grid; and for a power-law SED (is_PL_sed) the law is an XUV fit, so
      ! extending it into the FUV is not a spectrum.  The band form with a
      ! named incident flux is the precedent lyman_werner.f90 already sets
      ! for exactly this reason, and this module follows it.
      !
      ! The four bands are NOT free.  Their edges are fixed by the H2O
      ! branching ratios -- the three-body branch H2O + hv -> O + H + H
      ! opens only in the interval that contains Ly-alpha and is 12% there,
      ! so Ly-alpha needs its own band on branching grounds alone, quite
      ! apart from being a line rather than a continuum.
      !
      !   LW   912-1201 A   'Stellar LW flux'      (shared with H2)
      !   B2   Ly-alpha, 1215.67 A, from the existing 'Stellar Lya flux'
      !   B3  1231-1450 A   'Stellar FUV B3 flux'
      !   B4  1451-2304 A   'Stellar FUV B4 flux'
      !
      ! ONE INTERVAL, ONE FLUX, ONE BEAM.  Over 912-1201 A the same photons
      ! are absorbed by H2 in lines and by H2O and OH in a continuum.  That
      ! is physics; what is not physics is counting their ENERGY twice.  It
      ! would be counted twice in two independent ways, and both are closed
      ! here:
      !   * on the way IN, if the interval's flux were supplied both as
      !     'Stellar LW flux' and again inside a band reaching across it.
      !     Hence the first band IS the Lyman-Werner interval and carries
      !     that key's flux alone, and B2 starts where it ends.  The band
      !     was split at 1110 A between 2026-08 and 2026-09-06, with the
      !     1110-1201 A half carried by a separate 'Stellar FUV B1 flux'
      !     key; that key is retired and the parser stops a run that states
      !     it.  The split put the H2 line absorber astride the edge -- the
      !     self-shielding table's line list reaches 1201 A -- so the table
      !     rated 45 per cent more line absorptions than the 912-1110 A
      !     beam lost.  On the solar spectrum the HD 209458 b example uses,
      !     the 912-1201 A integral is 481.0 erg cm^-2 s^-1 and its two
      !     halves were 343.0 and 137.9.
      !   * on the way DOWN, if each absorber attenuated its own private
      !     copy of the beam.  Section 3 below makes the three absorbers
      !     share one beam, so the photons one removes are gone for the
      !     others and the sum of what they absorb is exactly what the beam
      !     loses.  Gate G4 of the design measures that identity, and
      !     output/FUV_bands.txt reports it band by band.
      !
      ! ---------------------------------------------------------------
      ! 2. CROSS SECTIONS, AND THE BAND AVERAGE
      !
      ! The photodissociation cross sections are the 'photodissociation'
      ! datasets of photochem/data/xsections/{H2O,OH}.h5, which Photochem
      ! concatenates from Huebner & Mukherjee (2015) below 6.3 nm, Heays,
      ! Bosman & van Dishoeck (2017), A&A 602, A105 (the Leiden database) to
      ! 192.056 nm and Ranjan et al. (2020) to 230.413 nm; those cut points
      ! are Photochem's, not any paper's.  Measured from the files: the H2O
      ! grid runs 0.1-230.413 nm and peaks at 2.75515e-17 cm^2 at 111.5 nm,
      ! the OH grid runs 0.06-264.9 nm and peaks at 1.36397e-17 cm^2 at
      ! 100 nm.
      !
      ! A band is entered with ONE number, the band-integrated energy flux
      ! at the planet, so the shape inside the band has to be assumed.  A
      ! flat F_lambda is assumed, as lyman_werner.f90 assumes for its own
      ! band, and the two moments that follow from it are tabulated below.
      ! With F_lambda constant the photon flux density is proportional to
      ! lambda, so
      !
      !   photon flux         N_b  = F_b / <hv>_b ,  <hv>_b = 2hc/(l1+l2)
      !   band cross section  s_b  = int s(l) l dl / int l dl
      !   mean photon energy  <E>_b = hc int s(l) dl / int s(l) l dl
      !
      ! and the dissociation rate is j_b = s_b N_b exp(-tau_b).
      !
      ! WHICH MEAN PHOTON ENERGY THE HEATING USES, AND WHY IT IS <hv>_b.
      ! The excess of the absorbed photon over the dissociation threshold
      ! goes to the fragments, so the heating needs the mean energy of the
      ! photons that are ABSORBED.  That mean is not one number: in an
      ! optically thin column the absorbed photons are drawn with weight
      ! s(l), so their mean is <E>_b; in a saturated column every photon of
      ! the band is absorbed, so their mean is the band's own <hv>_b.
      ! Measured, <E>_b/<hv>_b is 1.027 / 1.000 / 1.033 / 1.135 for
      ! H2O on LW..B4 and 0.970 / 1.000 / 1.014 / 1.182 for OH.
      !
      ! <hv>_b is adopted, for a reason that is a conservation law rather
      ! than a preference.  This model applies ONE transmission exp(-s_b N)
      ! to the whole band, so in the saturated limit it says that all N_b
      ! photons are absorbed; charging each of them <E>_b would then deposit
      ! N_b <E>_b = F_b <E>_b/<hv>_b, i.e. MORE ENERGY THAN THE BAND
      ! CARRIES -- 13.5% more on B4.  With <hv>_b the absorbed energy is
      ! F_b (1 - exp(-tau_b)), which cannot exceed the incident flux
      ! whatever the state of the gas; that bound is what gate G4 of the
      ! design asks of the band ledger, and output/FUV_bands.txt measures
      ! it.  The price is in the optically THIN limit, where the heating per
      ! dissociation is then wrong by the ratios above: low by 0 to 13.5%
      ! for H2O, and low for OH except on LW, where it is 3.0% high.  Both
      ! are far inside the band average's own spectral-shape uncertainty
      ! below, and <hv>_b is the accurate end of the two where it matters --
      ! the water layer this option exists for is optically thick in these
      ! bands.  <E>_b is kept and exported so the size of that error can be
      ! read off rather than recalled.
      !
      ! SPECTRAL-SHAPE SENSITIVITY, IN THE TERMS lyman_werner.f90 USES.
      ! s_b depends on the shape inside the band, and here that dependence
      ! is LARGE -- much larger than the +-25% lyman_werner.f90 quotes for
      ! its own band, because these cross sections fall by decades across a
      ! band while the Lyman-Werner lines do not.  Measured on the same
      ! files, replacing the flat F_lambda by a blackbody F_lambda changes
      ! s_b as follows.  The LW band is given as its two halves,
      ! 912-1110 A and
      ! 1110-1201 A, because the two move in opposite directions and a
      ! merged number would hide that:
      !
      !                    H2O                        OH
      !              3000K   5000K   9000K      3000K   5000K   9000K
      !   LW  short -60.0%  -43.4%  -22.9%     +89.3%  +70.6%  +41.0%
      !   LW  long   +1.0%   -2.1%   -1.9%     -23.4%  -14.4%   -6.8%
      !   B3        -67.3%  -45.8%  -21.7%     -17.2%  -15.0%   -8.3%
      !   B4        -98.5%  -85.4%  -42.7%     -99.2%  -90.5%  -52.5%
      !
      ! and, on a REAL stellar spectrum rather than a blackbody -- the
      ! HD 189733 b flux at the planet that the P1 Photochem reference arm
      ! was run with, weighted properly -- s_b/s_b(flat) is
      !
      !                H2O      OH
      !   LW  short   1.028   0.781
      !   LW  long    1.033   0.823
      !   B3          1.349   1.070
      !   B4          0.219   0.164
      !
      ! so both halves of LW are good to 3% for H2O and to 18-22% for OH on
      ! that star, and B3 to about 35%, while
      ! B4 IS OFF BY A FACTOR 4.6 TO 6.1.  B4's cross section falls three
      ! decades from 1451 A to 2304 A, so almost all of a cool star's B4
      ! energy sits where the molecule barely absorbs, and a single average
      ! cannot represent it.  A run that puts a large flux in B4 and finds
      ! it matters is outside this treatment and should say so.  B4 is
      ! carried because it is where the 1.00/0.00 branching applies and
      ! because leaving it out would silently drop photons, not because the
      ! average is good.
      !
      ! WHERE INSIDE THE LW BAND ITS SHAPE SENSITIVITY SITS.  The two LW
      ! rows above were measured while 1110 A was a band edge, and they are
      ! kept as halves rather than merged because they say where the
      ! sensitivity lives: nearly all of it is shortward of 1110 A (60% for
      ! H2O at 3000 K, and OH moves the OTHER WAY, +89%), while the long
      ! half is the best-behaved interval of the set for H2O (<= 2%).
      ! Averaged over the whole band the two trends partly cancel, so the
      ! single s_b below is less shape-sensitive than the short half and
      ! more so than the long one; a star strongly tilted across
      ! 912-1201 A cannot read its own error off one number, and the two
      ! half-band rows are what it should read instead.
      !
      ! THE FAILURE CASE IS A BAND DOMINATED BY ONE EMISSION LINE, and for
      ! an M dwarf, whose FUV is line-dominated (C II 1335, Si IV 1394/1403,
      ! C IV 1548/1551 all fall in B3 or B4), LW and B3 are exactly that
      ! case: a flat-continuum band average is then wrong by whatever the
      ! ratio of the line-position cross section to the band average happens
      ! to be.  That is the reason decision D1 of the design split the FUV
      ! into three keys rather than one -- it lets a line-dominated star be
      ! entered band by band -- and it is not a reason to trust the average
      ! WITHIN a band.  B2 has no averaging question of its own: the cross
      ! section is evaluated at 1215.67 A, and the band's flux really is the
      ! line -- on the HD 189733 b spectrum above, B2 carries 5300 erg
      ! cm^-2 s^-1 per nm against 43 for the 912-1201 A interval, so it is
      ! the line
      ! plus a negligible continuum.  That matters because s(H2O) varies by
      ! a factor 6 across the 28 A of B2 (5.57e-18 at 1202 A, 1.53e-17 at
      ! the line, 2.46e-18 at 1230 A), so the line-center value is right
      ! only as long as the band's energy is in the line.
      !
      ! ---------------------------------------------------------------
      ! 3. ATTENUATION
      !
      ! H2O and OH absorb the FUV as CONTINUA, so no self-shielding function
      ! is needed for them: the f_shield treatment of lyman_werner.f90
      ! exists because H2 absorbs in lines that saturate.  The star-ward
      ! attenuation is exp(-tau) with
      !
      !     tau_b = s_b(H2O) N(H2O) + s_b(OH) N(OH) ,
      !
      ! the columns built by calc_column_dens_one, the same radial
      ! integration and opa_pf weighting every other absorber column uses.
      ! The rate a cell sees is the MEAN over the cell, not the value at one
      ! of its faces: see water_photolysis_rate below, which is exact in tau
      ! at any grid spacing.  Applying exp(-tau) with the BAND-AVERAGED s_b
      ! is a separate approximation
      ! whose error grows with the spread of s(l) inside the band and with
      ! tau: the true band transmission is the s-weighted average of
      ! exp(-s(l) N), which exceeds exp(-s_b N) once the band is optically
      ! thick, because the transparent part of the band survives.  So the
      ! rates below are LOWER bounds deep in the water layer.
      !
      ! THE LW BAND HAS A THIRD ABSORBER: H2, AND THE BEAM IS SHARED.
      ! Over 912-1201 A the H2 Lyman and Werner lines take photons out of
      ! the same beam H2O and OH draw on.  Write the beam's transmission to
      ! the star-ward face of a cell as the product of the two,
      !
      !     T = (1 - A) exp(-tau_c) ,
      !
      ! with tau_c the H2O + OH continuum depth and A the fraction of the
      ! band the H2 lines have removed,
      !
      !     A(N_H2) = int_0^N sigma_pump(N') dN' ,
      !
      ! the column integral of the pump cross section of the SAME
      ! level-resolved table the dissociation rate comes from
      ! (h2_lw_band_photon_fraction_absorbed).  The beam therefore loses
      ! exactly the photons the rate spends: sigma_pump is sigma_diss/p_eff
      ! and p_eff is the dissociations per photon removed from the beam, so
      ! the two sides are one absorption with one normalization.  Then
      ! across
      ! one cell, with dA and dtau its own increments, the identity
      !
      !   T_out - T_in = e^{-tau_out} [ (1-A_out)(1 - e^{-dtau})
      !                                 + dA e^{-dtau} ]
      !
      ! is EXACT, and it splits the photons the cell removes into a
      ! continuum share and a line share with nothing left over.  The two
      ! shares are exactly what each absorber's rate must be, so the coupling
      ! costs one factor on each side and no approximation:
      !
      !   * H2O and OH in the LW band multiply their rate by (1 - A_out),
      !     the line transmission at the cell's star-ward face -- the
      !     tr_lines argument of water_photolysis_rate below;
      !   * H2 multiplies its rate by exp(-tau_c) at the cell's own depth,
      !     which is DB96 eq. (40)'s continuum term restored now that there
      !     is a continuum absorber to put in it (lyman_werner.f90).
      !
      ! Both factors are 1 for every band but LW and for every run without
      ! the other absorber, so nothing outside the shared interval changes.
      !
      ! WHERE THIS STOPS BEING TRUE.  A is the fraction of the band photons
      ! the lines take, so it cannot pass 1, and the table's own column
      ! integral respects that at the top of its column axis
      ! (h2_lw_band_photon_fraction_absorbed states the measured value).
      ! Above that column the table is clamped at its edge cross section,
      ! which goes on adding to A linearly, so A is capped at 1 there; a
      ! capped A means the lines have eaten the band and the continuum
      ! absorbers of the same interval see nothing, which is the physical
      ! limit rather than a defect.  The run reports the largest A it
      ! reached.
      !
      ! WHAT ELSE ABSORBS THESE BANDS, AND WHY IT IS NOT CARRIED.
      !  - Atomic H absorbs in its resonance LINES, never in its
      !    photoionization continuum: every band lies longward of the 912 A
      !    Lyman edge.  Band B2 IS one of those lines -- Ly-alpha, 1215.67 A
      !    -- and its absorber is the strongest in the atmosphere, so the
      !    transmission of the stellar Ly-alpha line through the atomic
      !    hydrogen column is applied to it as the tr_lines factor below
      !    (fuv_lw_photon_field in util_ion_eq.f90, the transmission from
      !    lya_rt.f90).  The higher Lyman-series lines inside the LW band
      !    (Ly-beta 1025.7, Ly-gamma 972.5, ...) are line absorbers in a
      !    continuum band and are not treated separately, exactly as
      !    lyman_werner.f90 leaves them; they can only reduce the rate.
      !  - Neutral low-IP metals photoionize in B3 and B4 (Mg I 1620 A,
      !    Fe I 1570 A, Ca I 2029 A, Na I 2412 A, K I 2857 A).  At solar
      !    abundance and s ~ 1e-18 cm^2 their optical depth is ~1e-21 N_H,
      !    i.e. 1e-3 at N_H = 1e18 cm^-2, where the H2O continuum is
      !    already thick.  Neglected, on the same measurement
      !    lyman_werner.f90 makes for its own band.
      !  - CO absorbs below 1118 A (in the LW band) in predissociating lines that
      !    self-shield.  It is not carried: CO is chemically frozen in this
      !    network (decision D4) so its photodissociation would have no
      !    consumer, and adding it as an absorber alone would remove
      !    photons from the water cycle without the compensating oxygen
      !    release.  A network that lets CO react must add both.
      !
      ! ---------------------------------------------------------------
      ! 4. HEATING
      !
      ! Each dissociation deposits hv - E_threshold as fragment kinetic
      ! energy, the precedent being e_lw_fragment_erg in lyman_werner.f90.
      ! The threshold energies are the 0 K dissociation energies formed
      ! from the Shomate table of oxygen_rates (photolysis_threshold_erg:
      ! the lowest interval's F coefficient is the 0 K enthalpy); the O(1D)
      ! channel adds the O(1D)
      ! excitation energy above O(3P) on top of the ground-state threshold,
      ! because that energy is carried away as electronic excitation and is
      ! returned to the gas later, through O6, not here.  The bond energy
      ! itself is paid by the photon and is NOT a thermal sink, exactly as
      ! for the Lyman-Werner channel.
      !
      ! Reference: Draine & Bertoldi (1996) is the precedent for the band
      ! form; the cross sections and yields are sourced in
      ! docs/a2_reaction_audit.md sections 4 and 7.

      use oxygen_rates, only: n_fuv_band, fuv_band_name,                  &
                              qy_H2O_OH_H, qy_H2O_H2_O1D, qy_H2O_O_H_H,   &
                              qy_OH_O_H,                                  &
                              photolysis_threshold_erg, e_excite_O1D_erg, &
                              ich_H2O_OH_H, ich_H2O_H2_O, ich_H2O_O_H_H,  &
                              ich_OH_O_H
      ! The Lyman-Werner band is shared with the H2 absorber, so its mean
      ! photon energy must be the SAME number lyman_werner.f90 converts its
      ! flux with: one beam cannot have two photon counts.
      use lyman_werner_photodissociation, only: e_lw_photon_erg
      ! (1 - exp(-d))/d, the fraction of the photons entering a cell that the
      ! cell absorbs, per unit of its own optical depth.  The XUV beam takes
      ! the same quantity for its cell means (utilities.f90,
      ! cell_mean_attenuation), and it is defined there once: one expression,
      ! one series branch and one switch point for both beams.
      use utils, only: absorbed_fraction_per_unit_depth

      implicit none
      private

      public :: ib_LW, ib_B2, ib_B3, ib_B4, n_fuv_band, fuv_band_name
      public :: sigma_H2O_band, sigma_OH_band
      public :: e_photon_H2O_band, e_photon_OH_band, e_photon_flat_band
      public :: water_photolysis_init
      public :: fuv_band_photon_flux, fuv_band_optical_depth
      public :: water_photolysis_rate, hydroxyl_photolysis_rate
      ! Re-exported so that a consumer of the band rates has the cell-mean
      ! factor they are built from without a second use statement, and so
      ! that there is no second definition of it to drift.
      public :: absorbed_fraction_per_unit_depth
      public :: heat_per_water_dissociation, heat_per_hydroxyl_dissociation

      integer, parameter :: dp = kind(1.0d0)

      ! Band identifiers, in the oxygen_rates band order.
      integer, parameter :: ib_LW = 1, ib_B2 = 2, ib_B3 = 3, ib_B4 = 4

      ! Band-averaged photodissociation cross sections [cm^2] over a flat
      ! F_lambda band, s_b = int s(l) l dl / int l dl, measured from the
      ! Photochem cross-section files named in the header (B2 is the value
      ! AT the Ly-alpha line, not a band average).
      real(dp), parameter :: sigma_H2O_band(n_fuv_band) =                 &
           (/ 8.51336791d-18, 1.52782477d-17,                             &
              4.40440310d-18, 1.24150238d-18 /)
      real(dp), parameter :: sigma_OH_band(n_fuv_band) =                  &
           (/ 4.74370570d-18, 4.57724155d-18,                             &
              1.40980849d-18, 6.99004135d-19 /)

      ! Cross-section-weighted mean photon energy of each band [erg],
      ! <E>_b = hc int s dl / int s l dl, i.e. the mean energy of the
      ! photons this species actually absorbs.  In eV:
      !   H2O  12.0526  10.1988   9.5511   7.4949
      !   OH   11.3829  10.1988   9.3766   7.8035
      real(dp), parameter :: e_photon_H2O_band(n_fuv_band) =              &
           (/ 1.93103584d-11, 1.63403379d-11,                            &
              1.53025281d-11, 1.20081719d-11 /)
      real(dp), parameter :: e_photon_OH_band(n_fuv_band) =               &
           (/ 1.82374916d-11, 1.63403379d-11,                            &
              1.50230218d-11, 1.25025409d-11 /)

      ! Mean photon energy of a flat F_lambda band, <hv>_b = 2hc/(l1+l2)
      ! [erg]; this is what converts the band ENERGY flux into a band PHOTON
      ! flux, and it is the same construction lyman_werner.f90 uses for its
      ! own band.  In eV: 11.7354, 10.1988, 9.2491, 6.6037.
      ! The LW entry is taken FROM lyman_werner.f90 rather than restated,
      ! because H2 and the two continuum absorbers share that beam and a
      ! shared beam has one photon count; the same integral measured here
      ! over 912-1201 A gives 1.88021378d-11, i.e. the two agree to the six
      ! digits that module carries.
      real(dp), parameter :: e_photon_flat_band(n_fuv_band) =             &
           (/ e_lw_photon_erg, 1.63403379d-11,                            &
              1.48186935d-11, 1.05802709d-11 /)

      ! Yield-weighted threshold energy of the H2O photolysis in each band,
      ! and the single OH threshold [erg per dissociation].  Filled once by
      ! water_photolysis_init (they come from the Shomate table, so they are
      ! not compile-time constants); read-only afterwards, which is what
      ! makes them safe to read from the OpenMP cell sweep.
      real(dp), save :: eth_H2O_band(n_fuv_band) = 0.0d0
      real(dp), save :: eth_OH_erg = 0.0d0
      logical,  save :: wp_ready = .false.

      contains

      ! Fill the threshold energies from the thermodynamic table.  Call once,
      ! from serial setup (input_read), before any cell sweep reads them.
      !
      ! The O(1D) channel's threshold is the ground-state H2 + O(3P)
      ! enthalpy PLUS the O(1D) excitation energy: that 1.96 eV leaves the
      ! photon as electronic excitation of the oxygen atom, not as fragment
      ! kinetic energy, so it is not deposited here.  It returns to the gas
      ! through O6 (O(1D) + H2 -> OH + H), whose exothermicity the network
      ! does not currently deposit either -- an omission recorded at the
      ! call site in ionization_equilibrium, not hidden here.
      subroutine water_photolysis_init
      integer :: ib
      do ib = 1, n_fuv_band
        eth_H2O_band(ib) =                                                &
             qy_H2O_OH_H(ib)  *photolysis_threshold_erg(ich_H2O_OH_H)     &
           + qy_H2O_H2_O1D(ib)*(photolysis_threshold_erg(ich_H2O_H2_O)    &
                                + e_excite_O1D_erg)                       &
           + qy_H2O_O_H_H(ib) *photolysis_threshold_erg(ich_H2O_O_H_H)
      end do
      eth_OH_erg = qy_OH_O_H(1)*photolysis_threshold_erg(ich_OH_O_H)
      wp_ready   = .true.
      end subroutine water_photolysis_init

      ! Photon flux of band ib [photons cm^-2 s^-1] for a band-integrated
      ! energy flux F_band [erg cm^-2 s^-1] at the planet.
      double precision function fuv_band_photon_flux(F_band, ib) result(Nph)
      real(dp), intent(in) :: F_band
      integer,  intent(in) :: ib
      if (F_band .le. 0.0d0) then
        Nph = 0.0d0
      else
        Nph = F_band/e_photon_flat_band(ib)
      end if
      end function fuv_band_photon_flux

      ! Star-ward continuum optical depth of band ib for H2O and OH columns
      ! N_H2O, N_OH [cm^-2].
      double precision function fuv_band_optical_depth(ib, N_H2O, N_OH)   &
                                result(tau)
      integer,  intent(in) :: ib
      real(dp), intent(in) :: N_H2O, N_OH
      tau = sigma_H2O_band(ib)*max(N_H2O, 0.0d0)                          &
          + sigma_OH_band(ib) *max(N_OH,  0.0d0)
      end function fuv_band_optical_depth

      ! Mean H2O photodissociation rate over ONE CELL in band ib [s^-1].
      ! Multiply by the quantum yields of oxygen_rates to split it into
      ! O3 / O4 / O5.
      !
      ! tau_out is the optical depth at the cell's star-ward face and
      ! dtau its own optical depth, so the photons that enter the cell are
      ! N_b tr_lines exp(-tau_out) and the ones the continuum absorbers
      ! take are N_b tr_lines exp(-tau_out) (1 - exp(-dtau)).  Those are
      ! shared between the two continuum absorbers in proportion to their
      ! contributions to dtau, which leaves this species' rate as
      !
      !     j = s N_b tr_lines exp(-tau_out) (1 - exp(-dtau))/dtau .
      !
      ! tr_lines is the transmission of the LINE absorber of the band down
      ! to the same face.  In the LW band it is the H2 Lyman-Werner lines,
      ! 1 - A (section 3), which is what keeps that band's photons from
      ! being absorbed twice; in the Ly-alpha band it is the fraction of the
      ! stellar line that penetrates the atomic hydrogen column to that
      ! face.  It is 1 in the three bands with no line absorber and in a run
      ! that carries neither H2 nor atomic H.
      !
      ! WHY NOT s N_b exp(-tau).  Evaluating the rate at one point of the
      ! cell and multiplying by the cell's absorbers is a rectangle rule in
      ! tau: it is exact only while dtau << 1, and it under-counts the
      ! absorption everywhere else, one-signed.  Measured on an HD 189733 b
      ! run before this form was used, the column sum of the absorbed
      ! photons fell 30% below the closed-form N_b(1 - exp(-tau)) in the
      ! Ly-alpha band, whose optical depth crosses unity inside a single
      ! cell.  The form above reduces to the rectangle rule as dtau -> 0
      ! and reproduces the closed form EXACTLY, cell by cell and summed
      ! over the column, whatever the grid -- which is what gate G4 of
      ! docs/a2_oxygen_option_design.md asks of the band ledger, and what
      ! output/FUV_bands.txt measures.
      double precision function water_photolysis_rate(F_band, ib, tau_out, &
                                dtau, tr_lines) result(j)
      real(dp), intent(in) :: F_band, tau_out, dtau, tr_lines
      integer,  intent(in) :: ib
      if (F_band .le. 0.0d0) then
        j = 0.0d0
        return
      end if
      j = sigma_H2O_band(ib)*fuv_band_photon_flux(F_band, ib)             &
          *max(tr_lines, 0.0d0)                                           &
          *exp(-max(tau_out, 0.0d0))                                      &
          *absorbed_fraction_per_unit_depth(max(dtau, 0.0d0))
      end function water_photolysis_rate

      ! The same for OH (O7, one merged channel of unit yield; see the
      ! caveat in oxygen_rates on what that merge is).
      double precision function hydroxyl_photolysis_rate(F_band, ib,      &
                                tau_out, dtau, tr_lines) result(j)
      real(dp), intent(in) :: F_band, tau_out, dtau, tr_lines
      integer,  intent(in) :: ib
      if (F_band .le. 0.0d0) then
        j = 0.0d0
        return
      end if
      j = sigma_OH_band(ib)*fuv_band_photon_flux(F_band, ib)              &
          *max(tr_lines, 0.0d0)                                           &
          *exp(-max(tau_out, 0.0d0))                                      &
          *absorbed_fraction_per_unit_depth(max(dtau, 0.0d0))
      end function hydroxyl_photolysis_rate

      ! Kinetic energy given to the fragments by one H2O dissociation in
      ! band ib [erg], averaged over the three channels with their yields.
      ! The mean absorbed photon is <hv>_b, the band's own flat-F_lambda
      ! mean; section 2 gives the conservation argument for that choice and
      ! the size of what it costs.  Non-negative by construction: a band
      ! whose mean photon sits below the yield-weighted threshold deposits
      ! nothing rather than cooling the gas (it would mean the band average
      ! has left the range where it represents the band).
      double precision function heat_per_water_dissociation(ib) result(e)
      integer, intent(in) :: ib
      ! THE THRESHOLD TABLE IS THE GUARD wp_ready NAMES.  Until
      ! water_photolysis_init has filled it, eth_H2O_band is zero and the
      ! expression below would return the WHOLE photon energy as fragment
      ! kinetic energy, i.e. deposit the bond energy a second time on top of
      ! the photon that already paid for it.  The setup fills the table
      ! exactly when the oxygen chemistry is on, and without that chemistry
      ! there is no H2O in the gas to dissociate, so the deposit there is
      ! zero.
      if (.not. wp_ready) then
        e = 0.0d0
        return
      end if
      e = max(e_photon_flat_band(ib) - eth_H2O_band(ib), 0.0d0)
      end function heat_per_water_dissociation

      ! The same for one OH dissociation in band ib [erg].
      double precision function heat_per_hydroxyl_dissociation(ib) result(e)
      integer, intent(in) :: ib
      ! Same guard as heat_per_water_dissociation: no threshold table means
      ! no OH in the gas either, and returning the whole photon energy would
      ! deposit the bond energy twice.
      if (.not. wp_ready) then
        e = 0.0d0
        return
      end if
      e = max(e_photon_flat_band(ib) - eth_OH_erg, 0.0d0)
      end function heat_per_hydroxyl_dissociation

      ! End of module
      end module water_photolysis
