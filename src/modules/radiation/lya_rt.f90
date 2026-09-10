   module lya_rt
   ! Ly-alpha radiative transfer by the escape-probability method,
   ! computed IN-LINE (every timestep) from the current state.
   !
   ! Because the line is closed with an analytic escape probability rather than a
   ! Monte Carlo, the field is cheap and is evaluated directly inside the run
   ! (like the H(n=2) feedback) instead of an offline pass. The result
   ! is the Voigt-profile-averaged Ly-alpha mean intensity J_lya(r) (Huang et al.
   ! 2023 Eq. 11), which excited_hydrogen uses to pump H(n=2) when jlya_mode = 2.
   ! The separate option of reading an externally-computed (e.g. real Monte Carlo)
   ! J_lya(r) profile is kept as jlya_mode = 1 (load_jlya_rt in excited_hydrogen).
   !
   ! Why an escape probability, not Monte Carlo: the WASP-121b column is extremely
   ! thick at line center (tau0 = 1.2e8 and a.tau0 = 1.2e5 at the base of the
   ! stored wasp_full regression output), so an acceleration-free MC is
   ! intractable and core-skipping is excluded.
   !
   ! METHOD (plane-parallel, substellar column):
   !   Dnu_D(r) = nu_lya sqrt(2 kT/mu)/c                       Doppler width
   !   a(r)     = (A_2p1s/4pi)/Dnu_D                           Voigt parameter
   !   tau(r)   = top-down rectangle sum of nhi*C_lya/Dnu_D    line-center depth
   !              (per cell, so a cell's own depth dtau and its two face
   !               depths are exact; see lya_line_center_optical_depth)
   !   beta(r)  = 1/(1 + <N>)                                  escape probability
   !              of a line photon EMITTED at the cell centre, with <N> the
   !              mean number of scatterings before escape of the published
   !              static plane-parallel damping-wing solution below.  tau_c
   !              = the CELL-CENTRE depth, since beta is a property of the
   !              material at a depth and not a beam to be averaged over a
   !              cell.
   !              WHY 1/(1 + <N>) AND NOT 1/<N>.  <N> counts the photons
   !              ABSORBED over the medium per photon generated (Harrington
   !              1973, MNRAS 162, 43, eq. 38), so a generated photon is
   !              EMITTED 1 + <N> times before it leaves.  The steady 2p
   !              budget of the closure below is n2p A beta = P, i.e. beta
   !              is the escape chance per emission, hence 1/(1 + <N>).
   !              This is exact at both ends -- beta -> 1 as tau -> 0 and
   !              beta -> 1/<N> in the damping-wing regime -- so no join
   !              between the thin and thick limits has to be chosen.
   !   <N>      = (4 sqrt(6)/pi^2) [(tau_up + tau_dn)/2] Phi(xi),
   !              xi = tau_dn/(tau_up + tau_dn),
   !              Phi(xi) = sum_{n odd} sin(n pi xi)/n^2, Phi(1/2) = u_2 =
   !              1 - 3^-2 + 5^-2 - ... = 0.9159656 (Catalan's constant).
   !              This is Neufeld (1990), ApJ 350, 216, eq. (3.27) with the
   !              continuum destruction opacity beta_N = 0 (his W[0,0] = 1):
   !              a source plane in a slab whose two faces lie at line-centre
   !              depths tau_up and tau_dn.  With the source in the mid-plane
   !              (tau_dn = tau_up = B) it reduces to Harrington (1973) eq.
   !              (40), <N> = (4 sqrt(6)/pi^2) u_2 B = 0.909316 B (a uniform
   !              source, his eq. 48, would give 0.664736 B).
   !              Redistribution R_II-A (Hummer 1962): Ly-alpha shifts by
   !              about one Doppler width per scattering and cannot be
   !              re-emitted directly in the far wing, so escape is a random
   !              walk in frequency and <N> is proportional to the optical
   !              depth and INDEPENDENT of the Voigt parameter a.  a leaves
   !              the escape probability entirely; it survives only in the
   !              stellar-beam penetration below, which is a one-flight
   !              transmission and a different question.
   !              WHAT WAS HERE BEFORE (removed 2026-09-07, item LYA-BETA):
   !              pi^(-1/4) sqrt(a/tau_c), the ONE-FLIGHT chance that a
   !              photon drawn afresh from the whole Voigt profile lands
   !              beyond the frequency at which the wing depth falls to
   !              unity.  That is the complete-redistribution value, an
   !              escape route partial redistribution does not give the
   !              line; measured on the stored wasp_full column it exceeded
   !              the published 1/<N> by 234x at the base and 14x at
   !              r = 1.20 R_p.
   !              VALIDITY.  The damping-wing solution requires
   !              (a tau)^(1/3) > 10 (Neufeld section Va).  Cells that fail
   !              it are counted by lya_wing_domain_record: there the closure
   !              is outside its published domain, and what carries it is the
   !              thin limit beta -> 1, which 1/(1 + <N>) reaches smoothly.
   !   THE TWO FACES.  tau_up = tau_c is the depth to the star-ward face.
   !              tau_dn is the depth to the planet-ward face: the rest of the
   !              column when the lower boundary is declared a pure absorber
   !              ("Lya absorbing bottom: True"), and tau_c itself for the
   !              default reflecting boundary, where the emitting plane is a
   !              symmetry plane and the mirror image is a slab of half
   !              thickness tau_c with a mid-plane source.  Huang et al.
   !              (2017) close their Monte Carlo domain with a pure absorber
   !              because the H2 layer beneath the base truly absorbs
   !              Ly-alpha through accidental resonances; the reflecting
   !              bottom keeps every downward photon.  The two faces are NOT
   !              two independent escape channels to be multiplied: they
   !              share ONE escaping population, which the slab solution
   !              above already counts once, and Neufeld eq. (2.25) splits
   !              between them in the ratio (tau_0 -+ tau_s)/(2 tau_0), i.e.
   !              tau_dn/(tau_up + tau_dn) leaves star-ward
   !              (lya_starward_escape_fraction).
   !   P(r)     = alpha_B ne nhii + C_1s2p ne nhi              Ly-alpha production
   !              (recombination cascade, 1 per case-B recomb; electron-impact
   !               1s->2p -- the two internal sources used by Huang)
   !   D(r)     = ne q_2p1s + Gamma_2 + ne C_2p2s P_2gamma     2p destruction
   !              (the channels that remove the atom from 2p without putting a
   !               photon into the line; hydrogen_n2_rates.f90)
   !   J_int(r) = (2 h nu^3/c^2)(g1s/g2p) P (1 - beta) / ((A_2p1s beta + D) nhi)
   !              (trapped source-function buildup: n2p = P/(A beta + D) and
   !               Jbar = S (1 - beta), the standard two-level escape-probability
   !               relation; Eq. 11 is the S <-> n2p half of it, and the
   !               destruction channels are retained in the denominator)
   !   J_star(r)= xi F_Lya_star/(4 pi Dnu_star) <T_s E>        stellar beam (dayside
   !              dilution xi as for ground-state EUV; Dnu_star is the BROAD
   !              stellar line width, T_s = erfc(x1/(sqrt(2) Xs)) the fraction of
   !              the stellar profile whose wings reach depth tau, and
   !              E = 1 + (lya_star_boost-1)(1-beta) T_s the bounded trapping
   !              buildup of the scattered beam.  <.> is the MEAN of the beam
   !              over the cell's own optical depth, not its value at a face:
   !              the pumping rate is what the cell's atoms undergo)
   !   J_lya = J_int + J_star.

   use global_parameters
   use hydrogen_n2_rates, only: lA_lya, nu_lya, A_2p1s, g1s, g2p, C_lya,     &
                                c1s2p_rate, alpha_B_hydrogen,                &
                                n2p_destruction_rate

   implicit none

   ! Coefficient of the static plane-parallel damping-wing solution for the
   ! mean number of scatterings before escape: 4 sqrt(6)/pi^2 (Neufeld 1990
   ! eq. 3.27 at zero destruction; Harrington 1973 eq. 40), and Catalan's
   ! constant u_2 = 1 - 3^-2 + 5^-2 - ... = Phi(1/2), the source-position
   ! factor of a mid-plane source, which Harrington gives below his eq. (40)
   ! as 0.9159656.  Their product 0.9093164 is his 0.909316.
   real*8, parameter :: slab_scatterings_c = 0.9927408002342284d0
   real*8, parameter :: catalan_u2         = 0.9159655941772190d0

   ! Domain record of that solution (Neufeld section Va: (a tau)^(1/3) > 10).
   ! Cell visits of jlya_escape_prob outside it, and the worst (smallest)
   ! (a tau)^(1/3) seen, read through lya_wing_domain_record.  Informational:
   ! outside the domain the escape probability is carried by the thin limit
   ! beta -> 1, which the closure reaches smoothly.
   real*8, parameter :: slab_wing_domain_min = 10.0d0
   integer :: lya_wing_cells_out  = 0
   integer :: lya_wing_cells_seen = 0
   real*8  :: lya_wing_worst      = huge(1.0d0)

   ! Diagnostic split of J_lya into internal (recomb+collisional) and stellar
   ! contributions, filled each call for output/Excited_H.txt.
   real*8, dimension(:), allocatable :: jint_arr, jstar_arr

   contains

   ! --------------------------------------------------------------- !

   subroutine lya_rt_allocate_arrays
   ! Allocate the grid-sized module arrays once the number of cells N is
   ! known; called from EXHALE_main right after input_read. The zeros are
   ! the initializers the declarations used to carry.

   allocate(jint_arr(1-Ng:N+Ng), jstar_arr(1-Ng:N+Ng))

   jint_arr  = 0.0d0
   jstar_arr = 0.0d0

   end subroutine lya_rt_allocate_arrays

   ! --------------------------------------------------------------- !

   subroutine lya_line_center_optical_depth(T_K, nhi, tau, dtau, tau_out)
   ! Ly-alpha line-centre optical depth of the stellar beam, cell by cell.
   ! T_K [K] and nhi [cm^-3] are physical cell-by-cell arrays.
   !
   !   tau(j)     depth from the star down to the INNER (planet-ward) face
   !              of cell j,
   !   dtau(j)    cell j's OWN line-centre depth (optional),
   !   tau_out(j) depth at cell j's star-ward face (optional), so that
   !              tau_out(j) + dtau(j) = tau(j) and tau_out(N+Ng) = 0.
   !
   ! DISCRETIZATION.  The rectangle rule of calc_column_dens, on the cell's
   ! own line-centre opacity nhi*C_lya/Dnu_D and its own width dr_j: the
   ! finite-volume state IS piecewise constant over a cell, so this is the
   ! exact depth of that state, and one cell's depth is exactly
   ! tau(j) - tau(j+1).  Every other column of the code (N1, N15, N2, NTR,
   ! NH2col, the H2O and OH columns) is accumulated by the same rule with
   ! the same widths, which is what lets the line depths and the continuum
   ! depths of the same cell be paired face by face
   ! (fuv_lw_photon_field in util_ion_eq.f90 pairs them).
   !
   ! WHICH CONSUMER READS WHICH.  A rate is what the cell's atoms undergo
   ! averaged over the cell, so the stellar BEAM is taken as the mean of its
   ! transmission across [tau_out, tau_out + dtau]
   ! (lya_stellar_beam_transmission_cell_mean below).  A point quantity is
   ! read at the cell centre, tau(j) - dtau(j)/2: the escape probability
   ! beta is one, since it is the probability that a photon created AT a
   ! depth leaves the line from that depth, so it is evaluated at the depth
   ! of the material in question and never averaged over a cell.
   real*8, dimension(1-Ng:N+Ng), intent(in)  :: T_K, nhi
   real*8, dimension(1-Ng:N+Ng), intent(out) :: tau
   real*8, dimension(1-Ng:N+Ng), intent(out), optional :: dtau, tau_out

   integer :: j
   real*8, dimension(1-Ng:N+Ng) :: DnuD, dt_cell

   do j = 1-Ng, N+Ng
      DnuD(j)    = nu_lya*sqrt(2.0d0*kb_erg*max(T_K(j),1.0d0)/mu)/c_light
      dt_cell(j) = max(nhi(j),0.0d0)*C_lya/max(DnuD(j),1.0d-30)             &
                 *dr_j(j)*R0
   enddo

   tau(N+Ng) = dt_cell(N+Ng)
   do j = N+Ng-1, 1-Ng, -1
      tau(j) = tau(j+1) + dt_cell(j)
   enddo

   if (present(dtau))    dtau = dt_cell
   if (present(tau_out)) tau_out = tau - dt_cell

   end subroutine lya_line_center_optical_depth

   ! --------------------------------------------------------------- !

   elemental double precision function                                      &
                  lya_slab_source_position_factor(xi_in) result(phi_s)
   ! Phi(xi) = sum_{n odd} sin(n pi xi)/n^2, the dependence of the trapping
   ! of a static plane-parallel slab on WHERE between its two faces the line
   ! photons are created: xi is the fractional position of the source plane,
   ! so xi = 1/2 is the mid-plane and xi -> 0 or 1 a source lying on a face.
   ! It is the sum of Neufeld (1990), ApJ 350, 216, eq. (3.27) at zero
   ! continuum destruction, where the two sines
   ! sin{n pi [1 -+ (tau_s/tau_0)]/2} add to 2 sin(n pi xi) for odd n and
   ! cancel for even n.  Phi(1/2) = 1 - 3^-2 + 5^-2 - ... = u_2, Catalan's
   ! constant, which returns Harrington (1973), MNRAS 162, 43, eq. (40).
   ! Phi(1 - xi) = Phi(xi): the slab does not know which face is which.
   !
   ! EVALUATION.  Term by term the series converges only as 1/n^2, so it is
   ! summed once and for all: Phi'(xi) = (pi/2) ln cot(pi xi/2) from
   ! sum_{n odd} cos(n theta)/n = ln cot(theta/2)/2, and integrating with
   ! the logarithmic end behaviour taken analytically,
   !   Phi(xi) = (pi/2) { xi [ln(2/(pi xi)) + 1]
   !                      + int_0^xi ln[(pi t/2) cot(pi t/2)] dt } .
   ! The remaining integrand is even and analytic on |t| < 1 and vanishes at
   ! t = 0, and the reduction to xi <= 1/2 by the symmetry keeps the path
   ! away from its singularity at t = 1, so ten-point Gauss-Legendre is at
   ! round-off (MEASURED 4e-14 relative at the worst point, xi = 1/2,
   ! against the series; src/tests/physics_probe/lya_escape_probability.f90
   ! carries the check).  xi = 1/2 is the reflecting-boundary default and is
   ! returned from the closed form directly.
   real*8, intent(in) :: xi_in

   integer, parameter :: n_gl = 10
   ! Ten-point Gauss-Legendre on [0,1]: the nodes (1 + x_i)/2 and weights
   ! w_i/2 of the standard rule on [-1,1].
   real*8, parameter :: gl_s(n_gl) = (/ 0.01304673574141414d0,              &
                                        0.06746831665550773d0,              &
                                        0.16029521585048778d0,              &
                                        0.28330230293537639d0,              &
                                        0.42556283050918439d0,              &
                                        0.57443716949081561d0,              &
                                        0.71669769706462361d0,              &
                                        0.83970478414951222d0,              &
                                        0.93253168334449227d0,              &
                                        0.98695326425858586d0 /)
   real*8, parameter :: gl_w(n_gl) = (/ 0.03333567215434407d0,              &
                                        0.07472567457529030d0,              &
                                        0.10954318125799102d0,              &
                                        0.13463335965499826d0,              &
                                        0.14776211235737644d0,              &
                                        0.14776211235737644d0,              &
                                        0.13463335965499826d0,              &
                                        0.10954318125799102d0,              &
                                        0.07472567457529030d0,              &
                                        0.03333567215434407d0 /)
   real*8  :: xi, t, x, acc
   integer :: i

   xi = min(max(xi_in, 0.0d0), 1.0d0)
   xi = min(xi, 1.0d0 - xi)
   if (xi .le. 0.0d0) then
      phi_s = 0.0d0
      return
   endif
   if (xi .eq. 0.5d0) then
      phi_s = catalan_u2
      return
   endif

   acc = 0.0d0
   do i = 1, n_gl
      t   = xi*gl_s(i)
      x   = 0.5d0*pi*t
      acc = acc + xi*gl_w(i)*log(x/tan(x))
   enddo
   phi_s = 0.5d0*pi*(xi*(log(2.0d0/(pi*xi)) + 1.0d0) + acc)

   end function lya_slab_source_position_factor

   ! --------------------------------------------------------------- !

   elemental double precision function                                      &
                  lya_slab_scatterings_before_escape(tau_up, tau_dn)        &
                  result(Nbar)
   ! Mean number of scatterings a Ly-alpha photon undergoes before it leaves
   ! a static plane-parallel slab whose star-ward face lies at line-centre
   ! optical depth tau_up from the emitting plane and whose planet-ward face
   ! lies at tau_dn:
   !
   !   <N> = (4 sqrt(6)/pi^2) [(tau_up + tau_dn)/2] Phi(tau_dn/(tau_up+tau_dn))
   !
   ! Neufeld (1990) eq. (3.27) with his continuum destruction opacity set to
   ! zero, reducing for tau_dn = tau_up to Harrington (1973) eq. (40),
   ! <N> = 0.909316 tau_up.  Damping wings, redistribution R_II-A: <N> is
   ! proportional to the optical depth and independent of the Voigt
   ! parameter.  Valid for (a tau)^(1/3) > 10 (Neufeld section Va); below
   ! that the slab is no longer in the damping-wing regime and only the
   ! thin limit <N> -> 0 that this expression also gives is meaningful.
   real*8, intent(in) :: tau_up, tau_dn

   real*8 :: t_up, t_dn, t_tot

   t_up  = max(tau_up, 0.0d0)
   t_dn  = max(tau_dn, 0.0d0)
   t_tot = t_up + t_dn
   if (t_tot .le. 0.0d0) then
      Nbar = 0.0d0
      return
   endif
   Nbar = slab_scatterings_c*0.5d0*t_tot                                   &
        *lya_slab_source_position_factor(t_dn/t_tot)

   end function lya_slab_scatterings_before_escape

   ! --------------------------------------------------------------- !

   elemental double precision function                                      &
                  lya_wing_escape_probability(tau_up, tau_dn) result(beta)
   ! Probability that a Ly-alpha photon EMITTED at the source plane of that
   ! slab leaves the line rather than being absorbed again:
   !
   !   beta = 1/(1 + <N>) .
   !
   ! <N> counts absorptions (Harrington 1973 eq. 38), so a photon generated
   ! once is emitted 1 + <N> times before it goes; the two-level budget the
   ! closure solves, n2p A beta = production, therefore needs the escape
   ! chance per EMISSION.  beta -> 1 as the slab thins and -> 1/<N> in the
   ! damping-wing regime, with no interpolation between the two.
   real*8, intent(in) :: tau_up, tau_dn

   beta = 1.0d0/(1.0d0 + lya_slab_scatterings_before_escape(tau_up, tau_dn))

   end function lya_wing_escape_probability

   ! --------------------------------------------------------------- !

   elemental double precision function                                      &
                  lya_starward_escape_fraction(tau_up, tau_dn) result(f_up)
   ! Share of the escaping photons of that slab that leave through the
   ! STAR-WARD face, (tau_0 - tau_s)/(2 tau_0) = tau_dn/(tau_up + tau_dn)
   ! (Neufeld 1990 eq. 2.25, the gambler's-ruin split of a random walk
   ! between two absorbing boundaries).  The two faces share ONE escaping
   ! population -- Neufeld's own remark under that equation is that the two
   ! fractions sum to one when there are no loss processes -- so they are a
   ! partition of the escape probability above and never two independent
   ! probabilities to be multiplied.  A source plane sitting on a face loses
   ! everything through the other one.
   real*8, intent(in) :: tau_up, tau_dn

   real*8 :: t_up, t_dn, t_tot

   t_up  = max(tau_up, 0.0d0)
   t_dn  = max(tau_dn, 0.0d0)
   t_tot = t_up + t_dn
   if (t_tot .le. 0.0d0) then
      f_up = 0.5d0
      return
   endif
   f_up = t_dn/t_tot

   end function lya_starward_escape_fraction

   ! --------------------------------------------------------------- !

   subroutine lya_wing_domain_reset
   ! Zero the domain record of the damping-wing solution.

   lya_wing_cells_out  = 0
   lya_wing_cells_seen = 0
   lya_wing_worst      = huge(1.0d0)

   end subroutine lya_wing_domain_reset

   ! --------------------------------------------------------------- !

   subroutine lya_wing_domain_record(ncells_out, ncells_seen, worst,        &
                                     threshold)
   ! The record accumulated by jlya_escape_prob: cell visits whose
   ! (a tau)^(1/3) fell below the stated limit of the damping-wing solution
   ! (Neufeld section Va), how many were seen, and the worst (smallest)
   ! (a tau)^(1/3) of the run.  Informational: outside the domain the
   ! closure is carried by its thin limit, beta -> 1.
   integer, intent(out) :: ncells_out, ncells_seen
   real*8,  intent(out) :: worst
   real*8,  intent(out), optional :: threshold

   ncells_out  = lya_wing_cells_out
   ncells_seen = lya_wing_cells_seen
   worst       = lya_wing_worst
   if (lya_wing_cells_seen .eq. 0) worst = 0.0d0
   if (present(threshold)) threshold = slab_wing_domain_min

   end subroutine lya_wing_domain_record

   ! --------------------------------------------------------------- !

   ! Fraction of the incident stellar Ly-alpha line that reaches line-centre
   ! optical depth tau through gas at temperature T_K [K].
   !
   ! Ly-alpha is resonantly SCATTERED, not destroyed, so the beam is NOT
   ! attenuated by exp(-tau): the broad stellar line penetrates through its
   ! own wings.  A stellar photon at Doppler offset x from line centre
   ! reaches depth tau once the wing optical depth there has fallen to
   ! unity, i.e. |x| > x1 = sqrt(a tau/sqrt(pi)) with a the Voigt parameter
   ! of the local gas.  The stellar profile is a Gaussian of Doppler width
   ! Xs = dv_star_lya/v_th in those units, so the surviving fraction of it
   ! is erfc(x1/(sqrt(2) Xs)).
   !
   ! Validity: the wing regime a tau >> 1.  It tends to erfc(0) = 1 as
   ! tau -> 0, so the optically thin limit is the free-streaming beam, and
   ! to 0 as tau -> infinity.
   !
   ! NOT VERIFIED AGAINST A PUBLISHED SOLUTION (item REF-NEUFELD,
   ! 2026-09-07).  x1 above is the depth at which a photon crosses the
   ! column in ONE flight; it is not the escape frequency of the published
   ! static-slab solutions, where a photon that has random-walked leaves at
   ! x_m = 0.88119 (a B)^(1/3) for a mid-plane source (Harrington 1973,
   ! MNRAS 162, 43, eq. 37), 0.71025 (a B)^(1/3) for a uniform source (his
   ! eq. 47), or x_s = 0.525 (a tau_0)^(1/3) in Neufeld (1990), ApJ 350,
   ! 216, eq. 4.32.  Those scale as (a tau)^(1/3) and x1 as (a tau)^(1/2):
   ! at the base of the stored wasp_full column x1 = 258 against
   ! 0.88119 (a tau)^(1/3) = 43, so this beam is required to reach 6 times
   ! further into the wing than a scattered photon would.  Whether the
   ! one-flight condition or the diffusive one is right for an INCIDENT
   ! beam is exactly the question Neufeld's illuminated slab answers
   ! (his section IIc and Fig. 5): he gives the transmitted fraction of a
   ! non-absorbing illuminated slab in closed form (his eq. 2.34, the same
   ! result Ambartsumian 1958 obtained for coherent scattering), and that
   ! is what this expression should be compared against.
   !
   ! It is NOT corrected to (a tau)^(1/3) (item LYA-BETA, 2026-09-07).  The
   ! escape frequency of the published slab solutions is where a photon that
   ! has ALREADY random-walked in frequency leaves; an incident stellar
   ! photon has not, so the frequency at which it crosses the column in one
   ! flight is the right question for this beam and the wrong one for the
   ! escape probability, which is why the Voigt parameter left the latter
   ! (see the header) and stays here.  What is open is the fate of the
   ! stellar photons INSIDE x1, which scatter rather than stop: they are
   ! neither transmitted by this expression nor counted anywhere else.
   !
   ! ONE definition of that transmission: the cell mean below takes the
   ! mean OF IT over a cell for the stellar term of jlya_escape_prob (the
   ! n = 2 pumping field), and the Ly-alpha photolysis band of the
   ! molecular layer (fuv_lw_photon_field in utils_ion_eq) reads it at a
   ! face.  The same stellar beam reaching the same depth cannot carry two
   ! transmissions.
   elemental double precision function                                      &
                  lya_stellar_beam_transmission(T_K, tau) result(Tstar)
   real*8, intent(in) :: T_K, tau

   Tstar = erfc(lya_wing_penetration_coefficient(T_K)                       &
                *sqrt(max(tau,0.0d0)))

   end function lya_stellar_beam_transmission

   ! --------------------------------------------------------------- !

   ! The erfc argument of that transmission per sqrt(line-centre depth):
   ! x1/(sqrt(2) Xs) = c sqrt(tau) with
   !   c = sqrt(a/sqrt(pi))/(sqrt(2) Xs) ,
   ! a the Voigt parameter of the gas and Xs = dv_star_lya/v_th the stellar
   ! line width in Doppler units.  ONE definition of the depth dependence,
   ! read by the transmission itself and by the quadrature that takes its
   ! mean over a cell, which is what makes the mean a mean OF the production
   ! transmission.
   elemental double precision function                                      &
                  lya_wing_penetration_coefficient(T_K) result(c)
   real*8, intent(in) :: T_K

   real*8 :: Tl, DnuD, avoigt, vth, Xs

   Tl     = max(T_K, 1.0d0)
   DnuD   = nu_lya*sqrt(2.0d0*kb_erg*Tl/mu)/c_light      ! Doppler width [Hz]
   avoigt = (A_2p1s/(4.0d0*pi))/max(DnuD,1.0d-30)        ! Voigt parameter
   vth    = sqrt(2.0d0*kb_erg*Tl/mu)                     ! thermal speed [cm/s]
   Xs     = (dv_star_lya*1.0d5)/max(vth,1.0d-30)         ! stellar width [Doppler]
   c      = sqrt(avoigt/sqrt(pi))/(sqrt(2.0d0)*max(Xs,1.0d-30))

   end function lya_wing_penetration_coefficient

   ! --------------------------------------------------------------- !

   elemental subroutine lya_stellar_beam_transmission_cell_mean            &
                           (T_K, tau_out, dtau, tbar, tbar2)
   ! MEAN over one cell of the stellar Ly-alpha beam transmission, and of
   ! its square, for a cell whose star-ward face sits at line-centre depth
   ! tau_out and whose own depth is dtau:
   !
   !   tbar  = (1/dtau) int_{tau_out}^{tau_out+dtau} T_s(tau) dtau ,
   !   tbar2 = (1/dtau) int                          T_s(tau)^2 dtau .
   !
   ! WHY A MEAN AND NOT A FACE VALUE.  The pumping rate of a cell is what
   ! its atoms undergo, averaged over the cell, and the beam falls across
   ! the cell by exactly that cell's own optical depth.  A face value is
   ! one-sided at every cell in the same direction, and the sum over the
   ! column of (cell field) x (cell depth) then misses the depth integral of
   ! the beam, which is the number of scatterings the column performs.  This
   ! is the same discretization the XUV photoionization sums take
   ! (cell_mean_attenuation in utilities.f90) and the FUV photolysis rates
   ! take (water_photolysis.f90).
   !
   ! WHY A QUADRATURE.  The Ly-alpha beam is resonantly SCATTERED, so its
   ! transmission is not exponential in the column: it is
   ! T_s = erfc(c sqrt(tau)), whose mean has no closed form.  With the
   ! substitution tau = u^2 the integrand becomes 2 u erfc(c u), entire in
   ! u, which removes the square-root behaviour T_s ~ 1 - (2c/sqrt(pi))
   ! sqrt(tau) that a rule in tau itself would have to resolve at a cell
   ! whose star-ward face is at zero depth.  Three-point Gauss-Legendre in u
   ! is exact through fifth order per segment, with a remainder
   ! h^6 f^(6)/2016000; the sixth derivative of erfc(x) carries the Hermite
   ! factor (2x)^6 of exp(-x^2), so the segments are cut at
   ! h = 0.25/max(x_in,1) in the erfc argument x = c u, which holds the
   ! relative remainder at (2 x h)^6/2016000 = 8e-9 wherever the beam is
   ! not already extinguished.  The driver
   ! src/tests/physics_probe/lya_beam_cell_mean.f90 measures what is left
   ! against a fine reference.
   !
   ! tbar2 exists because the trapping buildup E of the penetrating beam is
   ! itself linear in T_s, so the field the cell sees, T_s E, is quadratic
   ! in it: the exact cell mean of T_s E at fixed escape probability is
   ! tbar + (lya_star_boost - 1)(1 - beta) tbar2.
   real*8, intent(in)  :: T_K, tau_out, dtau
   real*8, intent(out) :: tbar, tbar2

   ! Three-point Gauss-Legendre on [0,1]: nodes (1 -+ sqrt(3/5))/2 and 1/2,
   ! weights 5/18, 8/18, 5/18.
   integer, parameter :: n_gl = 3
   real*8, parameter :: gl_s(n_gl) = (/ 0.1127016653792583d0,               &
                                        0.5d0,                              &
                                        0.8872983346207417d0 /)
   real*8, parameter :: gl_w(n_gl) = (/ 5.0d0/18.0d0, 8.0d0/18.0d0,         &
                                        5.0d0/18.0d0 /)
   real*8, parameter  :: seg_arg   = 0.25d0
   integer, parameter :: nseg_max  = 256
   real*8  :: c, t_lo, d, u0, u1, du, u, e, w, acc, acc2
   integer :: nseg, i, g

   t_lo = max(tau_out, 0.0d0)
   d    = max(dtau,    0.0d0)
   c    = lya_wing_penetration_coefficient(T_K)

   if (d .le. 0.0d0) then
      tbar  = erfc(c*sqrt(t_lo))
      tbar2 = tbar*tbar
      return
   endif

   u0   = sqrt(t_lo)
   u1   = sqrt(t_lo + d)
   nseg = min(max(int(c*(u1 - u0)*max(c*u1,1.0d0)/seg_arg) + 1, 1),        &
              nseg_max)
   du   = (u1 - u0)/dble(nseg)
   acc  = 0.0d0
   acc2 = 0.0d0
   do i = 0,nseg-1
      do g = 1,n_gl
         u = u0 + (dble(i) + gl_s(g))*du
         e = erfc(c*u)
         ! dtau = 2 u du, so the weights sum to u1^2 - u0^2 = dtau exactly.
         w = gl_w(g)*du*2.0d0*u
         acc  = acc  + w*e
         acc2 = acc2 + w*e*e
      enddo
   enddo
   tbar  = acc /d
   tbar2 = acc2/d

   end subroutine lya_stellar_beam_transmission_cell_mean

   ! --------------------------------------------------------------- !

   elemental double precision function                                      &
                  lya_photosphere_attenuation_cell_mean(tau_out, dtau)      &
                  result(f)
   ! MEAN over one cell of the parameterized Ly-alpha photosphere factor
   ! 1/(1 + tau) of Huang et al. (2017) Eq. (6), which the jlya_mode = 0
   ! field is attenuated by:
   !
   !   (1/dtau) int_{tau_out}^{tau_out+dtau} dtau/(1 + tau)
   !     = ln[(1 + tau_out + dtau)/(1 + tau_out)]/dtau ,
   !
   ! exact for the piecewise-constant absorber the finite-volume state is.
   !
   ! The logarithm is taken of 1 + x with x = dtau/(1 + tau_out), which is
   ! below 1e-7 in most of the column: forming 1 + x in double precision
   ! there throws away the low bits of x and leaves the result relatively
   ! wrong by eps/x, 2e-9 at x = 5e-8.  Kahan's identity
   ! ln(1 + x) = x ln(y)/(y - 1) with y = fl(1 + x) restores those bits
   ! from the SAME rounded y, and is exact to a few ulp at every x, with no
   ! threshold to choose.
   real*8, intent(in) :: tau_out, dtau

   real*8 :: t_lo, d, x, y

   t_lo = max(tau_out, 0.0d0)
   d    = max(dtau,    0.0d0)
   x    = d/(1.0d0 + t_lo)
   y    = 1.0d0 + x
   if (y .eq. 1.0d0) then
      f = 1.0d0/(1.0d0 + t_lo)
   else
      f = x*log(y)/(y - 1.0d0)/d
   endif

   end function lya_photosphere_attenuation_cell_mean

   ! --------------------------------------------------------------- !

   subroutine jlya_escape_prob(T_K, nhi, nhii, ne, v_in, Jlya, tau)
   ! Voigt-averaged Ly-alpha mean intensity J_lya(r) [erg s^-1 cm^-2 Hz^-1 sr^-1]
   ! by the escape-probability method, and the line-centre optical depth
   ! tau(r) at each cell's inner face.
   ! T_K/nhi/nhii/ne are physical (cgs) cell-by-cell arrays; v_in is the dimensionless
   ! radial velocity (v_in*v0 = cm/s). Called in-line from excited_H_update.
   real*8, dimension(1-Ng:N+Ng), intent(in)  :: T_K, nhi, nhii, ne, v_in
   real*8, dimension(1-Ng:N+Ng), intent(out) :: Jlya, tau

   integer :: j, jm, jp
   real*8, dimension(1-Ng:N+Ng) :: DnuD, dtau, tau_out
   real*8 :: avoigt, beta_esc, Tl, aB, C1s2p, Prec, Pcol, Jint, Jstar, Jpref, xi
   real*8 :: Tstar, Tstar2, tau_c, Dnu_star, D2p
   real*8 :: C_sob, dvdr, tau_sob, beta_sob, beta_tot
   real*8 :: tau_dn, atau13

   ! Doppler width (needed for both tau and the source function).
   do j = 1-Ng, N+Ng
      DnuD(j) = nu_lya*sqrt(2.0d0*kb_erg*max(T_K(j),1.0d0)/mu)/c_light
   enddo

   ! Line-centre depth at each cell's inner face, that cell's own depth, and
   ! the depth at its star-ward face.
   call lya_line_center_optical_depth(T_K, nhi, tau, dtau, tau_out)

   ! Dayside dilution for the stellar beam: global_parameters'
   ! dayside_dilution(), the same factor every other stellar band uses.
   xi = dayside_dilution()

   Jpref = 2.0d0*hp_erg*nu_lya**3.0d0/c_light**2.0d0     ! 2 h nu^3 / c^2
   ! Stellar Ly-alpha line width [Hz]: the BROAD stellar profile sets the
   ! profile-averaged mean intensity (Jbar = J_nu(nu0) for a stellar field flat
   ! across the narrow planetary line), NOT the planetary Doppler width Dnu_D.
   ! Dividing by Dnu_D would make J_star rise spuriously as the wind cools.
   Dnu_star = nu_lya*(dv_star_lya*1.0d5)/c_light
   ! Sobolev (velocity-gradient) escape coefficient: tau_S = C_sob n1s/|dv/dr|,
   ! C_sob = (pi e^2/m_e c) f_lya lambda = sqrt(pi) C_lya lambda.
   C_sob = sqrt(pi)*C_lya*lA_lya
   do j = 1-Ng, N+Ng
      avoigt   = (A_2p1s/(4.0d0*pi))/max(DnuD(j),1.0d-30)
      ! Depth of the cell's CENTRE.  beta is the chance that a photon
      ! created at a depth random-walks out through the wings from there:
      ! a point quantity of the material, evaluated at the depth of the
      ! material, and not a beam to be averaged over the cell.
      tau_c    = max(tau(j) - 0.5d0*dtau(j), 0.0d0)
      ! Depth to the planet-ward face of the slab the photon is trapped in.
      ! Pure absorber: the rest of the column, tau(1-Ng) being the depth at
      ! the inner face of the innermost cell (max() only guards roundoff).
      ! Reflecting default: the emitting plane is a symmetry plane, and the
      ! mirror image is a slab of half thickness tau_c with a mid-plane
      ! source, i.e. Harrington's B = tau_c.
      if (lya_bottom_absorber) then
         tau_dn = max(tau(1-Ng) - tau_c, 0.0d0)
      else
         tau_dn = tau_c
      endif
      ! ONE escaping population for the two faces (see the header): the
      ! static damping-wing slab solution, not a product of two channels.
      beta_esc = lya_wing_escape_probability(tau_c, tau_dn)
      ! Domain record of that solution, (a tau)^(1/3) > 10 (Neufeld section
      ! Va), on the full slab depth the photon random-walks across.
      atau13 = (avoigt*(tau_c + tau_dn))**(1.0d0/3.0d0)
      lya_wing_cells_seen = lya_wing_cells_seen + 1
      if (atau13 .lt. slab_wing_domain_min)                                &
         lya_wing_cells_out = lya_wing_cells_out + 1
      lya_wing_worst = min(lya_wing_worst, atau13)

      ! Sobolev escape via the wind velocity gradient, a channel parallel to
      ! the static wing escape: the photon Doppler-shifts out of resonance
      ! over the Sobolev length. |dv/dr| in cgs = (v0/R0)|d v_in/d r|.
      ! (1 - exp(-tau_S))/tau_S is the Sobolev (1960) escape probability on
      ! the RADIAL ray; the angle average over the direction dependence of
      ! the Sobolev depth in a spherical flow is not taken. Neufeld (1990)
      ! is a static solution and supplies nothing for this channel: his
      ! title says static, and his Appendix C states the assumption as
      ! "a real system with no velocity gradients".
      jm = max(j-1, 1-Ng); jp = min(j+1, N+Ng)
      dvdr     = (v0/R0)*abs(v_in(jp) - v_in(jm))/max(abs(r(jp) - r(jm)),1.0d-30)
      tau_sob  = C_sob*max(nhi(j),0.0d0)/max(dvdr,1.0d-30)
      beta_sob = min(1.0d0,(1.0d0 - exp(-min(tau_sob,200.0d0)))/max(tau_sob,1.0d-30))
      ! Total escape: out of the static slab OR out of resonance through the
      ! velocity gradient.  Written as the union of two independent chances
      ! rather than as 1 - (1-beta_esc)(1-beta_sob): the latter forms 1 - x
      ! for x of order 1e-8 at the base and recovers it by subtraction, which
      ! throws away eight digits of a probability the 2p budget divides by.
      beta_tot = beta_esc + beta_sob*(1.0d0 - beta_esc)
      ! beta_tot below feeds the (1-beta_tot) trapping factor of Jint, the
      ! A_2p1s*beta_tot term in its denominator, and the buildup factor E of
      ! Jstar, so a declared absorbing bottom removes photons consistently
      ! from the internal field and from the trapped stellar beam -- no
      ! further change is needed.

      Tl    = max(T_K(j), 1.0d0)
      aB    = alpha_B_hydrogen(Tl)                             ! case-B [cm^3/s]
      C1s2p = c1s2p_rate(Tl)                                   ! 1s->2p [cm^3/s]
      Prec  = aB*ne(j)*nhii(j)
      Pcol  = C1s2p*ne(j)*nhi(j)

      ! Trapped internal field, escape-probability closure Jbar = S(1-beta):
      ! S = (2hv3/c2)(g1s/g2p) n2p/n1s with the trapped pile-up
      ! n2p = P/(A beta + D). D collects the channels that empty 2p without
      ! returning a photon to the line -- collisional de-excitation, n=2
      ! photoionization, and l-mixing followed by two-photon decay. Where the
      ! line is thick enough that A beta drops to D (the base), leaving D out
      ! over-estimates n2p and hence the pumping field.
      ! The (1-beta) factor makes Jbar -> S in the thick limit and -> 0 when thin
      ! (beta -> 1), instead of the source function S diverging as n1s -> 0 in the
      ! ionized outer wind (which spuriously raised Jbar outward).
      D2p   = n2p_destruction_rate(Tl, ne(j), gamma2_bal, gamma2_bal)
      Jint  = Jpref*(g1s/g2p)*(Prec + Pcol)*(1.0d0 - beta_tot)             &
            /(max(A_2p1s*beta_tot + D2p, 1.0d-30)*max(nhi(j),1.0d-30))

      ! Stellar beam: resonantly SCATTERED, not destroyed -- so it is not killed by
      ! exp(-tau); instead the broad stellar line penetrates through its wings and
      ! FILLS the region it reaches with the free-streaming mean intensity
      ! F/(4 pi Dnu_D). A stellar photon at Doppler offset x reaches depth r if
      ! |x| > x1 = sqrt(a tau/sqrt(pi)) (planetary wing optical depth 1); the
      ! fraction of a Gaussian stellar profile of Doppler width Xs = dv_star/v_th
      ! beyond x1 is T_star = erfc(x1/(sqrt(2) Xs)) -- smooth, so the penetration
      ! photosphere is not a sharp edge (a flat profile gives a kink). (No 1/beta
      ! buildup: that applies to a VOLUMETRIC source, not an incident beam.)
      ! The stellar term still dominates the internal sources where it penetrates.
      ! lya_stellar_beam_transmission above holds that expression, and
      ! lya_stellar_beam_transmission_cell_mean its mean over a cell; the
      ! FUV photolysis band of the molecular layer reads the same one.
      ! The MEAN of the beam over this cell's own optical depth, and the
      ! mean of its square, so that the quadratic buildup below is the exact
      ! cell mean of T_s E at this cell's escape probability.
      call lya_stellar_beam_transmission_cell_mean(T_K(j), tau_out(j),      &
                                                   dtau(j), Tstar, Tstar2)
      ! Bounded trapping buildup of the penetrating stellar beam: the scattered
      ! photons raise the mean intensity above free-streaming by
      !   E = 1 + (lya_star_boost - 1)*(1 - beta)*T_star,
      ! which is ~lya_star_boost where the beam fully penetrates and is trapped
      ! (T_star~1, beta<<1, the wind) and -> 1 both in the thin outer wind
      ! (beta -> 1) and at the deep penetration edge (T_star -> 0), so it does not
      ! over-count the base. (The full 1/beta blows up; a flat cap over-enhances
      ! the deep edge.) lya_star_boost tuned to Huang+2023 Fig. 11.
      Jstar = xi*F_Lya_star/(4.0d0*pi*Dnu_star)                              &
            * (Tstar + (lya_star_boost - 1.0d0)*(1.0d0 - beta_tot)*Tstar2)

      jint_arr(j)  = Jint
      jstar_arr(j) = Jstar
      Jlya(j) = Jint + Jstar
   enddo

   end subroutine jlya_escape_prob

   ! --------------------------------------------------------------- !

   end module lya_rt
