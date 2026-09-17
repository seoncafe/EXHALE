      module Numerical_Fluxes
      ! Collection of numerical flux functions
      
      use global_parameters
      use Conversion
      ! S_estimate_HLLC (Toro's adaptive p_star estimate) used to be `use`d
      ! here and never called: the HLLC branch below builds its own
      ! Davis/Einfeldt speeds inline.  The module has been removed.
      use S_estimate_ROE
      use caloric_eos, only: energy_density_from_pressure,             &
                            adiabatic_index_from_state
  
      implicit none

      ! Faces at which the Roe flux was not used because the estimated star
      ! state does not exist (two rarefactions separating fast enough leave
      ! vacuum) or is not admissible; those faces take the HLLE flux, which
      ! is positively conservative (Einfeldt et al. 1991, J. Comput. Phys.
      ! 92, 273).  Summed over a whole run and reported in the run summary.
      ! Only "Numerical flux: ROE" can raise it.
      integer :: n_faces_roe_hlle = 0

      ! Faces at which the Lax-Friedrichs flux was the flux of the scheme,
      ! that is the faces evaluated under "Numerical flux: LLF", summed over
      ! a whole run and reported in the run summary.  The same routine is
      ! also called by the positivity repair of RK_integration; those faces
      ! are counted there (n_faces_flux_positivity_limited) and not here, so
      ! that this counter reports the discretization the run selected and not
      ! the repairs it needed.
      integer :: n_faces_llf = 0

      ! ---- low-Mach scaling of the normal velocity jump -------------- !
      ! "Low Mach velocity jump: True|False" (default False), acting on the
      ! ROE branch only.  The normal velocity jump that enters the acoustic
      ! expansion coefficients of the Roe dissipation is replaced by
      !     dvel -> min(Ma_Roe, 1) dvel,   Ma_Roe = |U_Roe| / a_Roe
      ! (Rieper 2011, J. Comput. Phys. 230, 5263, his eq. 3.15 with the
      ! local Mach number 3.16, Ma = (|U|+|V|)/a from the Roe-averaged
      ! velocities; in one dimension the tangential component V is absent,
      ! so the factor is |U_Roe|/a_Roe).  Without it the upwind velocity
      ! jump carries an artificial viscosity of the momentum that is one
      ! order in the Mach number LARGER than the momentum update it damps
      ! (his eq. 2.14), and the discrete solution does not follow the
      ! incompressible limit of the continuous equations; with it the
      ! viscosity is returned to the order of the update (his eq. 3.17).
      ! The eigenvalues, the eigenvectors, the central flux, the
      ! well-balanced pressure departure and the admissibility test are
      ! untouched, and for Ma_Roe >= 1 the factor is exactly 1, so a
      ! supersonic face is the unmodified Roe flux.
      !
      ! Validity.  The asymptotic sorting of his section 3.2 assumes the
      ! local Mach number 3.16 and the global Mach number used as the
      ! expansion parameter do not deviate substantially; where they do,
      ! the factor is applied outside the regime it was derived in.  It is
      ! an ACCURACY correction, not a stiffness one: his section 5 states
      ! that the stiffness of the equations is not removed, explicit steps
      ! keeping dt = O(Ma) and implicit solves remaining slowly convergent.
      ! The same section says HLL, Rusanov and van Leer flux splitting are
      ! not suited to the correction (too much diffusion on the shear and
      ! entropy waves) and that an HLLC adaptation is possible but needs
      ! its own derivation, which is why this option exists on the ROE
      ! branch alone.
      !
      ! Default .false.: with it off every arithmetic operation of the
      ! branch is the one the goldens were taken with.
      logical :: low_mach_velocity_jump = .false.

      contains
      
      ! Subroutine for the numerical flux.
      !
      ! jL and jR name the cells the two face states were reconstructed
      ! from.  Each side's total energy and sound speed are then built with
      ! ITS OWN heat capacity, which is what the caloric EOS makes different
      ! across a face inside the H2 front.
      !
      ! The one place a SINGLE gamma is still needed is the Roe average:
      ! a_avg = sqrt((gamma-1)(H_avg - v_avg^2/2)) and the p_star estimate of
      ! speed_estimate_ROE are both derived for one ideal gas on both sides
      ! (a two-gamma Roe average needs the Vinokur-Montagne/Glaister
      ! extension).  Feeding those two gammas would produce a wave speed
      ! with no bound property, so they are given the arithmetic mean of the
      ! two.  This is a WAVE-SPEED ESTIMATE, not a definition of the steady
      ! state: the converged solution is the root of a residual built from
      ! the conservative fluxes, and those carry each side's own EOS above.
      ! Where the two sides agree -- every face of an atomic gas -- the mean
      ! is that value exactly, (a + a)/2 = a being exact in binary
      ! floating point, so nothing moves off the constant-gamma arithmetic.
      ! The HLLC branch needs no such average: its Davis/Einfeldt speeds
      ! min(v-c) / max(v+c) already take one sound speed per side.
      ! WELL-BALANCED OPTION ("Well balanced:").  With the optional arguments
      ! present the PRESSURE JUMP of the Riemann problem is not formed as
      ! pR - pL, the difference of two O(1) face pressures, but as
      !     dp = dp_eq + dev_R - dev_L
      ! from the mismatch of the two neighboring hydrostatic equilibria at
      ! this face and the two face states' departures from their own
      ! equilibrium, all three of the size of the departure.  On a discrete
      ! hydrostatic equilibrium the three vanish together, the dissipation of
      ! every row vanishes with them, and the face is the stationary contact
      ! discontinuity the well-balanced property requires the flux to resolve
      ! exactly (Kaeppeli and Mishra 2016, A&A 587, A94, their eq. 10; ROE and
      ! HLLC have that property, the Lax-Friedrichs flux does not).
      ! The two optional outputs are the face pressure measured from each
      ! side's own equilibrium, q_up = p_out - P_up(jL) and
      ! q_dn = p_out - P_dn(jR), which is how RK_rhs forms the pressure force
      ! without ever subtracting two O(1) numbers.
      subroutine Num_flux(WL,WR,NF,p_out,jL,jR,                          &
                          dev_L,dev_R,dp_eq,q_up,q_dn)

      real*8,intent(in) :: WL(3),WR(3)
      integer,intent(in) :: jL,jR
      real*8, intent(in),  optional :: dev_L,dev_R,dp_eq
      real*8, intent(out), optional :: q_up,q_dn
      logical :: wb, hllc_left
      real*8  :: dp_wb
      real*8 :: gam_L,gam_R,gam_face
      real*8 :: uL(3),uR(3)
      real*8 :: FL(3),FR(3)
      real*8 :: rhoL,vL,pL,aL,EL,HL,SL,phiL
      real*8 :: rhoR,vR,pR,aR,ER,HR,SR,phiR
      real*8 :: S_star
      real*8 :: usL(3),usR(3)
      real*8 :: s,rho_avg,v_avg,H_avg,a_avg          
      real*8 :: l1,l2,l3
      real*8 :: l1L,l1R,l3L,l3R 
      real*8 :: v_star,aL_star,aR_star
      type(roe_star_state) :: star
      integer :: roe_status
      real*8 :: drho,dvel,dp
      real*8 :: a1,a2,a3   
      real*8, dimension(3) :: K1,K2,K3
      real*8, intent(out) :: NF(3),p_out   
      
      wb = present(dev_L) .and. present(dev_R) .and. present(dp_eq)
      hllc_left = .true.
      dp_wb = 0.0d0
      if (wb) dp_wb = dp_eq + dev_R - dev_L

      ! Exctract left state
      rhoL = WL(1)
      vL   = WL(2)
      pL   = WL(3)
      EL   = 0.5*rhoL*vL*vL + energy_density_from_pressure(jL,rhoL,pL)
      HL   = (EL+pL)/rhoL
      gam_L = adiabatic_index_from_state(jL,rhoL,pL)
      aL   = sqrt(gam_L*pL/rhoL)   
      
      ! Exctract right state
      rhoR = WR(1)
      vR   = WR(2)
      pR   = WR(3)
      ER   = 0.5*rhoR*vR*vR + energy_density_from_pressure(jR,rhoR,pR)
      HR   = (ER+pR)/rhoR
      gam_R = adiabatic_index_from_state(jR,rhoR,pR)
      aR   = sqrt(gam_R*pR/rhoR)
      gam_face = 0.5*(gam_L + gam_R)
            
      ! Evaluate numerical flux      
      select case(flux)
      
      case ('LLF') ! Local Lax Friedrichs

         call lax_friedrichs_flux(WL,WR,NF,p_out,jL,jR)

         ! p_out is the average of the two face pressures
         if (wb) then
            q_up = 0.5d0*(dev_L + dev_R + dp_eq)
            q_dn = q_up - dp_eq
         endif

         !$omp atomic
         n_faces_llf = n_faces_llf + 1

      !----------------------------------------------!
      
      case('HLLC') ! HLLC solver
            
         call W_to_U_comp(WL,uL,jL)
         call W_to_U_comp(WR,uR,jR)
   
         ! Speed estimates
         SL = min(0.0,min(vL-aL, vR-aR))
         SR = max(0.0,max(vL+aL, vR+aR))
         phiL = rhoL*(SL-vL)
         phiR = rhoR*(SR-vR)
         if (wb) then
            S_star = (dp_wb + phiL*vL - phiR*vR)/(phiL - phiR)
         else
            S_star = (pR - pL + phiL*vL - phiR*vR)/(phiL - phiR)
         endif
         
         ! HLLC state approximation
         
         ! Left state
         usL(1) = 1.0
         usL(2) = S_star
         usL(3) = EL/rhoL + (S_star-vL)*(S_star + pL/phiL)            
         usL = rhoL*(SL-vL)/(SL-S_star)*usL

         ! Right state
         usR(1) = 1.0
         usR(2) = S_star
         usR(3) = ER/rhoR + (S_star-vR)*(S_star + pR/phiR)           
         usR = rhoR*(SR-vR)/(SR-S_star)*usR
         
         ! Evaluate HLLC flux 
         
         if(SL.ge.(0.0)) then
               
            call Phys_flux(WL,NF,jL)
            p_out = pL
            hllc_left = .true.
               
         elseif(SL.lt.(0.0).and.S_star.ge.(0.0)) then
         
            call Phys_flux(WL,NF,jL)
            NF = NF + SL * (usL-uL)
            p_out = pL
            hllc_left = .true.
               
         elseif(S_star.lt.(0.0).and.SR.ge.(0.0)) then
         
            call Phys_flux(WR,NF,jR)
            NF = NF + SR * (usR-uR)
            p_out = pR
            hllc_left = .false.
         else
            call Phys_flux(WR,NF,jR)
            p_out = pR
            hllc_left = .false.
               
         endif

         ! p_out is one side's own face pressure here, so the departure it
         ! carries is that side's.  hllc_left says which side the branch
         ! above took.
         if (wb) then
            if (hllc_left) then
               q_up = dev_L
               q_dn = dev_L - dp_eq
            else
               q_dn = dev_R
               q_up = dev_R + dp_eq
            endif
         endif

      !----------------------------------------------!
      
      case ('ROE')
      
         ! Construct averaged states
         
         ! Auxiliary parameter
         s = sqrt(rhoL)/(sqrt(rhoL)+sqrt(rhoR))
         
         ! Density
         rho_avg = sqrt(rhoL*rhoR)
         
         ! Velocity
         v_avg = s*vL + (1.0-s)*vR
         
         ! Entalpy
         H_avg = s*HL + (1.0-s)*HR
         
         ! Sound speed
         a_avg = sqrt((gam_face-1.0)*(H_avg-0.5*v_avg*v_avg))
      
      	!---------------------------------!
      	
      	! Average eigenvalues
      
      	l1 = v_avg - a_avg
      	l2 = v_avg
      	l3 = v_avg + a_avg
      	
      	!---------------------------------!
      
         ! Star state of the face Riemann problem, Toro (2009) chapter 9.
         ! It supplies the entropy fix below, and it says whether the
         ! linearized Roe flux is entitled to this face at all.
         call speed_estimate_ROE(WL,WR,gam_face,star,roe_status)

         if (roe_status .ne. ROE_STAR_OK) then

            ! No admissible star state: the flow separates into vacuum, or
            ! no branch of the estimate returned a positive star pressure
            ! and positive star densities.  The Roe flux is a linearization
            ! about the average state and is not positivity preserving
            ! there: one conservative update with it leaves a negative
            ! thermal energy for a vacuum-producing expansion.  The face
            ! takes the HLLE flux instead, which is positively conservative
            ! whenever its two speeds bound the physical waves (Einfeldt
            ! 1988, SIAM J. Numer. Anal. 25, 294; Einfeldt, Munz, Roe &
            ! Sjogreen 1991, J. Comput. Phys. 92, 273).  Each bound pairs
            ! one side's own signal speed with the Roe average, so the
            ! bounds enclose the waves of the linearized problem as well.
            SL = min(vL - aL, v_avg - a_avg)
            SR = max(vR + aR, v_avg + a_avg)

            call Phys_flux(WL,FL,jL)
            call Phys_flux(WR,FR,jR)
            call hlle_flux(WL,WR,FL,FR,SL,SR,jL,jR,NF)

            !$omp atomic
            n_faces_roe_hlle = n_faces_roe_hlle + 1

         else

            ! Entropy correction

            ! Approximate velocities of the star state
            v_star  = star%u_star
            aL_star = star%c_L_star
            aR_star = star%c_R_star

            ! Intermediate eigenvalues
            l1L = vL - aL
            l1R = v_star - aL_star
            l3L = v_star + aR_star
            l3R = vR + aR

            ! Modify the eigenvalues if a rarefaction is present

            ! Left rarefaction
            if (l1L.lt.(0.0).and.l1R.gt.(0.0)) then
                  l1 = l1L*(l1R-l1)/(l1R-l1L)
            endif

            ! Right rarefaction
            if (l3L.lt.(0.0).and.l3R.gt.(0.0)) then
                  l3 = l3R*(l3-l3L)/(l3R-l3L)
            endif

            !---------------------------------!

            ! Averaged eigenevectors

            ! K1
            K1(1) = 1.0
            K1(2) = v_avg - a_avg
            K1(3) = H_avg - v_avg*a_avg

            ! K2
            K2(1) = 1.0
            K2(2) = v_avg
            K2(3) = 0.5*v_avg*v_avg

            ! K3
            K3(1) = 1.0
            K3(2) = v_avg + a_avg
            K3(3) = H_avg + v_avg*a_avg

            !---------------------------------!

            ! Evaluate the conserved-variables differences.  Under the
            ! well-balanced option the pressure jump is the one built from the
            ! departures (see the header), not the difference of the two
            ! O(1) face pressures.
            drho = rhoR - rhoL
            dvel = vR - vL
            dp   = pR - PL
            if (wb) dp = dp_wb

            ! Low-Mach scaling of the normal velocity jump, Rieper (2011)
            ! eq. 3.15 with the local Mach number 3.16 of the Roe averages
            ! (one dimension: no tangential component, so Ma = |U|/a).  It
            ! enters a1 and a3, the two acoustic expansion coefficients, and
            ! nothing else: a2 carries drho and dp only.  See the module
            ! header for the validity of the correction.
            if (low_mach_velocity_jump)                                  &
               dvel = min(abs(v_avg)/a_avg, 1.0d0)*dvel

            ! Evaluate the expansion coefficients
            a1 = 0.5/a_avg**2.0*(dp - rho_avg*a_avg*dvel)
            a2 = drho - dp/a_avg**2.0
            a3 = 0.5/a_avg**2.0*(dp + rho_avg*a_avg*dvel)

            !---------------------------------!

            ! Evaluate the left and right fluxes
            call Phys_flux(WL,FL,jL)
            call Phys_flux(WR,FR,jR)

            ! Evaluate the flux at interface ( eq.[11.29] Toro )
            NF = 0.5*(FR+FL)   &
               - 0.5*(a1*abs(l1)*K1 + a2*abs(l2)*K2 + a3*abs(l3)*K3)

         endif

         ! Output pressure
         p_out = 0.5*(pR + pL)

         ! p_out is the average of the two face pressures, on the Roe branch
         ! and on the HLLE fallback alike
         if (wb) then
            q_up = 0.5d0*(dev_L + dev_R + dp_eq)
            q_dn = q_up - dp_eq
         endif

      case default

         write(*,*) 'ERROR: unknown numerical flux: ', trim(flux)
         write(*,*) '  allowed: LLF, HLLC, ROE'
         error stop 1

      end select

      ! End of subroutine
      end subroutine Num_flux

      !-----------------------------------------------------------!

      ! Local Lax-Friedrichs (Rusanov) flux: the arithmetic average of the
      ! two physical fluxes plus a jump term carrying the largest signal
      ! speed of the pair.
      !
      ! The viscosity coefficient is the spectral radius of the flux
      ! Jacobian on either side.  The eigenvalues of the Euler system are
      ! v-a, v and v+a, so that radius is |v| + a and the coefficient is
      !     alpha = max(|v_L| + a_L, |v_R| + a_R)
      ! (Rusanov 1961, J. Comput. Math. Phys. USSR 1, 267; it is also the
      ! HLLE flux of the symmetric speed pair -alpha, +alpha).  The
      ! bound is what makes the flux monotone for the scalar problem and
      ! what the positivity lemma below is proved under; it is also
      ! symmetric under (rho,v,p) -> (rho,-v,p) with the sides exchanged,
      ! which is an exact symmetry of the Euler equations, so the numerical
      ! flux inherits it.
      !
      ! It is the numerical flux selected by "Numerical flux: LLF", and it is
      ! also the flux the positivity repair in RK_integration substitutes at a
      ! single interface: with first-order (cell-average) input states the
      ! Lax-Friedrichs update maps positive density and positive internal
      ! energy to positive density and positive internal energy whenever
      ! dt(|v|+c)/dr <= 1 (Perthame & Shu 1996, Numer. Math. 73, 119; the LF
      ! lemma restated in Zhang & Shu 2010, J. Comput. Phys. 229, 3091), a
      ! bound the CFL number of a run (default 0.6) respects.  Both proofs
      ! assume the coefficient above; with it, the lemma applies to this
      ! flux as written.
      subroutine lax_friedrichs_flux(WL,WR,NF,p_out,jL,jR)

      real*8, intent(in) :: WL(3),WR(3)
      integer, intent(in) :: jL,jR
      real*8 :: uL(3),uR(3)
      real*8 :: FL(3),FR(3)
      real*8 :: rhoL,vL,pL,aL
      real*8 :: rhoR,vR,pR,aR
      real*8 :: a1
      real*8, intent(out) :: NF(3),p_out

      ! Exctract left state
      rhoL = WL(1)
      vL   = WL(2)
      pL   = WL(3)
      aL   = sqrt(adiabatic_index_from_state(jL,rhoL,pL)*pL/rhoL)

      ! Exctract right state
      rhoR = WR(1)
      vR   = WR(2)
      pR   = WR(3)
      aR   = sqrt(adiabatic_index_from_state(jR,rhoR,pR)*pR/rhoR)

      ! Evaluate physical left and right flux
      call Phys_flux(WL,FL,jL)
      call Phys_flux(WR,FR,jR)

      ! Largest signal speed of the pair: the spectral radius |v| + a of the
      ! flux Jacobian, taken on whichever side is larger (see the header).
      a1 = max(abs(vL)+aL,abs(vR)+aR)

      ! Get vector of conservative variables
      call W_to_U_comp(WL,uL,jL)
      call W_to_U_comp(WR,uR,jR)

      ! Evaluate numerical flux
      NF = 0.5*(FL + FR - a1*(uR-uL))

      ! Output pressure
      p_out = 0.5*(pR + pL)

      ! End of subroutine
      end subroutine lax_friedrichs_flux

      !-----------------------------------------------------------!

      ! Harten-Lax-van Leer flux with Einfeldt's wave bounds (HLLE): the
      ! exact flux of the two-wave approximate Riemann problem whose signal
      ! speeds SL and SR bound the physical waves.  Under that bound the
      ! intermediate state is an average of admissible states, so the flux
      ! is positively conservative -- it maps positive density and positive
      ! internal energy to positive density and positive internal energy
      ! under the CFL bound (Einfeldt 1988, SIAM J. Numer. Anal. 25, 294;
      ! Einfeldt, Munz, Roe & Sjogreen 1991, J. Comput. Phys. 92, 273).
      ! It resolves no contact wave, which is why it is used only at the
      ! faces where the Roe flux has no admissible star state.
      subroutine hlle_flux(WL,WR,FL,FR,SL,SR,jL,jR,NF)

      real*8, intent(in) :: WL(3),WR(3),FL(3),FR(3)
      real*8, intent(in) :: SL,SR
      integer, intent(in) :: jL,jR
      real*8 :: uL(3),uR(3)
      real*8, intent(out) :: NF(3)

      if (SL .ge. 0.0d0) then

         NF = FL

      elseif (SR .le. 0.0d0) then

         NF = FR

      else

         ! Get vector of conservative variables
         call W_to_U_comp(WL,uL,jL)
         call W_to_U_comp(WR,uR,jR)

         NF = (SR*FL - SL*FR + SL*SR*(uR-uL))/(SR-SL)

      endif

      ! End of subroutine
      end subroutine hlle_flux

      !-----------------------------------------------------------!

      ! Subroutine to compute the physical flux function
      subroutine Phys_flux(W,PF,jcell)
      
      real*8, intent(in)  :: W(3)
      integer, intent(in) :: jcell
      real*8 :: rho,v,p,E
      real*8, intent(out) :: PF(3)
      
      ! Get physical variables
      rho = W(1)
      v   = W(2)
      p   = W(3)
      E   = 0.5*rho*v*v + energy_density_from_pressure(jcell,rho,p)
            
      ! Output exact flux vector
      PF(1) = rho*v
      PF(2) = rho*v*v
      PF(3) = v*(E+p)
      
      ! Add pressure for PLM discretization.  Not under the well-balanced
      ! option: there the pressure force of both discretizations is assembled in
      ! RK_rhs from the departure of the face pressure from the cell's own
      ! hydrostatic equilibrium, so the momentum flux carried through the
      ! face is the kinetic one alone.
      if (use_plm .and. .not. well_balanced) PF(2) = PF(2) + p
      
      ! End of subroutine
      end subroutine Phys_flux     
      
      ! End of module
      end module Numerical_Fluxes
