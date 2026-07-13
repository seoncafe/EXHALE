	module initial_conditions
	! Set general initial conditions (isothermal atmosphere or transonic wind)

	use global_parameters
	use species_table, only: isp_HI, isp_HII, isp_HeI, isp_HeII,      &
	                         isp_HeIII, isp_HeTR,                      &
	                         n_melem, melem_i0, melem_top, mion_fsp
	use grav_func

	implicit none

	contains

	subroutine set_IC(W,T,f_sp)
	! Subroutine to set IC for the model

	integer :: j,i_rhalf
	integer :: e,k,c,i0
	real*8 :: r_half, minrho
	real*8 :: b0_eff
	real*8 :: c2, cs, xi
	logical :: wind_ok
	! Hot-Parker warm-seed IC: Parker velocity head-start kept in temp arrays
	! while the density is the (stable, balanced) cold hydrostatic profile.
	real*8, dimension(1-Ng:N+Ng) :: rho_p, v_p
	real*8 :: swin
	real*8, dimension(3,1-Ng:N+Ng), intent(out) :: W
	real*8, dimension(1-Ng:N+Ng),   intent(out) :: T
	real*8, dimension(1-Ng:N+Ng,n_species), intent(out) :: f_sp


	!--- Set initial conditions for thermodynamic variables ---!

	! Automatic IC selection ("IC mode: auto"): pick the IC family from the
	! cold sonic-point topology of the actual potential before anything is
	! built (v2 decision tree, docs/auto_ic_design.md). Explicit legacy keys
	! (Transonic IC / Hot Parker IC) take precedence and skip the selector.
	if (ic_mode .eq. 3) call select_IC_auto

	! Density + velocity.  Two IC families (see transonic_ic in parameters.f90):
	!   (1) transonic isothermal-wind IC (transonic_ic=.true.): solve the steady
	!       isothermal-wind profile from the ATES potential (Bernoulli integral)
	!       and set rho from mass conservation.  Needed for deep-RLOF bases whose
	!       sonic point sits near L1 and which launch no wind from the
	!       hydrostatic IC.
	!   (2) isothermal hydrostatic atmosphere + small v seed (default): the
	!       original ATES IC.
	wind_ok = .false.
	if (hot_parker_ic) then
		! Hot-Parker warm-seed IC. The benchmark + root-cause analysis
		! (docs/initial_condition_benchmark) showed that seeding the FULL
		! isothermal Parker wind (its density AND velocity) destabilizes the run:
		! the hot, low-density Parker column is globally inconsistent with ATES's
		! NON-isothermal equilibrium and drains onto the base, driving a self-
		! amplifying inflow (the cold hydrostatic IC breathes too, but self-limits
		! and converges). The robust warm seed therefore keeps the STABLE cold
		! hydrostatic density and overlays only (i) a warm T + ionization ramp
		! (which shortcuts the stiff, dominant thermal/ionization relaxation) and
		! (ii) the Parker VELOCITY as a transonic head-start, both anchored to the
		! cold static base by a spatial smoothstep window (cold below r=1, warming
		! to the wind value at r>=hp_base_rtr). The Parker velocity is obtained for
		! the head-start only (density discarded).
		c2 = (ntot_bc + dp_bc)/rho_bc * (T_wind_ic/T0) * 2.0d0
		call wind_profile(rho_p, v_p, c2, wind_ok)
		if (wind_ok) then
			write(*,'(A,F8.1,A)') '    (set_IC.f90) Using hot-Parker warm-seed IC, T_wind =', &
			                      T_wind_ic, ' K'
		else
			v_p = 0.0d0    ! no interior sonic point: fall back to the small linear seed
			write(*,*) '   (set_IC.f90) Hot-Parker warm-seed IC (no Parker ' // &
			           'velocity head-start; using small v seed)'
		endif
		write(*,*)
	else if (transonic_ic) then
		c2 = (ntot_bc + dp_bc)/rho_bc
		call wind_profile(W(1,:), W(2,:), c2, wind_ok)
		if (wind_ok) then
			write(*,*) '   (set_IC.f90) Using transonic isothermal-wind IC'
			write(*,*)
		else
			write(*,*) '   (set_IC.f90) Transonic IC requested but no ' // &
			           'interior sonic point; falling back to hydrostatic IC'
			write(*,*)
		endif
	endif

	! Hydrostatic density for the default IC AND the hot-Parker warm-seed IC
	! (hot-Parker keeps this stable density; only T/ionization/velocity are warmed).
	if (hot_parker_ic .or. .not. wind_ok) then
		! Density (isothermal hydrostatic)
		b0_eff = 1.0d0	! Change if the planet b0 is too low - only for IC
		do
			W(1,:) = (/ (rho_bc*exp(b0_eff*(-Gphi_c(j) + Gphi_c(0))), 	&
					j = 1-Ng,N+Ng) /)

			! Calculate minimum of density profile
			r_half = 0.5e0*(r_max + 1.0e0)
			i_rhalf = minloc(abs(r-r_half), dim = 1)
			minrho = W(1,i_rhalf)

			if (minrho .gt. 1.0e-8) then
				b0_eff = b0_eff + 0.2
			else
				exit
			endif

		enddo

		! Write b0_eff to output
		write(*,'(A31,F5.1,A7)') '    (set_IC.f90) Using b0_eff =',b0_eff,' for IC'
		write(*,*)

		! Velocity: small linear seed (default), or the Parker head-start (hot-Parker)
		W(2,:) = 0.5*(r-r(0))
	endif

	! Fix density in outer layers
	where (W(1,:).lt.(1.0e-8)) W(1,:) = 1.0e-8

	if (hot_parker_ic) then
		! Warm-seed overlay on the hydrostatic density. xi is a spatial smoothstep
		! window (0 at the cold static base -> 1 in the warm ionized wind at
		! r>=hp_base_rtr); it ramps T, the H/He ionization, the pressure, and the
		! Parker velocity head-start. At the base xi=0 so T->T0, p->1+dp_bc, v->0,
		! matching the lower BC exactly (no inverted gradient); above hp_base_rtr
		! the gas is warm (T_wind), ionized, and carries the transonic velocity.
		! Initialize ALL species columns (f_sp is intent(out); metals and
	! molecular columns stay zero unless set below or loaded).
	f_sp = 0.0d0

	do j = 1-Ng, N+Ng
			swin = min(max((r(j) - 1.0d0)/(hp_base_rtr - 1.0d0), 0.0d0), 1.0d0)
			xi   = swin*swin*(3.0d0 - 2.0d0*swin)               ! smoothstep window
			W(2,j) = v_p(j) * swin                              ! Parker velocity head-start, 0 at base
			T(j)   = 1.0d0 + (T_wind_ic/T0 - 1.0d0)*xi
			W(3,j) = W(1,j)*((ntot_bc + dp_bc)/rho_bc)*T(j)*(1.0d0 + xi)
			f_sp(j,isp_HI)    = (1.0d0 - (dp_bc + (1.0d0-dp_bc)*xi))/mass_per_H
			f_sp(j,isp_HII)   = (dp_bc + (1.0d0-dp_bc)*xi)/mass_per_H
			f_sp(j,isp_HeI)   = HeH*(1.0d0 - xi)/mass_per_H
			f_sp(j,isp_HeII)  = HeH*xi/mass_per_H
			f_sp(j,isp_HeIII) = 0.0d0
			f_sp(j,isp_HeTR)  = 0.0d0
		enddo
	else
		! Pressure (= rho * c_iso^2, isothermal; consistent with the base BC,
		!  where p_base = ntot_bc + dp_bc and rho_base = rho_bc)
		W(3,:) = (ntot_bc + dp_bc)*W(1,:)/rho_bc
		! Temperature
		T = 1.0
		! Ionized fractions (mostly neutral)
		f_sp(:,isp_HI)    = (1.0 - dp_bc)/mass_per_H
		f_sp(:,isp_HII)   = dp_bc/mass_per_H
		f_sp(:,isp_HeI)   = HeH*(1.0 - dp_bc)/mass_per_H
		f_sp(:,isp_HeII)  = 1.0d-10*HeH/mass_per_H
		f_sp(:,isp_HeIII) = 0.0
		f_sp(:,isp_HeTR)  = 0.0
	endif

   ! Metals: start mostly neutral, element by element from the abundance
   ! array (f_X = melem_ab/mass_per_H, so n_X = melem_ab*n_H; higher stages zero).
   ! Adding an element is a species_table + abundance change only.
   do e = 1, n_melem
      i0 = melem_i0(e)
      do k = 0, melem_top(e)
         c = mion_fsp(i0+k)
         if (k .eq. 0) then
            f_sp(:,c) = melem_ab(e)/mass_per_H
         else
            f_sp(:,c) = 0.0
         endif
      enddo
   enddo

	! Dump the initial condition (dimensional) for inspection / IC benchmarking.
	open(unit = 77, file = 'output/IC_dump.txt', status = 'replace')
	write(77,'(A)') '# r[R_p]        n[cm-3]        v[cm/s]         p[cgs]          T[K]'
	do j = 1-Ng, N+Ng
		write(77,'(5(ES15.6,1X))') r(j), W(1,j)*n0, W(2,j)*v0, W(3,j)*p0, T(j)*T0
	enddo
	close(77)

	! End of subroutine
	end subroutine set_IC

	!-------------------------------------------------------!

	subroutine select_IC_auto
	! Automatic IC-family selection ("IC mode: auto"); v2 decision tree of
	! docs/auto_ic_design.md. A single exact probe -- does the COLD base
	! sound speed admit an interior sonic point in the actual (Roche or
	! spherical) potential? -- separates the regimes with no tunable
	! threshold:
	!   interior sonic point -> transonic IC. Catches (a) deep-RLOF bases,
	!     where the Roche dphi/dr -> 0 toward L1 always crosses the
	!     critical condition (e.g. WASP-121b), and (b) low-gravity
	!     boil-off planets, where r_c ~ b0/(2 c0^2) falls inside the
	!     domain; both launch no wind (or breathe) from a cold start.
	!   no interior sonic point -> cold hydrostatic IC (classic EUV-heated
	!     wind, e.g. HD 209458 b: r_c ~ 54 R_p >> r_max; launches fine).
	! The hot-Parker warm seed is deliberately NOT auto-selected: the IC
	! benchmark showed no speedup where the cold start works, and a warm
	! sonic point exists for essentially every hot Jupiter. It remains a
	! manual option. b0 (= Kubyshkina's escape parameter Lambda) is logged
	! as the regime diagnostic.
	real*8  :: c2_cold, rc
	logical :: have_rc

	! Explicit legacy keys take precedence over auto.
	if (transonic_ic .or. hot_parker_ic) then
		write(*,*) '   (select_IC_auto) explicit IC key present; ' // &
		           'auto selection skipped'
		return
	endif

	! Cold-base isothermal sound speed^2 (same c2 the transonic IC uses).
	c2_cold = (ntot_bc + dp_bc)/rho_bc
	call find_sonic(c2_cold, rc, have_rc)

	if (have_rc) then
		transonic_ic = .true.
		write(*,'(A,F7.3,A,F7.1,A)')                                    &
			'    (select_IC_auto) interior cold sonic point at r_c =',  &
			rc, ' R_p (b0 =', b0, '): transonic IC'
	else
		write(*,'(A,F7.1,A)')                                           &
			'    (select_IC_auto) no interior cold sonic point (b0 =',  &
			b0, '): cold hydrostatic IC'
	endif

	! Re-evaluate force_start: the assignment in input_read.f90 ran before
	! auto could flip the flags (mirrors that logic; read in the main loop,
	! which starts after set_IC).
	force_start = (transonic_ic .or. hot_parker_ic) .and. .not. do_only_pp

	! End of subroutine
	end subroutine select_IC_auto

	!-------------------------------------------------------!

	subroutine wind_profile(rho_w, v_w, c2, ok)
	! Transonic isothermal-wind initial profile via the algebraic (Bernoulli)
	! integral of the steady isothermal-wind equation in the ATES potential.
	!
	! Steady, isothermal (T=1), spherical (area ~ r^2) flow in the dimensionless
	! potential Gphi = phi(r) obeys
	!     1/2 v^2 - c2*ln v - 2*c2*ln r + phi(r) = B,                       (*)
	! with c2 the dimensionless isothermal sound speed^2 and B fixed by the
	! sonic point (v = cs at r = rc, where Dphi(rc) = 2*c2/rc).  At each radius
	! (*) has a subsonic root (v<cs) and a supersonic root (v>cs); the transonic
	! wind takes the subsonic root for r<rc and the supersonic root for r>rc.
	! rho then follows from steady mass conservation rho*v*r^2 = Mdot, with Mdot
	! pinned to the base density BC (rho = rho_bc at r(0)).
	!
	! ok = .false. signals no interior sonic point (a subsonic breeze / bound
	! case), so set_IC keeps the hydrostatic IC.
	real*8, dimension(1-Ng:N+Ng), intent(out) :: rho_w, v_w
	real*8, intent(in) :: c2
	logical, intent(out) :: ok
	real*8 :: cs, rc, Bc, Kr, hmin, Mdot, vbase, rj
	integer :: j

	! c2 = dimensionless isothermal sound speed^2 (passed in: cold base value for
	! the transonic IC, or the warm c2_hot for the hot-Parker IC).
	cs = sqrt(c2)

	! Locate the sonic point
	call find_sonic(c2, rc, ok)
	if (.not. ok) return

	! Bernoulli constant from the sonic point (v = cs at r = rc); hmin is the
	! minimum of h(v) = 1/2 v^2 - c2 ln v, attained at v = cs.
	Bc   = 0.5d0*c2 - c2*log(cs) - 2.0d0*c2*log(rc) + phi(rc)
	hmin = 0.5d0*c2 - c2*log(cs)

	! Solve (*) for v on the appropriate branch at every grid point
	do j = 1-Ng, N+Ng
		rj = r(j)
		Kr = Bc + 2.0d0*c2*log(rj) - phi(rj)   ! solve h(v) = Kr
		if (Kr .le. hmin) then
			v_w(j) = cs                        ! at (or numerically at) rc
		else if (rj .lt. rc) then
			v_w(j) = wind_root(Kr, c2, cs, .false.)   ! subsonic root in (0,cs)
		else
			v_w(j) = wind_root(Kr, c2, cs, .true.)    ! supersonic root in (cs,inf)
		endif
	enddo

	! Density from mass conservation, normalized to the base BC (rho_bc at r(0))
	vbase = v_w(0)
	Mdot  = rho_bc*vbase*r(0)*r(0)
	do j = 1-Ng, N+Ng
		rho_w(j) = Mdot/(v_w(j)*r(j)*r(j))
	enddo

	! End of subroutine
	end subroutine wind_profile

	!-------------------------------------------------------!

	subroutine find_sonic(c2, rc, have_rc)
	! Locate the transonic (sonic) point r_c in [r(0), r_max] where the
	! isothermal-wind critical condition Dphi(r_c) = 2*c2/r_c holds.
	! f(r) = Dphi(r) - 2*c2/r is > 0 at the base (gravity dominated) and turns
	! negative toward L1 (Dphi -> 0 there in Roche mode).  Bisection on the
	! sign change.  have_rc = .false. signals no interior crossing (a subsonic
	! breeze / bound case), so the caller falls back to the hydrostatic IC.
	! If the base itself is already at/above the sonic point (f(base) <= 0),
	! the flow is supersonic from the base and rc is set to the base radius.
	real*8, intent(in)  :: c2
	real*8, intent(out) :: rc
	logical, intent(out) :: have_rc
	real*8 :: rlo, rhi, flo, fhi, rmid, fmid
	integer :: it

	rlo = r(0)
	rhi = r_max
	flo = Dphi(rlo) - 2.0d0*c2/rlo
	fhi = Dphi(rhi) - 2.0d0*c2/rhi

	have_rc = .false.
	rc = rlo

	if (flo .le. 0.0d0) then
		! Base already at/above the sonic point: supersonic from the base.
		have_rc = .true.
		rc = rlo
		return
	endif
	if (fhi .ge. 0.0d0) then
		! Gravity dominates throughout the (truncated) domain: no interior
		! sonic point -> subsonic breeze.  Caller keeps the hydrostatic IC.
		have_rc = .false.
		return
	endif

	! Bisect f(r) = 0 (f decreasing through the crossing)
	do it = 1, 200
		rmid = 0.5d0*(rlo + rhi)
		fmid = Dphi(rmid) - 2.0d0*c2/rmid
		if (fmid .gt. 0.0d0) then
			rlo = rmid
		else
			rhi = rmid
		endif
		if ((rhi - rlo) .lt. 1.0d-12*rmid) exit
	enddo
	rc = 0.5d0*(rlo + rhi)
	have_rc = .true.

	! End of subroutine
	end subroutine find_sonic

	!-------------------------------------------------------!

	double precision function wind_root(Kr, c2, cs, supersonic)
	! Solve h(v) = 1/2 v^2 - c2*ln v = Kr for v on the requested branch.
	! h decreases on (0,cs) and increases on (cs,inf), so each branch has a
	! single root; bisection is unconditionally robust.
	real*8, intent(in)  :: Kr, c2, cs
	logical, intent(in) :: supersonic
	real*8 :: vlo, vhi, vmid, hmid
	integer :: it

	if (supersonic) then
		vlo = cs
		vhi = max(2.0d0*cs, sqrt(2.0d0*max(Kr, c2)))
		! grow the upper bracket until h(vhi) > Kr
		do it = 1, 200
			if (0.5d0*vhi*vhi - c2*log(vhi) .gt. Kr) exit
			vhi = 2.0d0*vhi
		enddo
	else
		vhi = cs
		vlo = cs*1.0d-12                  ! h -> +inf as v -> 0, so a root exists
	endif

	do it = 1, 200
		vmid = 0.5d0*(vlo + vhi)
		hmid = 0.5d0*vmid*vmid - c2*log(vmid)
		if (supersonic) then
			! h increasing: if h<Kr the root is to the right
			if (hmid .lt. Kr) then
				vlo = vmid
			else
				vhi = vmid
			endif
		else
			! h decreasing: if h<Kr the root is to the left (smaller v)
			if (hmid .lt. Kr) then
				vhi = vmid
			else
				vlo = vmid
			endif
		endif
		if (abs(vhi - vlo) .lt. 1.0d-12*max(vmid,1.0d-30)) exit
	enddo
	wind_root = 0.5d0*(vlo + vhi)

	! End of function
	end function wind_root

	! End of module
	end module initial_conditions
