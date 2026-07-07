      module wae_continuation
      ! static-BC continuation driver: ramp the system parameters
      ! (Ftot, gravity, star) from a shipped seed windsoln to a target,
      ! re-relaxing at each adaptive step. Mirrors the Python wrapper's
      ! ramp_var/ramp_to with static_bcs=True (the base BCs stay at the
      ! seed values). In-memory guess hand-off replaces the C
      ! windsoln<->guess.inp file round-trip.
      use wae_config,  only: nsp => wae_nspecies, m => wae_m, ne => wae_ne, &
                             nb => wae_nb, itmax => wae_itmax,             &
                             conv => wae_conv, slowc => wae_slowc,         &
                             VSCALE=>wae_vscale, ZSCALE=>wae_zscale,       &
                             TEMPSCALE=>wae_tempscale, FPSCALE=>wae_fpscale, &
                             NCOLSCALE=>wae_ncolscale
      use wae_params,  only: par => wae_par, wae_G, wae_cs0_val, wae_K, wae_MH
      use wae_config,  only: T0 => wae_T0
      use wae_spectrum,only: wae_HX, wae_atomic_mass, wae_Ftot,             &
                             wae_npts, wae_hc_over_wl, wae_wPhi_wl,         &
                             wae_sigma_wl, wae_ion_pot
      use wae_grid,    only: wae_x
      use wae_nr,      only: wae_solvde
      use wae_status,  only: wae_err
      use wae_config,  only: RHO0 => wae_RHO0, addpts => wae_addpts,        &
                             NCOL0 => wae_NCOL0
      use wae_intode,  only: wae_integrate_ode
      use wae_glq_rates, only: wae_glq_rates_eval
      implicit none

      ! relaxation work arrays + indices/scales (set once)
      real*8  :: cy(ne, m)            ! current relaxation solution (the guess)
      integer :: cv(ne)               ! indexv
      real*8  :: sv(ne)               ! scalv
      real*8  :: rhoscale = 100.0d0
      real*8, parameter :: WAE_PI     = 3.14159265358979323846d0
      real*8, parameter :: WAE_SIG_SB = 5.670508538505d-05  ! wind-ae const.sig_SB
      contains

      subroutine setup_indices_scales()
      integer :: i, j
      do i = 1, ne
         cv(i) = i
      end do
      cv(1) = 3 + 2*nsp; cv(2) = 4 + 2*nsp; cv(3) = 1; cv(4) = 2 + nsp
      do j = 0, nsp-1
         cv(5+j) = 2 + j; cv(5+j+nsp) = 3 + nsp + j
      end do
      sv(1) = VSCALE; sv(2) = ZSCALE; sv(3) = rhoscale; sv(4) = TEMPSCALE
      do j = 0, nsp-1
         sv(5+j) = FPSCALE; sv(5+nsp+j) = NCOLSCALE
      end do
      end subroutine setup_indices_scales

      integer function solve()
      ! relax cy in place from itself as the guess; return wae_err.
      real*8, allocatable :: c(:,:,:)
      real*8 :: s(ne, 2*ne+1)
      allocate(c(ne, ne-nb+1, m+1)); c = 0.0d0; s = 0.0d0
      wae_err = 0
      call wae_solvde(itmax, conv, slowc, sv, cv, ne, nb, m, cy, c, s)
      deallocate(c)
      solve = wae_err
      end function solve

      subroutine apply_var(which, val)
      ! set parameter `which` to val and refresh derived quantities
      integer, intent(in) :: which
      real*8,  intent(in) :: val
      select case (which)
      case (1); par%Ftot = val; wae_Ftot = val      ! glq uses wae_Ftot copy
      case (2); par%Mp = val
      case (3); par%Rp = val
      case (4); par%Mstar = val
      case (5); par%semimajor = val
      case (6); par%Lstar = val
      end select
      ! H0 = CS0^2 Rp^2 /(G Mp) depends on Mp, Rp
      par%H0 = wae_cs0_val**2 * par%Rp**2 / wae_G / par%Mp
      end subroutine apply_var

      real*8 function get_var(which)
      integer, intent(in) :: which
      select case (which)
      case (1); get_var = par%Ftot
      case (2); get_var = par%Mp
      case (3); get_var = par%Rp
      case (4); get_var = par%Mstar
      case (5); get_var = par%semimajor
      case (6); get_var = par%Lstar
      end select
      end function get_var

      integer function ramp_var(which, target, label, static_bcs)
      ! Adaptive multiplicative ramp of one system parameter. With static_bcs
      ! absent or .true., the seed base BCs are held fixed (static-BC mode). With
      ! static_bcs=.false. (self-consistent-BC mode), after 5 consecutive failed steps the base BCs
      ! and the sonic-point column density are re-converged (ramp_base_bcs +
      ! converge_Ncol_sp) before retrying -- this is what lets far-from-seed
      ! planets keep moving. Returns 0 on success, 101 on failure.
      integer, intent(in) :: which
      real*8,  intent(in) :: target
      character(len=*), intent(in) :: label
      logical, intent(in), optional :: static_bcs
      real*8  :: cur, trial, delta, flip
      real*8  :: ysave(ne, m)
      integer :: failed, ierr, junk, conv_cntr
      logical :: static
      static = .true.
      if (present(static_bcs)) static = static_bcs
      cur = get_var(which)
      ramp_var = 0
      if (target .eq. 0.0d0 .or. abs(cur-target)/abs(target) .lt. 1.0d-10) then
         write(*,'(A,A,A)') '  ', trim(label), ' already at target.'
         return
      end if
      flip = 1.0d0
      if (cur .gt. target) flip = -1.0d0
      delta = 0.02d0*flip
      failed = 0
      conv_cntr = 0
      do while (abs(cur-target)/abs(target) .gt. 1.0d-10)
         trial = cur*(1.0d0 + delta)
         if (flip*trial .gt. flip*target) trial = target
         ysave = cy
         call apply_var(which, trial)
         ierr = solve()
         if (ierr .ne. 0) then
            cy = ysave                 ! restore guess
            call apply_var(which, cur) ! restore param
            delta = delta/2.0d0
            failed = failed + 1
            if ((.not. static) .and. failed .eq. 5) then
               write(*,'(A,A)') '   ...intermediate BC reconvergence during ', &
                  trim(label)
               junk = ramp_base_bcs(.false.)
               junk = converge_mol_atomic()
               junk = converge_Ncol_sp()
               failed = 0
               delta = 0.02d0*flip     ! reset the step after moving the BCs
            else if (failed .gt. 60) then
               write(*,'(A,A,A,ES12.5)') '  RAMP FAIL ', trim(label),    &
                  ' stuck at ', cur
               ramp_var = 101; return
            end if
         else
            cur = trial
            if (.not. (flip .lt. 0.0d0 .and. delta .le. -0.5d0)) delta = delta*2.0d0
            failed = 0
            ! self-consistent-BC (proactive, cheap): periodically check the molecular/atomic
            ! transition and turn the bolometric layer off once the base enters
            ! the wind. converge_mol_atomic does NO solve unless it actually
            ! turns bolo off (drop_index<=10), so for a seed-adjacent planet
            ! (bolo stays on) it is just a glq scan -- HD209458b stays fast,
            ! while HD189733b turns bolo off early instead of waiting for a
            ! 5-fail stall. The expensive base-BC/Ncol re-convergence stays
            ! reactive only (the failed==5 branch above).
            if (.not. static) then
               conv_cntr = conv_cntr + 1
               if (conv_cntr .ge. 10) then
                  junk = converge_mol_atomic()
                  conv_cntr = 0
               end if
            end if
         end if
      end do
      write(*,'(A,A,A,ES12.5)') '  ramped ', trim(label), ' -> ', cur
      end function ramp_var

      integer function ramp_to(Ftot_t, Mp_t, Rp_t, Mstar_t, a_t, Lstar_t,    &
                               static_bcs)
      ! ramp order matches the wrapper: Ftot -> gravity (Rp,Mp) ->
      ! star (Mstar, semimajor, Lstar). With static_bcs absent/.true. the seed
      ! BCs are held fixed (static-BC mode); with static_bcs=.false. (self-consistent-BC mode) the base BCs
      ! re-converge when a ramp stalls and once more at the end.
      real*8, intent(in) :: Ftot_t, Mp_t, Rp_t, Mstar_t, a_t, Lstar_t
      logical, intent(in), optional :: static_bcs
      integer :: r
      logical :: static
      static = .true.
      if (present(static_bcs)) static = static_bcs
      ramp_to = 0
      r = ramp_var(1, Ftot_t,  'Ftot',      static); if (r.ne.0) then; ramp_to=r; return; end if
      r = ramp_var(3, Rp_t,    'Rp',        static); if (r.ne.0) then; ramp_to=r; return; end if
      r = ramp_var(2, Mp_t,    'Mp',        static); if (r.ne.0) then; ramp_to=r; return; end if
      r = ramp_var(4, Mstar_t, 'Mstar',     static); if (r.ne.0) then; ramp_to=r; return; end if
      r = ramp_var(5, a_t,     'semimajor', static); if (r.ne.0) then; ramp_to=r; return; end if
      if (Lstar_t .gt. 0.0d0) then
         r = ramp_var(6, Lstar_t, 'Lstar',  static); if (r.ne.0) then; ramp_to=r; return; end if
      end if
      ! final BC self-consistency pass (self-consistent-BC mode only)
      if (.not. static) then
         r = ramp_base_bcs(.false.)
         r = converge_mol_atomic()
         r = converge_Ncol_sp()
      end if
      end function ramp_to

      ! ====================================================================
      !  self-consistent-BC machinery
      ! ====================================================================

      subroutine base_bcs(R_out, Rmax_out, rho_out, T_out)
      ! Port of relax_wrapper.base_bcs (defaults Kappa_opt=4e-3, Kappa_IR=1e-2,
      ! adiabat=False; base_press auto from the current base pressure). Computes
      ! the self-consistent base BCs for the CURRENT planet params (par) and base
      ! solution (cy): assumes Rp is the optical slant-path tau=1 surface, sets
      ! the simulation base at the vertical-IR tau=1 radius R_IR if R_IR>Rp, else
      ! at the microbar-pressure radius R_mbar. Returns R in units of Rp and
      ! rho/T normalized to (RHO0, T0) -- directly comparable to par%{Rmin,
      ! rho_rmin, T_rmin}.
      real*8, intent(out) :: R_out, Rmax_out, rho_out, T_out
      real*8, parameter :: KAP_OPT = 4.0d-3, KAP_IR = 1.0d-2
      real*8 :: bolo, mu, mu_denom, P, base_press
      real*8 :: F_opt, T_skin, T_eff, rho_Rp, P_Rp, cs2, R_mbar, rho_mbar
      real*8 :: Hsc, rho_IR, R_IR
      integer :: j
      bolo = par%bolo_heat_cool
      ! mu at base = calc_mu()[0]: molecular (molec_adjust*mH, weighted by the
      ! bolo smoothing erf which is 1 at the base) blended with the atomic
      ! mH/mu_denom for the (1-bolo) remainder.
      mu_denom = 0.0d0
      do j = 1, nsp
         mu_denom = mu_denom + (wae_MH/par%atomic_mass(j))*par%HX(j)        &
                               *(2.0d0 - cy(4+j,1))
      end do
      mu = (wae_MH/mu_denom)*(1.0d0 - bolo) + par%molec_adjust*wae_MH*bolo
      ! base pressure [microbar] from the current base state, rounded to nearest
      ! 10, clamped >= 1 (matches np.round(P/10)*10 with a >=1 floor).
      P = (cy(3,1)*RHO0)*wae_K*(cy(4,1)*T0)/mu
      base_press = dnint(P/10.0d0)*10.0d0
      if (base_press .lt. 1.0d0) base_press = 1.0d0
      ! skin / effective temperature from the stellar bolometric flux
      F_opt  = par%Lstar/(4.0d0*WAE_PI*par%semimajor**2)
      T_skin = (F_opt*(KAP_OPT+KAP_IR/4.0d0)/(2.0d0*WAE_SIG_SB*KAP_IR))**0.25d0
      T_eff  = (F_opt/(4.0d0*WAE_SIG_SB))**0.25d0
      ! slant-path tau=1 base density and the microbar-pressure radius
      rho_Rp = sqrt(mu*wae_G*par%Mp/(8.0d0*par%Rp**3*wae_K*T_skin))/KAP_OPT
      P_Rp   = rho_Rp*wae_K*T_skin/mu
      cs2    = wae_K*T_skin/mu
      R_mbar = 1.0d0/((cs2/(wae_G*par%Mp))*log(base_press/P_Rp)             &
                      + 1.0d0/par%Rp)
      rho_mbar = base_press*mu/(wae_K*T_skin)
      ! vertical IR tau=1 radius (adiabat=False branch)
      Hsc    = wae_K*T_eff*par%Rp**2/(mu*wae_G*par%Mp)
      rho_IR = 1.0d0/(KAP_IR*Hsc)
      if (par%Rp/Hsc .gt. 690.0d0) then
         R_IR = 0.0d0
      else
         R_IR = par%Rp**2/(Hsc*log(rho_IR/rho_Rp*exp(par%Rp/Hsc)))
      end if
      Rmax_out = par%Rmax
      if (R_IR .gt. par%Rp) then
         R_out = R_IR/par%Rp; rho_out = rho_IR/RHO0; T_out = T_eff/T0
      else
         if (bolo .eq. 0.0d0) T_skin = (sum(cy(4,1:10))/10.0d0)*T0
         R_out = R_mbar/par%Rp; rho_out = rho_mbar/RHO0; T_out = T_skin/T0
      end if
      end subroutine base_bcs

      real*8 function get_bc(field)
      ! base BC field: 1=Rmin (units Rp), 2=rho_rmin (RHO0), 3=T_rmin (T0)
      integer, intent(in) :: field
      select case (field)
      case (1); get_bc = par%Rmin
      case (2); get_bc = par%rho_rmin
      case (3); get_bc = par%T_rmin
      end select
      end function get_bc

      subroutine set_bc(field, val)
      integer, intent(in) :: field
      real*8,  intent(in) :: val
      select case (field)
      case (1); par%Rmin = val
      case (2); par%rho_rmin = val
      case (3); par%T_rmin = val
      end select
      end subroutine set_bc

      integer function ramp_bc(field, goal, label)
      ! adaptive stepper for one base BC (mirrors ramp_Rmin/ramp_T_rmin/
      ! ramp_rho_rmin): set the field, re-solve; on failure restore the guess
      ! and halve the step. rhoscale tracks the rho order of magnitude so the
      ! relaxation scale stays valid (the C code recompiles for this; here it is
      ! just a runtime scalv entry).
      integer, intent(in) :: field
      real*8,  intent(in) :: goal
      character(len=*), intent(in) :: label
      real*8  :: cur, trial, delta, flip, ysave(ne, m)
      integer :: fails, ierr
      cur = get_bc(field); ramp_bc = 0
      if (goal .eq. 0.0d0 .or. abs(cur-goal)/abs(goal) .lt. 2.0d-3) return
      flip = 1.0d0; if (cur .gt. goal) flip = -1.0d0
      delta = 0.05d0*flip; fails = 0
      do while (abs(cur-goal)/abs(goal) .gt. 2.0d-3)
         trial = cur*(1.0d0 + delta)
         if (flip*trial .gt. flip*goal) trial = goal
         ysave = cy
         call set_bc(field, trial)
         if (field .eq. 2) sv(3) = 10.0d0**floor(log10(trial*0.01d0)) ! rho scale
         ierr = solve()
         if (ierr .ne. 0) then
            cy = ysave; call set_bc(field, cur)
            delta = delta/2.0d0; fails = fails + 1
            if (fails .gt. 40) then
               write(*,'(A,A,A,ES12.5)') '   ramp_bc give up ', trim(label), &
                  ' at ', cur
               ramp_bc = 1; return
            end if
         else
            cur = trial
            if (.not. (flip .lt. 0.0d0 .and. delta .le. -0.5d0))             &
               delta = delta*2.0d0
            fails = 0
         end if
      end do
      end function ramp_bc

      integer function ramp_base_bcs(static_bcs)
      ! Port of relax_wrapper.ramp_base_bcs: ramp Rmin/rho_rmin/T_rmin to the
      ! self-consistent base_bcs() target. Try a combined jump; on failure ramp
      ! each field individually. static_bcs=.true. -> no-op (keep seed BCs).
      logical, intent(in) :: static_bcs
      real*8  :: gR, gRmax, grho, gT, cR, crho, cT, ysave(ne, m)
      integer :: r1, r2, r3
      ramp_base_bcs = 0
      if (static_bcs) return
      call base_bcs(gR, gRmax, grho, gT)
      cR = par%Rmin; crho = par%rho_rmin; cT = par%T_rmin
      if (abs((gR-cR)/gR) + abs((grho-crho)/grho) + abs((gT-cT)/gT)          &
          .lt. 1.0d-2) return
      write(*,'(A,F6.3,A,F6.3,A,ES10.3,A,ES10.3,A,F6.0,A,F6.0,A)')           &
         '   reconverging base BCs: Rmin ', cR, '->', gR, ', rho ',          &
         crho*RHO0, '->', grho*RHO0, ', T ', cT*T0, '->', gT*T0, ' K'
      ysave = cy
      par%Rmin = gR; par%rho_rmin = grho; par%T_rmin = gT
      sv(3) = 10.0d0**floor(log10(grho*0.01d0))
      if (solve() .ne. 0) then
         cy = ysave; par%Rmin = cR; par%rho_rmin = crho; par%T_rmin = cT
         sv(3) = rhoscale
         r1 = ramp_bc(1, gR,   'Rmin')
         r3 = ramp_bc(3, gT,   'T_rmin')
         r2 = ramp_bc(2, grho, 'rho_rmin')
         if (r1 .ne. 0 .or. r2 .ne. 0 .or. r3 .ne. 0) ramp_base_bcs = 1
      end if
      end function ramp_base_bcs

      subroutine self_consistent_Ncol(goal)
      ! Port of relax_wrapper.self_consistent_Ncol (method=1): integrate the
      ! current relaxed solution outward, then Ncol_sp[j] = integral of the
      ! neutral number density of species j from the sonic point outward,
      ! normalized by NCOL0. Needs the outward integration (wae_integrate_ode).
      real*8, intent(out) :: goal(nsp)
      real*8, allocatable :: sr(:),srho(:),svv(:),sT(:),sq(:),sz(:)
      real*8, allocatable :: sYs(:,:),sNcol(:,:)
      real*8  :: dr, nneut, csum(nsp)
      integer :: tp, ntot, row, j
      tp = m + addpts
      allocate(sr(tp),srho(tp),svv(tp),sT(tp),sq(tp),sz(tp))
      allocate(sYs(tp,nsp),sNcol(tp,nsp))
      do row = 1, m
         sq(row)=wae_x(row); sz(row)=cy(2,row)
         sr(row)=wae_x(row)*cy(2,row)+par%Rmin
         srho(row)=cy(3,row); svv(row)=cy(1,row); sT(row)=cy(4,row)
         do j = 1, nsp
            sYs(row,j)=cy(4+j,row); sNcol(row,j)=cy(4+nsp+j,row)
         end do
      end do
      call wae_integrate_ode(sr, srho, svv, sT, sYs, sNcol, sq, sz, ntot)
      ! cumulative neutral column from the top down to the sonic point (row m)
      csum = 0.0d0
      do row = ntot, m, -1
         if (row .gt. 1) then
            dr = (sr(row)-sr(row-1))*par%Rp
         else
            dr = 0.0d0
         end if
         do j = 1, nsp
            nneut = par%HX(j)*(srho(row)*RHO0)*sYs(row,j)/par%atomic_mass(j)
            csum(j) = csum(j) + nneut*dr
         end do
      end do
      do j = 1, nsp
         goal(j) = csum(j)/NCOL0
      end do
      deallocate(sr,srho,svv,sT,sq,sz,sYs,sNcol)
      end subroutine self_consistent_Ncol

      integer function converge_Ncol_sp()
      ! Port of relax_wrapper.converge_Ncol_sp: fixed-point the sonic-point
      ! column-density BC to its self-consistent value (mean rel. diff < 8%),
      ! re-solving each iteration. Best-effort for a warm-start IC: on a solve
      ! failure it reverts and returns nonzero rather than the C code's
      ! Ncol-raising retry ladder.
      real*8  :: goal(nsp), cur(nsp), reldiff, ysave(ne, m)
      integer :: it, j, ierr
      converge_Ncol_sp = 0
      do it = 1, 10
         call self_consistent_Ncol(goal)
         do j = 1, nsp
            if (goal(j) .eq. 0.0d0) goal(j) = 1.0d-10*goal(1)
         end do
         cur = par%Ncol_sp
         reldiff = 0.0d0
         do j = 1, nsp
            reldiff = reldiff + (goal(j)-cur(j))/goal(j)
         end do
         reldiff = abs(reldiff/dble(nsp))
         if (reldiff .lt. 0.08d0) return
         ysave = cy
         par%Ncol_sp = goal
         ierr = solve()
         if (ierr .ne. 0) then
            cy = ysave; par%Ncol_sp = cur
            converge_Ncol_sp = 1; return
         end if
      end do
      end function converge_Ncol_sp

      integer function erf_drop_index()
      ! Port of relax_wrapper.erf_velocity's transition finder: the first grid
      ! point (from the base) where photoionization heating exceeds the
      ! magnitude of PdV (adiabatic-expansion) cooling -- i.e. where the wind
      ! "launches". heat = sum_s n_neutral,s * heatarr_s (wae_glq_rates_eval,
      ! per volume); cool = (k_B T / (molec_adjust m_H)) v drho/dr (< 0).
      real*8  :: heat, cool, rho_p, T_p, v_p, drhodr, csmol2, molad
      real*8  :: rp(m), ionarr(nsp*(nsp+1)), heatarr(nsp), nneut
      integer :: i, j
      molad = par%molec_adjust
      if (molad .le. 0.0d0) molad = 1.0d0
      do i = 1, m
         rp(i) = (wae_x(i)*cy(2,i) + par%Rmin)*par%Rp
      end do
      erf_drop_index = 0
      do i = 1, m
         rho_p = cy(3,i)*RHO0; T_p = cy(4,i)*T0; v_p = cy(1,i)*wae_cs0_val
         call wae_glq_rates_eval(cy(4+nsp+1:4+2*nsp,i), cy(5:4+nsp,i),       &
                                 ionarr, heatarr)
         heat = 0.0d0
         do j = 1, nsp
            nneut = cy(4+j,i)*par%HX(j)*rho_p/par%atomic_mass(j)
            heat = heat + nneut*heatarr(j)
         end do
         if (i .eq. 1) then
            drhodr = (cy(3,2)-cy(3,1))*RHO0/(rp(2)-rp(1))
         else if (i .eq. m) then
            drhodr = (cy(3,m)-cy(3,m-1))*RHO0/(rp(m)-rp(m-1))
         else
            drhodr = (cy(3,i+1)-cy(3,i-1))*RHO0/(rp(i+1)-rp(i-1))
         end if
         csmol2 = wae_K*T_p/(molad*wae_MH)
         cool   = csmol2*v_p*drhodr
         if (heat .gt. -cool) then
            erf_drop_index = i; return
         end if
      end do
      end function erf_drop_index

      integer function turn_off_bolo()
      ! Port of relax_wrapper.turn_off_bolo: ramp the bolometric-heating/cooling
      ! flag bolo_heat_cool 1 -> 0 and re-solve. Because the molecular-region
      ! smoothing erf is scaled by bolo_heat_cool (wae_get_mu/spQ), turning it
      ! off also drops the mean molecular weight to its atomic value -- the
      ! seed (bolo=1) -> HD189733b-class (bolo=0) transition.
      real*8  :: ysave(ne, m), bcur, delta
      integer :: fails
      turn_off_bolo = 0
      if (par%bolo_heat_cool .le. 0.0d0) return
      ysave = cy; bcur = par%bolo_heat_cool
      par%bolo_heat_cool = 0.0d0                 ! try a direct jump first
      if (solve() .eq. 0) then
         write(*,'(A)') '   ...bolometric heating/cooling turned off (jump).'
         return
      end if
      cy = ysave; par%bolo_heat_cool = bcur      ! ramp it down adaptively
      delta = bcur; fails = 0
      do while (par%bolo_heat_cool .gt. 1.0d-6)
         ysave = cy
         par%bolo_heat_cool = max(0.0d0, bcur - delta)
         if (solve() .eq. 0) then
            bcur = par%bolo_heat_cool; fails = 0; delta = delta*1.5d0
         else
            cy = ysave; par%bolo_heat_cool = bcur
            delta = delta/2.0d0; fails = fails + 1
            if (fails .gt. 12) then
               write(*,'(A)') '   ...turn_off_bolo failed.'
               turn_off_bolo = 1; return
            end if
         end if
      end do
      par%bolo_heat_cool = 0.0d0
      if (solve() .ne. 0) turn_off_bolo = 1
      write(*,'(A)') '   ...bolometric heating/cooling turned off (ramped).'
      end function turn_off_bolo

      integer function converge_mol_atomic()
      ! Port of relax_wrapper.converge_mol_atomic_transition (the path used in
      ! the ramp). With the molecular layer already off (bolo<=0) it is a no-op.
      ! Otherwise, if the molecular/atomic transition has reached the base (the
      ! wind launches within ~10 cells of R_min, i.e. R_min is inside the wind),
      ! turn the bolometric layer off; this is what lets a far-from-seed,
      ! strongly-bound planet (whose solution is atomic, bolo=0) keep ramping.
      integer :: di
      converge_mol_atomic = 0
      if (par%bolo_heat_cool .le. 0.0d0) return
      di = erf_drop_index()
      if (di .le. 10) then
         converge_mol_atomic = turn_off_bolo()
      end if
      end function converge_mol_atomic

      integer function ramp_tidal(goal)
      ! Ramp par%tidalforce from its current value (the seed's, normally 1) to
      ! `goal` (0 for a spherical EXHALE domain), re-solving in small steps --
      ! a port of the reference turn_off_tidal_grav. A discrete 1->0 jump would
      ! disrupt the tidally-converged seed and stall the first relaxation, so it
      ! must be ramped. Call after setup_indices_scales, before ramp_to.
      real*8, intent(in) :: goal
      real*8  :: cur, step, ysave(ne, m)
      integer :: fails
      ramp_tidal = 0
      cur = par%tidalforce
      if (abs(cur-goal) .lt. 1.0d-9) return
      step = 0.1d0; fails = 0
      do while (abs(par%tidalforce - goal) .gt. 1.0d-9)
         ysave = cy
         if (goal .lt. cur) then
            par%tidalforce = max(goal, cur - step)
         else
            par%tidalforce = min(goal, cur + step)
         end if
         if (solve() .eq. 0) then
            cur = par%tidalforce; fails = 0
         else
            cy = ysave; par%tidalforce = cur
            step = step/2.0d0; fails = fails + 1
            if (fails .gt. 14) then
               write(*,'(A)') '   ...ramp_tidal failed.'
               ramp_tidal = 1; return
            end if
         end if
      end do
      write(*,'(A,F4.1)') '   ...ramped tidalforce -> ', goal
      end function ramp_tidal

      subroutine load_seed(wsfile)
      ! parse a windsoln header (#plnt_prms/#bcs/#phys_prms/#tech/#flags)
      ! into wae_par + wae_spectrum copies, and load the first M data rows
      ! into cy / wae_x. BCs stay fixed at these seed values (static_bcs).
      character(len=*), intent(in) :: wsfile
      character(len=8192) :: line
      integer :: u, ios, row, j
      real*8  :: pp(6), bb(4+2*nsp+2), tc(4), fl(4), cval(10)
      wae_cs0_val = sqrt(wae_K*T0/wae_MH)   ! CS0 (needed for H0)
      open(newunit=u, file=wsfile, status='old', action='read')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (index(line,'#plnt_prms:') .gt. 0) then
            call csv_after(line, pp, 6)
            par%Mp=pp(1); par%Rp=pp(2); par%Mstar=pp(3)
            par%semimajor=pp(4); par%Ftot=pp(5); par%Lstar=pp(6)
         else if (index(line,'#bcs:') .gt. 0) then
            call csv_after(line, bb, 4+2*nsp+2)
            par%Rmin=bb(1); par%Rmax=bb(2); par%rho_rmin=bb(3); par%T_rmin=bb(4)
            do j = 1, nsp
               par%Ys_rmin(j) = bb(4+j); par%Ncol_sp(j) = bb(4+nsp+j)
            end do
            par%erf_drop(1)=bb(4+2*nsp+1); par%erf_drop(2)=bb(4+2*nsp+2)
         else if (index(line,'#phys_prms:') .gt. 0) then
            call phys_after(line)
         else if (index(line,'#tech:') .gt. 0) then
            call csv_after(line, tc, 4)
            par%breezeparam=tc(1); par%rapidity=tc(2); par%erfn=tc(3); par%mach_limit=tc(4)
         else if (index(line,'#flags:') .gt. 0) then
            call csv_after(line, fl, 4)
            par%lyacool=nint(fl(1)); par%tidalforce=fl(2)
            par%bolo_heat_cool=fl(3); par%integrate_outward=nint(fl(4))
         end if
      end do
      ! derived
      par%H0 = wae_cs0_val**2 * par%Rp**2 / wae_G / par%Mp
      wae_Ftot = par%Ftot
      wae_HX = par%HX
      wae_atomic_mass = par%atomic_mass
      ! data rows -> cy, wae_x (first M)
      rewind(u); row = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         if (row .ge. m) exit
         call csv_after_pos(line, cval, 10, 0)
         row = row + 1
         wae_x(row)=cval(9)
         cy(1,row)=cval(3); cy(2,row)=cval(10); cy(3,row)=cval(2); cy(4,row)=cval(4)
         cy(5,row)=cval(5); cy(6,row)=cval(6); cy(7,row)=cval(7); cy(8,row)=cval(8)
      end do
      close(u)
      end subroutine load_seed

      subroutine dump_seed(wsfile)
      ! Inverse of load_seed: write the current converged relaxation solution
      ! (par + cy + wae_x + the loaded spectrum) to a windsoln/seed CSV in the
      ! same format as inputdata/windae_seed.csv and the Broome+2025 grid files,
      ! so it can be re-loaded by load_seed and banked as a new starting seed
      ! (e.g. the first spherical, tidalforce=0 solution, or one close to a hard
      ! target). The #scales/#spec_prms/## spectrum block is reconstructed for
      ! file compatibility; load_seed ignores it (the spectrum is read
      ! separately from windae_spectrum.inp).
      character(len=*), intent(in) :: wsfile
      integer :: u, row, j
      character(len=4096) :: vars
      real*8 :: r_phys
      ! fixed spectral-window constants, identical across all wind-ae seeds
      real*8, parameter :: WL(6) = (/ 5.99999999999999946d-08,             &
         9.11999999999999908d-06, 5.99999999999999946d-08,                 &
         9.11999999999999908d-06, 5.99999999999999946d-08,                 &
         9.11999999999999908d-06 /)
      open(newunit=u, file=wsfile, status='replace', action='write')
      ! ---- header ----
      write(u,'(A,I0)') '#nspecies: ', nsp
      vars = '#vars: r,rho,v,T'
      do j = 1, nsp
         vars = trim(vars)//',Ys_'//trim(par%species(j))
      end do
      do j = 1, nsp
         vars = trim(vars)//',Ncol_'//trim(par%species(j))
      end do
      vars = trim(vars)//',q,z'
      write(u,'(A)') trim(vars)
      write(u,'(A)', advance='no') '#scales: '
      write(u,'(*(ES24.17E2,:,","))') par%Rp, RHO0, wae_cs0_val, T0,        &
         1.0d0, 1.0d0, NCOL0, NCOL0, 1.0d0, par%Rp
      write(u,'(A)', advance='no') '#plnt_prms: '
      write(u,'(*(ES24.17E2,:,","))') par%Mp, par%Rp, par%Mstar,           &
         par%semimajor, par%Ftot, par%Lstar
      write(u,'(A)', advance='no') '#phys_prms: '
      do j = 1, nsp
         write(u,'(ES24.17E2,A)', advance='no') par%HX(j), ','
      end do
      do j = 1, nsp
         write(u,'(A,A)', advance='no') trim(par%species(j)), ','
      end do
      do j = 1, nsp
         write(u,'(ES24.17E2,A)', advance='no') par%atomic_mass(j), ','
      end do
      write(u,'(ES24.17E2)') par%molec_adjust
      write(u,'(A)', advance='no') '#bcs: '
      write(u,'(*(ES24.17E2,:,","))') par%Rmin, par%Rmax, par%rho_rmin,    &
         par%T_rmin, (par%Ys_rmin(j), j=1,nsp), (par%Ncol_sp(j), j=1,nsp), &
         par%erf_drop(1), par%erf_drop(2)
      write(u,'(A)', advance='no') '#tech: '
      write(u,'(*(ES24.17E2,:,","))') par%breezeparam, par%rapidity,       &
         par%erfn, par%mach_limit
      write(u,'(A,I0,A,F7.5,A,F7.5,A,I0)') '#flags: ', par%lyacool, ',',   &
         par%tidalforce, ',', par%bolo_heat_cool, ',', par%integrate_outward
      write(u,'(A)') '#add_prms: 0'
      ! spectrum block (reconstructed; load_seed reads windae_spectrum.inp)
      write(u,'(A,I0,A,I0,A)', advance='no') '#spec_prms: ', wae_npts, ',', &
         nsp, ',2009-01-01,scaled-solar,full,'
      do j = 1, 6
         write(u,'(ES24.17E2,A)', advance='no') WL(j), ','
      end do
      write(u,'(*(ES24.17E2,:,","))') (wae_ion_pot(j), j=1,nsp)
      write(u,'(A)') '# $hc/\lambda_i$, $w_i\Phi_{\lambda_i}/F_{tot}$, '//  &
         '$\sigma_{\lambda_i,HI}$, $\sigma_{\lambda_i,HeI}$'
      do row = 1, wae_npts
         write(u,'(A)', advance='no') '## '
         write(u,'(*(ES24.17E2,:,","))') wae_hc_over_wl(row),              &
            wae_wPhi_wl(row), (wae_sigma_wl(row,j), j=1,nsp)
      end do
      ! ---- data rows: first m relaxation points (inverse of load_seed) ----
      do row = 1, m
         r_phys = wae_x(row)*cy(2,row) + par%Rmin
         write(u,'(*(ES24.17E2,:,","))') r_phys, cy(3,row), cy(1,row),     &
            cy(4,row), (cy(4+j,row), j=1,nsp), (cy(4+nsp+j,row), j=1,nsp),  &
            wae_x(row), cy(2,row)
      end do
      close(u)
      write(*,'(A,A)') '   (dump_seed) wrote converged seed -> ', trim(wsfile)
      end subroutine dump_seed

      subroutine pick_nearest_seed(Mp_t, Rp_t, Ftot_t, griddir, path_out)
      ! Scan <griddir>/manifest.csv (one "path,Mp,Rp,Ftot" row per converged
      ! grid solution) and return the path whose (log Mp, log Rp, log Ftot) is
      ! nearest the target planet, so the static-/self-consistent-BC ramp starts from a close
      ! solution instead of the single distant fiducial seed. Falls back to
      ! inputdata/windae_seed.csv (with a warning) if the manifest is missing.
      real*8, intent(in) :: Mp_t, Rp_t, Ftot_t
      character(len=*), intent(in)  :: griddir
      character(len=*), intent(out) :: path_out
      character(len=4096) :: mfile, line, bestpath, p
      integer :: u, ios, c1, n
      real*8 :: Mp_g, Rp_g, Ft_g, d, dbest
      ! (Mp,Rp) set the wind's gravity/structure and are hard to continue, so
      ! they dominate the pick; Ftot ramps easily, so it is only a mild
      ! tiebreaker (e.g. choosing between the hi/lo-flux copy of one Mp,Rp node).
      real*8, parameter :: WF = 0.02d0
      mfile = trim(griddir)//'/manifest.csv'
      open(newunit=u, file=mfile, status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
         write(*,'(A)') '   (pick_nearest_seed) no manifest at '//trim(mfile)// &
                        '; falling back to inputdata/windae_seed.csv'
         path_out = 'inputdata/windae_seed.csv'
         return
      end if
      dbest = huge(1.0d0); bestpath = ''; n = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (len_trim(line) .eq. 0) cycle
         c1 = index(line, ',')
         if (c1 .le. 1) cycle
         p = line(1:c1-1)
         read(line(c1+1:),*,iostat=ios) Mp_g, Rp_g, Ft_g
         if (ios .ne. 0) cycle
         if (Mp_g .le. 0.0d0 .or. Rp_g .le. 0.0d0) cycle
         d = (log(Mp_g/Mp_t))**2 + (log(Rp_g/Rp_t))**2
         if (Ftot_t .gt. 0.0d0 .and. Ft_g .gt. 0.0d0)                       &
            d = d + WF*(log(Ft_g/Ftot_t))**2
         if (d .lt. dbest) then
            dbest = d; bestpath = p
         end if
         n = n + 1
      end do
      close(u)
      if (len_trim(bestpath) .eq. 0) then
         write(*,'(A)') '   (pick_nearest_seed) empty manifest; falling back '// &
                        'to inputdata/windae_seed.csv'
         path_out = 'inputdata/windae_seed.csv'
      else
         path_out = trim(bestpath)
         write(*,'(A,I0,A)')  '   (pick_nearest_seed) scanned ', n,            &
                        ' grid solutions'
         write(*,'(A,A)')     '   (pick_nearest_seed) nearest seed -> ',       &
                        trim(path_out)
         write(*,'(A,ES9.2)') '   (pick_nearest_seed) log-distance^2 = ', dbest
      end if
      end subroutine pick_nearest_seed

      ! ---- parsing helpers ----
      subroutine csv_after(line, vals, n)
      character(len=*), intent(in) :: line
      integer, intent(in) :: n
      real*8, intent(out) :: vals(n)
      call csv_after_pos(line, vals, n, index(line,':'))
      end subroutine csv_after

      subroutine csv_after_pos(line, vals, n, after)
      character(len=*), intent(in) :: line
      integer, intent(in) :: n, after
      real*8, intent(out) :: vals(n)
      integer :: p, jj, b, L
      character(len=64) :: tok
      p = after + 1; L = len_trim(line)
      do jj = 1, n
         do while (p .le. L .and. line(p:p) .eq. ' ')
            p = p + 1
         end do
         b = index(line(p:L), ',')
         if (b .eq. 0) then
            tok = adjustl(line(p:L)); p = L + 1
         else
            tok = adjustl(line(p:p+b-2)); p = p + b
         end if
         read(tok,*) vals(jj)
      end do
      end subroutine csv_after_pos

      subroutine phys_after(line)
      ! #phys_prms: HX(1..nsp), name(1..nsp), amass(1..nsp), molec_adjust
      character(len=*), intent(in) :: line
      integer :: p, jj, b, L
      character(len=64) :: tok
      p = index(line,':') + 1; L = len_trim(line)
      do jj = 1, nsp
         call next_tok(line, p, L, tok); read(tok,*) par%HX(jj)
      end do
      do jj = 1, nsp
         call next_tok(line, p, L, tok); par%species(jj) = trim(tok)
      end do
      do jj = 1, nsp
         call next_tok(line, p, L, tok); read(tok,*) par%atomic_mass(jj)
      end do
      call next_tok(line, p, L, tok); read(tok,*) par%molec_adjust
      ! Z/N_e are not in the windsoln header; assume H/He
      par%Z(1)=1; par%Z(2)=2; par%N_e(1)=1; par%N_e(2)=2
      end subroutine phys_after

      subroutine next_tok(line, p, L, tok)
      character(len=*), intent(in) :: line
      integer, intent(inout) :: p
      integer, intent(in) :: L
      character(len=*), intent(out) :: tok
      integer :: b
      do while (p .le. L .and. line(p:p) .eq. ' ')
         p = p + 1
      end do
      b = index(line(p:L), ',')
      if (b .eq. 0) then
         tok = adjustl(line(p:L)); p = L + 1
      else
         tok = adjustl(line(p:p+b-2)); p = p + b
      end if
      end subroutine next_tok

      end module wae_continuation
