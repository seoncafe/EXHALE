module energy_semi_implicit
   use global_parameters
   use species_table, only: n_mion, mion_fsp,                        &
                            isp_HI, isp_HII, isp_HeI, isp_HeII,       &
                            isp_HeIII, isp_HeTR,                      &
                            isp_H2, isp_H2p, isp_H3p, isp_HeHp,       &
                            isp_OH, isp_H2O, isp_CO
   use utils
   use utils_ion_eq
   use caloric_eos, only: caloric_mixture_active,                      &
                          internal_energy_per_particle,                &
                          heat_capacity_per_particle,                  &
                          energy_density_from_pressure

   implicit none

   ! Activations of the temperature floor of the semi-implicit energy update
   ! (the clamp in the Newton loop below). A cell that lands on the floor has
   ! NOT converged to a physical temperature: the update overshot and the
   ! clamp absorbed it, and at 0.01 T0 the chemical network is frozen, so
   ! almost any composition satisfies its reaction balance there
   ! (docs/Update_EXHALE.md section 113). Left silent, a run resting on the
   ! floor looks exactly like a run that resolved a cold layer, which is why
   ! these are counted and reported at the end of every run, as the
   ! non-root acceptances of section 113 are.
   !   hits        total activations over the run
   !   first/last  first and last step on which the floor was reached
   !   cell_hits   activations of each cell, so the number of DISTINCT cells
   !               that ever touched the floor can be reported
   integer, save :: n_energy_floor_hits = 0
   integer, save :: energy_floor_first_step = -1
   integer, save :: energy_floor_last_step  = -1
   integer, allocatable, save :: energy_floor_cell_hits(:)

contains

   integer function n_energy_floor_cells()
   ! Number of DISTINCT cells that have reached the temperature floor of the
   ! energy update at least once. Written as a module function because the
   ! marching loop of EXHALE_main carries a local integer named "count",
   ! which shadows the Fortran intrinsic of the same name at that call site.
   integer :: j
   n_energy_floor_cells = 0
   if (.not. allocated(energy_floor_cell_hits)) return
   do j = lbound(energy_floor_cell_hits,1), ubound(energy_floor_cell_hits,1)
      if (energy_floor_cell_hits(j) .gt. 0)                                &
         n_energy_floor_cells = n_energy_floor_cells + 1
   enddo
   end function n_energy_floor_cells

   subroutine solve_energy_semi_implicit(u, W, dt, heat, cool, f_sp, step)
      ! Input/Output variables
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: W
      ! Cell-by-cell pseudo-time steps (uniform = global dt unless
      ! "Time stepping: Local"; the cell solve below is local anyway).
      real*8, dimension(1-Ng:N+Ng), intent(in) :: dt
      ! Heating of the state this step starts from, from ioniz_eq. It is
      ! held fixed across the temperature update (the photoheating rates do
      ! not follow T within one step).
      real*8, dimension(1-Ng:N+Ng), intent(in) :: heat
      ! On entry: the cooling ioniz_eq evaluated at T_old and at the
      ! post-sweep composition -- the same composition and the same
      ! temperature this routine would evaluate it at, so it is taken as the
      ! first cooling of the Newton iteration instead of being recomputed.
      ! On exit: the cooling at the temperature this step lands on.
      real*8, dimension(1-Ng:N+Ng), intent(inout) :: cool
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      ! Marching step index, for the temperature-floor report only.
      integer, intent(in) :: step

      ! Local arrays for dimensional calculations
      real*8, dimension(1-Ng:N+Ng) :: T_K, T_trial, T_perturbed, T_old
      real*8, dimension(1-Ng:N+Ng) :: nhi, nhii, nhei, nheii, nheiii, nheiTR
      ! Metal ion densities (canonical species_table order)
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
      ! Molecular densities: H2, H2+, H3+, HeH+ -- adimensional (nmol_l, the
      ! units calc_ne is called with here) and cgs (nmol_dim, the units
      ! eval_cool takes alongside nhi..nm)
      real*8, dimension(1-Ng:N+Ng,4) :: nmol_l, nmol_dim
      ! Oxygen-chemistry carriers OH, H2O, CO -- adimensional; zero without
      ! the option, so the particle count is then unchanged.
      ! Oxygen carriers OH, H2O, CO -- adimensional (nox_l, the units
      ! calc_ntot is called with here) and cgs (nox_dim, the units eval_cool
      ! needs for the H2O and CO infrared bands)
      real*8, dimension(1-Ng:N+Ng,3) :: nox_l, nox_dim
      real*8, dimension(1-Ng:N+Ng) :: rchiiB, rcheiiB, rcheiiiB
      real*8, dimension(1-Ng:N+Ng) :: a_ion_HI, a_ion_HeI, a_ion_HeII
      ! Metal rates for each ion returned by eval_cool but unused here
      real*8, dimension(1-Ng:N+Ng,n_mion) :: rec_m, aion_m
      real*8, dimension(1-Ng:N+Ng) :: cool_dim, cool_trial, cool_perturbed
      real*8, dimension(1-Ng:N+Ng) :: F, dF_dT, dC_dT, delta_T
      real*8, dimension(1-Ng:N+Ng) :: rho, v, p
      real*8, dimension(1-Ng:N+Ng) :: zero_arr

      ! Adimensional number densities
      real*8, dimension(1-Ng:N+Ng) :: ne_ad, n_tot_ad
      real*8 :: c_factor, c_energy
      integer :: j, iter, im

      ! Zero array for the calls that take no metal ions
      zero_arr = 0.0d0

      ! Extract primitive variables
      rho = W(1,:)
      v   = W(2,:)
      p   = W(3,:)

      ! Dimensional species densities (n0 is the density normalization)
      nhi  = rho*f_sp(:,isp_HI)*n0
      nhii = rho*f_sp(:,isp_HII)*n0
      if (thereis_He) then
         nhei   = rho*f_sp(:,isp_HeI)*n0
         nheii  = rho*f_sp(:,isp_HeII)*n0
         nheiii = rho*f_sp(:,isp_HeIII)*n0
         if (thereis_HeITR) then
            nheiTR = rho*f_sp(:,isp_HeTR)*n0
         else
            nheiTR = 0.0d0
         endif
      else
         nhei   = 0.0d0
         nheii  = 0.0d0
         nheiii = 0.0d0
         nheiTR = 0.0d0
      endif
      ! Metal ion densities in canonical species_table order
      ! (mion_fsp = [7..18], so nm(:,im) reproduces each ion's
      ! rho*f_sp(:,col)*n0 expressions bit-for-bit).
      do im = 1,n_mion
         nm(:,im) = rho*f_sp(:,mion_fsp(im))*n0
      enddo
      ! Molecular densities (adimensional, matching the nm/n0 units used below).
      ! Zero for non-molecular runs (f_sp molecular columns are zero), so
      ! calc_ne/calc_ntot add exactly zero and the atomic result is bitwise
      ! unchanged; for molecular runs they restore the neutral-H2 particle count
      ! and the molecular-ion electrons to n_tot/ne (first-order near a
      ! molecular base, where H2 dominates the particle budget).
      nmol_l(:,1) = rho*f_sp(:,isp_H2)
      nmol_l(:,2) = rho*f_sp(:,isp_H2p)
      nmol_l(:,3) = rho*f_sp(:,isp_H3p)
      nmol_l(:,4) = rho*f_sp(:,isp_HeHp)
      nmol_dim    = nmol_l*n0

      ! Compute adimensional total and electron densities for T calculation
      ! (nm/n0 = adimensional metal densities; adds the metal electrons and
      ! nuclei under the eos_metals policy)
      call calc_ne(rho*f_sp(:,isp_HII), rho*f_sp(:,isp_HeII), rho*f_sp(:,isp_HeIII), ne_ad, nm/n0, nmol_l)
      ! The oxygen-chemistry carriers are gas particles as well, and their
      ! oxygen and carbon nuclei have been removed from the metal ion
      ! columns by the ionization solve; leaving them out would make this
      ! energy solve and ioniz_eq disagree about the particle count of the
      ! same state. They are neutral, so calc_ne is unaffected.
      nox_l = 0.0d0
      if (thereis_oxychem) then
         nox_l(:,1) = rho*f_sp(:,isp_OH)
         nox_l(:,2) = rho*f_sp(:,isp_H2O)
         nox_l(:,3) = rho*f_sp(:,isp_CO)
      endif
      nox_dim = nox_l*n0
      ! The He 2^3S column is inside the HeI column (bsp_is_excited_level),
      ! so it is not passed and the triplet branch disappears with it.
      if (thereis_He) then
         if (thereis_oxychem) then
            call calc_ntot(rho*f_sp(:,isp_HI), rho*f_sp(:,isp_HII), rho*f_sp(:,isp_HeI), &
                        rho*f_sp(:,isp_HeII), rho*f_sp(:,isp_HeIII), n_tot_ad, nm/n0, nmol_l, nox_l)
         else
         call calc_ntot(rho*f_sp(:,isp_HI), rho*f_sp(:,isp_HII), rho*f_sp(:,isp_HeI), &
                        rho*f_sp(:,isp_HeII), rho*f_sp(:,isp_HeIII), n_tot_ad, nm/n0, nmol_l)
         endif
      else
         call calc_ntot(rho*f_sp(:,isp_HI), rho*f_sp(:,isp_HII), zero_arr, &
                        zero_arr, zero_arr, n_tot_ad, nm/n0, nmol_l)
      endif

      ! Old temperature (adimensional). The pressure was formed as
      ! (n_tot + n_e) T from the state ioniz_eq returned and the same T it
      ! was evaluated at, and the particle counts here come from that same
      ! state through the same calc_ne/calc_ntot, so this reproduces the
      ! ioniz_eq temperature to round-off.
      T_old = p / (n_tot_ad + ne_ad)

      ! Initial guess for Newton-Raphson is the old temperature
      T_trial = T_old

      ! 1. Cooling at T_old. This is exactly what ioniz_eq returned: it
      ! evaluates the cooling after its cell sweep, at the composition it
      ! hands over -- the composition f_sp carries here -- and at T_old.
      ! Re-evaluating it would repeat that traversal for the same number.
      cool_trial = cool

      ! 2. Evaluate cooling at perturbed temperature to get derivative.
      ! nheiTR adds the He 2^3S channels (collisional ionization, 10830 A
      ! excitation, 2^3S -> 2^1S/2^1P conversion) to the cooling that acts on
      ! the temperature update; the array is zero when the triplet is off.
      ! nmol_dim gives eval_cool the same electron density the equilibrium
      ! solver uses (molecular ions included); zero for an atomic run.
      delta_T = max(1.0d-5, 1.0d-5 * T_trial)
      T_perturbed = T_trial + delta_T
      T_K = T_perturbed * T0
      call eval_cool(T_K, nhi, nhii, nhei, nheii, nheiii, nm, &
                     rchiiB, rcheiiB, rcheiiiB, rec_m, &
                     a_ion_HI, a_ion_HeI, a_ion_HeII, aion_m, &
                     cool_dim, nheiTR = nheiTR, nmol = nmol_dim,   &
                     nox = nox_dim)
      cool_perturbed = cool_dim / q0

      dC_dT = (cool_perturbed - cool_trial) / delta_T

      ! Newton-Raphson iteration loop (2 iterations, keeping dC_dT constant)
      do iter = 1, 2
         ! 3. Compute residual F and derivative dF/dT, then update T_trial
         do j = 1-Ng, N+Ng
            ! The equation being solved is the ENERGY balance,
            !    e(T_new) - e(T_old) = dt (heat - cool(T_new)) ,
            ! per (n_tot + n_e).  With a constant heat capacity that is the
            ! temperature form below, rescaled by (gamma - 1); with the
            ! caloric EOS the internal energy is not proportional to T and
            ! only the energy form is the balance.  The atomic branch keeps
            ! the temperature form verbatim so its arithmetic is unchanged.
            if (caloric_mixture_active) then
               c_energy = dt(j) / (n_tot_ad(j) + ne_ad(j))
               F(j) = internal_energy_per_particle(j, T_trial(j))       &
                    - internal_energy_per_particle(j, T_old(j))         &
                    - c_energy * (heat(j) - cool_trial(j))
            else
            c_factor = dt(j) * (gamma_ad - 1.0d0) / (n_tot_ad(j) + ne_ad(j))
            F(j) = T_trial(j) - T_old(j) - c_factor * (heat(j) - cool_trial(j))
            endif
            ! Both-branch-stable damping derivative. The exact Newton derivative
            ! is dF/dT = 1 + c_factor*dC_dT; on the FALLING cooling branch
            ! (dC_dT < 0, i.e. T past the ~2e4 K cooling peak) that can drop
            ! below 1 (or negative), so the previous safeguard max(0,dC_dT)
            ! clamped it to dF_dT = 1 -- a fully explicit, UNDAMPED step that
            ! lets a hot cell run away/oscillate (the source of the
            ! metal-rich-wind base/front limit cycle that floors du). Using
            ! |dC_dT| keeps dF_dT >= 1 AND adds damping (dF_dT > 1) on the
            ! falling branch, stabilizing it. On the rising branch dC_dT > 0 so
            ! |dC_dT| = dC_dT = max(0,dC_dT): IDENTICAL to before, hence a no-op
            ! for any model whose T never crosses the cooling peak. Damping only
            ! affects the iteration path, not the fixed point F = 0, so the
            ! converged steady state is unchanged.
            if (caloric_mixture_active) then
               dF_dT(j) = heat_capacity_per_particle(j, T_trial(j))      &
                        + c_energy * abs(dC_dT(j))
            else
            dF_dT(j) = 1.0d0 + c_factor * abs(dC_dT(j))
            endif

            ! Newton-Raphson step
            T_trial(j) = T_trial(j) - F(j) / dF_dT(j)

            ! Apply temperature floor (guarding against negative values).
            ! Counted: reaching it means this cell's energy update did not
            ! land on a physical temperature (see the counters above).
            if (T_trial(j) .lt. 0.01d0) then
               T_trial(j) = 0.01d0
               if (.not. allocated(energy_floor_cell_hits)) then
                  allocate(energy_floor_cell_hits(1-Ng:N+Ng))
                  energy_floor_cell_hits = 0
               endif
               n_energy_floor_hits = n_energy_floor_hits + 1
               energy_floor_cell_hits(j) = energy_floor_cell_hits(j) + 1
               if (energy_floor_first_step .lt. 0)                        &
                  energy_floor_first_step = step
               energy_floor_last_step = step
            endif
         end do

         ! Update cool_trial for the second iteration
         if (iter .eq. 1) then
            T_K = T_trial * T0
            call eval_cool(T_K, nhi, nhii, nhei, nheii, nheiii, nm, &
                           rchiiB, rcheiiB, rcheiiiB, rec_m, &
                           a_ion_HI, a_ion_HeI, a_ion_HeII, aion_m, &
                           cool_dim, nheiTR = nheiTR, nmol = nmol_dim,   &
                     nox = nox_dim)
            cool_trial = cool_dim / q0
         end if
      end do

      ! Update the primitive pressure and conservative energy density
      p = (n_tot_ad + ne_ad) * T_trial
      W(3,:) = p
      if (caloric_mixture_active) then
         do j = 1-Ng, N+Ng
            u(3,j) = 0.5d0*rho(j)*v(j)**2.0                             &
                   + energy_density_from_pressure(j, rho(j), p(j))
         enddo
      else
      u(3,:) = 0.5d0 * rho * v**2.0 + p / (gamma_ad - 1.0d0)
      endif

      ! Set the output cool array (adimensional)
      cool = cool_trial

   end subroutine solve_energy_semi_implicit

end module energy_semi_implicit
