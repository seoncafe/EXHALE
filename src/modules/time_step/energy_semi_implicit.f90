module energy_semi_implicit
   use global_parameters
   use species_table, only: n_mion, mion_fsp,                        &
                            isp_HI, isp_HII, isp_HeI, isp_HeII,       &
                            isp_HeIII, isp_HeTR
   use utils
   use utils_ion_eq

   implicit none

contains

   subroutine solve_energy_semi_implicit(u, W, dt, heat, cool, f_sp)
      ! Input/Output variables
      real*8, dimension(1-Ng:N+Ng,3), intent(inout) :: u
      real*8, dimension(1-Ng:N+Ng,3), intent(inout) :: W
      ! Per-cell pseudo-time steps (uniform = global dt unless
      ! "Time stepping: Local"; the cell solve below is local anyway).
      real*8, dimension(1-Ng:N+Ng), intent(in) :: dt
      real*8, dimension(1-Ng:N+Ng), intent(in) :: heat
      real*8, dimension(1-Ng:N+Ng), intent(inout) :: cool
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp

      ! Local arrays for dimensional calculations
      real*8, dimension(1-Ng:N+Ng) :: T_K, T_trial, T_perturbed, T_old
      real*8, dimension(1-Ng:N+Ng) :: nhi, nhii, nhei, nheii, nheiii, nheiTR
      ! Metal ion densities (canonical species_table order)
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
      real*8, dimension(1-Ng:N+Ng) :: rchiiB, rcheiiB, rcheiiiB
      real*8, dimension(1-Ng:N+Ng) :: a_ion_HI, a_ion_HeI, a_ion_HeII
      ! Per-ion metal rates returned by eval_cool but unused here
      real*8, dimension(1-Ng:N+Ng,n_mion) :: rec_m, aion_m
      real*8, dimension(1-Ng:N+Ng) :: cool_dim, cool_trial, cool_perturbed
      real*8, dimension(1-Ng:N+Ng) :: F, dF_dT, dC_dT, delta_T
      real*8, dimension(1-Ng:N+Ng) :: rho, v, p
      real*8, dimension(1-Ng:N+Ng) :: zero_arr

      ! Adimensional number densities
      real*8, dimension(1-Ng:N+Ng) :: ne_ad, n_tot_ad
      real*8 :: c_factor
      integer :: j, iter, im

      ! Initialize helper zero array
      zero_arr = 0.0d0

      ! Extract primitive variables
      rho = W(:,1)
      v   = W(:,2)
      p   = W(:,3)

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
      ! (mion_fsp = [7..18], so nm(:,im) reproduces the per-ion
      ! rho*f_sp(:,col)*n0 expressions bit-for-bit).
      do im = 1,n_mion
         nm(:,im) = rho*f_sp(:,mion_fsp(im))*n0
      enddo

      ! Compute adimensional total and electron densities for T calculation
      ! (nm/n0 = adimensional metal densities; adds the metal electrons and
      ! nuclei under the eos_metals policy)
      call calc_ne(rho*f_sp(:,isp_HII), rho*f_sp(:,isp_HeII), rho*f_sp(:,isp_HeIII), ne_ad, nm/n0)
      if (thereis_He) then
         if (thereis_HeITR) then
            call calc_ntot(rho*f_sp(:,isp_HI), rho*f_sp(:,isp_HII), rho*f_sp(:,isp_HeI), &
                           rho*f_sp(:,isp_HeII), rho*f_sp(:,isp_HeIII), rho*f_sp(:,isp_HeTR), n_tot_ad, nm/n0)
         else
            call calc_ntot(rho*f_sp(:,isp_HI), rho*f_sp(:,isp_HII), rho*f_sp(:,isp_HeI), &
                           rho*f_sp(:,isp_HeII), rho*f_sp(:,isp_HeIII), zero_arr, n_tot_ad, nm/n0)
         endif
      else
         call calc_ntot(rho*f_sp(:,isp_HI), rho*f_sp(:,isp_HII), zero_arr, &
                        zero_arr, zero_arr, zero_arr, n_tot_ad, nm/n0)
      endif

      ! Old temperature (adimensional)
      T_old = p / (n_tot_ad + ne_ad)

      ! Initial guess for Newton-Raphson is the old temperature
      T_trial = T_old

      ! 1. Evaluate cooling at T_old (initial guess T_trial is T_old)
      T_K = T_trial * T0
      call eval_cool(T_K, nhi, nhii, nhei, nheii, nheiii, nm, &
                     rchiiB, rcheiiB, rcheiiiB, rec_m, &
                     a_ion_HI, a_ion_HeI, a_ion_HeII, aion_m, &
                     cool_dim)
      cool_trial = cool_dim / q0

      ! 2. Evaluate cooling at perturbed temperature to get derivative
      delta_T = max(1.0d-5, 1.0d-5 * T_trial)
      T_perturbed = T_trial + delta_T
      T_K = T_perturbed * T0
      call eval_cool(T_K, nhi, nhii, nhei, nheii, nheiii, nm, &
                     rchiiB, rcheiiB, rcheiiiB, rec_m, &
                     a_ion_HI, a_ion_HeI, a_ion_HeII, aion_m, &
                     cool_dim)
      cool_perturbed = cool_dim / q0

      dC_dT = (cool_perturbed - cool_trial) / delta_T

      ! Newton-Raphson iteration loop (2 iterations, keeping dC_dT constant)
      do iter = 1, 2
         ! 3. Compute residual F and derivative dF/dT, then update T_trial
         do j = 1-Ng, N+Ng
            c_factor = dt(j) * (g - 1.0d0) / (n_tot_ad(j) + ne_ad(j))
            F(j) = T_trial(j) - T_old(j) - c_factor * (heat(j) - cool_trial(j))
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
            dF_dT(j) = 1.0d0 + c_factor * abs(dC_dT(j))

            ! Newton-Raphson step
            T_trial(j) = T_trial(j) - F(j) / dF_dT(j)

            ! Apply temperature floor (guarding against negative values)
            if (T_trial(j) .lt. 0.01d0) T_trial(j) = 0.01d0
         end do

         ! Update cool_trial for the second iteration
         if (iter .eq. 1) then
            T_K = T_trial * T0
            call eval_cool(T_K, nhi, nhii, nhei, nheii, nheiii, nm, &
                           rchiiB, rcheiiB, rcheiiiB, rec_m, &
                           a_ion_HI, a_ion_HeI, a_ion_HeII, aion_m, &
                           cool_dim)
            cool_trial = cool_dim / q0
         end if
      end do

      ! Update the primitive pressure and conservative energy density
      p = (n_tot_ad + ne_ad) * T_trial
      W(:,3) = p
      u(:,3) = 0.5d0 * rho * v**2.0 + p / (g - 1.0d0)

      ! Set the output cool array (adimensional)
      cool = cool_trial

   end subroutine solve_energy_semi_implicit

end module energy_semi_implicit
