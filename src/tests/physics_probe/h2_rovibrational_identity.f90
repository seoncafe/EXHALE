      program h2_rovibrational_identity
      ! The thermodynamic identities that tie the H2 partition function used
      ! by the chemistry to the H2 internal energy used by the caloric
      ! equation of state, and the equilibrium constant that follows from
      ! both.
      !
      ! Origin: the H2 block of docs/audit_20260905/audit_probe.f90, which
      ! prints the two energies and their difference.  This program asserts
      ! the identities and exits nonzero when one is violated.
      !
      ! Production routines exercised:
      !   h2_partition_function and h2_rovibrational_sum of
      !     src/modules/states/caloric_eos.f90 -- the one Boltzmann sum over
      !     the one H2 level ladder, which is what the equilibrium constant
      !     is built from and what the internal energy is a moment of;
      !   h2_rovibrational_energy_and_heat_capacity of the same module --
      !     the tabulated accessor the energy equation actually calls;
      !   keq_H_H_to_H2 and h2_thermochemistry_init of
      !     src/modules/lower_atmosphere/mol_rates.f90 -- the tabulated
      !     H + H <-> H2 equilibrium constant, whose reverse is R12;
      !   equilibrium_constant_conc of
      !     src/modules/lower_atmosphere/oxygen_rates.f90 -- the same
      !     equilibrium constant from the NIST-JANAF Shomate enthalpies and
      !     entropies, an independent source.
      !
      ! REFERENCE.  For any internal partition function Q(T) the mean
      ! internal energy of the same level set is
      !     u/k = T^2 d ln Q / dT = T d ln Q / d ln T,
      ! exactly, with no free constant.  Two routines that claim to describe
      ! the same molecule must satisfy it; a violation means they carry
      ! different level sets or different degeneracies, and the equilibrium
      ! constant and the energy equation then disagree about how much energy
      ! an H2 molecule holds.  The logarithmic derivative is taken with a
      ! centered difference of relative step 1e-5, whose own truncation error
      ! is of order 1e-10.
      !
      ! One further exact statement, used by the third group: with
      !     ln K_eq = ln Lam(2m_H) + ln Q + D0 hc/(kT) - 2 ln(4 Lam(m_H)),
      !     Lam(m) = (2 pi m kB T/h^2)^(3/2),
      ! the equilibrium constant of the tabulated production path must
      ! reproduce that closed form to the tabulation error mol_rates states
      ! (4.7e-4 in ln K_eq at the cold end, less above).
      use global_parameters, only: T0, caloric_eos_monatomic,             &
                                   kb_erg, hp_erg, c_light, mu, pi
      use caloric_eos, only: h2_rovibrational_energy_and_heat_capacity,   &
                             h2_rovibrational_sum, h2_partition_function
      use mol_rates,   only: keq_H_H_to_H2, h2_thermochemistry_init
      use oxygen_rates, only: equilibrium_constant_conc, ith_H, ith_H2
      use assertion_report
      implicit none

      real*8, parameter :: temps(5) =                                     &
           (/ 300.0d0, 1000.0d0, 2000.0d0, 4000.0d0, 8000.0d0 /)
      ! The three temperatures at which the equilibrium constant is checked.
      ! They lie inside the 100-20000 K tabulation of mol_rates and inside
      ! the 6000 K ceiling of the Shomate table of oxygen_rates.
      real*8, parameter :: temps_keq(3) = (/ 1000.0d0, 2000.0d0, 4000.0d0 /)
      ! MEASURED 2026-09-06 on the ladder of Roueff et al. (2019) against the
      ! Huber and Herzberg (1979) model ladder this module set retired: the
      ! ratio Q(observed)/Q(model) at the three temperatures above.  This is
      ! the movement B1 decision 9 adopts, recorded so that a silent return
      ! to the model ladder fails here.  It is a transition record, not an
      ! ownership statement, and goes when the retired model does.
      real*8, parameter :: q_ratio_ref(3) =                               &
           (/ 1.017599d0, 1.036501d0, 1.077285d0 /)
      real*8, parameter :: eps = 1.0d-5
      real*8 :: tk, u_eos, c_eos, u_chem, u_sum, c_sum, q
      real*8 :: hc_k, lam_h, lam_h2, keq_closed, keq_code, keq_janaf
      character(len=48) :: label
      integer :: i

      T0 = 1000.0d0
      caloric_eos_monatomic = .false.
      call h2_thermochemistry_init
      hc_k = hp_erg*c_light/kb_erg          ! K per cm^-1

      ! (1) The chemistry's partition function and the energy the equation of
      ! state carries describe the same molecule.  The energy side is the
      ! tabulated accessor, so this also covers the cubic Hermite table.
      do i = 1, size(temps)
         tk = temps(i)
         call h2_rovibrational_energy_and_heat_capacity(tk, u_eos, c_eos)
         u_chem = tk*(log(h2_partition_function(tk*exp(eps)))             &
                    - log(h2_partition_function(tk*exp(-eps))))/(2.0d0*eps)
         write(label,'(a,i0)') 'h2_rovibrational_energy_at_', int(tk)
         call check_relative(label, u_chem, u_eos, 1.0d-3)
      enddo

      ! (2) The identity itself, against the exact level sum rather than its
      ! table, at the truncation error of the centered difference.
      do i = 1, size(temps)
         tk = temps(i)
         call h2_rovibrational_sum(tk, u_sum, c_sum, q)
         u_chem = tk*(log(h2_partition_function(tk*exp(eps)))             &
                    - log(h2_partition_function(tk*exp(-eps))))/(2.0d0*eps)
         write(label,'(a,i0)') 'h2_partition_function_identity_at_', int(tk)
         call check_relative(label, u_chem, u_sum, 1.0d-8)
      enddo

      ! (3) The equilibrium constant is the free energy of that same Q.
      do i = 1, size(temps_keq)
         tk     = temps_keq(i)
         lam_h  = (2.0d0*pi*mu*kb_erg*tk/hp_erg**2)**1.5d0
         lam_h2 = (2.0d0*pi*2.0d0*mu*kb_erg*tk/hp_erg**2)**1.5d0
         keq_closed = exp(log(lam_h2) + log(h2_partition_function(tk))    &
                        + 36118.11d0*hc_k/tk - 2.0d0*log(4.0d0*lam_h))
         keq_code   = keq_H_H_to_H2(tk)
         write(label,'(a,i0)') 'h2_equilibrium_constant_table_at_', int(tk)
         call check_relative(label, keq_code, keq_closed, 5.0d-4)

         ! Independent source: the NIST-JANAF Shomate enthalpies and
         ! entropies of H and H2 that oxygen_rates carries.  Nothing of the
         ! level ladder enters it.
         keq_janaf = equilibrium_constant_conc((/ ith_H, ith_H /),        &
                                               (/ ith_H2 /), tk)
         ! MEASURED 2026-09-06: the two agree to 1.7e-4 at 1000 K, 4.2e-5
         ! at 2000 K and 2.3e-4 at 4000 K, so 1e-3 is a gate on the level
         ! ladder and not on the Shomate fit's own resolution.
         write(label,'(a,i0)') 'h2_equilibrium_constant_janaf_at_', int(tk)
         call check_relative(label, keq_code, keq_janaf, 1.0d-3)

         ! The movement of decision 9, against the retired model ladder.
         write(label,'(a,i0)') 'h2_partition_function_movement_at_', int(tk)
         call check_relative(label,                                       &
              h2_partition_function(tk)/q_huber_herzberg(tk),             &
              q_ratio_ref(i), 1.0d-4)
      enddo

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'h2_rovibrational_identity: ',               &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'h2_rovibrational_identity: all assertions passed'

      contains

      double precision function q_huber_herzberg(T) result(q)
      ! The H2 internal partition function as mol_rates built it until
      ! 2026-09-06: a rigid rotor with a vibration-rotation coupling and no
      ! centrifugal distortion, from the Huber and Herzberg (1979) constants,
      ! summed to the level where a term energy passes D0.  It exists ONLY so
      ! that the movement adopted by B1 decision 9 is recorded by a running
      ! assertion instead of by a comment; no production path reaches it.
      real*8, intent(in) :: T
      real*8, parameter :: D0_cm = 36118.11d0, we = 4401.213d0,           &
                           wexe = 121.336d0, be = 60.853d0, ae = 3.062d0
      real*8  :: hck, g_v0, g_v, b_v, v_half, e_vj, w_ns
      integer :: iv, jr
      hck  = hp_erg*c_light/kb_erg
      g_v0 = 0.5d0*we - 0.25d0*wexe
      q    = 0.0d0
      iv   = 0
      do
         v_half = dble(iv) + 0.5d0
         b_v    = be - ae*v_half
         g_v    = we*v_half - wexe*v_half*v_half - g_v0
         if (b_v .le. 0.0d0 .or. g_v .gt. D0_cm) exit
         jr = 0
         do
            e_vj = g_v + b_v*dble(jr)*dble(jr+1)
            if (e_vj .gt. D0_cm) exit
            if (mod(jr,2) .eq. 1) then
               w_ns = 3.0d0
            else
               w_ns = 1.0d0
            endif
            q  = q + w_ns*dble(2*jr+1)*exp(-e_vj*hck/T)
            jr = jr + 1
         enddo
         iv = iv + 1
      enddo
      end function q_huber_herzberg

      end program h2_rovibrational_identity
