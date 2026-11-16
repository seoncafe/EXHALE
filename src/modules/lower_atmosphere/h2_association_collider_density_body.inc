!  THE H2-EQUIVALENT COLLIDER DENSITY, ONE SOURCE TEXT FOR TWO ARITHMETICS:
!  included by h2_association_collider_density_real (real*8 densities) and
!  h2_association_collider_density_dual (dual densities) of mol_rates.  The
!  physics is stated above h2_association_collider_density_real.
      k_h2  = k3b_H_H_to_H2(T)
      n_eff = positive_part(n_H2)                                            &
            + positive_part(n_HI)*k3b_H_H_to_H2_atomic_H()/k_h2              &
            + positive_part(n_He)*k3b_H_H_to_H2_monatomic(T)/k_h2
