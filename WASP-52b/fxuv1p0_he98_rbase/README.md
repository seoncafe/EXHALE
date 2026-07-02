# fxuv1p0_he98_rbase — corrected-base-radius rerun (2026-07-02)

Same as `../fxuv1p0_he98/` except the 1-ubar base radius corrected from the
transit-radius shortcut (1.27 R_J) to the Tier-1 equilibrium-column value
(1.43718 R_J; see docs/lower_atmosphere_coupling.pdf, Examples).

Newton-converged result vs the old base (both rerun with the same binary):
  log Mdot: 11.97 -> 12.14 (x1.5)
  He 10830 line center: 27.11% -> 27.56% (x1.02); EW 0.713 -> 0.722 A (x1.01)
=> the base misplacement biases Mdot, not the He 10830 observable (the tau=1
   line-forming surface stays nearly fixed in physical units).
