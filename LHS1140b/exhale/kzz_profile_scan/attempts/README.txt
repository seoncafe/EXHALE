Attempts that were superseded, kept because they are the reason the scan is
shaped the way it is.

WHY THERE ARE HALF-DECADE RUNGS (3.16e7, 1.78e7, 1.33e7)

ladder_kzz.log   1e9 -> 1e8 at reservoir He/H = 11.1 went through.  1e8 ->
                 1e7 in one decade left the wind composition drifting at
                 1.7e-1 after the first outer pass and the next steady solve
                 returned info = 2.  Run directory:
                 ../kzz1e7/heh11p1_attempt_decade_step.
ladder_kzz2.log  the same step retried through the half decade 3.16e7, which
                 went through; 3.16e7 -> 1e7 then failed the same way.  That
                 is what put the EXHALE_DIFF_OMEGA retry into run_arm.sh.
ladder_kzz3_stopped_at_1e7.log
ladder_kzz3b.log
ladder_kzz3_stopped_at_1e7_cycle.log
                 three further 3.16e7 -> 1e7 attempts, with the residual
                 target loosened to 1.0e-3 and then 3.0e-3 and with
                 EXHALE_DIFF_OMEGA at 0.25 and 0.125.  With omega = 0.125 and
                 a 3.0e-3 target the steady solve does succeed
                 (info = 0, ||R|| = 1.77e-3) but the composition outer loop
                 does not: the drift oscillates between 1.2e-1 and 2.3e-1
                 over all 20 passes without a trend.  Run directory:
                 ../kzz1e7/heh11p1_attempt_from_3p16e7.
kzz1e7_heh11p1_attempt_from_3p16e7.log
                 the log of that last attempt.

The step that finally worked was 3.16e7 -> 1.78e7 -> 1.33e7 -> 1e7, and the
last two rungs converged at the STANDARD settings (resid_tol 4.0e-4, no
omega override, drift 6e-4).  They converged onto the hot branch, which is
how the two states of the wind were found: the cool state that the ladder
was following does not exist at 1.33e7 and below at that reservoir ratio.

WHY THE RESERVOIR STEP IS CAPPED AT 1.25

search_1p78e7_jumped_branch.log
                 the secant search at K_zz = 1.78e7 took its first step from
                 He/H = 11.1 to 15.383, a factor 1.42, and the new arm came
                 back with EW = 4.106 %A against 0.979 at the seed -- it had
                 crossed to the hot branch (log10 Mdot 7.83 against 7.42).
                 The cap in search_heh.sh was 1.5 at the time and is 1.25
                 now.  The arm itself is kept as ../kzz1p78e7/heh15p383.
