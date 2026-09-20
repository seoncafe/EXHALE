! NOTE (2026-09-18, after item D7a): this probe calls cx_add_to_fvec with the
! argument list of the entry text it reviewed. Item D7a gave that routine a
! helium-row orientation argument and a cell-temperature guard, so this file
! no longer compiles against the tree; it is kept as the dated record of the
! review's counterexample, which src/tests/charge_exchange_rows/ now carries
! (the same Si I + He II numbers, RED on the entry text and GREEN after).
program plan_revision_source_probe
   ! Targeted calls to the unchanged production charge-exchange module.
   ! This is not an execution of the full atmosphere or ionization solver.
   use species_table, only: n_melem, iel_Si
   use charge_exchange, only: cx_full, cx_init, cx_set_cell, &
                             cx_metal_base, cx_add_to_fvec, he_h_cx_fvec
   implicit none
   real*8 :: nm0(n_melem), nm1(n_melem), nm2(n_melem), fv(24), rate
   character(len=24) :: narrow_positive, narrow_negative
   character(len=26) :: wide_negative

   write(narrow_positive,'(ES24.17E3)') 2.0d-4
   write(narrow_negative,'(ES24.17E3)') -2.0d-4
   write(wide_negative,'(ES26.17E3)') -2.0d-4
   print '(A,A,A)', 'D1 ES24.17E3 positive: [', narrow_positive, ']'
   print '(A,A,A)', 'D1 ES24.17E3 negative: [', narrow_negative, ']'
   print '(A,A,A)', 'D1 ES26.17E3 negative: [', wide_negative, ']'

   ! With only Si I and He II present, the only active event is
   ! Si I + He II -> Si II + He I. The TR row 2 must gain +rate.
   cx_full = .true.
   call cx_init()
   call cx_set_cell(1.0d4)
   cx_metal_base = 5
   nm0 = 0.0d0; nm1 = 0.0d0; nm2 = 0.0d0
   nm0(iel_Si) = 1.0d0
   fv = 0.0d0
   call cx_add_to_fvec(24, fv, nm0, nm1, nm2, &
                       0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0)
   rate = fv(cx_metal_base + 2*(iel_Si - 1))
   print '(A,ES24.16)', 'D7 Si II production from production module: ', rate
   print '(A,ES24.16)', 'D7 generic row 2 increment: ', fv(2)
   print '(A,ES24.16)', 'D7 expected TR He I row increment: ', rate
   print '(A,ES24.16)', 'D7 He II source with TR mapping -f2-f3: ', -fv(2)-fv(3)
   print '(A,ES24.16)', 'D7 expected He II source: ', -rate
   print '(A,ES24.16)', 'D7 total charge source with TR interpretation: ', &
                       rate-fv(2)+fv(3)
   if (.not. (rate > 0.0d0 .and. fv(2) < 0.0d0)) error stop 1

   ! carrier_source adds this correction after mol_heh_rows has returned
   ! pHp/lHp. This checks the correction, not the full carrier routine.
   fv = 0.0d0
   call he_h_cx_fvec(fv, 1.0d0, 0.0d0, &
                    0.0d0, 1.0d0, 1.0d0, 0.0d0, 1.0d0)
   print '(A,ES24.16)', 'D3 H II correction absent from base pHp/lHp: ', fv(1)
end program plan_revision_source_probe
