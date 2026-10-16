! Scratch driver: EXHALE's photoabsorption cross sections [1e-18 cm^2] on a
! list of photon energies [eV] (stdin: n, then n energies).  Columns:
! E, HI, HeI, HeII, HeI(2^3S), H2 total, then the 17 metal stages in
! Cross_sections column order.
program dump_xs
use global_parameters
use Cross_sections
implicit none
integer :: nE_in, iE, kk
real*8 :: Eph
real*8 :: row(22)
read(*,*) nE_in
do iE = 1, nE_in
   read(*,*) Eph
   row(1) = sigma(Eph, ih, e_th_HI)
   row(2) = sigma_HeI(Eph)
   row(3) = sigma(Eph, ihe, e_th_HeII)
   row(4) = sigma_HeI23S(Eph)
   row(5) = sigma_H2(Eph)
   do kk = 1, 17
      row(5+kk) = metal_photoion_sigma(kk, Eph)
   enddo
   write(*,'(es16.8,22es16.8)') Eph, row
enddo
end program dump_xs
