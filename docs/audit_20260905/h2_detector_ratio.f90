program h2_detector_ratio
  use Cross_sections, only: frac_H2_dissociative_ionization
  use h2_photo_channels
  implicit none
  real*8 :: f, r, sig(4), ion_events, branch(3)
  f=frac_H2_dissociative_ionization(80.d0)
  r=f/(1.d0-f)
  h2_double_ionization_model='chung80'
  call h2_channel_cross_sections(80.d0,1.d0,sig)
  ion_events=r/(1.d0+r)
  branch=[1.d0/(1.d0+r),0.8d0*ion_events,0.2d0*ion_events]
  write(*,'(a,es22.12)') 'tabulated signal ratio r: ',r
  write(*,'(a,3es22.12)') 'current M, S, D fractions: ',sig(1:3)
  write(*,'(a,3es22.12)') 'event interpretation M, S, D: ',branch
  write(*,'(a,2es22.12)') 'current detector ratio, current proton ratio: ', &
     (sig(2)+sig(3))/sig(1),(sig(2)+2.d0*sig(3))/sig(1)
  write(*,'(a,2es22.12)') 'sum and D/(S+D): ',sum(sig),sig(3)/(sig(2)+sig(3))
end program h2_detector_ratio
