   module hydrogen_n2_rates
   ! Atomic data and rate coefficients of the hydrogen n=2 manifold and of the
   ! Ly-alpha line, in one place.
   !
   ! Both the H(n=2) statistical-equilibrium solver (excited_hydrogen.f90) and
   ! the Ly-alpha escape-probability transfer (lya_rt.f90) need the same
   ! constants (A_2p1s, A_2s1s, the statistical weights, the line opacity
   ! constant) and the same collisional rate coefficients. They used to carry
   ! private copies of both. The copies drifted: the reverse-rate expressions
   ! written out inside n2_populations evaluated their statistical weights
   ! against shadowing dummy arguments and so disagreed with the identical
   ! expressions a few lines further down in the same file. Keeping one
   ! definition here removes that failure mode.
   !
   ! Collisional data: Christie, Arras & Li (2013, ApJ 772, 144) Table 2.
   ! Recombination: Draine (2011) case B and its 2s/2p split.
   ! All quantities cgs.

   use global_parameters, only: hp_erg, c_light

   implicit none

   ! ----- Ly-alpha line and n=2 level data ----- !
   real*8, parameter :: lA_lya   = 1215.6701d-8        ! Ly-alpha wavelength [cm]
   real*8, parameter :: nu_lya   = c_light/lA_lya      ! Ly-alpha frequency [s^-1]
   real*8, parameter :: A_2p1s   = 6.3d8               ! A(2p->1s) [s^-1]
   real*8, parameter :: A_2s1s   = 8.26d0              ! A(2s->1s) two-photon [s^-1]
   real*8, parameter :: g1s = 2.0d0, g2s = 2.0d0, g2p = 6.0d0
   real*8, parameter :: E21_erg  = 1.634d-11           ! 1s-2s/2p gap, 10.2 eV [erg]
   ! Einstein-B in the J_nu (mean-intensity) convention: B*J gives s^-1.
   real*8, parameter :: B21_lya = A_2p1s*c_light**2.0/(2.0d0*hp_erg*nu_lya**3.0)
   real*8, parameter :: B12_lya = (g2p/g1s)*B21_lya
   ! Line-opacity constant: sqrt(pi) e^2/(m_e c) f_lya [cm^2 Hz], so that the
   ! line-center cross section is C_lya/Dnu_D.
   real*8, parameter :: f_lya    = 0.4162d0            ! Ly-alpha oscillator strength
   real*8, parameter :: C_lya    = 1.49736d-2*f_lya    ! [cm^2 Hz]

   contains

   ! --------------------------------------------------------------- !
   ! Collisional rate coefficients [cm^3 s^-1]. T in K.

   real*8 function c1s2s_rate(T)
   ! 1s->2s collisional excitation.
   real*8, intent(in) :: T
   real*8 :: Tl, t4
   Tl = max(T,1.0d0)
   t4 = Tl/1.0d4
   c1s2s_rate = 1.21d-8*(1.0d0/t4)**0.455d0*exp(-118400.0d0/Tl)
   end function c1s2s_rate

   real*8 function c1s2p_rate(T)
   ! 1s->2p collisional excitation.
   real*8, intent(in) :: T
   real*8 :: Tl, t4
   Tl = max(T,1.0d0)
   t4 = Tl/1.0d4
   c1s2p_rate = 1.71d-8*(1.0d0/t4)**0.077d0*exp(-118400.0d0/Tl)
   end function c1s2p_rate

   real*8 function c2s1s_rate(T)
   ! 2s->1s collisional de-excitation, by detailed balance from c1s2s_rate in
   ! the analytically-cancelled form (the Boltzmann factor cancels, avoiding
   ! 0*inf at very low T).
   real*8, intent(in) :: T
   real*8 :: t4
   t4 = max(T,1.0d0)/1.0d4
   c2s1s_rate = 1.21d-8*(1.0d0/t4)**0.455d0*(g1s/g2s)
   end function c2s1s_rate

   real*8 function c2p1s_rate(T)
   ! 2p->1s collisional de-excitation (detailed balance from c1s2p_rate).
   real*8, intent(in) :: T
   real*8 :: t4
   t4 = max(T,1.0d0)/1.0d4
   c2p1s_rate = 1.71d-8*(1.0d0/t4)**0.077d0*(g1s/g2p)
   end function c2p1s_rate

   real*8 function c2s2p_rate(T)
   ! 2s->2p l-mixing by electron impact.
   real*8, intent(in) :: T
   real*8 :: Tl
   Tl = max(T,1.0d0)
   c2s2p_rate = 6.21d-5*(log(Tl/1.02d0) - 0.57721d0)/sqrt(Tl)
   end function c2s2p_rate

   real*8 function c2p2s_rate(T)
   ! 2p->2s l-mixing (detailed balance; the levels are degenerate, so only the
   ! statistical weights enter).
   real*8, intent(in) :: T
   c2p2s_rate = c2s2p_rate(T)*(g2s/g2p)
   end function c2p2s_rate

   ! --------------------------------------------------------------- !
   ! Recombination into the n=2 levels [cm^3 s^-1]. Draine (2011), Table 2
   ! entries R2 (case B), R8 (2s) and R9 (2p).

   real*8 function alpha_B_hydrogen(T)
   real*8, intent(in) :: T
   real*8 :: t4
   t4 = max(T,1.0d0)/1.0d4
   alpha_B_hydrogen = 2.54d-13*t4**(-0.8163d0 - 0.0208d0*log(t4))
   end function alpha_B_hydrogen

   real*8 function alpha_2s_hydrogen(T)
   real*8, intent(in) :: T
   real*8 :: t4
   t4 = max(T,1.0d0)/1.0d4
   alpha_2s_hydrogen = (0.282d0 + 0.047d0*t4 - 0.006d0*t4**2.0)             &
                     *alpha_B_hydrogen(T)
   end function alpha_2s_hydrogen

   real*8 function alpha_2p_hydrogen(T)
   real*8, intent(in) :: T
   alpha_2p_hydrogen = alpha_B_hydrogen(T) - alpha_2s_hydrogen(T)
   end function alpha_2p_hydrogen

   ! --------------------------------------------------------------- !

   real*8 function n2p_destruction_rate(T, ne_l, gam_ion_2s, gam_ion_2p)
   ! Rate [s^-1] at which a 2p atom is removed WITHOUT emitting a Ly-alpha
   ! photon into the line. Three channels:
   !
   !   ne q(2p->1s)                 collisional de-excitation (10.2 eV back to
   !                                the electron gas)
   !   Gamma_2p                     photoionization of n=2 by the stellar
   !                                Balmer continuum
   !   ne C(2p->2s) * P_2gamma      l-mixing into 2s followed by two-photon
   !                                decay, with the 2s branching ratio
   !                                P_2gamma = A_2s1s/(A_2s1s + ne q(2s->1s)
   !                                                   + ne C(2s->2p) + Gamma_2s)
   !
   ! Together with the escape rate beta*A_2p1s these close the 2p budget:
   ! n2p = P/(beta A_2p1s + n2p_destruction_rate). Dropping the destruction
   ! term over-estimates the trapped population wherever beta is small enough
   ! that beta*A_2p1s falls to the destruction rate, i.e. at the base.
   real*8, intent(in) :: T, ne_l, gam_ion_2s, gam_ion_2p
   real*8 :: L2s, P_2gamma

   L2s = A_2s1s + ne_l*(c2s1s_rate(T) + c2s2p_rate(T)) + gam_ion_2s
   P_2gamma = A_2s1s/max(L2s, 1.0d-30)
   n2p_destruction_rate = ne_l*c2p1s_rate(T) + gam_ion_2p                   &
                        + ne_l*c2p2s_rate(T)*P_2gamma

   end function n2p_destruction_rate

   ! End of module
   end module hydrogen_n2_rates
