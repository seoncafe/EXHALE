# heh100 under the integral-matched SED: residual 1.181e-3

Under the integral-matched SED this composition is the one case the JFNK
finish does not take below the 1e-3 target: from every seed tried (its own
fiducialA-converged state, heh1000, heh10, a 20k-step re-relaxation) the
solver stagnates (info=2), and from the best seed it stalls at
||R|| = 1.181e-3 -- 18 % above target -- unchanged across four Newton
restarts. The adopted solution is that stagnation point. The other five
compositions sit at 4.6e-4 to 9.8e-4, so heh100's extra residual moves the
observables in the third digit at most, and its Mdot (7.53) is continuous
with its neighbors (heh10: 7.54, heh1000: 7.54).
