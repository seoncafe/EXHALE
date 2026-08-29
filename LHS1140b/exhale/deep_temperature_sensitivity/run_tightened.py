#!/usr/bin/env python3
"""The same handoff run, with the photochemical stopping test tightened.

`photochem_to_lower_profile.py` is called unchanged; the only thing this
wrapper does is replace `EvoAtmosphereGasGiant` with a subclass that sets

    var.conv_longdy      the largest relative change of any mixing ratio
                         between t and t/2 that still counts as steady
                         (`gasgiants.py` sets 1e-2)
    var.conv_longdydt    the same change per unit time (1e-6 by default)
    var.equilibrium_time the integration time past which the run is declared
                         steady whatever the state is doing (1e17 s)

so that the exit is decided by the convergence test alone and at a stated
tolerance.  Everything else -- climate solve, mechanism, grid, boundary
conditions, the written file -- is the production path.

usage: run_tightened.py <longdy> <longdydt> <equilibrium_time> -- <adapter args>
"""
import os
import sys

EX = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00'
sys.path.insert(0, os.path.join(EX, 'src', 'utils'))

LONGDY = float(sys.argv[1])
LONGDYDT = float(sys.argv[2])
EQTIME = float(sys.argv[3])
assert sys.argv[4] == '--'
sys.argv = [sys.argv[0]] + sys.argv[5:]

import photochem_to_lower_profile as ad                      # noqa: E402
from photochem.extensions import gasgiants                   # noqa: E402

_Base = gasgiants.EvoAtmosphereGasGiant


class StatedTolerance(_Base):
    def __init__(self, *a, **k):
        _Base.__init__(self, *a, **k)
        self.var.conv_longdy = LONGDY
        self.var.conv_longdydt = LONGDYDT
        self.var.equilibrium_time = EQTIME


gasgiants.EvoAtmosphereGasGiant = StatedTolerance
print('# stopping test: conv_longdy = %.1e, conv_longdydt = %.1e, '
      'equilibrium_time = %.1e s' % (LONGDY, LONGDYDT, EQTIME), flush=True)
ad.main()
