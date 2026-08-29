#!/usr/bin/env python3
"""`<reservoir He/H> <red EW [%A]> <last k dir> <meas_F_H> <meas_F_He>` for one
converged arm, so the shell search can chain arms without parsing files."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import measure as M

case = sys.argv[1]
r = M.last_row(case)
kdir = os.path.join(case, 'k%02d' % int(r[0]))
_, resv = M.case_kzz_heh(case)
print('%.6f %.6f %s %s %s'
      % (resv, M.red_ew(os.path.join(kdir, 'tpm_He10830.txt')),
         os.path.join(kdir, 'output'), r[3], r[4]))
