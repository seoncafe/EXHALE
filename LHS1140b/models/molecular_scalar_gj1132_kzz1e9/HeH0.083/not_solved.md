# not solved: one scalar movement bound over a column holding a slow H2 front and a far wind

Continued from `molecular_scalar_gj1132_wellmixed/HeH0.083`; this is the
ONE continuation of forty passes the classification allows, and it was not
repeated. The refusing row is the H2 carrier balance at **cell 500,
r = 29.031 R_p** (the domain edge), measure **2.972e-03** against 1.0e-05.

The trajectory does not descend: 2.07e-03 at pass 11, 3.29e-03 at pass 24,
2.972e-03 at pass 40 -- a band, not a fall. Every relaxation ended on the
movement bound.

What differs from the wellmixed pair is WHERE the bound is attained: here it
is the outer wind itself (cell 466, r = 16.50 R_p) and not the front. One
number over a column holds the front in one case and the outer wind in the
other, which is the statement item L22 is about.

The measurement is `docs/lhs1140b_stationary_L7e_20260915.md` section
23: the row decomposition of the refusing cell with the two faces separated,
and the pass trajectory beside it. A regional or coupled relaxation is
opened as plan item L22, a PROPOSAL only.

`MODELS.md` sections 7 and 8 are written by `models/status.py` from this
directory and carry the verdict alone, which is why the reason is here.
