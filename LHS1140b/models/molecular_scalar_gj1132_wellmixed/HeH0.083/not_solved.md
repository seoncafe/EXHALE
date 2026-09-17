# not solved: one scalar movement bound over a column holding a slow H2 front and a far wind

Refused at outer pass 24 by the stagnation rule. The refusing row is the H2
carrier balance at **cell 500, r = 29.031 R_p** (the domain edge), measure
**2.462e-02** against a tolerance of 1.0e-05.

What the row is made of there: the chemistry is 3 per cent of the residual
(net source 1.96e-05 against a residual of 6.75e-04); the two advective face
fluxes are each 1.3e-02 and cancel to -3.23e-04; the diffusive inflow
9.79e-04 has no outward diffusive face, and advection removes a third of it.
The mass flux is outward at both faces.

What holds it: every carrier relaxation of every pass ended on the
composition movement bound and none on its own residual. The bound was cut
5.0e-03 -> 2.5e-03 -> 1.25e-03 -> floor 1.0e-03 and the movement of a pass
fell in proportion. The cell that ATTAINS the bound is the H2 front (cell
231, r = 1.2565 R_p, its own row measure 6.6e-03), four times inside the
cell that refuses. The measure has been flat between 2.41e-02 and 2.46e-02
since pass 15 and the front has not moved more than one cell since pass 19.

The measurement is `docs/lhs1140b_stationary_L7e_20260915.md` section
23: the row decomposition of the refusing cell with the two faces separated,
and the pass trajectory beside it. A regional or coupled relaxation is
opened as plan item L22, a PROPOSAL only.

`MODELS.md` sections 7 and 8 are written by `models/status.py` from this
directory and carry the verdict alone, which is why the reason is here.
