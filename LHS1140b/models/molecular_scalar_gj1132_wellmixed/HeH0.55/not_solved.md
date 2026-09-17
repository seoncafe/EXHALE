# not solved: one scalar movement bound over a column holding a slow H2 front and a far wind

Ran the forty passes out. The refusing row is the H2 carrier balance at
**cell 306, r = 1.9517 R_p**, measure **2.477e-02** against 1.0e-05.

What the row is made of there (T = 4462 K, x(H2) = 2.66e-06): diffusive
supply -2.29e-02 against chemical destruction -1.91e-02, with advection
-1.74e-03 (7 per cent) and radiation 3.95e-05, 0.09 per cent of the loss.
The destruction is the He+ + H2 channel, 97 per cent of it; the production
is H2+ + H.

What holds it: as for HeH0.083, every relaxation ended on the movement
bound. The state has not arrived rather than stalled -- the refusing cell
migrates outward pass by pass, 233 at pass 1 to 306 at pass 40, at the pace
of the x2 = 1e-02 radius of the front (1.4816 to 1.5075 R_p over the last
five passes, about one cell a pass) -- and the cell that attains the bound
is that front (cell 225, r = 1.2309 R_p, its own row measure 5.1e-03).

The measurement is `docs/lhs1140b_stationary_L7e_20260915.md` section
23: the row decomposition of the refusing cell with the two faces separated,
and the pass trajectory beside it. A regional or coupled relaxation is
opened as plan item L22, a PROPOSAL only.

`MODELS.md` sections 7 and 8 are written by `models/status.py` from this
directory and carry the verdict alone, which is why the reason is here.
