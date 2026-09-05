# Named cases the G23 build does not converge, and what blocks them

| case | info | `||R||` | flux spread | its | how it ended | blocking row | worst cell |
|---|---|---|---|---|---|---|---|
| armD_D2_LW_newton | 2 | 1.221E+00 | 2.132E+01 | 0 | no descent | energy (mass 1.422E-01, mom 9.999E-01, energy 1.000E+00) | j=384 r=1.9095 |
| armD_D2_newton | 2 | 1.283E+00 | 2.802E+01 | 8 | no descent | mass 8/8 | j=491 r=4.341 (8) |
| armD_D2_newton_bigstack | 2 | 1.283E+00 | 2.802E+01 | 8 | no descent | mass 8/8 | j=491 r=4.341 (8) |

# Named cases the A build does not converge, and what blocks them

| case | info | `||R||` | flux spread | its | how it ended | blocking row | worst cell |
|---|---|---|---|---|---|---|---|
| armD_D2_LW_newton | 2 | 1.482E-01 | 2.132E+01 | 0 | no descent | energy (mass 1.482E-01, mom 9.795E-02, energy 6.716E-02) | j=371 r=1.7761 |
| armD_D2_newton | 2 | 3.043E-01 | 2.802E+01 | 8 | no descent | mass 8/8 | j=491 r=4.341 (8) |
| armD_D2_newton_bigstack | 2 | 3.043E-01 | 2.802E+01 | 8 | no descent | mass 8/8 | j=491 r=4.341 (8) |
| arm_heh1_x2matched | 2 | 1.123E-04 | 4.534E+00 | 40 | aborted | mass 22/40, energy 18/40 | j=274 r=1.235 (16); j=276 r=1.241 (9); j=275 r=1.238 (7) |
