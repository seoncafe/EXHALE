# The element face-flux table of the `K_zz` cases

`src/modules/functions/binary_element_diffusion.f90` writes
`output/element_flux_profile.txt` -- the face fluxes, the stage-resolved
binary coefficient `D_eff` and the eddy term -- only when
`EXHALE_DIFFUSION_CHECK=1` is set or a lower-atmosphere profile is in use.
The prescribed-composition cases of `models/` are solved without it, so they
do not carry that table, and the homopause panel of
`docs/figures/lhs1140b_kzz_profiles.pdf` -- the radius at which `D_eff`
equals the run's own `K_zz` -- has nothing to read.

Each directory here is the post-processing pass of one solved case repeated
with the flag on: the case's own `input.inp` with `Load IC? True`,
`Do only PP: True` and `CFL: 1.0e-12`, its solved state as the initial
condition, nothing else changed. The pass takes no time step that moves the
state (MEASURED against the case's stored `Hydro_ioniz_adv.txt`: largest
relative difference 1.4e-4), so the table describes the solved case and the
case directory itself is left as its own run wrote it.

```bash
d=diffusion_check/atomic_scalar_gj1132_kzz1e9_HeH0.55
c=atomic_scalar_gj1132_kzz1e9/HeH0.55
mkdir -p $d/output
cp -f $c/input.inp $d/
cp -f $c/output/Hydro_ioniz.txt $d/output/Hydro_ioniz_IC.txt
cp -f $c/output/Ion_species.txt $d/output/Ion_species_IC.txt
sed -i -e 's/^Load IC?.*/Load IC? True/' -e 's/^Do only PP:.*/Do only PP: True/' \
       -e '/^Restart intent:/d' -e '/^Solver:/d' -e '/^du_th /d' $d/input.inp
echo 'CFL: 1.0e-12' >> $d/input.inp
(cd $d && OMP_NUM_THREADS=2 EXHALE_DIFFUSION_CHECK=1 ../../../../EXHALE.x > run.log 2>&1)
```
