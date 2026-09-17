set -e
L=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L22
R=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
cd $L/i3_alt
for f in Hydro_ioniz Ion_species; do \cp -f output/$f.txt output/${f}_IC.txt; done
\cp -f input.inp input.inp.solved
sed -i -e 's/^Load IC?.*/Load IC? True/' -e 's/^Do only PP:.*/Do only PP: True/' \
       -e '/^Restart intent:/d' -e '/^Solver:/d' -e '/^du_th /d' input.inp
echo 'CFL: 1.0e-12' >> input.inp
OMP_NUM_THREADS=8 $L/EXHALE_L22e.x > pp.log 2>&1
\mv -f input.inp.solved input.inp
. $R/LHS1140b/winered_hires_y.sh
MPLBACKEND=Agg PYTHONPATH=$R python3 $R/EXHALE_transit.py > transit.log 2>&1
echo PP_AND_TRANSIT_DONE
