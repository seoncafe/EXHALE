#!/bin/bash
# Rebuild results.txt from the run directories.  Everything below the header
# is read out of the arms themselves by measure.py and crossing.py, so the
# table cannot drift from what was actually solved.
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
SC=$EX/LHS1140b/exhale/kzz_profile_scan
FC=$EX/LHS1140b/exhale/flux_closure
REF="$FC/heh3 $FC/heh5 $FC/heh8 $FC/heh9p7 $FC/heh10p3 $FC/heh10p7 $FC/heh11p1 $FC/heh12_ctl"

{
cat $SC/results_header.txt
echo
echo '=========================================================================='
echo 'ALL ARMS'
echo '=========================================================================='
echo '(the K_zz = 1e9 rows named ../flux_closure/... are the arms that were'
echo ' already on disk; every kzz*/ row was solved by this scan)'
echo
python3 $SC/measure.py $REF $(ls -d $SC/kzz*/heh* | grep -v attempt)
echo
echo '=========================================================================='
echo 'CROSSING WITH THE MEASURED EQUIVALENT WIDTH -- COOL BRANCH'
echo '=========================================================================='
python3 $SC/crossing.py --cool $REF $SC/kzz1e9 $SC/kzz1e8 $SC/kzz3p16e7 \
        $SC/kzz1p78e7 $SC/kzz1p33e7 $SC/kzz1e7 2>&1 | grep -v '(skipped'
echo
echo '=========================================================================='
echo 'CROSSING WITH THE MEASURED EQUIVALENT WIDTH -- HOT BRANCH'
echo '=========================================================================='
python3 $SC/crossing.py --hot $SC/kzz1e9 $SC/kzz1e8 $SC/kzz3p16e7 \
        $SC/kzz1p78e7 $SC/kzz1p33e7 $SC/kzz1e7 2>&1 | grep -v '(skipped'
echo
echo '=========================================================================='
echo 'CHEMISTRY-ONLY PROBES (clima_probe/)'
echo '=========================================================================='
echo 'Does the chemistry+climate handoff exist at all?  Trial fluxes fixed at'
echo 'F_H = 4.0e6, F_He = 2.2e7 g/s, so these bound the handoff and not the'
echo 'closure.  "atomic-H share" is the share of the hydrogen nuclei that the'
echo 'CONVERGED model top carries as atomic H; above 0.01 the adapter refuses,'
echo 'because the elemental hydrogen flux was imposed entirely on H2.'
echo
for d in $SC/clima_probe/heh*_k*/; do
  n=$(basename $d)
  s=$(grep -h "atomic-H share" $d/probe.log 2>/dev/null | sed 's/.*imposition, //;s/ converged//')
  w=$(grep -hoE "REFUSED: [^.]{0,90}|ClimaException: .*" $d/probe.log 2>/dev/null | head -1)
  printf '%-24s atomic-H share %-10s %s\n' "$n" "${s:-n/a}" "$w"
done
} > $SC/results.txt
echo "wrote $SC/results.txt"
