#!/bin/bash
cd /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L25
printf '%s\n' roe_fid_HeH2.13 roe_s030_HeH2.13 roe_s010_HeH2.13 roe_pfid_HeH9.7 g1000_roe g500_hllc g1000_hllc   | xargs -P 3 -I{} bash -c 'bash /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L25/launch_{}.sh > /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L25/{}.launch.log 2>&1; echo "$(date +%T) {} done"'
