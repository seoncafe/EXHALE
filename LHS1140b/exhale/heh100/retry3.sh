#!/bin/bash
cd "$(dirname "$0")/.."          # exhale/
bash finish_case.sh heh100 heh10 && { echo CHAIN_OK; exit 0; }
for i in 1 2 3; do
  echo "== own-state restart $i =="
  bash finish_case.sh heh100 && { echo CHAIN_OK; exit 0; }
done
echo CHAIN_FAILED
