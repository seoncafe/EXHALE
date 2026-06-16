#!/bin/bash
# Run TPM.py for each converged WASP-52b case, saving model He/Halpha curves.
ROOT=/nfs/mocafe/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE
export MPLBACKEND=Agg
export TPM_RSTAR_RSUN=0.79      # WASP-52 R_star [R_sun]
export TPM_TSTAR=5000.0         # WASP-52 Teff [K] (K2V)
export TPM_ROTP=1.7497          # tidally-locked = orbital period [days]
export TPM_HE_LMIN=10825.0      # widen He window: fast wind broadens the line
export TPM_HE_LMAX=10837.0      #   past the default +/-1.5 A and matches obs span
export TPM_HE_N=481
for d in fxuv1p0_solar fxuv0p5_solar fxuv1p0_he98 fxuv0p5_he98 fxuv1p0_solar_met fxuv0p5_solar_met; do
  CASE="$ROOT/WASP-52b/$d"
  if [ -f "$CASE/output/Hydro_ioniz_adv.txt" ]; then
    TPM_PATH="$CASE" TPM_SAVE_PREFIX="$CASE/" python3 "$ROOT/TPM.py" > "$CASE/tpm.log" 2>&1
    echo "$d  rc=$?  curves=$(ls "$CASE"/tpm_*.txt 2>/dev/null | wc -l)"
  else
    echo "$d: no converged _adv output yet"
  fi
done
echo "TPM ALL DONE"
