# L22 step 2, 2026-09-16: identities of the measurement

The frozen states reloaded (md5 of the files copied into the run directories as
`output/Hydro_ioniz_IC.txt` and `output/Ion_species_IC.txt`):

| case | file | md5 |
|---|---|---|
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | `output/Hydro_ioniz.txt` | `c35719db84ecc19926f38a4db5718fa2` |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | `output/Ion_species.txt` | `f91db8eecc98743b9c82cec59549b182` |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | `output/Hydro_ioniz.txt` | `5a9de6ae680246c6badbec424d60b8dd` |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | `output/Ion_species.txt` | `b0635219ad30ce5551a066e3b890c90e` |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | `output/Hydro_ioniz.txt` | `1229b961323d7ac710146650da72454a` |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | `output/Ion_species.txt` | `5c5bfcde1ea939adf39941bdd29002c6` |

| item | value |
|---|---|
| measuring binary | `EXHALE_L22.x`, md5 `f23da8d84fafc2f0587b610add56e714` (private build `build_L22`, deleted at the end of the item) |
| tree binary, untouched | `EXHALE.x`, md5 `c2e9c9990b9f14f1be8cd77abca68945` |
| edited source | `src/modules/lower_atmosphere/diffusive_photochemistry.f90`, md5 `2ac1efac64119da315a4fdccd9c13ff9` after the edit |
| threads | `OMP_NUM_THREADS=1` in every run |
| run controls | `EXHALE_PTC_DTAU0=1.0e8`, `EXHALE_OUTER_PASSES=3`, `EXHALE_CARRIER_DEBUG=1`, `EXHALE_CARRIER_ROW_TERMS=1` |

Mask files in this directory: `wm083_narrow.msk`, `wm055_narrow.msk`,
`kz083_narrow.msk`, `farwind.msk`; each carries its own definition in its header.
