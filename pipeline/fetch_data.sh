#!/usr/bin/env bash
# Download every raw input used by build_dataset.py into data/raw/ (≈2.1 GB, mostly LMR SST fields).
set -euo pipefail
cd "$(dirname "$0")/../data" && mkdir -p raw/large && cd raw

CPC=https://www.cpc.ncep.noaa.gov/data/indices
curl -fsSO $CPC/oni.ascii.txt
curl -fsSO $CPC/RONI.ascii.txt
curl -fsSO $CPC/Rnino34.ascii.txt

T=https://www.ncei.noaa.gov/pub/data/paleo/icecore/trop
for f in quelccaya/quelccaya2013-noaa.txt \
         quelccaya/mayewski2023/mayewski2023-qu-18_isotopes.txt quelccaya/mayewski2023/mayewski2023-qu-18_ic.txt \
         huascaran/weber2023/weber2023-hs2019.txt huascaran/weber2026/weber2026-dust.txt \
         huascaran/thompson1995-annualt-noaa.txt \
         illimani/illimani2010nh4ann-noaa.txt \
         illimani/lindau2021/lindau2021-il2017_isotopes.txt illimani/lindau2021/lindau2021-il2017_ic.txt; do
  curl -fsSO "$T/$f"
done

curl -fsS -o ne_50m_admin_0_countries.geojson \
  https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_0_countries.geojson
curl -fsS -o etopo_andes.csv \
  "https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.csv?altitude%5B(-24):6:(4)%5D%5B(-84):6:(-60)%5D"

cd large
curl -fsS -o ersst.v6.sst.mnmean.nc https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v6/sst.mnmean.nc
L=https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/tardif2019lmr/v2_1
for f in posterior_climate_indices_MCruns_ensemble_full_LMRv2.1.nc \
         sst_MCruns_ensemble_mean_LMRv2.1.nc sst_MCruns_ensemble_spread_LMRv2.1.nc; do
  curl -fsSO "$L/$f"
done
echo "done — now run: python3 pipeline/build_dataset.py"
