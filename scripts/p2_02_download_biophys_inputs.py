"""
p2_02_download_biophys_inputs.py - Scarica cio' che serve per stimare LAI, copertura e FAPAR.

1. Angoli di Sentinel-2 (sole e vista) a 20 m, via openEO.
   La rete neurale di SNAP li usa: la stessa vegetazione riflette in modo diverso
   a seconda di come e' illuminata e da dove e' osservata.
2. I coefficienti ufficiali delle reti neurali di ESA SNAP (versione 2.1), dal repository
   GitHub del toolbox Sentinel-2 (s2tbx). Sono piccoli file di testo: pesi, normalizzazioni,
   intervalli validi. Li scarichiamo invece di copiarli a mano: fonte tracciabile.
3. Il LAI del servizio BIOPAR (openEO, versione a 3 bande), per confronto.

Output (esclusi da Git):
  data/raw/p1_20220719/s2_angles.nc
  data/raw/p1_20220719/biopar_lai.nc
  data/ancillary/snap_biophys_2_1/<VARIABILE>/<VARIABILE>_<parte>

Come eseguirlo:
    python scripts/p2_02_download_biophys_inputs.py
"""

import sys
import urllib.error
import urllib.request
from pathlib import Path

import openeo
import pandas as pd
import yaml

from config_utils import load_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"
OPENEO_URL = "openeo.dataspace.copernicus.eu"
# Parti di ogni rete neurale SNAP (nomi dei file nel repository s2tbx)
SNAP_PARTS = ["Normalisation", "Denormalisation", "Weights_Layer1_Neurons", "Weights_Layer1_Bias",
              "Weights_Layer2_Neurons", "Weights_Layer2_Bias", "ExtremeCases", "DefinitionDomain_MinMax"]


def download_snap_coefficients(cfg):
    bio = cfg["biophysical"]
    out_root = PROJECT_ROOT / cfg["paths"]["aux"] / f"snap_biophys_{bio['snap_version']}"
    for var in bio["snap_variables"]:
        (out_root / var).mkdir(parents=True, exist_ok=True)
        for part in SNAP_PARTS:
            out = out_root / var / f"{var}_{part}"
            if out.exists():
                continue
            url = f"{bio['snap_url']}/{bio['snap_version']}/{var}/{var}_{part}"
            try:
                urllib.request.urlretrieve(url, out)
            except urllib.error.HTTPError as error:
                if part == "DefinitionDomain_MinMax" and error.code == 404:
                    # Per FCOVER il file non esiste nel repository SNAP: il dominio degli input
                    # e' comune a tutte le reti (stesso insieme di simulazioni) e verra' preso da LAI.
                    print(f"  {var}_{part} non presente in SNAP: usero' quello comune (LAI)")
                    continue
                sys.exit(f"Impossibile scaricare {url}: {error}")
            except Exception as error:
                sys.exit(f"Impossibile scaricare {url}: {error}")
        print(f"Coefficienti SNAP {var}: OK")


def run_job(cube, out_file, title):
    if out_file.exists():
        print(f"{out_file.name} esiste gia': salto.")
        return
    print(f"\nAvvio elaborazione: {title} ...")
    cube.execute_batch(outputfile=str(out_file), out_format="netCDF", title=title)
    print(f"Salvato {out_file.name} ({out_file.stat().st_size / 1e6:.1f} MB)")


def main():
    cfg = load_config()
    bio, d = cfg["biophysical"], cfg["download"]
    aoi = cfg["training_aoi"]
    extent = {k: aoi[k] for k in ("west", "south", "east", "north")}
    day = cfg["reference_date"]
    temporal = [day, (pd.Timestamp(day) + pd.Timedelta(days=1)).strftime("%Y-%m-%d")]
    raw = PROJECT_ROOT / cfg["paths"]["raw"]

    download_snap_coefficients(cfg)

    print(f"\nConnessione a {OPENEO_URL} ...")
    connection = openeo.connect(OPENEO_URL)
    connection.authenticate_oidc()

    # 1. Angoli: variano lentamente nello spazio, quindi interpolazione bilineare
    angles = connection.load_collection(d["s2_collection"], spatial_extent=extent,
                                        temporal_extent=temporal, bands=bio["angle_bands"])
    angles = angles.resample_spatial(resolution=d["s2_resolution_m"], projection=d["output_crs"],
                                     method="bilinear")
    run_job(angles, raw / "s2_angles.nc", f"eo_water p2 {day} - angoli Sentinel-2")

    # 3. BIOPAR: un "processo definito dall'utente" (UDP) pubblicato da VITO/ESA.
    # datacube_from_process lo richiama per nome e indirizzo (namespace), come una funzione.
    biopar = connection.datacube_from_process(
        "biopar", namespace=bio["biopar_namespace"],
        temporal_extent=temporal, spatial_extent=extent, biopar_type=bio["biopar_type"])
    # BIOPAR produce 10 m: lo portiamo a 20 m con la media, sulla stessa griglia degli altri dati
    biopar = biopar.resample_spatial(resolution=d["s2_resolution_m"], projection=d["output_crs"],
                                     method="average")
    run_job(biopar, raw / "biopar_lai.nc", f"eo_water p2 {day} - BIOPAR {bio['biopar_type']}")

    print("\nFatto. Prossimo comando: python scripts/p2_03_biophysical.py")


if __name__ == "__main__":
    main()
