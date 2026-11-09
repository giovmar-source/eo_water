"""
p2_06_download_cropmap.py - Ritaglio della mappa colturale europea EUCROPMAP 2022 sul Delta del Po.

Perche': TSEB dipende dall'altezza della chioma (rugosita' -> resistenza aerodinamica -> H).
Finora l'abbiamo stimata dal LAI, ma un campo di mais e una risaia con lo stesso LAI hanno
altezze molto diverse (circa 2,5 m contro 1 m). La mappa colturale ci dice COSA c'e' in ogni campo.

Dato: EUCROPMAP 2022 (JRC; d'Andrimont et al., Scientific Data 2024), 10 m, 19 colture
+ classi di copertura del suolo, ottenuto con foresta casuale su Sentinel-1 e Sentinel-2
e dati a terra LUCAS 2022. Accuratezza dichiarata ~71% sulle 19 colture (UE).
Licenza CC BY 4.0.

Il file europeo pesa 6,8 GB, ma e' un GeoTIFF a blocchi: GDAL legge via HTTP solo i blocchi
che servono (come per Landsat in p1_02). Scarichiamo quindi pochi MB.

ATTENZIONE alle classi: e' un dato CATEGORICO. Per cambiare griglia si usa la MODA (la classe
piu' frequente), mai la media: la media tra "mais" (216) e "riso" (217) non e' una coltura!

Output:
  data/ancillary/eucropmap_2022_training_aoi.tif   (10 m, CRS originale; escluso da Git)
  data/ancillary/eucropmap_2022_20m.tif            (20 m sulla griglia della LST, moda)
  outputs/tables/p2_06_crop_classes.csv      (superficie per classe nell'area di analisi)

Come eseguirlo:
    python scripts/p2_06_download_cropmap.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import rioxarray
import xarray as xr
import yaml

from config_utils import load_config
from rasterio.enums import Resampling

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"


def main():
    cfg_all = load_config()
    cm = cfg_all["cropmap"]
    t, a = cfg_all["training_aoi"], cfg_all["analysis_aoi"]
    aux = PROJECT_ROOT / cfg_all["paths"]["aux"]
    proc = PROJECT_ROOT / cfg_all["paths"]["processed"]
    aux.mkdir(parents=True, exist_ok=True)

    # 1. Lettura parziale via HTTP: rioxarray passa l'URL a GDAL (/vsicurl/)
    print(f"Apro {cm['url']} (lettura parziale)...")
    full = rioxarray.open_rasterio(cm["url"], masked=False).squeeze("band", drop=True)
    print(f"  CRS {full.rio.crs}, pixel {full.rio.resolution()}, "
          f"dimensioni europee {full.shape}")
    clip = full.rio.clip_box(t["west"], t["south"], t["east"], t["north"], crs="EPSG:4326")
    clip = clip.load()                                     # qui avviene il download vero
    clip.rio.to_raster(aux / "eucropmap_2022_training_aoi.tif", compress="deflate")
    print(f"  ritaglio area di addestramento: {clip.shape}, salvato")

    # 2. Sulla griglia a 20 m della LST affinata, con la MODA
    lst = rioxarray.open_rasterio(proc / cfg_all["tseb"]["lst_file"]).squeeze("band", drop=True)
    crop20 = clip.rio.reproject_match(lst, resampling=Resampling.mode)
    crop20.rio.to_raster(aux / "eucropmap_2022_20m.tif", compress="deflate")
    print(f"  griglia 20 m: {crop20.shape}, salvata")

    # 3. Che colture ci sono nell'area di analisi?
    sub = crop20.rio.clip_box(a["west"], a["south"], a["east"], a["north"], crs="EPSG:4326")
    vals, counts = np.unique(sub.values, return_counts=True)
    names = {int(k): v for k, v in cm["legend"].items()}
    df = pd.DataFrame({"codice": vals, "classe": [names.get(int(v), "?") for v in vals],
                       "ettari": counts * 400 / 1e4})        # pixel 20x20 m = 400 m2
    df["percentuale"] = 100 * df["ettari"] / df["ettari"].sum()
    df = df.sort_values("ettari", ascending=False).round(1)
    out = PROJECT_ROOT / cfg_all["paths"]["tables"] / "p2_06_crop_classes.csv"
    df.to_csv(out, index=False)
    print("\nClassi nell'area di analisi:")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
