"""
p1_01_find_dates.py - Trova i giorni con Sentinel-2, Sentinel-3 e Landsat 8/9 sereni sul Delta.

Perche': per validare lo sharpening del termico a livello di campo serve un giorno con
  - Sentinel-3 (termico a ~1 km, da affinare),
  - Sentinel-2 (ottico a 20 m, i predittori),
  - Landsat 8 o 9 (termico a 100 m, il riferimento indipendente),
tutti lo stesso giorno e senza nuvole SULLA NOSTRA AREA.

Due fasi:
  1. Ricerca nei cataloghi STAC (veloce, nessun download): giorni con tutti e tre i satelliti
     e copertura nuvolosa "di scena" non eccessiva. Soglie larghe, perche' la copertura di
     scena riguarda 110 km (S2) o 185 km (Landsat), non i nostri 40 km.
  2. Controllo sull'area: per ogni giorno candidato si leggono SOLO i ritagli delle maschere
     di qualita' (Landsat QA_PIXEL, Sentinel-2 SCL) via HTTP (Cloud Optimized GeoTIFF), e si
     calcola la frazione di terraferma limpida dentro l'area di addestramento.

Fonti: Sentinel-2 e Landsat da Microsoft Planetary Computer (STAC + COG firmati),
       Sentinel-3 dal catalogo STAC di Copernicus Data Space (con ripetizioni se occupato).

Output:
  outputs/tables/p1_coincident_dates.csv  (tutti i candidati con la frazione limpida sull'area)

Come eseguirlo:
    python scripts/p1_01_find_dates.py
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd
import planetary_computer
import rioxarray
from pystac_client import Client
from pystac_client.stac_api_io import StacApiIO
from urllib3 import Retry

from config_utils import load_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Ripetizioni automatiche se il server risponde "troppe richieste" (429) o e' occupato (5xx)
RETRY = Retry(total=6, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504],
              allowed_methods=None, respect_retry_after_header=True)

LANDSAT_CLOUD_BITS = (1 << 1) | (1 << 3) | (1 << 4)   # nuvola dilatata, nuvola, ombra
LANDSAT_WATER_BIT = 1 << 7
LANDSAT_FILL_BIT = 1 << 0
S2_CLOUD_SCL = [3, 8, 9, 10]      # ombra, nuvola media, nuvola alta, cirri
S2_IGNORE_SCL = [0, 1, 6]         # nessun dato, saturo, acqua (non contano ne' come limpido ne' nuvola)


def search(client, collection, bbox, start, end, query=None):
    """Tutte le scene di una collezione che toccano l'area nel periodo, come DataFrame."""
    items = list(client.search(collections=[collection], bbox=bbox,
                               datetime=f"{start}/{end}", query=query, limit=100).items())
    rows = []
    for it in items:
        t = pd.Timestamp(it.datetime).tz_convert(None)
        rows.append({"item": it, "datetime_utc": t, "date": t.normalize(),
                     "platform": it.properties.get("platform"),
                     "cloud": it.properties.get("eo:cloud_cover")})
    return pd.DataFrame(rows)


def clip_counts(href, bbox):
    """Legge solo il ritaglio sull'area di una maschera COG e restituisce l'array."""
    with rioxarray.open_rasterio(href) as src:
        return src.rio.clip_box(*bbox, crs="EPSG:4326").squeeze("band", drop=True).values


def landsat_clear_fraction(items, bbox):
    """Frazione di terraferma limpida sull'area, dalle maschere QA_PIXEL di Landsat."""
    clear = valid = 0
    for it in items:
        qa = clip_counts(it.assets["qa_pixel"].href, bbox).astype("int64")
        land = ((qa & LANDSAT_FILL_BIT) == 0) & ((qa & LANDSAT_WATER_BIT) == 0)
        valid += int(land.sum())
        clear += int((land & ((qa & LANDSAT_CLOUD_BITS) == 0)).sum())
    return clear / valid if valid else 0.0


def s2_clear_fraction(items, bbox):
    """Frazione di terraferma limpida sull'area, dalla classificazione SCL di Sentinel-2."""
    clear = valid = 0
    for it in items:
        scl = clip_counts(it.assets["SCL"].href, bbox)
        counted = ~np.isin(scl, S2_IGNORE_SCL)
        valid += int(counted.sum())
        clear += int((counted & ~np.isin(scl, S2_CLOUD_SCL)).sum())
    return clear / valid if valid else 0.0


def main():
    cfg = load_config(verbose=False)
    ds = cfg["date_search"]
    aoi = cfg["training_aoi"]
    bbox = [aoi["west"], aoi["south"], aoi["east"], aoi["north"]]
    cat = cfg["catalogs"]
    m0, m1 = ds["months"]

    pc = Client.open(cat["landsat_stac"], modifier=planetary_computer.sign_inplace,
                     stac_io=StacApiIO(max_retries=RETRY))
    cdse = Client.open(cat["cdse_stac"], stac_io=StacApiIO(max_retries=RETRY))

    # --- Fase 1: cataloghi, un anno alla volta (richieste piu' piccole) ---
    cand = []
    for year in ds["years"]:
        start = f"{year}-{m0:02d}-01"
        end = (pd.Timestamp(f"{year}-{m1:02d}-01") + pd.offsets.MonthEnd(1)).strftime("%Y-%m-%d")
        print(f"\n{year}: cerco Landsat, Sentinel-2, Sentinel-3 ...", flush=True)
        ls = search(pc, cat["landsat_collection"], bbox, start, end)
        if ls.empty:
            continue
        # Solo Landsat 8 e 9 (il 7 ha un difetto di acquisizione dal 2003)
        ls = ls[ls["platform"].isin(cat["landsat_platforms"])]
        ls = ls[ls["cloud"] <= ds["max_cloud_landsat_scene"]]
        s2 = search(pc, cat["s2_collection_pc"], bbox, start, end)
        s2 = s2[s2["cloud"] <= ds["max_cloud_s2_scene"]]
        time.sleep(2)                                     # gentilezza verso il server CDSE
        s3 = search(cdse, cat["s3_collection"], bbox, start, end)
        h0, h1 = ds["s3_day_hours_utc"]
        s3 = s3[(s3["datetime_utc"].dt.hour >= h0) & (s3["datetime_utc"].dt.hour < h1)]
        days = sorted(set(ls["date"]) & set(s2["date"]) & set(s3["date"]))
        print(f"  Landsat {ls['date'].nunique()} giorni, S2 {s2['date'].nunique()}, "
              f"S3 diurno {s3['date'].nunique()} -> in comune {len(days)}")
        for d in days:
            cand.append((d, ls[ls["date"] == d], s2[s2["date"] == d], s3[s3["date"] == d]))

    # --- Fase 2: nuvole sulla nostra area (lettura dei soli ritagli delle maschere) ---
    print(f"\nControllo le nuvole sull'area per {len(cand)} giorni candidati ...")
    rows = []
    for d, ls_d, s2_d, s3_d in cand:
        try:
            f_ls = landsat_clear_fraction(ls_d["item"], bbox)
            f_s2 = s2_clear_fraction(s2_d["item"], bbox) if f_ls >= ds["min_clear_fraction"] else np.nan
        except Exception as error:                        # un file irraggiungibile non ferma tutto
            print(f"  {d:%Y-%m-%d}: errore di lettura ({error.__class__.__name__}), salto")
            continue
        ok = f_ls >= ds["min_clear_fraction"] and f_s2 >= ds["min_clear_fraction"]
        rows.append({
            "data": f"{d:%Y-%m-%d}",
            "landsat": "+".join(sorted(ls_d["platform"].unique())),
            "landsat_ora": ls_d["datetime_utc"].min().strftime("%H:%M"),
            "s3_ore": ", ".join(sorted(s3_d["datetime_utc"].dt.strftime("%H:%M").unique())),
            "limpido_landsat": round(f_ls, 3),
            "limpido_s2": round(f_s2, 3) if np.isfinite(f_s2) else np.nan,
            "utilizzabile": ok,
        })
        print(f"  {d:%Y-%m-%d}: Landsat limpido {f_ls:.0%}, S2 limpido "
              f"{'-' if not np.isfinite(f_s2) else f'{f_s2:.0%}'}  {'OK' if ok else ''}")

    table = pd.DataFrame(rows)
    out = PROJECT_ROOT / cfg["paths"]["tables_common"] / "p1_coincident_dates.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, index=False)
    good = table[table["utilizzabile"]] if not table.empty else table
    print(f"\nGiorni utilizzabili (>= {ds['min_clear_fraction']:.0%} limpido per Landsat e S2): {len(good)}")
    if not good.empty:
        print(good.to_string(index=False))
    print(f"Tabella completa: {out.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
