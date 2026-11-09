"""
p1_02_download.py - Scarica i dati del giorno di riferimento per lo sharpening del termico.

Per il giorno indicato in config/phase1_sharpening.yaml (reference_date) e
sull'AREA DI ADDESTRAMENTO (~40 x 40 km) scarica:
  1. Sentinel-2 L2A: tutte le bande a 20 m + SCL          (openEO, Copernicus Data Space)
  2. Sentinel-3 SLSTR LST + flag di qualita' a ~1 km       (openEO, Copernicus Data Space)
  3. Landsat 9: temperatura superficiale + qualita' pixel  (STAC, Microsoft Planetary Computer)
     -> il riferimento indipendente per la validazione.

Output in data/raw/p1_20220719 (escluso da Git):
  s2_l2a.nc, s3_lst.nc, landsat_lwir11.tif, landsat_qa_pixel.tif, landsat_metadata.json

Come eseguirlo:
    python scripts/p1_02_download.py

Se un file esiste gia' non viene riscaricato (per rifarlo, cancellalo).
"""

import json
import sys
import time
from pathlib import Path

import openeo
import pandas as pd
import planetary_computer
import rioxarray
import yaml

from config_utils import load_config
from pystac_client import Client
from pystac_client.stac_api_io import StacApiIO
from urllib3.util.retry import Retry

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"
OPENEO_URL = "openeo.dataspace.copernicus.eu"

RETRY = Retry(total=6, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504],
              allowed_methods=None, respect_retry_after_header=True)


# ----------------------------------------------------------------------------
# Parte 1 e 2: Sentinel-2 e Sentinel-3 con openEO
# ----------------------------------------------------------------------------

def openeo_connect():
    """Login a openEO (se il token e' gia' salvato non chiede nulla)."""
    print(f"Connessione a {OPENEO_URL} ...")
    connection = openeo.connect(OPENEO_URL)
    connection.authenticate_oidc()
    print("Login riuscito.\n")
    return connection


def check_bands(connection, collection_id, wanted):
    """Ferma lo script subito se una banda richiesta non esiste nella collezione."""
    meta = connection.describe_collection(collection_id)
    available = meta.get("cube:dimensions", {}).get("bands", {}).get("values", [])
    missing = [b for b in wanted if b not in available]
    if missing:
        sys.exit(f"{collection_id}: bande non trovate {missing}. Disponibili: {available}")
    print(f"{collection_id}: bande OK {wanted}")


def run_job(cube, out_file, title):
    """Esegue un batch job su CDSE e salva il NetCDF (salta se il file esiste)."""
    if out_file.exists():
        print(f"{out_file.name} esiste gia': salto.")
        return
    print(f"\nAvvio elaborazione: {title} ...")
    cube.execute_batch(outputfile=str(out_file), out_format="netCDF", title=title)
    print(f"Salvato {out_file.name} ({out_file.stat().st_size / 1e6:.1f} MB)")


def download_sentinel(cfg, out_dir):
    d = cfg["download"]
    aoi = cfg["training_aoi"]
    extent = {k: aoi[k] for k in ("west", "south", "east", "north")}
    day = cfg["reference_date"]
    # Intervallo di un giorno: openEO considera la data finale ESCLUSA
    next_day = (pd.Timestamp(day) + pd.Timedelta(days=1)).strftime("%Y-%m-%d")
    temporal = [day, next_day]

    connection = openeo_connect()
    check_bands(connection, d["s2_collection"], d["s2_bands"])
    check_bands(connection, d["s3_collection"], d["s3_bands"])

    # Sentinel-3: scarichiamo tutti i passaggi del giorno; la scelta del passaggio si fa dopo,
    # in locale (config_utils.pick_s3_pass, oppure l'ora fissata nel config per quella data)
    s3 = connection.load_collection(d["s3_collection"], spatial_extent=extent,
                                    temporal_extent=temporal, bands=d["s3_bands"])
    s3 = s3.resample_spatial(resolution=d["s3_resolution_m"], projection=d["output_crs"], method="near")
    run_job(s3, out_dir / "s3_lst.nc", f"eo_water p1 {day} - Sentinel-3 LST")

    # Sentinel-2: tutte le bande a 20 m. "near" non mescola valori, importante per SCL.
    s2 = connection.load_collection(d["s2_collection"], spatial_extent=extent,
                                    temporal_extent=temporal, bands=d["s2_bands"])
    s2 = s2.resample_spatial(resolution=d["s2_resolution_m"], projection=d["output_crs"], method="near")
    run_job(s2, out_dir / "s2_l2a.nc", f"eo_water p1 {day} - Sentinel-2 L2A")


# ----------------------------------------------------------------------------
# Parte 3: Landsat 9 con STAC (Microsoft Planetary Computer)
# ----------------------------------------------------------------------------

def overlap_fraction(item_bbox, aoi_bbox):
    """Frazione dell'area di interesse coperta dal riquadro della scena (in gradi, approssimato)."""
    w = max(item_bbox[0], aoi_bbox[0]); s = max(item_bbox[1], aoi_bbox[1])
    e = min(item_bbox[2], aoi_bbox[2]); n = min(item_bbox[3], aoi_bbox[3])
    if e <= w or n <= s:
        return 0.0
    return ((e - w) * (n - s)) / ((aoi_bbox[2] - aoi_bbox[0]) * (aoi_bbox[3] - aoi_bbox[1]))


def download_landsat(cfg, out_dir):
    d = cfg["download"]
    cat = cfg["catalogs"]
    aoi = cfg["training_aoi"]
    bbox = [aoi["west"], aoi["south"], aoi["east"], aoi["north"]]
    day = cfg["reference_date"]

    outputs = {key: out_dir / f"landsat_{key}.tif" for key in d["landsat_assets"]}
    meta_file = out_dir / "landsat_metadata.json"
    if all(p.exists() for p in outputs.values()) and meta_file.exists():
        print("File Landsat gia' presenti: salto.")
        return

    print("\nCerco la scena Landsat ...")
    # modifier=sign_inplace: aggiunge ai link dei file una "firma" temporanea, necessaria
    # per scaricarli da Planetary Computer (la ricerca invece e' libera).
    catalog = Client.open(cat["landsat_stac"], modifier=planetary_computer.sign_inplace,
                          stac_io=StacApiIO(max_retries=RETRY))
    items = list(catalog.search(
        collections=[cat["landsat_collection"]],
        bbox=bbox,
        datetime=f"{day}T00:00:00Z/{day}T23:59:59Z",
        query={"platform": {"in": d["landsat_platforms"]}},
    ).items())
    if not items:
        # Landsat passa ogni 8 giorni (8 e 9 insieme): spesso non c'e'. Non e' un errore:
        # la validazione si fara' con ECOSTRESS, se disponibile (vedi p1_05).
        print(f"Nessuna scena Landsat ({', '.join(d['landsat_platforms'])}) il {day}: "
              "salto Landsat.")
        return

    # Se piu' scene toccano l'area, prendiamo quella che la copre di piu'
    item = max(items, key=lambda it: overlap_fraction(it.bbox, bbox))
    cover = overlap_fraction(item.bbox, bbox)
    print(f"Scena: {item.id}  ({item.datetime:%Y-%m-%d %H:%M} UTC, copertura area ~{cover:.0%})")
    if cover < 0.99:
        print("ATTENZIONE: la scena non copre tutta l'area di addestramento; il bordo restera' vuoto.")

    for key, out in outputs.items():
        if out.exists():
            print(f"{out.name} esiste gia': salto.")
            continue
        href = item.assets[key].href
        print(f"Leggo {key} (solo il ritaglio sull'area, non la scena intera) ...")
        # Il file e' un Cloud Optimized GeoTIFF: GDAL legge via internet solo i blocchi
        # che servono. clip_box accetta il riquadro in gradi (crs=EPSG:4326) e lo
        # converte nella proiezione della scena.
        with rioxarray.open_rasterio(href) as src:
            clip = src.rio.clip_box(*bbox, crs="EPSG:4326").squeeze("band", drop=True).load()
        clip.rio.to_raster(out)
        print(f"Salvato {out.name}: {clip.shape[0]} x {clip.shape[1]} pixel, CRS {clip.rio.crs}")

    # Metadati: servono per convertire i valori grezzi in Kelvin (scala e offset)
    band_info = item.assets["lwir11"].extra_fields.get("raster:bands", [{}])[0]
    meta = {
        "item_id": item.id,
        "datetime_utc": item.datetime.isoformat(),
        "platform": item.properties.get("platform"),
        "cloud_cover": item.properties.get("eo:cloud_cover"),
        # Valori ufficiali USGS Collection 2 Level-2, usati se il catalogo non li riporta
        "st_scale": band_info.get("scale", 0.00341802),
        "st_offset": band_info.get("offset", 149.0),
        "st_fill_value": 0,
    }
    meta_file.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"Salvato {meta_file.name}")


def main():
    cfg = load_config()
    out_dir = PROJECT_ROOT / cfg["paths"]["raw"]
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Giorno di riferimento: {cfg['reference_date']}")
    print(f"Area di addestramento: {cfg['training_aoi']['name']}\n")

    download_sentinel(cfg, out_dir)
    time.sleep(2)
    download_landsat(cfg, out_dir)
    print("\nDownload completato.")


if __name__ == "__main__":
    main()
