"""
step02_download.py - Scarica Sentinel-3 LST e Sentinel-2 L2A sulla stessa area e periodo.

Cosa fa:
  1. Si collega a openEO su Copernicus Data Space (CDSE) e fa il login
     con il tuo account Copernicus (si apre una pagina nel browser).
  2. Controlla che le bande richieste esistano davvero nelle collezioni.
  3. Lancia due elaborazioni sul cloud CDSE (batch job):
       - Sentinel-3 SLSTR: temperatura superficiale (LST) a ~1 km + flag di qualita'
       - Sentinel-2 L2A: rosso, vicino infrarosso e classificazione scena a 20 m
     entrambe riproiettate nella stessa proiezione metrica (UTM 33N).
  4. Salva i risultati in NetCDF in data/raw/step02 (esclusa da Git).

Come eseguirlo (dalla cartella del progetto, ambiente eo_water attivo):
    python scripts/step02_download.py

Note:
  - Usa pochi crediti openEO (area piccola, 15 giorni).
  - Se un file esiste gia' non viene riscaricato: per riscaricarlo cancellalo.
"""

import logging
import sys
from pathlib import Path

import openeo
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "step02_first_look.yaml"
BACKEND_URL = "openeo.dataspace.copernicus.eu"


def load_config(path):
    """Legge il file YAML di configurazione."""
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def connect():
    """Connessione e login a openEO su CDSE.

    authenticate_oidc() usa il 'device flow': stampa un link (e spesso apre il
    browser) dove fai login con l'account Copernicus Data Space. Il token viene
    salvato dalla libreria nella tua cartella utente, NON nel progetto: le volte
    successive il login sara' automatico.
    """
    print(f"Connessione a {BACKEND_URL} ...")
    connection = openeo.connect(BACKEND_URL)
    connection.authenticate_oidc()
    print("Login riuscito.\n")
    return connection


def check_bands(connection, collection_id, wanted_bands):
    """Verifica che le bande richieste esistano nella collezione (evita errori dopo minuti di attesa)."""
    meta = connection.describe_collection(collection_id)
    available = meta.get("cube:dimensions", {}).get("bands", {}).get("values", [])
    missing = [b for b in wanted_bands if b not in available]
    print(f"{collection_id}: {len(available)} bande disponibili")
    if missing:
        print(f"  ERRORE: bande non trovate: {missing}")
        print(f"  Bande disponibili: {available}")
        sys.exit(1)
    print(f"  OK: {wanted_bands}")


def run_job(cube, out_file, title):
    """Esegue un batch job sul cloud CDSE e scarica il risultato in NetCDF."""
    if out_file.exists():
        print(f"{out_file.name} esiste gia': salto il download.")
        return
    print(f"\nAvvio elaborazione: {title} (puo' richiedere qualche minuto)...")
    cube.execute_batch(outputfile=str(out_file), out_format="netCDF", title=title)
    size_mb = out_file.stat().st_size / 1e6
    print(f"Salvato {out_file.relative_to(PROJECT_ROOT)} ({size_mb:.1f} MB)")


class Tee:
    """Scrive l'output sia a schermo sia in un file di log (cosi' si puo' rileggere dopo)."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for s in self.streams:
            s.write(text)
            s.flush()

    def flush(self):
        for s in self.streams:
            s.flush()


def setup_logging():
    """Duplica stampa ed errori in outputs/logs/step02_download.log."""
    log_dir = PROJECT_ROOT / "outputs" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = open(log_dir / "step02_download.log", "w", encoding="utf-8")
    sys.stdout = Tee(sys.__stdout__, log_file)
    sys.stderr = Tee(sys.__stderr__, log_file)
    # Mostra anche i messaggi di avanzamento di openEO (stato dei batch job)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s: %(message)s", stream=sys.stdout)


def main():
    setup_logging()
    cfg = load_config(CONFIG_FILE)
    aoi = cfg["aoi"]
    spatial_extent = {k: aoi[k] for k in ("west", "south", "east", "north")}
    temporal_extent = [cfg["period"]["start"], cfg["period"]["end"]]
    crs = cfg["output_crs"]

    out_dir = PROJECT_ROOT / cfg["paths"]["raw"]
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Area: {aoi['description']}")
    print(f"Periodo: {temporal_extent[0]} -> {temporal_extent[1]}\n")

    connection = connect()

    s3 = cfg["sentinel3"]
    s2 = cfg["sentinel2"]
    check_bands(connection, s3["collection"], s3["bands"])
    check_bands(connection, s2["collection"], s2["bands"])

    # --- Sentinel-3 SLSTR: temperatura superficiale a ~1 km ---
    lst_cube = connection.load_collection(
        s3["collection"],
        spatial_extent=spatial_extent,
        temporal_extent=temporal_extent,
        bands=s3["bands"],
    )
    # "near" (vicino piu' prossimo): non mescola i valori e non altera i flag di qualita'
    lst_cube = lst_cube.resample_spatial(resolution=s3["resolution_m"], projection=crs, method="near")
    run_job(lst_cube, out_dir / "s3_lst.nc", "eo_water step02 - Sentinel-3 LST")

    # --- Sentinel-2 L2A: rosso, NIR e classificazione della scena a 20 m ---
    s2_cube = connection.load_collection(
        s2["collection"],
        spatial_extent=spatial_extent,
        temporal_extent=temporal_extent,
        bands=s2["bands"],
        max_cloud_cover=s2["max_cloud_cover"],
    )
    s2_cube = s2_cube.resample_spatial(resolution=s2["resolution_m"], projection=crs, method="near")
    run_job(s2_cube, out_dir / "s2_l2a.nc", "eo_water step02 - Sentinel-2 L2A")

    print("\nDownload completato. Prossimo comando: python scripts/step02_first_look.py")


if __name__ == "__main__":
    main()
