"""
check_env.py - Verifica dell'ambiente di lavoro (Fase 0, Step 1)

Cosa fa:
  1. Controlla che tutte le librerie si importino e ne stampa la versione.
  2. Esegue un piccolo test reale: crea un raster georeferenziato,
     lo salva come GeoTIFF e lo rilegge (verifica GDAL/rasterio/rioxarray).
  3. Controlla che il file di configurazione del Climate Data Store esista
     (NON ne legge il contenuto: le credenziali restano private).
  4. Salva tutto l'output anche in check_env_output.txt, cosi' posso leggerlo.

Come eseguirlo (terminale con l'ambiente eo_water attivo):
    python check_env.py
"""

import importlib
import sys
import tempfile
from pathlib import Path

# Librerie da controllare: nome usato in "import" -> nome leggibile
PACKAGES = {
    "numpy": "NumPy",
    "pandas": "Pandas",
    "scipy": "SciPy",
    "statsmodels": "Statsmodels",
    "sklearn": "scikit-learn",
    "xarray": "Xarray",
    "dask": "Dask",
    "netCDF4": "netCDF4",
    "rasterio": "Rasterio",
    "rioxarray": "rioxarray",
    "geopandas": "GeoPandas",
    "shapely": "Shapely",
    "pyproj": "pyproj",
    "matplotlib": "Matplotlib",
    "cartopy": "Cartopy",
    "pystac_client": "pystac-client",
    "odc.stac": "odc-stac",
    "openeo": "openEO",
    "cdsapi": "cdsapi",
    "earthaccess": "earthaccess",
    "yaml": "PyYAML",
}

LINES = []  # raccoglie l'output per salvarlo su file


def log(text=""):
    """Stampa a schermo e memorizza la riga."""
    print(text)
    LINES.append(text)


def check_imports():
    """Prova a importare ogni libreria e restituisce l'elenco di quelle mancanti."""
    log("=" * 60)
    log(f"Python {sys.version.split()[0]}  ({sys.executable})")
    log("=" * 60)
    missing = []
    for module_name, label in PACKAGES.items():
        try:
            module = importlib.import_module(module_name)
            version = getattr(module, "__version__", "versione non indicata")
            log(f"  OK     {label:<15} {version}")
        except Exception as error:  # ImportError o errori di DLL su Windows
            log(f"  ERRORE {label:<15} -> {error}")
            missing.append(label)
    return missing


def check_raster_io():
    """Crea un piccolo raster georeferenziato, lo scrive e lo rilegge."""
    import numpy as np
    import rioxarray  # noqa: F401  (attiva l'accessor .rio su xarray)
    import xarray as xr

    log("\nTest lettura/scrittura GeoTIFF...")
    # Griglia 10x10 con coordinate in metri (EPSG:3035, proiezione europea LAEA)
    data = np.arange(100, dtype="float32").reshape(10, 10)
    x = 4_300_000 + np.arange(10) * 20.0   # passo 20 m
    y = 2_400_000 - np.arange(10) * 20.0
    da = xr.DataArray(data, coords={"y": y, "x": x}, dims=("y", "x"))
    da = da.rio.write_crs("EPSG:3035")

    # ignore_cleanup_errors: su Windows un file ancora aperto non si puo' cancellare;
    # cosi' un eventuale residuo nella cartella temporanea non blocca il test.
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        path = Path(tmp) / "test.tif"
        da.rio.to_raster(path)
        # "with" garantisce che il file venga chiuso; .load() legge i dati in
        # memoria prima della chiusura (rioxarray altrimenti legge "pigramente").
        with rioxarray.open_rasterio(path) as src:
            back = src.squeeze().load()
        same_values = bool(np.allclose(back.values, data))
        same_crs = back.rio.crs.to_epsg() == 3035

    if same_values and same_crs:
        log("  OK     GeoTIFF scritto e riletto correttamente (valori e CRS)")
        return True
    log("  ERRORE il raster riletto non coincide con l'originale")
    return False


def check_cds_config():
    """Controlla solo che il file .cdsapirc esista, senza leggerne il contenuto."""
    log("\nConfigurazione Climate Data Store...")
    cfg = Path.home() / ".cdsapirc"
    if cfg.is_file():
        log("  OK     file .cdsapirc trovato nella cartella utente")
        return True
    log("  MANCA  file .cdsapirc nella cartella utente (Parte 5 delle istruzioni)")
    return False


if __name__ == "__main__":
    missing = check_imports()
    raster_ok = check_raster_io() if not missing else False
    cds_ok = check_cds_config()

    log("\n" + "=" * 60)
    if not missing and raster_ok and cds_ok:
        log("TUTTO PRONTO: l'ambiente funziona.")
    else:
        if missing:
            log("Librerie con problemi: " + ", ".join(missing))
        if not raster_ok:
            log("Test raster non superato o non eseguito.")
        if not cds_ok:
            log("File .cdsapirc da creare.")
    log("=" * 60)

    out = Path(__file__).with_name("check_env_output.txt")
    out.write_text("\n".join(LINES), encoding="utf-8")
    print(f"\nOutput salvato in {out.name}")
