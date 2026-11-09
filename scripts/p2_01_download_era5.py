"""
p2_01_download_era5.py - Scarica il meteo ERA5-Land del giorno di riferimento e lo controlla.

Perche': il bilancio energetico (TSEB) ha bisogno, oltre alla temperatura superficiale,
delle condizioni dell'aria al momento del passaggio del satellite:
  - temperatura dell'aria e umidita' (punto di rugiada -> pressione di vapore),
  - vento (per la resistenza aerodinamica),
  - pressione,
  - radiazione solare e infrarossa in arrivo (energia disponibile).
Servono anche i TOTALI GIORNALIERI di radiazione, per passare da ET istantanea a giornaliera.

ERA5-Land: rianalisi ECMWF, ~9 km (0,1°), oraria. Si scarica dal Climate Data Store con
la libreria cdsapi, usando il token salvato in .cdsapirc (cartella utente).

ATTENZIONE alle variabili di radiazione di ERA5-Land: sono CUMULATE dalle 00 UTC
(J/m^2 sommati dall'inizio del giorno). Per avere il flusso medio di un'ora si fa la
differenza tra due ore consecutive e si divide per 3600 secondi. Il valore alle 00 UTC
del giorno dopo contiene il totale dell'intera giornata.

ATTENZIONE 2: ERA5-Land copre SOLO LA TERRA. Le celle classificate "mare" (anche gran
parte del Delta del Po, sottile lingua di terra tra mare e lagune) non hanno valori.
Le riempiamo con il valore della CELLA DI TERRA PIU' VICINA (le celle costiere sono le
piu' rappresentative del Delta). E' un'approssimazione da dichiarare e verificare con
stazioni meteo locali.

Output:
  data/raw/p1_20220719/era5land.nc               (dato originale)
  data/processed/p1_20220719/era5land_filled.nc  (celle mare riempite)
  outputs/p2_01_summary.txt

Come eseguirlo:
    python scripts/p2_01_download_era5.py
"""

import zipfile
from pathlib import Path

import cdsapi
import numpy as np
import pandas as pd
import xarray as xr
import yaml

from config_utils import load_config
from scipy.ndimage import distance_transform_edt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"

LINES = []


def log(text=""):
    print(text)
    LINES.append(text)


def download(cfg, out_file):
    """Richiesta al Climate Data Store: il giorno di riferimento (tutte le ore) + le 00 UTC del giorno dopo."""
    met = cfg["meteo"]
    aoi = cfg["training_aoi"]
    mg = met["margin_deg"]
    # Il CDS vuole l'area come [nord, ovest, sud, est]
    area = [aoi["north"] + mg, aoi["west"] - mg, aoi["south"] - mg, aoi["east"] + mg]
    day = pd.Timestamp(cfg["reference_date"])
    nxt = day + pd.Timedelta(days=1)
    # Due giorni interi (pochi MB): ci serve anche il cumulato delle 00 UTC del giorno dopo
    request = {
        "variable": met["variables"],
        "year": [str(day.year)],
        "month": sorted({f"{day.month:02d}", f"{nxt.month:02d}"}),
        "day": sorted({f"{day.day:02d}", f"{nxt.day:02d}"}),
        "time": [f"{h:02d}:00" for h in range(24)],
        "area": area,
        "data_format": "netcdf",
        "download_format": "unarchived",
    }
    print("Richiesta al Climate Data Store (puo' restare in coda qualche minuto) ...")
    client = cdsapi.Client()
    client.retrieve(met["dataset"], request).download(str(out_file))

    # A volte il CDS consegna uno zip con piu' NetCDF dentro: in quel caso li unisce
    if zipfile.is_zipfile(out_file):
        print("Il CDS ha consegnato uno zip: estraggo e unisco i file ...")
        tmp_dir = out_file.with_suffix("")
        with zipfile.ZipFile(out_file) as z:
            z.extractall(tmp_dir)
        parts = [xr.open_dataset(f).load() for f in sorted(tmp_dir.glob("*.nc"))]
        merged = xr.merge(parts, compat="override")
        out_file.unlink()
        merged.to_netcdf(out_file)


def open_era5(path):
    """Apre il NetCDF e uniforma il nome della dimensione temporale a 'time'."""
    ds = xr.open_dataset(path).load()
    ds.close()
    for name in ("valid_time", "t"):
        if name in ds.dims:
            ds = ds.rename({name: "time"})
    return ds


def fill_from_nearest_land(ds):
    """Riempie le celle senza dati (mare) con il valore della cella di terra piu' vicina.

    distance_transform_edt, per ogni cella "mare", restituisce gli indici (riga, colonna)
    della cella "terra" piu' vicina: basta poi copiare i valori da quegli indici.
    La maschera di terra e' la stessa per tutte le ore e le variabili.
    """
    land = np.isfinite(ds["t2m"].isel(time=0).values)
    _, (rows, cols) = distance_transform_edt(~land, return_indices=True)
    filled = ds.copy()
    for name, da in ds.data_vars.items():
        if set(("latitude", "longitude")) <= set(da.dims):
            v = da.transpose(..., "latitude", "longitude").values
            filled[name] = da.transpose(..., "latitude", "longitude").copy(data=v[..., rows, cols])
    filled["land_mask"] = (("latitude", "longitude"), land.astype("int8"))
    return filled, land


def main():
    cfg = load_config()
    out_dir = PROJECT_ROOT / cfg["paths"]["raw"]
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "era5land.nc"
    if out_file.exists():
        print(f"{out_file.name} esiste gia': salto il download.")
    else:
        download(cfg, out_file)

    raw_ds = open_era5(out_file)
    ds, land = fill_from_nearest_land(raw_ds)
    log(f"Celle ERA5-Land con dati (terra): {int(land.sum())} su {land.size}; "
        f"le altre riempite con la cella di terra piu' vicina")
    proc = PROJECT_ROOT / cfg["paths"]["processed"]
    proc.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(proc / "era5land_filled.nc")
    log(f"Salvato {cfg['paths']['processed']}/era5land_filled.nc")
    log(f"Variabili: {list(raw_ds.data_vars)}")
    log(f"Ore: da {pd.Timestamp(ds.time.values[0])} a {pd.Timestamp(ds.time.values[-1])}, "
        f"griglia {ds.sizes.get('latitude')} x {ds.sizes.get('longitude')}")

    # --- Controllo: valori al centro dell'area di analisi ---
    a = cfg["analysis_aoi"]
    lat_c, lon_c = (a["south"] + a["north"]) / 2, (a["west"] + a["east"]) / 2
    pt = ds.sel(latitude=lat_c, longitude=lon_c, method="nearest")
    raw_pt = raw_ds.sel(latitude=lat_c, longitude=lon_c, method="nearest")
    if np.isnan(float(raw_pt["t2m"].isel(time=0))):
        log(f"Il centro dell'area ({lat_c:.2f}N, {lon_c:.2f}E) cade in una cella 'mare': "
            f"uso i valori della cella di terra piu' vicina.")
    day = pd.Timestamp(cfg["reference_date"])
    h = cfg["meteo"]["overpass_hour_utc"]
    t_now = day + pd.Timedelta(hours=h)
    t_prev = t_now - pd.Timedelta(hours=1)
    now, prev = pt.sel(time=t_now), pt.sel(time=t_prev)

    ta_c = float(now["t2m"]) - 273.15
    td_c = float(now["d2m"]) - 273.15
    # Pressione di vapore dal punto di rugiada (formula di Tetens), in hPa
    ea = 6.112 * np.exp(17.67 * td_c / (td_c + 243.5))
    es = 6.112 * np.exp(17.67 * ta_c / (ta_c + 243.5))
    wind = float(np.hypot(now["u10"], now["v10"]))
    # Radiazione: differenza dei cumulati tra due ore / 3600 s = flusso medio in W/m^2
    sw = (float(now["ssrd"]) - float(prev["ssrd"])) / 3600
    lw = (float(now["strd"]) - float(prev["strd"])) / 3600
    # Totale del giorno = cumulato alle 00 UTC del giorno dopo
    sw_day = float(pt["ssrd"].sel(time=day + pd.Timedelta(days=1))) / 1e6

    log(f"\nCentro dell'area di analisi, ore {t_prev:%H}-{t_now:%H} UTC del {day:%d/%m/%Y}:")
    log(f"  Temperatura aria 2 m:        {ta_c:5.1f} °C")
    log(f"  Punto di rugiada:            {td_c:5.1f} °C  -> umidita' relativa {100 * ea / es:4.0f}%")
    log(f"  Pressione di vapore:         {ea:5.1f} hPa")
    log(f"  Vento 10 m:                  {wind:5.1f} m/s")
    log(f"  Pressione al suolo:          {float(now['sp']) / 100:6.1f} hPa")
    log(f"  Radiazione solare (media):   {sw:5.0f} W/m²")
    log(f"  Radiazione infrarossa:       {lw:5.0f} W/m²")
    log(f"  Radiazione solare del giorno: {sw_day:4.1f} MJ/m²")
    log("\nControlli di plausibilita' (luglio, pianura padana, cielo sereno, meta' mattina):")
    log("  aria 27-33 °C, radiazione solare 650-850 W/m², totale giornaliero 25-30 MJ/m².")

    # Gradiente costa-entroterra della temperatura dell'aria alle 10 UTC (solo celle di terra)
    t10 = raw_ds["t2m"].sel(time=t_now) - 273.15
    log(f"\nTemperatura dell'aria alle {h} UTC sulle celle di terra: "
        f"da {float(t10.min()):.1f} °C (costa) a {float(t10.max()):.1f} °C (entroterra)")

    (PROJECT_ROOT / cfg["paths"]["summaries"] / "p2_01_summary.txt").write_text("\n".join(LINES), encoding="utf-8")


if __name__ == "__main__":
    main()
