"""
p1_03_training_table.py - Costruisce la tabella di addestramento per lo sharpening.

Idea: il modello impara a 1 km la relazione  bande Sentinel-2  ->  LST Sentinel-3.
Per farlo servono, per ogni pixel da 1 km:
  - i PREDITTORI: le riflettanze Sentinel-2 mediate sui ~2500 pixel da 20 m che contiene;
  - il BERSAGLIO: la LST misurata da Sentinel-3 in quel pixel.

Passaggi:
  1. Sentinel-2 (20 m): riflettanze, NDVI e NDMI; maschera di "terra limpida" da SCL.
  2. Sentinel-3 (1 km): passaggio delle 09:58 UTC, maschera nuvole (bit 14).
  3. Sfocatura di Sentinel-2 come la "vede" Sentinel-3 (PSF gaussiana, sigma da config),
     poi aggregazione 20 m -> 1 km: media di ogni banda + frazione di pixel validi
     + variabilita' interna dell'NDVI (omogeneita').
  4. Filtro: tengo i pixel da 1 km con frazione valida >= soglia (config).
  5. Salvo la tabella e una figura di controllo.

Output:
  data/processed/p1_20220719/training_table.csv   (esclusa da Git)
  outputs/figures/p1_03_training_overview.png
  outputs/p1_03_summary.txt

Come eseguirlo:
    python scripts/p1_03_training_table.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray  # noqa: F401  (attiva .rio su xarray)
import xarray as xr
import yaml
from rasterio.enums import Resampling
from scipy.ndimage import gaussian_filter

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"
CLOUD_BIT_S3 = 14
REFLECTANCE_BANDS = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]

LINES = []


def log(text=""):
    print(text)
    LINES.append(text)


def open_nc(path):
    """Apre un NetCDF di openEO, lo legge in memoria e uniforma la dimensione temporale a 'time'."""
    ds = xr.open_dataset(path).load()
    ds.close()
    if "t" in ds.dims:
        ds = ds.rename({"t": "time"})
    return ds


def as_raster(da, crs):
    """Dichiara proiezione e 'nessun dato' = NaN: indispensabile prima di riproiettare/aggregare."""
    return da.rio.write_crs(crs).rio.write_nodata(np.nan, encoded=False)


def load_sentinel2(path, valid_scl, crs):
    """Riflettanze Sentinel-2 (0-1) con NaN fuori dalla terra limpida, piu' NDVI e NDMI."""
    s2 = open_nc(path)
    if s2.sizes["time"] > 1:
        # Piu' tasselli dello stesso giorno: li uniamo prendendo il primo valore valido
        log(f"Sentinel-2: {s2.sizes['time']} acquisizioni nel giorno, le unisco.")
    s2 = s2.max(dim="time", skipna=True) if s2.sizes["time"] > 1 else s2.isel(time=0)

    scl = s2["SCL"]
    valid = scl.isin(valid_scl)

    bands = {}
    for b in REFLECTANCE_BANDS:
        refl = s2[b].astype("float32") / 10000          # da valori x10000 a riflettanza 0-1
        bands[b] = refl.where(valid & (refl > 0))       # solo terra limpida e valori fisici
    ds = xr.Dataset(bands)
    ds["NDVI"] = (ds["B08"] - ds["B04"]) / (ds["B08"] + ds["B04"])
    ds["NDMI"] = (ds["B08"] - ds["B11"]) / (ds["B08"] + ds["B11"])
    ds["valid"] = valid.astype("float32")                # 1 = terra limpida, 0 = altro
    for v in ds.data_vars:
        ds[v] = as_raster(ds[v], crs)
    log(f"Sentinel-2: griglia {ds.sizes['y']} x {ds.sizes['x']} pixel da 20 m, "
        f"terra limpida {float(valid.mean()):.0%}")
    return ds


def load_sentinel3(path, time_utc, crs):
    """LST Sentinel-3 (Kelvin) del passaggio scelto, con le nuvole (bit 14) a NaN."""
    s3 = open_nc(path)
    times = pd.to_datetime(s3["time"].values)
    target = pd.Timestamp(f"{times[0].date()} {time_utc}")
    i = int(np.argmin(np.abs(times - target)))
    log(f"Sentinel-3: uso il passaggio delle {times[i]:%H:%M} UTC (disponibili: "
        f"{', '.join(t.strftime('%H:%M') for t in times)})")
    scene = s3.isel(time=i)
    flags = scene["confidence_in"].fillna(0).astype("int64")
    cloud = (flags & (1 << CLOUD_BIT_S3)) != 0
    lst = scene["LST"].where(~cloud)
    log(f"Sentinel-3: {int(scene['LST'].notnull().sum())} pixel con LST, "
        f"{int(lst.notnull().sum())} senza nuvole")
    return as_raster(lst, crs)


def blur_like_sensor(da, sigma_m, pixel_m):
    """Sfoca un'immagine con un filtro gaussiano, ignorando i pixel NaN.

    Tecnica della "convoluzione normalizzata": si sfocano separatamente
    (valori con NaN -> 0) e (maschera dei pixel validi), poi si divide.
    Cosi' i NaN non "contaminano" i vicini e i bordi non vengono scuriti.
    """
    if not sigma_m:
        return da
    sigma_px = sigma_m / pixel_m
    values = da.values
    valid = np.isfinite(values)
    num = gaussian_filter(np.where(valid, values, 0.0), sigma_px, mode="constant")
    den = gaussian_filter(valid.astype("float64"), sigma_px, mode="constant")
    with np.errstate(invalid="ignore", divide="ignore"):
        blurred = num / den
    blurred[den < 0.05] = np.nan          # troppo pochi pixel validi intorno: non affidabile
    return da.copy(data=blurred.astype("float32"))


def aggregate_to(fine, coarse):
    """Media dei pixel fini dentro ogni pixel grossolano (stessa proiezione, griglie diverse)."""
    return fine.rio.reproject_match(coarse, resampling=Resampling.average, nodata=np.nan)


def main():
    cfg = yaml.safe_load(open(CONFIG_FILE, encoding="utf-8"))
    tr = cfg["training"]
    crs = cfg["download"]["output_crs"]
    raw = PROJECT_ROOT / cfg["paths"]["raw"]

    s2 = load_sentinel2(raw / "s2_l2a.nc", tr["s2_valid_scl"], crs)
    lst = load_sentinel3(raw / "s3_lst.nc", tr["s3_time_utc"], crs)

    # --- Aggregazione 20 m -> 1 km ---
    sigma = tr.get("psf_sigma_m", 0)
    pixel_m = cfg["download"]["s2_resolution_m"]
    log(f"\nSfocatura PSF di Sentinel-2: sigma {sigma} m; aggrego sulla griglia di Sentinel-3 ...")
    agg = {}
    for name in REFLECTANCE_BANDS + ["NDVI", "NDMI"]:
        blurred = as_raster(blur_like_sensor(s2[name], sigma, pixel_m), crs)
        agg[name] = aggregate_to(blurred, lst)
    # Frazione valida: la media di una maschera 0/1 e' la percentuale di pixel validi
    # (fillna crea un nuovo array e perde la proiezione: va ridichiarata con as_raster)
    valid_fraction = aggregate_to(as_raster(s2["valid"].fillna(0), crs), lst)
    # Omogeneita': deviazione standard dell'NDVI dentro il pixel = sqrt(media(x^2) - media(x)^2)
    # (calcolata senza sfocatura: misura l'eterogeneita' reale dentro il chilometro)
    ndvi_mean = aggregate_to(s2["NDVI"], lst)
    ndvi_sq = aggregate_to(as_raster(s2["NDVI"] ** 2, crs), lst)
    ndvi_std = np.sqrt(np.maximum(ndvi_sq - ndvi_mean ** 2, 0))

    # --- Tabella: una riga per pixel da 1 km ---
    table = pd.DataFrame({name: arr.values.ravel() for name, arr in agg.items()})
    table["NDVI_std"] = ndvi_std.values.ravel()
    table["valid_fraction"] = valid_fraction.values.ravel()
    xx, yy = np.meshgrid(lst["x"].values, lst["y"].values)
    table["x"] = xx.ravel()
    table["y"] = yy.ravel()
    table["LST_K"] = lst.values.ravel()

    has_lst = table["LST_K"].notna()
    log(f"\nPixel da 1 km nella griglia: {len(table)}, con LST valida: {int(has_lst.sum())}")
    log("Pixel utilizzabili al variare della soglia di frazione valida:")
    for f in tr["sensitivity_fractions"]:
        n = int((has_lst & (table["valid_fraction"] >= f)).sum())
        log(f"  >= {f:.0%}: {n}")

    keep = has_lst & (table["valid_fraction"] >= tr["min_valid_fraction"]) & table["NDVI"].notna()
    training = table[keep].reset_index(drop=True)
    log(f"\nSoglia scelta {tr['min_valid_fraction']:.0%} -> {len(training)} pixel di addestramento")
    log(f"Predittori per pixel: {len(REFLECTANCE_BANDS) + 2} -> "
        f"{len(training) / (len(REFLECTANCE_BANDS) + 2):.0f} esempi per predittore")

    # --- Primo sguardo alla fisica: come varia la LST con ciascun predittore? ---
    log("\nCorrelazione di Spearman tra LST e predittori (segno atteso: NDVI e NDMI negativi):")
    corr = training[REFLECTANCE_BANDS + ["NDVI", "NDMI"]].corrwith(training["LST_K"], method="spearman")
    for name, r in corr.sort_values().items():
        log(f"  {name:5s} {r:+.2f}")

    out_dir = PROJECT_ROOT / cfg["paths"]["processed"]
    out_dir.mkdir(parents=True, exist_ok=True)
    training.to_csv(out_dir / "training_table.csv", index=False)
    log(f"\nTabella salvata: {(out_dir / 'training_table.csv').relative_to(PROJECT_ROOT)}")

    # --- Figura di controllo ---
    lst_c = training["LST_K"] - 273.15
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    im = axes[0].imshow(valid_fraction.values, cmap="viridis", vmin=0, vmax=1,
                        extent=[lst.x.min() - 500, lst.x.max() + 500, lst.y.min() - 500, lst.y.max() + 500])
    axes[0].set_title("Frazione di terra limpida per pixel da 1 km")
    fig.colorbar(im, ax=axes[0], shrink=0.8)
    axes[0].ticklabel_format(style="plain")
    axes[0].tick_params(axis="x", labelrotation=30)
    for ax, var in zip(axes[1:], ["NDVI", "NDMI"]):
        sc = ax.scatter(training[var], lst_c, c=training["NDVI_std"], s=10, cmap="magma_r")
        ax.set_xlabel(f"{var} medio nel pixel da 1 km")
        ax.set_ylabel("LST Sentinel-3 (°C)")
        ax.set_title(f"LST contro {var} (Spearman {corr[var]:+.2f})")
        ax.grid(alpha=0.3)
    fig.colorbar(sc, ax=axes[2], shrink=0.8, label="Deviazione std NDVI nel pixel (eterogeneita')")
    fig.suptitle(f"Tabella di addestramento - {cfg['reference_date']}, {len(training)} pixel da 1 km")
    fig_dir = PROJECT_ROOT / cfg["paths"]["figures"]
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / "p1_03_training_overview.png", dpi=150)
    plt.close(fig)
    log(f"Figura salvata: {cfg['paths']['figures']}/p1_03_training_overview.png")

    (PROJECT_ROOT / "outputs" / "p1_03_summary.txt").write_text("\n".join(LINES), encoding="utf-8")


if __name__ == "__main__":
    main()
