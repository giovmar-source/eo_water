"""
p1_05_sharpen.py - Sharpening della LST Sentinel-3 (1 km) a 20 m e validazione con Landsat 9.

Passaggi:
  1. Addestra i modelli sulla tabella a 1 km (predittori Sentinel-2 sfocati con la PSF).
  2. Li applica ai pixel Sentinel-2 da 20 m (non sfocati) -> LST "prevista" a 20 m.
  3. CORREZIONE DEI RESIDUI: sfoca la previsione con la stessa PSF, la aggrega a 1 km
     e la confronta con la LST misurata da Sentinel-3. La differenza (residuo) viene
     interpolata in modo morbido a 20 m e sommata: cosi' la media su ogni chilometro
     torna a coincidere con la misura del satellite (conservazione dell'energia osservata).
  4. VALIDAZIONE con Landsat 9 (termico misurato a 100 m), alla scala di 100 m:
       - confronto con la LST Sentinel-3 semplicemente interpolata (nessuno sharpening):
         e' il riferimento da battere;
       - metriche assolute (bias, RMSE) e "di forma" (correlazione, RMSE senza bias,
         variabilita' ricostruita dentro ogni chilometro).

Riusa le funzioni degli script p1_03 e p1_04 (stessa cartella).

Output:
  data/processed/p1_20220719/lst_sharpened_<modello>.tif   (20 m, Kelvin, esclusi da Git)
  outputs/tables/p1_05_validation.csv
  outputs/figures/p1_05_sharpening_maps.png
  outputs/figures/p1_05_validation_scatter.png
  outputs/p1_05_summary.txt

Come eseguirlo:
    python scripts/p1_05_sharpen.py
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray
import xarray as xr
import yaml
from rasterio.enums import Resampling

from p1_03_training_table import (REFLECTANCE_BANDS, aggregate_to, as_raster,
                                  blur_like_sensor, load_sentinel2, load_sentinel3)
from p1_04_models import build_models

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"
# Bit di QA_PIXEL di Landsat da escludere: 0 riempimento, 1 nuvola dilatata, 3 nuvola, 4 ombra, 7 acqua
LANDSAT_BAD_BITS = (1 << 0) | (1 << 1) | (1 << 3) | (1 << 4) | (1 << 7)

LINES = []


def log(text=""):
    print(text)
    LINES.append(text)


# ----------------------------------------------------------------------------
# 1-2. Addestramento e previsione a 20 m
# ----------------------------------------------------------------------------

def predict_fine(model, s2, predictors, chunk):
    """Applica il modello a ogni pixel Sentinel-2 valido (a blocchi, per non saturare la memoria)."""
    stack = np.stack([s2[p].values.ravel() for p in predictors], axis=1)
    ok = np.all(np.isfinite(stack), axis=1)
    out = np.full(stack.shape[0], np.nan, dtype="float32")
    idx = np.flatnonzero(ok)
    for start in range(0, idx.size, chunk):
        sel = idx[start:start + chunk]
        out[sel] = model.predict(pd.DataFrame(stack[sel], columns=predictors))
    template = s2[predictors[0]]
    return template.copy(data=out.reshape(template.shape))


# ----------------------------------------------------------------------------
# 3. Correzione dei residui
# ----------------------------------------------------------------------------

def residual_correction(pred20, lst1km, sigma_m, pixel_m, crs):
    """Residuo = LST misurata (1 km) - previsione sfocata e aggregata (1 km), riportato a 20 m."""
    pred_blur = as_raster(blur_like_sensor(pred20, sigma_m, pixel_m), crs)
    pred_1km = aggregate_to(pred_blur, lst1km)
    residual = (lst1km - pred_1km)
    # Dove manca il residuo (nuvole, mare) non correggiamo: residuo = 0
    residual = as_raster(residual.fillna(0.0), crs)
    # Interpolazione bilineare: il residuo varia in modo morbido, senza "quadrettoni" da 1 km
    residual20 = residual.rio.reproject_match(pred20, resampling=Resampling.bilinear)
    corrected = pred20 + residual20.values
    return as_raster(corrected, crs), residual


# ----------------------------------------------------------------------------
# 4. Validazione con Landsat
# ----------------------------------------------------------------------------

def load_landsat_kelvin(raw):
    """Temperatura Landsat in Kelvin, con nuvole/ombre/acqua/riempimento a NaN."""
    meta = json.loads((raw / "landsat_metadata.json").read_text(encoding="utf-8"))
    with rioxarray.open_rasterio(raw / "landsat_lwir11.tif") as src:
        st = src.squeeze("band", drop=True).load()
    with rioxarray.open_rasterio(raw / "landsat_qa_pixel.tif") as src:
        qa = src.squeeze("band", drop=True).load()
    kelvin = st.astype("float64") * meta["st_scale"] + meta["st_offset"]
    bad = (st == meta["st_fill_value"]) | ((qa.astype("int64") & LANDSAT_BAD_BITS) != 0)
    kelvin = kelvin.where(~bad).rio.write_crs(st.rio.crs).rio.write_nodata(np.nan, encoded=False)
    log(f"Landsat: {meta['item_id']}, {meta['datetime_utc'][:16]} UTC, "
        f"pixel validi {float(kelvin.notnull().mean()):.0%}")
    return kelvin


def to_grid(da, grid):
    """Media sulla griglia di validazione + frazione di pixel validi in ogni cella."""
    da = da.rio.write_nodata(np.nan, encoded=False)
    mean = da.rio.reproject_match(grid, resampling=Resampling.average, nodata=np.nan)
    frac = (da.notnull().astype("float32").rio.write_crs(da.rio.crs).rio.write_nodata(-1, encoded=False)
            .rio.reproject_match(grid, resampling=Resampling.average, nodata=-1))
    return mean, frac.where(frac >= 0, 0)


def metrics(est, ref):
    err = est - ref
    bias = float(np.mean(err))
    return {
        "n": int(err.size),
        "bias_K": bias,
        "RMSE_K": float(np.sqrt(np.mean(err ** 2))),
        "RMSE_senza_bias_K": float(np.sqrt(np.mean((err - bias) ** 2))),  # errore "di forma"
        "r": float(np.corrcoef(est, ref)[0, 1]),
    }


def within_km_std(values100, grid100, lst1km):
    """Deviazione standard dei valori a 100 m dentro ogni pixel da 1 km (variabilita' interna)."""
    v = as_raster(grid100.copy(data=values100), grid100.rio.crs)
    mean = v.rio.reproject_match(lst1km, resampling=Resampling.average, nodata=np.nan)
    sq = as_raster(v ** 2, grid100.rio.crs).rio.reproject_match(lst1km, resampling=Resampling.average, nodata=np.nan)
    return np.sqrt(np.maximum(sq - mean ** 2, 0)).values.ravel()


def main():
    cfg = yaml.safe_load(open(CONFIG_FILE, encoding="utf-8"))
    tr, sh, m = cfg["training"], cfg["sharpening"], cfg["models"]
    crs = cfg["download"]["output_crs"]
    pixel_m = cfg["download"]["s2_resolution_m"]
    sigma = tr.get("psf_sigma_m", 0)
    raw = PROJECT_ROOT / cfg["paths"]["raw"]
    proc = PROJECT_ROOT / cfg["paths"]["processed"]

    # --- Dati ---
    s2 = load_sentinel2(raw / "s2_l2a.nc", tr["s2_valid_scl"], crs)
    lst = load_sentinel3(raw / "s3_lst.nc", tr["s3_time_utc"], crs)
    table = pd.read_csv(proc / "training_table.csv")
    X, y = table[m["predictors"]], table["LST_K"]
    landsat = load_landsat_kelvin(raw)

    # --- Griglia di validazione a 100 m, allineata a Sentinel-2 ---
    grid100 = as_raster(s2["B04"], crs).rio.reproject(
        f"EPSG:{crs}", resolution=sh["validation_resolution_m"], resampling=Resampling.average, nodata=np.nan)
    land100, land_frac = to_grid(landsat, grid100)

    # Riferimento da battere: Sentinel-3 interpolato (bilineare) senza sharpening
    s3_interp = lst.rio.reproject_match(grid100, resampling=Resampling.bilinear, nodata=np.nan)

    results, maps = [], {}
    models = build_models(m, len(m["predictors"]))
    for name in sh["models"]:
        model, cols = models[name]
        log(f"\n=== {name} ===")
        model.fit(X[cols], y)
        pred20 = predict_fine(model, s2, cols, sh["predict_chunk"])
        log(f"Previsione a 20 m: {int(np.isfinite(pred20.values).sum())} pixel")
        sharp20, residual = residual_correction(as_raster(pred20, crs), lst, sigma, pixel_m, crs)
        log(f"Residui a 1 km: media {float(residual.where(residual != 0).mean()):+.2f} K, "
            f"deviazione std {float(residual.where(residual != 0).std()):.2f} K")
        maps[name] = (as_raster(pred20, crs), sharp20)
        sharp20.rio.to_raster(proc / f"lst_sharpened_{name.replace(' ', '_').lower()}.tif")

        for label, fine in [(f"{name} (senza correzione)", as_raster(pred20, crs)),
                            (f"{name} (con correzione)", sharp20)]:
            est100, est_frac = to_grid(fine, grid100)
            mask = ((est_frac >= sh["min_valid_fraction"]) & (land_frac >= sh["min_valid_fraction"])
                    & est100.notnull() & land100.notnull() & s3_interp.notnull())
            results.append({"metodo": label, **metrics(est100.values[mask.values], land100.values[mask.values])})
            if label.endswith("(con correzione)"):
                maps[name] += (est100, mask)

    # Riferimento senza sharpening, sulle stesse celle del primo modello
    mask = maps[sh["models"][0]][3]
    base = {"metodo": "Sentinel-3 interpolato (nessuno sharpening)",
            **metrics(s3_interp.values[mask.values], land100.values[mask.values])}
    results.insert(0, base)

    res = pd.DataFrame(results)
    log("\nValidazione con Landsat 9 a 100 m (stesse celle per tutti i metodi):")
    with pd.option_context("display.width", 160):
        log(res.round(3).to_string(index=False))

    # --- Variabilita' ricostruita dentro ogni chilometro ---
    log("\nVariabilita' dentro ogni pixel da 1 km (deviazione std dei valori a 100 m, mediana):")
    common = mask.values
    ref_std = within_km_std(np.where(common, land100.values, np.nan), grid100, lst)
    log(f"  Landsat (vera):            {np.nanmedian(ref_std):.2f} K")
    for name in sh["models"]:
        est100 = maps[name][2]
        est_std = within_km_std(np.where(common, est100.values, np.nan), grid100, lst)
        log(f"  {name:25s} {np.nanmedian(est_std):.2f} K  "
            f"(rapporto ricostruita/vera {np.nanmedian(est_std) / np.nanmedian(ref_std):.0%})")

    tab_dir = PROJECT_ROOT / cfg["paths"]["tables"]
    res.to_csv(tab_dir / "p1_05_validation.csv", index=False)

    # --- Figure ---
    fig_dir = PROJECT_ROOT / cfg["paths"]["figures"]
    best = sh["models"][0]
    a = cfg["analysis_aoi"]
    from pyproj import Transformer
    tf = Transformer.from_crs(4326, crs, always_xy=True)
    x0, y0 = tf.transform(a["west"], a["south"])
    x1, y1 = tf.transform(a["east"], a["north"])

    def crop(da):
        return da.sel(x=slice(x0, x1), y=slice(y1, y0)) if da.y[0] > da.y[-1] else da.sel(x=slice(x0, x1), y=slice(y0, y1))

    land_utm = landsat.rio.reproject(f"EPSG:{crs}", resolution=30, nodata=np.nan)
    panels = [
        ("Sentinel-3 misurata, ~1 km", crop(lst)),
        (f"Sharpening 20 m ({best}, con correzione)", crop(maps[best][1])),
        ("Landsat 9 (riferimento), 30 m", crop(land_utm)),
        ("Sentinel-2 NDMI, 20 m", crop(s2["NDMI"])),
    ]
    vals = np.concatenate([np.ravel(p[1].values - 273.15) for p in panels[:3]])
    vmin, vmax = np.nanpercentile(vals, [2, 98])
    fig, axes = plt.subplots(2, 2, figsize=(13, 12), constrained_layout=True)
    for ax, (title, da) in zip(axes.ravel(), panels):
        ext = [float(da.x.min()), float(da.x.max()), float(da.y.min()), float(da.y.max())]
        if "NDMI" in title:
            im = ax.imshow(da.values, extent=ext, cmap="BrBG", vmin=-0.3, vmax=0.4)
            fig.colorbar(im, ax=ax, shrink=0.8, label="NDMI")
        else:
            im = ax.imshow(da.values - 273.15, extent=ext, cmap="inferno", vmin=vmin, vmax=vmax)
            fig.colorbar(im, ax=ax, shrink=0.8, label="Temperatura superficiale (°C)")
        ax.set_title(title)
        ax.ticklabel_format(style="plain")
        ax.tick_params(axis="x", labelrotation=30)
    fig.suptitle(f"Delta del Po, {cfg['reference_date']} ~09:58 UTC - area di analisi")
    fig.savefig(fig_dir / "p1_05_sharpening_maps.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    ref = land100.values[mask.values] - 273.15
    rows = [res.iloc[0], res[res.metodo == f"{best} (con correzione)"].iloc[0]]
    for ax, row, (label, est) in zip(axes, rows, [
            ("Sentinel-3 interpolato", s3_interp.values[mask.values] - 273.15),
            (f"Sharpening ({best})", maps[best][2].values[mask.values] - 273.15)]):
        ax.hexbin(ref, est, gridsize=60, mincnt=1, cmap="viridis", bins="log")
        lims = [np.nanpercentile(ref, 0.5), np.nanpercentile(ref, 99.5)]
        ax.plot(lims, lims, "w--", lw=1)
        ax.set_xlim(lims); ax.set_ylim(lims)
        ax.set_xlabel("Landsat 9 a 100 m (°C)")
        ax.set_ylabel(f"{label} a 100 m (°C)")
        ax.set_title(f"{label}\nr {row.r:.2f}, RMSE {row.RMSE_K:.2f} K, senza bias {row.RMSE_senza_bias_K:.2f} K")
    fig.savefig(fig_dir / "p1_05_validation_scatter.png", dpi=150)
    plt.close(fig)

    log(f"\nMappe: {cfg['paths']['figures']}/p1_05_sharpening_maps.png")
    log(f"Confronto: {cfg['paths']['figures']}/p1_05_validation_scatter.png")
    (PROJECT_ROOT / "outputs" / "p1_05_summary.txt").write_text("\n".join(LINES), encoding="utf-8")


if __name__ == "__main__":
    main()
