"""
p2_03_biophysical.py - LAI, copertura (FCOVER) e FAPAR a 20 m con il processore biofisico di SNAP.

Il processore di ESA SNAP (Weiss & Baret) e' una piccola rete neurale, addestrata su
milioni di simulazioni del modello fisico PROSAIL (come la luce interagisce con foglie e suolo).
Per ogni pixel:
  input  : 8 riflettanze (B03, B04, B05, B06, B07, B8A, B11, B12)
           + cos(zenit di vista), cos(zenit solare), cos(azimut relativo)      -> 11 valori
  1. normalizzazione di ogni input tra -1 e +1 (con i minimi/massimi del file Normalisation)
  2. strato nascosto: 5 neuroni, h = tanh(W1 . x + b1)
  3. uscita: y = W2 . h + b2
  4. de-normalizzazione nell'unita' fisica (es. LAI in m2/m2)
  5. controlli: input fuori dal dominio di addestramento -> NaN; uscite leggermente fuori
     dai limiti fisici -> riportate al limite; molto fuori -> NaN.

Poi confronta il LAI di SNAP (8 bande) con il LAI del servizio BIOPAR (3 bande):
ci aspettiamo accordo sulla vegetazione rada e differenze sulla vegetazione fitta
(saturazione delle bande visibili/NIR).

Output:
  data/processed/p1_20220719/{lai,fcover,fapar}_snap.tif   (20 m, esclusi da Git)
  outputs/figures/p2_03_biophysical.png
  outputs/p2_03_summary.txt

Come eseguirlo:
    python scripts/p2_03_biophysical.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import yaml

from p1_03_training_table import as_raster, open_nc

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"

LINES = []


def log(text=""):
    print(text)
    LINES.append(text)


# ----------------------------------------------------------------------------
# Rete neurale SNAP
# ----------------------------------------------------------------------------

def load_snap_network(folder, var):
    """Legge i coefficienti di una rete SNAP dai file di testo (valori separati da virgola)."""
    def read(part):
        path = folder / var / f"{var}_{part}"
        if not path.exists() and part == "DefinitionDomain_MinMax":
            # FCOVER non ha questo file in SNAP: il dominio degli input e' lo stesso per tutte
            # le reti (verificato: LAI, FAPAR e LAI_Cab hanno file identici)
            path = folder / "LAI" / "LAI_DefinitionDomain_MinMax"
        return np.loadtxt(path, delimiter=",", ndmin=2)
    return {
        "norm": read("Normalisation"),              # 11 righe: [min, max] di ogni input
        "denorm": read("Denormalisation").ravel(),   # [min, max] dell'uscita
        "w1": read("Weights_Layer1_Neurons"),        # 5 neuroni x 11 input
        "b1": read("Weights_Layer1_Bias").ravel(),   # 5 valori
        "w2": read("Weights_Layer2_Neurons").ravel(),  # 5 valori
        "b2": float(read("Weights_Layer2_Bias").ravel()[0]),
        "extreme": read("ExtremeCases").ravel(),     # [tolleranza, min, max]
        "domain": read("DefinitionDomain_MinMax"),   # 2 righe x 8 bande: [min; max]
    }


def run_snap_network(net, X):
    """Applica la rete a una matrice X (pixel x 11 input). Restituisce un vettore di uscite."""
    xmin, xmax = net["norm"][:, 0], net["norm"][:, 1]
    xn = 2 * (X - xmin) / (xmax - xmin) - 1                  # 1. normalizzazione in [-1, 1]
    hidden = np.tanh(xn @ net["w1"].T + net["b1"])           # 2. strato nascosto (tanh)
    y = hidden @ net["w2"] + net["b2"]                       # 3. uscita lineare
    dmin, dmax = net["denorm"]
    out = 0.5 * (y + 1) * (dmax - dmin) + dmin               # 4. unita' fisiche

    # 5a. Input fuori dal dominio su cui la rete e' stata addestrata -> non affidabile
    bands = X[:, :8]
    outside = np.any((bands < net["domain"][0]) | (bands > net["domain"][1]), axis=1)
    out[outside] = np.nan
    # 5b. Uscite fuori dai limiti fisici: entro la tolleranza si riportano al limite, oltre -> NaN
    tol, lo, hi = net["extreme"]
    tol = abs(tol)
    out[(out < lo - tol) | (out > hi + tol)] = np.nan
    out = np.clip(out, lo, hi)
    return out, outside


def main():
    cfg = yaml.safe_load(open(CONFIG_FILE, encoding="utf-8"))
    bio, tr = cfg["biophysical"], cfg["training"]
    crs = cfg["download"]["output_crs"]
    raw = PROJECT_ROOT / cfg["paths"]["raw"]
    proc = PROJECT_ROOT / cfg["paths"]["processed"]
    coef_dir = PROJECT_ROOT / cfg["paths"]["aux"] / f"snap_biophys_{bio['snap_version']}"

    # --- Sentinel-2: riflettanze 0-1 sui pixel di terra limpida ---
    s2 = open_nc(raw / "s2_l2a.nc")
    s2 = s2.max(dim="time", skipna=True) if s2.sizes["time"] > 1 else s2.isel(time=0)
    valid = s2["SCL"].isin(tr["s2_valid_scl"])
    refl = {b: (s2[b].astype("float32") / 10000).where(valid) for b in bio["snap_bands"] + ["B08"]}
    ndvi = (refl["B08"] - refl["B04"]) / (refl["B08"] + refl["B04"])

    # --- Angoli (gradi) sulla stessa griglia ---
    ang = open_nc(raw / "s2_angles.nc")
    ang = ang.max(dim="time", skipna=True) if ang.sizes["time"] > 1 else ang.isel(time=0)
    ang = ang.interp_like(s2["B04"], method="nearest")
    sza, saa = ang["sunZenithAngles"], ang["sunAzimuthAngles"]
    vza, vaa = ang["viewZenithMean"], ang["viewAzimuthMean"]
    log(f"Angoli medi: zenit solare {float(sza.mean()):.1f}°, azimut solare {float(saa.mean()):.1f}°, "
        f"zenit di vista {float(vza.mean()):.1f}°, azimut di vista {float(vaa.mean()):.1f}°")

    # --- Matrice degli input: una riga per pixel, 11 colonne ---
    cols = [refl[b].values.ravel() for b in bio["snap_bands"]]
    cols += [np.cos(np.radians(vza.values)).ravel(),
             np.cos(np.radians(sza.values)).ravel(),
             np.cos(np.radians(saa.values - vaa.values)).ravel()]
    X = np.stack(cols, axis=1)
    ok = np.all(np.isfinite(X), axis=1)
    log(f"Pixel di terra limpida con tutti gli input: {int(ok.sum())}")

    template = s2["B04"]
    results = {}
    for var in bio["snap_variables"]:
        net = load_snap_network(coef_dir, var)
        out = np.full(X.shape[0], np.nan, dtype="float32")
        vals, outside = run_snap_network(net, X[ok])
        out[ok] = vals
        da = template.copy(data=out.reshape(template.shape))
        # Il "modello" B04 si porta dietro attributi del NetCDF (_FillValue, scala...) che non
        # valgono per il LAI: li togliamo, altrimenti la scrittura del GeoTIFF va in conflitto.
        da.attrs, da.encoding, da.name = {}, {}, var
        da = as_raster(da, crs)
        results[var] = da
        da.rio.to_raster(proc / f"{var.lower()}_snap.tif")
        log(f"{var:6s}: mediana {float(da.median()):.2f}, 5°-95° percentile "
            f"{float(da.quantile(0.05)):.2f}-{float(da.quantile(0.95)):.2f}, "
            f"fuori dominio di addestramento {outside.mean():.1%}")

    lai = results["LAI"]

    # --- Confronto con BIOPAR (se disponibile) ---
    biopar_file = raw / "biopar_lai.nc"
    biopar = None
    if biopar_file.exists():
        b = open_nc(biopar_file)
        if "time" in b.dims:
            b = b.max(dim="time", skipna=True) if b.sizes["time"] > 1 else b.isel(time=0)
        name = bio["biopar_type"] if bio["biopar_type"] in b.data_vars else list(b.data_vars)[-1]
        biopar = b[name].interp_like(template, method="nearest").where(valid)
        log(f"\nBIOPAR {name}: mediana {float(biopar.median()):.2f}")
        df = pd.DataFrame({"NDVI": ndvi.values.ravel(), "SNAP": lai.values.ravel(),
                           "BIOPAR": biopar.values.ravel()}).dropna()
        df["classe NDVI"] = pd.cut(df["NDVI"], [-1, 0.3, 0.5, 0.7, 0.8, 0.85, 1.0])
        summary = df.groupby("classe NDVI", observed=True).agg(
            pixel=("SNAP", "size"), LAI_SNAP=("SNAP", "median"), LAI_BIOPAR=("BIOPAR", "median"))
        summary["differenza"] = summary["LAI_SNAP"] - summary["LAI_BIOPAR"]
        log("\nLAI mediano per classe di NDVI (SNAP 8 bande vs BIOPAR 3 bande):")
        log(summary.round(2).to_string())
    else:
        log("\nbiopar_lai.nc non trovato: salto il confronto con BIOPAR.")

    # --- Figura ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)
    ext = [float(template.x.min()), float(template.x.max()), float(template.y.min()), float(template.y.max())]
    im = axes[0, 0].imshow(lai.values, extent=ext, cmap="YlGn", vmin=0, vmax=6)
    axes[0, 0].set_title("LAI - SNAP (8 bande + angoli), 20 m")
    fig.colorbar(im, ax=axes[0, 0], shrink=0.8, label="LAI (m²/m²)")
    im = axes[0, 1].imshow(results["FCOVER"].values, extent=ext, cmap="YlGn", vmin=0, vmax=1)
    axes[0, 1].set_title("Frazione di copertura (FCOVER) - SNAP")
    fig.colorbar(im, ax=axes[0, 1], shrink=0.8, label="FCOVER")
    for ax in axes[0]:
        ax.ticklabel_format(style="plain"); ax.tick_params(axis="x", labelrotation=30)

    sample = np.random.default_rng(0).choice(int(ok.sum()), size=min(200000, int(ok.sum())), replace=False)
    nd = ndvi.values.ravel()[ok][sample]
    axes[1, 0].hexbin(nd, lai.values.ravel()[ok][sample], gridsize=80, mincnt=1, bins="log", cmap="viridis")
    axes[1, 0].set_xlabel("NDVI"); axes[1, 0].set_ylabel("LAI SNAP")
    axes[1, 0].set_title("NDVI contro LAI (SNAP): l'NDVI 'satura' sopra ~0,85")
    axes[1, 0].set_xlim(0, 1); axes[1, 0].set_ylim(0, 7)

    if biopar is not None:
        bp = biopar.values.ravel()[ok][sample]
        m = np.isfinite(bp)
        axes[1, 1].hexbin(lai.values.ravel()[ok][sample][m], bp[m], gridsize=80, mincnt=1, bins="log", cmap="viridis")
        axes[1, 1].plot([0, 7], [0, 7], "w--", lw=1)
        axes[1, 1].set_xlabel("LAI SNAP (8 bande)"); axes[1, 1].set_ylabel("LAI BIOPAR (3 bande)")
        axes[1, 1].set_title("SNAP contro BIOPAR")
        axes[1, 1].set_xlim(0, 7); axes[1, 1].set_ylim(0, 7)
    else:
        axes[1, 1].axis("off")
    fig.suptitle(f"Parametri biofisici da Sentinel-2 - {cfg['reference_date']}")
    fig.savefig(PROJECT_ROOT / cfg["paths"]["figures"] / "p2_03_biophysical.png", dpi=150)
    plt.close(fig)
    log(f"\nFigura: {cfg['paths']['figures']}/p2_03_biophysical.png")
    (PROJECT_ROOT / "outputs" / "p2_03_summary.txt").write_text("\n".join(LINES), encoding="utf-8")


if __name__ == "__main__":
    main()
