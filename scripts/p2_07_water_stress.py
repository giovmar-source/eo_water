"""
p2_07_water_stress.py - Stress idrico della vegetazione a 20 m (Delta del Po, 19/07/2022).

Stress idrico = la pianta traspira MENO di quanto potrebbe con tutta l'acqua che vuole.
Quando l'acqua manca, gli stomi si chiudono, la foglia "suda" meno e si scalda.

Due indici, costruiti in modo indipendente, per non fidarsi di uno solo:

1. CWSI_TSEB (assoluto, dal modello fisico)
      CWSI = 1 - LE_C / LE_C,pot
   LE_C     = traspirazione stimata da TSEB (Step 2c)
   LE_C,pot = traspirazione potenziale di Priestley-Taylor (alpha 1,26, stessa frazione verde)
   0 = chioma al potenziale, 1 = chioma che non traspira.
   Pregi: ha un significato fisico (frazione di traspirazione persa).
   Difetti: eredita gli errori di Ta, rugosita', LST. Inoltre TSEB riduce alpha a passi di 0,1,
   quindi il CWSI cambia a scatti di ~0,08.

2. CWSI contestuale (relativo, dalla sola temperatura)
   Dentro ogni classe di LAI (campi "ugualmente verdi") si prende come riferimento
     T_wet = 2° percentile della LST  (i campi piu' freschi: ben irrigati)
     T_dry = 98° percentile           (i piu' caldi: senz'acqua)
   e  CWSI_ctx = (T - T_wet) / (T_dry - T_wet), limitato tra 0 e 1.
   Pregi: non usa ne' Ta ne' TSEB -> immune all'errore sull'aria.
   Difetti: e' RELATIVO alla scena (se tutti i campi fossero stressati, il piu' fresco
   risulterebbe comunque "0"); confronta campi diversi come se avessero la stessa coltura.

Se i due indici concordano su QUALI campi sono stressati (correlazione di rango), la mappa
e' robusta anche se il livello assoluto e' incerto.

Output:
  data/processed/p1_20220719/cwsi_tseb.tif, cwsi_context.tif   (esclusi da Git)
  outputs/figures/p2_07_water_stress.png
  outputs/tables/p2_07_stress_by_crop.csv, p2_07_stress_by_lai.csv
  outputs/p2_07_summary.txt

Come eseguirlo (dopo p2_04_tseb.py):
    python scripts/p2_07_water_stress.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray  # noqa: F401
import xarray as xr
import yaml

from config_utils import load_config
from rasterio.enums import Resampling
from scipy.stats import spearmanr

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"
LINES = []


def log(text=""):
    print(text)
    LINES.append(text)


def read(path, like=None, resampling=Resampling.nearest):
    """Legge un GeoTIFF a una banda; se 'like' e' dato, lo porta sulla stessa griglia."""
    da = xr.open_dataarray(path, engine="rasterio").squeeze("band", drop=True)
    if like is not None:
        da = da.rio.reproject_match(like, resampling=resampling)
    return da


def contextual_cwsi(T, lai, mask, bins, p_wet, p_dry):
    """
    CWSI relativo: per ogni classe di LAI, posizione della temperatura tra il pixel
    piu' fresco (T_wet) e il piu' caldo (T_dry) della stessa classe.
    Le classi di LAI servono a confrontare campi con la stessa quota di suolo visibile:
    un campo rado e' piu' caldo anche se ben irrigato, perche' il satellite vede piu' suolo.
    """
    out = np.full(T.shape, np.nan)
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        sel = mask & (lai >= lo) & (lai < hi)
        if sel.sum() < 500:
            continue
        t_wet, t_dry = np.percentile(T[sel], [p_wet, p_dry])
        out[sel] = np.clip((T[sel] - t_wet) / (t_dry - t_wet), 0, 1)
        rows.append({"LAI": f"{lo:g}-{hi:g}", "pixel": int(sel.sum()),
                     "T_wet_C": round(t_wet - 273.15, 2), "T_dry_C": round(t_dry - 273.15, 2)})
    return out, pd.DataFrame(rows)


def main():
    cfg_all = load_config()
    cfg, scfg = cfg_all["tseb"], cfg_all["stress"]
    proc = PROJECT_ROOT / cfg_all["paths"]["processed"]
    anc = PROJECT_ROOT / cfg_all["paths"]["aux"]
    fig_dir = PROJECT_ROOT / cfg_all["paths"]["figures"]
    tab_dir = PROJECT_ROOT / cfg_all["paths"]["tables"]

    # --- Input: tutti sulla griglia dei risultati TSEB (area di analisi, 20 m) ---
    ratio = read(proc / "tseb_t_ratio.tif")
    grid = ratio
    lai = read(proc / cfg["lai_file"], grid).values
    Tr = read(proc / cfg["lst_file"], grid).values
    fg = read(proc / "tseb_f_g.tif", grid).values
    crop = read(anc / cfg["cropmap_file"], grid).values
    ratio = ratio.values

    # --- Quali pixel sono "vegetazione" su cui ha senso parlare di stress ---
    veg = (np.isfinite(ratio) & np.isfinite(Tr) & (lai >= scfg["lai_min"])
           & (fg >= scfg["f_g_min"]) & ~np.isin(crop, cfg["exclude_classes"]))
    log(f"Pixel vegetati: {veg.sum():,} ({veg.sum() * 400 / 1e4:,.0f} ha) "
        f"[LAI >= {scfg['lai_min']}, f_g >= {scfg['f_g_min']}, niente acqua/edificato]")

    # --- Indice 1: dal modello TSEB ---
    cwsi_t = np.where(veg, 1 - ratio, np.nan)

    # --- Indice 2: contestuale, dalla sola temperatura ---
    cwsi_c, refs = contextual_cwsi(Tr, lai, veg, scfg["lai_bins"],
                                   scfg["p_wet"], scfg["p_dry"])
    log("\nRiferimenti caldo/freddo per classe di LAI (indice contestuale):")
    log(refs.to_string(index=False))

    # --- Accordo tra i due indici ---
    both = veg & np.isfinite(cwsi_c)
    rho_all = spearmanr(cwsi_t[both], cwsi_c[both]).statistic
    log(f"\nAccordo sul RANGO (Spearman) tra CWSI TSEB e contestuale: {rho_all:.2f}")

    thr = scfg["stress_threshold"]
    names = {int(k): v for k, v in cfg_all["cropmap"]["legend"].items()}
    df = pd.DataFrame({"cwsi_tseb": cwsi_t[both], "cwsi_ctx": cwsi_c[both],
                       "LAI": lai[both], "coltura": pd.Series(crop[both]).map(names).fillna("altro")})
    df["classe_LAI"] = pd.cut(df["LAI"], scfg["lai_bins"], right=False)

    def summary(g):
        return pd.Series({
            "ettari": len(g) * 0.04,                      # pixel 20 x 20 m = 0,04 ha
            "CWSI_TSEB_mediana": g["cwsi_tseb"].median(),
            f"TSEB_stress>{thr}_%": 100 * (g["cwsi_tseb"] > thr).mean(),
            "CWSI_ctx_mediana": g["cwsi_ctx"].median(),
            f"ctx_stress>{thr}_%": 100 * (g["cwsi_ctx"] > thr).mean(),
            "Spearman": spearmanr(g["cwsi_tseb"], g["cwsi_ctx"]).statistic if len(g) > 50 else np.nan,
        })

    by_lai = df.groupby("classe_LAI", observed=True).apply(summary).round(2)
    by_crop = df.groupby("coltura").apply(summary)
    by_crop = by_crop[by_crop["ettari"] > 50].sort_values("ettari", ascending=False).round(2)
    by_lai.to_csv(tab_dir / "p2_07_stress_by_lai.csv")
    by_crop.to_csv(tab_dir / "p2_07_stress_by_crop.csv")
    log("\nPer classe di LAI:")
    log(by_lai.to_string())
    log("\nPer coltura (classi > 50 ha):")
    log(by_crop.to_string())

    # --- Salvataggio mappe ---
    for name, arr in (("cwsi_tseb", cwsi_t), ("cwsi_context", cwsi_c)):
        da = grid.copy(data=arr.astype("float32"))
        da.attrs, da.encoding, da.name = {}, {}, name
        da.rio.write_nodata(np.nan, encoded=False).rio.to_raster(
            proc / f"{name}.tif", compress="deflate")

    # --- Figura ---
    fig, axs = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)
    ext = [grid.x.min() / 1e3, grid.x.max() / 1e3, grid.y.min() / 1e3, grid.y.max() / 1e3]
    for ax, arr, title in ((axs[0, 0], cwsi_t, "CWSI da TSEB (1 - T/T_pot)"),
                           (axs[0, 1], cwsi_c, "CWSI contestuale (stessa classe di LAI)")):
        ax.set_facecolor("0.85")                       # grigio = non vegetato / escluso
        im = ax.imshow(arr, extent=ext, cmap="RdYlBu_r", vmin=0, vmax=1)
        ax.set_title(title)
        ax.set_xlabel("UTM 33N x (km)")
        ax.set_ylabel("y (km)")
        fig.colorbar(im, ax=ax, shrink=0.8, label="0 = nessuno stress, 1 = massimo")
    ax = axs[1, 0]
    sample = df.sample(min(len(df), 40000), random_state=0)
    ax.hexbin(sample["cwsi_ctx"], sample["cwsi_tseb"], gridsize=40, bins="log", cmap="viridis")
    ax.set_xlabel("CWSI contestuale")
    ax.set_ylabel("CWSI TSEB")
    ax.set_title(f"Accordo tra i due indici (Spearman {rho_all:.2f})")
    ax = axs[1, 1]
    top = by_crop.head(6)
    x = np.arange(len(top))
    # Solo TSEB: l'indice contestuale e' RELATIVO alla scena, la sua "percentuale stressata"
    # dipende solo dai percentili scelti e non e' un livello di stress confrontabile.
    ax.bar(x, top[f"TSEB_stress>{thr}_%"], 0.6, color="tab:red")
    for xi, (ha, rho) in enumerate(zip(top["ettari"], top["Spearman"])):
        ax.text(xi, 2, f"{ha:,.0f} ha\nrho {rho:.2f}", ha="center", va="bottom", fontsize=8,
                color="white")
    ax.set_xticks(x, [t[:22] for t in top.index], rotation=30, ha="right")
    ax.set_ylabel(f"% della superficie con CWSI TSEB > {thr}")
    ax.set_title("Superficie stressata per coltura (TSEB; rho = accordo con il contestuale)")
    fig.suptitle(f"Stress idrico - Delta del Po - {cfg_all['reference_date']}", fontsize=14)
    fig.savefig(fig_dir / "p2_07_water_stress.png", dpi=150)
    log(f"Figura: {cfg_all['paths']['figures']}/p2_07_water_stress.png")

    (PROJECT_ROOT / cfg_all["paths"]["summaries"] / "p2_07_summary.txt").write_text("\n".join(LINES) + "\n",
                                                                 encoding="utf-8")


if __name__ == "__main__":
    main()
