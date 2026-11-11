"""
p3_02_compare_ecostress.py - Confronto della nostra catena con ECOSTRESS (LST, ET, stress).

Cosa confrontiamo (stessa data, stessa area):
  LST     nostra LST affinata (S3 + S2)       vs  ECOSTRESS LST
  ET      ET giornaliera di TSEB              vs  ECOSTRESS ETdaily (modello PT-JPL)
  stress  CWSI di TSEB (0 = nessuno stress)   vs  1 - ESI di ECOSTRESS (stesso verso)

ATTENZIONE alla scala. ECOSTRESS e' nominalmente a 70 m, ma la nitidezza reale delle scene
estive 2022 sul Delta e' molto minore (campi sfumati, strisce): verificato confrontando le
immagini e la correlazione con l'NDMI di Sentinel-2. Per questo il confronto si fa su una
griglia piu' grossolana (compare_resolution_m, default 500 m), dove entrambi sono affidabili.
Si riporta anche il confronto a 100 m, per trasparenza.

Le metriche:
  r (Pearson) e rho (Spearman, sui ranghi), bias (nostro - ECOSTRESS), RMSE.
  Due modelli di ET indipendenti non devono coincidere: conta che concordino su DOVE si
  evapora di piu' o di meno (rho) e che lo scarto medio sia ragionevole.

Output:
  outputs/<data>/tables/p3_02_ecostress_comparison.csv
  outputs/<data>/figures/p3_02_ecostress_comparison.png
  outputs/<data>/p3_02_summary.txt

Come eseguirlo (dopo run_chain per quella data e p3_01):
    python scripts/p3_02_compare_ecostress.py --date 2022-08-11
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray
import yaml
from rasterio.enums import Resampling
from scipy.stats import spearmanr

from config_utils import load_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]
LINES = []


def log(text=""):
    print(text)
    LINES.append(text)


def read(path, masked=True):
    return rioxarray.open_rasterio(path, masked=masked).squeeze("band", drop=True)


def find_ecostress(cfg):
    """Cartella ECOSTRESS dello stesso giorno, la piu' vicina all'ora di Sentinel-3."""
    root = PROJECT_ROOT / "data" / "raw" / "ecostress"
    day = cfg["reference_date"].replace("-", "")
    folders = sorted(root.glob(f"{day}T*"))
    if not folders:
        raise SystemExit(f"Nessun dato ECOSTRESS per il {cfg['reference_date']}: lancia prima p3_01.")
    t_s3 = pd.Timestamp(f"{cfg['reference_date']} {cfg['training']['s3_time_utc']}")
    folder = min(folders, key=lambda f: abs(pd.Timestamp(f.name) - t_s3))
    tile = sorted(folder.glob("*_LST.tif"))[0].name.split("_")[0]
    return folder, tile


def metrics(ours, ref):
    m = np.isfinite(ours) & np.isfinite(ref)
    a, b = ours[m], ref[m]
    return {"n": int(m.sum()), "r": np.corrcoef(a, b)[0, 1], "rho": spearmanr(a, b).statistic,
            "bias": float(np.mean(a - b)), "RMSE": float(np.sqrt(np.mean((a - b) ** 2))),
            "media_nostra": float(np.mean(a)), "media_ECOSTRESS": float(np.mean(b))}


def main():
    cfg = load_config()
    proc = PROJECT_ROOT / cfg["paths"]["processed"]
    res_m = cfg["ecostress"].get("compare_resolution_m", 500)
    folder, tile = find_ecostress(cfg)
    log(f"ECOSTRESS {folder.name} UTC (tassello {tile}) contro la nostra catena "
        f"(Sentinel-3 {cfg['training']['s3_time_utc']} UTC)")

    def eco(name):
        return read(folder / f"{tile}_{name}.tif")

    bad = (eco("cloud") == 1) | (eco("water") == 1)       # nuvole e acqua di ECOSTRESS
    eco_vars = {"LST": eco("LST").where(~bad),
                "ET": eco("ETdaily").where(~bad),
                "stress": (1 - eco("ESI")).where(~bad)}      # 1 - ESI: 0 = nessuno stress

    lst = read(proc / cfg["tseb"]["lst_file"])
    ours = {"LST": lst,
            "ET": read(proc / f"tseb_et_day{cfg['tseb']['run_suffix']}.tif"),
            "stress": read(proc / f"cwsi_tseb{cfg['tseb']['run_suffix']}.tif")}           # solo pixel vegetati (Step 2d)
    units = {"LST": "K", "ET": "mm/giorno", "stress": "0-1"}

    rows, pairs = [], {}
    for res in (100, res_m):
        # Griglia comune in UTM 33N alla risoluzione del confronto, dentro l'area di analisi
        grid = ours["ET"].rio.reproject(ours["ET"].rio.crs, resolution=res,
                                        resampling=Resampling.average, nodata=np.nan)
        for var in ("LST", "ET", "stress"):
            o = ours[var].rio.reproject_match(grid, resampling=Resampling.average).values.ravel()
            e = eco_vars[var].rio.reproject_match(grid, resampling=Resampling.average).values.ravel()
            row = {"variabile": var, "risoluzione_m": res, **metrics(o, e)}
            rows.append(row)
            if res == res_m:
                pairs[var] = (o, e)

    df = pd.DataFrame(rows).round(3)
    df.to_csv(PROJECT_ROOT / cfg["paths"]["tables"] / "p3_02_ecostress_comparison.csv", index=False)
    log("\nConfronto (bias = nostro - ECOSTRESS):")
    log(df.to_string(index=False))

    fig, axs = plt.subplots(1, 3, figsize=(17, 5.5), constrained_layout=True)
    for ax, var in zip(axs, ("LST", "ET", "stress")):
        o, e = pairs[var]
        m = np.isfinite(o) & np.isfinite(e)
        x, y = e[m], o[m]
        if var == "LST":
            x, y = x - 273.15, y - 273.15
        ax.hexbin(x, y, gridsize=35, bins="log", cmap="viridis", mincnt=1)
        lims = [np.nanpercentile(np.r_[x, y], 1), np.nanpercentile(np.r_[x, y], 99)]
        ax.plot(lims, lims, "w--", lw=1)
        r = df[(df.variabile == var) & (df.risoluzione_m == res_m)].iloc[0]
        unit = "°C" if var == "LST" else units[var]
        ax.set_xlabel(f"ECOSTRESS ({unit})")
        ax.set_ylabel(f"Nostra catena ({unit})")
        ax.set_title(f"{var} a {res_m} m: r {r.r:.2f}, rho {r.rho:.2f}, bias {r.bias:+.2f}")
    fig.suptitle(f"Confronto con ECOSTRESS - {cfg['reference_date']} "
                 f"(stress: CWSI TSEB contro 1 - ESI, solo vegetazione)")
    fig.savefig(PROJECT_ROOT / cfg["paths"]["figures"] / "p3_02_ecostress_comparison.png", dpi=150)
    log(f"Figura: {cfg['paths']['figures']}/p3_02_ecostress_comparison.png")

    (PROJECT_ROOT / cfg["paths"]["summaries"] / "p3_02_summary.txt").write_text(
        "\n".join(LINES) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
