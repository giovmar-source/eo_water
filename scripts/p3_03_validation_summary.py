"""
p3_03_validation_summary.py - Riassume la validazione dello sharpening su tutte le date elaborate.

Legge, per ogni cartella outputs/<data>/, la tabella p1_05_validation.csv e il riassunto
p1_05_summary.txt, e costruisce:
  - una tabella con, per ogni data: riferimento usato (Landsat o ECOSTRESS), ora di S3,
    r / RMSE senza bias di Sentinel-3 interpolato e dello sharpening (foresta casuale con
    correzione dei residui), e la variabilita' sotto il km ricostruita;
  - una figura "a manubrio": per ogni data, r prima (S3 interpolato) e dopo (sharpening).

Perche' conta: una validazione su un giorno solo puo' essere fortuna. Su molte date, stagioni
e anni diversi si vede se il metodo funziona IN GENERALE, e quanto varia la sua qualita'.

Output:
  outputs/tables/p3_03_validation_all_dates.csv
  outputs/figures/p3_03_validation_all_dates.png

Come eseguirlo (dopo run_chain su piu' date):
    python scripts/p3_03_validation_summary.py
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config_utils import load_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BASE = "Sentinel-3 interpolato (nessuno sharpening)"
SHARP = "Foresta casuale (con correzione)"


def read_summary(folder):
    """Estrae dal riassunto di p1_05 il riferimento usato e la variabilita' sotto il km."""
    text = (folder / "p1_05_summary.txt").read_text(encoding="utf-8")
    ref = "Landsat" if text.startswith("Landsat") else ("ECOSTRESS" if "ECOSTRESS:" in text else "?")
    ratio = re.search(r"Foresta casuale\s+[\d.]+ K\s+\(rapporto ricostruita/vera (\d+)%\)", text)
    s3 = None
    p13 = folder / "p1_03_summary.txt"
    if p13.exists():
        m = re.search(r"uso il passaggio delle (\d\d:\d\d)", p13.read_text(encoding="utf-8"))
        s3 = m.group(1) if m else None
    return ref, (int(ratio.group(1)) if ratio else np.nan), s3


def main():
    cfg = load_config(verbose=False)
    rows = []
    for folder in sorted((PROJECT_ROOT / "outputs").glob("20*")):
        table = folder / "tables" / "p1_05_validation.csv"
        if not table.exists():
            continue
        v = pd.read_csv(table).set_index("metodo")
        if BASE not in v.index or SHARP not in v.index:
            continue
        ref, ratio, s3 = read_summary(folder)
        d = pd.Timestamp(folder.name)
        rows.append({
            "data": f"{d:%Y-%m-%d}", "mese": d.month, "riferimento": ref, "ora_S3": s3,
            "celle_100m": int(v.loc[SHARP, "n"]),
            "r_S3": v.loc[BASE, "r"], "r_sharpening": v.loc[SHARP, "r"],
            "RMSEub_S3_K": v.loc[BASE, "RMSE_senza_bias_K"],
            "RMSEub_sharpening_K": v.loc[SHARP, "RMSE_senza_bias_K"],
            "bias_sharpening_K": v.loc[SHARP, "bias_K"],
            "variabilita_ricostruita_%": ratio,
        })
    if not rows:
        raise SystemExit("Nessuna validazione trovata in outputs/<data>/tables/p1_05_validation.csv")

    df = pd.DataFrame(rows)
    df["guadagno_r"] = df["r_sharpening"] - df["r_S3"]
    out_tab = PROJECT_ROOT / cfg["paths"]["tables_common"] / "p3_03_validation_all_dates.csv"
    df.round(3).to_csv(out_tab, index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(df.round(2).to_string(index=False))

    for ref in sorted(df["riferimento"].unique()):
        sub = df[df["riferimento"] == ref]
        print(f"\n{ref} ({len(sub)} date): r mediano S3 {sub['r_S3'].median():.2f} -> sharpening "
              f"{sub['r_sharpening'].median():.2f}; RMSE senza bias mediano "
              f"{sub['RMSEub_S3_K'].median():.2f} -> {sub['RMSEub_sharpening_K'].median():.2f} K; "
              f"sharpening migliore in {int((sub['guadagno_r'] > 0).sum())} date su {len(sub)}")

    # Figura a manubrio: un segmento per data, dal pallino grigio (S3) al pallino colorato
    fig, ax = plt.subplots(figsize=(9, 0.45 * len(df) + 2), constrained_layout=True)
    y = np.arange(len(df))[::-1]
    colors = df["riferimento"].map({"Landsat": "tab:blue", "ECOSTRESS": "tab:orange"}).fillna("k")
    ax.hlines(y, df["r_S3"], df["r_sharpening"], color="0.7", lw=2)
    ax.scatter(df["r_S3"], y, color="0.5", zorder=3, label="Sentinel-3 interpolato")
    ax.scatter(df["r_sharpening"], y, color=colors, zorder=3, label="sharpening (colore = riferimento)")
    ax.set_yticks(y, [f"{d}  ({r})" for d, r in zip(df["data"], df["riferimento"])])
    ax.set_xlabel("correlazione con il riferimento a 100 m (r)")
    ax.set_title("Validazione dello sharpening su tutte le date")
    ax.grid(axis="x", alpha=0.3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=2, frameon=False)
    out_fig = PROJECT_ROOT / cfg["paths"]["figures_common"] / "p3_03_validation_all_dates.png"
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, dpi=150)
    print(f"\nTabella: {out_tab.relative_to(PROJECT_ROOT)}\nFigura: {out_fig.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
