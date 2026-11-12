"""
p4_03_wapor.py - Confronto con WaPOR (FAO): anche un prodotto indipendente vede il calo del 2022?

WaPOR versione 3, Livello 1: evapotraspirazione reale + intercettazione (AETI), a 300 m,
ogni decade (10 giorni), su tutto il globo dal 2018. Modello diverso dal nostro (ETLook,
FAO/IHE Delft), dati diversi: e' un confronto INDIPENDENTE dell'anomalia.

Cosa fa:
  1. Scarica, per ogni estate 2018-2025, le 9 decadi da giugno ad agosto. Di ogni file globale
     (122880 x 61440 pixel) legge SOLO il ritaglio sull'area: poche decine di KB.
  2. Calcola l'ET estiva (somma di mm/giorno x giorni della decade) di ogni anno e l'anomalia
     2022 rispetto agli anni di riferimento (gli stessi di p4_02) e a tutti gli altri anni.
  3. Confronta, pixel per pixel a 300 m, la NOSTRA anomalia (5 date a 20 m, da p4_02) con quella
     di WaPOR, e controlla se anche WaPOR vede il gradiente costa-entroterra.

Perche' un processo separato per ogni file: su Windows la lettura via internet con GDAL puo'
lasciare il processo "appeso" a fine lettura (vedi p1_02). Ogni ritaglio lo legge un piccolo
processo "operaio"; il processo principale aspetta che il file sia scritto e poi lo chiude d'ufficio.

Output:
  data/raw/wapor/AETI_<decade>.tif                       (ritagli, esclusi da Git)
  outputs/tables/p4_03_wapor_by_year.csv, p4_03_wapor_vs_ours.csv
  outputs/figures/p4_03_wapor.png
  outputs/p4_03_summary.txt

Come eseguirlo (dopo p4_02):
    python scripts/p4_03_wapor.py --altezza lai
"""

import calendar
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

BASE_URL = ("https://storage.googleapis.com/fao-gismgr-wapor-3-data/DATA/WAPOR-3/MAPSET/"
            "L1-AETI-D/WAPOR-3.L1-AETI-D.{code}.tif")
SCALE = 0.1          # i valori sono interi: mm/giorno = valore x 0,1 (dal catalogo WaPOR)
NODATA = -9999


def worker(url, out, bbox):
    """Processo operaio: legge il ritaglio di un file WaPOR e lo salva. Poi segnala 'fatto'."""
    os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")   # non elencare la cartella remota
    import rioxarray
    with rioxarray.open_rasterio(url, masked=True) as src:
        clip = src.rio.clip_box(*bbox, crs="EPSG:4326").squeeze("band", drop=True).load()
    clip.rio.to_raster(out)
    Path(str(out) + ".done").write_text("ok")


def dekads(year, months):
    """Codici delle decadi WaPOR (es. 2022-07-D1) e i loro giorni (la terza decade va a fine mese)."""
    out = []
    for m in range(months[0], months[1] + 1):
        last = calendar.monthrange(year, m)[1]
        for d, days in ((1, 10), (2, 10), (3, last - 20)):
            out.append((f"{year}-{m:02d}-D{d}", days))
    return out


def fetch(code, out_dir, bbox, timeout_s=180):
    """Scarica il ritaglio di una decade con un processo operaio. Restituisce il percorso o None."""
    out = out_dir / f"AETI_{code}.tif"
    if out.exists():
        return out
    tmp = out_dir / f"_tmp_{code}.tif"
    done = Path(str(tmp) + ".done")
    for f in (tmp, done):
        if f.exists():
            f.unlink()
    cmd = [sys.executable, __file__, "--worker", BASE_URL.format(code=code), str(tmp),
           *[str(b) for b in bbox]]
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    t0 = time.time()
    while time.time() - t0 < timeout_s and not done.exists() and proc.poll() is None:
        time.sleep(1)
    if proc.poll() is None:
        proc.kill()                       # chiude d'ufficio l'operaio (anche se e' rimasto appeso)
    if done.exists():
        done.unlink()
        tmp.replace(out)
        return out
    err = proc.stderr.read().decode(errors="ignore")[-300:] if proc.poll() is not None else "tempo scaduto"
    print(f"    {code}: non riuscito ({err.strip().splitlines()[-1] if err.strip() else 'errore'})")
    return None


def main():
    import matplotlib.pyplot as plt
    import pandas as pd
    import rioxarray
    from rasterio.enums import Resampling
    from scipy.stats import spearmanr

    from config_utils import PROJECT_ROOT, load_config

    cfg = load_config(verbose=False)
    wc, dc = cfg["wapor"], cfg["drought"]
    sfx = cfg["tseb"]["run_suffix"]
    a = cfg["analysis_aoi"]
    m = wc["margin_deg"]
    bbox = (a["west"] - m, a["south"] - m, a["east"] + m, a["north"] + m)
    out_dir = PROJECT_ROOT / "data" / "raw" / "wapor"
    out_dir.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(t=""):
        print(t)
        lines.append(t)

    # --- 1. Download dei ritagli ---
    log(f"WaPOR L1 AETI, decadi {wc['months'][0]}-{wc['months'][1]}, anni {wc['years'][0]}-{wc['years'][-1]}")
    season = {}
    for year in wc["years"]:
        arrs, days_tot, grid = [], 0, None
        for code, days in dekads(year, wc["months"]):
            f = fetch(code, out_dir, bbox)
            if f is None:
                continue
            with rioxarray.open_rasterio(f, masked=True) as src:
                da = src.squeeze("band", drop=True).load()
            da = da.where(da != NODATA) * SCALE                       # mm/giorno
            grid = da if grid is None else grid
            arrs.append(da.rio.reproject_match(grid).values * days)   # mm nella decade
            days_tot += days
        if arrs:
            total = np.nansum(np.stack(arrs), axis=0)
            total[np.all(~np.isfinite(np.stack(arrs)), axis=0)] = np.nan
            season[year] = (grid.copy(data=total), len(arrs), days_tot)
            log(f"  {year}: {len(arrs)} decadi")
    if dc["target_year"] not in season:
        raise SystemExit("Mancano le decadi del 2022: impossibile il confronto.")

    # --- 2. ET estiva per anno sull'area di analisi e anomalia 2022 ---
    ref_grid = season[dc["target_year"]][0]
    inner = ref_grid.rio.clip_box(a["west"], a["south"], a["east"], a["north"], crs="EPSG:4326")
    rows = []
    for y, (da, n, d) in season.items():
        v = da.rio.reproject_match(inner).values
        rows.append({"anno": y, "decadi": n, "giorni": d, "ET_estiva_mm": np.nanmean(v),
                     "ET_media_mm_giorno": np.nanmean(v) / d})
    by_year = pd.DataFrame(rows).round(2)
    by_year.to_csv(PROJECT_ROOT / cfg["paths"]["tables_common"] / "p4_03_wapor_by_year.csv", index=False)
    log("\nET estiva WaPOR (giugno-agosto), media dell'area di analisi:")
    log(by_year.to_string(index=False))
    t = by_year.set_index("anno")["ET_media_mm_giorno"]
    ref_med = t.loc[[y for y in dc["reference_years"] if y in t.index]].median()
    others = t.drop(index=dc["target_year"]).median()
    log(f"\nAnomalia 2022 WaPOR: contro {dc['reference_years']} {100 * (t[dc['target_year']] / ref_med - 1):+.0f}%, "
        f"contro tutti gli altri anni {100 * (t[dc['target_year']] / others - 1):+.0f}%")

    # --- 3. Confronto pixel per pixel con la nostra anomalia (aggregata a 300 m) ---
    drought = PROJECT_ROOT / "data" / "processed" / "drought"
    ours = {}
    for name in ("et_2022", "et_rif", "distanza_mare_km"):
        with rioxarray.open_rasterio(drought / f"{name}{sfx}.tif", masked=True) as src:
            da = src.squeeze("band", drop=True).load()
        ours[name] = da.rio.reproject_match(inner, resampling=Resampling.average).values
    w22 = season[dc["target_year"]][0].rio.reproject_match(inner).values / season[dc["target_year"]][2]
    w_ref = np.nanmedian(np.stack([season[y][0].rio.reproject_match(inner).values / season[y][2]
                                   for y in dc["reference_years"] if y in season]), axis=0)
    our_pct = 100 * (ours["et_2022"] / ours["et_rif"] - 1)
    wap_pct = 100 * (w22 / w_ref - 1)
    ok = np.isfinite(our_pct) & np.isfinite(wap_pct)
    rho = spearmanr(our_pct[ok], wap_pct[ok]).statistic
    r = np.corrcoef(our_pct[ok], wap_pct[ok])[0, 1]
    log(f"\nAnomalia 2022 pixel per pixel a 300 m ({ok.sum()} pixel): "
        f"nostra mediana {np.median(our_pct[ok]):+.0f}%, WaPOR {np.median(wap_pct[ok]):+.0f}%, "
        f"r {r:.2f}, rho {rho:.2f}")

    bins = dc["distance_bins_km"]
    df = pd.DataFrame({"fascia_km": pd.cut(ours["distanza_mare_km"][ok], bins, include_lowest=True),
                       "nostra_%": our_pct[ok], "WaPOR_%": wap_pct[ok]})
    by_dist = df.groupby("fascia_km", observed=True).agg(pixel=("WaPOR_%", "size"),
                                                         nostra_pct=("nostra_%", "median"),
                                                         wapor_pct=("WaPOR_%", "median")).round(1)
    by_dist.to_csv(PROJECT_ROOT / cfg["paths"]["tables_common"] / f"p4_03_wapor_vs_ours{sfx}.csv")
    log("\nAnomalia 2022 per distanza dal mare (mediane a 300 m):")
    log(by_dist.to_string())

    # --- Figura ---
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.2), constrained_layout=True)
    axs[0].bar(by_year["anno"], by_year["ET_media_mm_giorno"],
               color=["tab:red" if y == dc["target_year"] else "0.6" for y in by_year["anno"]])
    axs[0].set_ylabel("ET media giugno-agosto (mm/giorno)")
    axs[0].set_title("WaPOR L1: ET estiva del Delta per anno")
    axs[1].hexbin(wap_pct[ok], our_pct[ok], gridsize=30, bins="log", cmap="viridis", mincnt=1)
    axs[1].plot([-80, 20], [-80, 20], "w--", lw=1)
    axs[1].set_xlabel("anomalia 2022 WaPOR (%)")
    axs[1].set_ylabel("anomalia 2022 nostra (%)")
    axs[1].set_title(f"Pixel da 300 m: r {r:.2f}, rho {rho:.2f}")
    x = np.arange(len(by_dist))
    axs[2].bar(x - 0.2, by_dist["nostra_pct"], 0.4, label="nostra (5 date, 20 m)")
    axs[2].bar(x + 0.2, by_dist["wapor_pct"], 0.4, label="WaPOR (9 decadi, 300 m)")
    axs[2].set_xticks(x, [str(i) for i in by_dist.index], rotation=30)
    axs[2].set_xlabel("distanza dal mare (km)")
    axs[2].set_ylabel("anomalia 2022 (%)")
    axs[2].set_title("Gradiente costa-entroterra")
    axs[2].legend()
    fig.savefig(PROJECT_ROOT / cfg["paths"]["figures_common"] / f"p4_03_wapor{sfx}.png", dpi=150)
    log(f"Figura: outputs/figures/p4_03_wapor{sfx}.png")
    (PROJECT_ROOT / "outputs" / f"p4_03_summary{sfx}.txt").write_text("\n".join(lines) + "\n",
                                                                       encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        worker(sys.argv[2], sys.argv[3], tuple(float(v) for v in sys.argv[4:8]))
    else:
        main()
