"""
p4_02_drought_comparison.py - Quanto la siccita' del 2022 ha cambiato l'acqua consumata nel Delta?

Confronto, pixel per pixel (20 m), tra l'anno di siccita' (2022) e gli anni di riferimento
(2023, 2024, 2025, scelti sulla pioggia misurata a terra, vedi p4_01).

Per ogni anno si usano 5 date nelle STESSE finestre stagionali (inizio giugno, fine giugno,
meta' luglio, fine luglio / inizio agosto, meta'-fine agosto), elaborate con la catena completa
(sharpening + TSEB + stress) e la stessa variante dell'altezza della chioma (--altezza).

Tre grandezze:
  ET   evapotraspirazione giornaliera (mm/giorno): l'acqua consumata.
  EF   frazione evaporativa = LE / (Rn - G): la quota dell'energia disponibile usata per
       evaporare. Dipende poco dal meteo del giorno (piu' o meno Sole), quindi isola meglio
       l'effetto della disponibilita' d'acqua.
  CWSI indice di stress della vegetazione (0 = nessuno stress), solo pixel vegetati.

Come si costruisce il confronto:
  1. per ogni anno, la MEDIA delle date disponibili in ogni pixel (almeno min_dates date valide);
  2. il riferimento e' la MEDIANA delle medie dei tre anni di riferimento (robusta a un anno anomalo);
  3. anomalia = 2022 - riferimento (in mm/giorno e in %).

Analisi spaziale: l'anomalia in funzione della distanza dal mare (classe "nessun dato" della
mappa EUCROPMAP dell'area di addestramento, che sul Delta corrisponde a mare e lagune). Ipotesi:
nel 2022 il cuneo salino ha reso l'acqua irrigua piu' salata vicino alla costa.

Date escluse automaticamente: quelle senza correzione dell'aria dalla stazione (p2_05),
perche' TSEB senza correzione e' fortemente distorto (vedi 14/07/2023).

Output:
  outputs/tables/p4_02_by_date.csv, p4_02_by_year.csv, p4_02_by_distance.csv, p4_02_by_class.csv
  outputs/figures/p4_02_timeseries.png, p4_02_anomaly_maps.png
  outputs/p4_02_summary.txt

Come eseguirlo (dopo run_chain sulle date del config):
    python scripts/p4_02_drought_comparison.py --altezza lai
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray
import yaml
from rasterio.enums import Resampling
from scipy.ndimage import distance_transform_edt

from config_utils import PROJECT_ROOT, load_config

LINES = []


def log(text=""):
    print(text)
    LINES.append(text)


def read(path, like=None):
    """Legge un GeoTIFF a una banda (NaN dove manca il dato); opzionalmente sulla griglia 'like'."""
    with rioxarray.open_rasterio(path, masked=True) as src:
        da = src.squeeze("band", drop=True).load()
    if like is not None and (da.shape != like.shape or not np.allclose(da.x, like.x)):
        da = da.rio.reproject_match(like, resampling=Resampling.nearest)
    return da


def distance_from_sea(cfg, grid):
    """
    Distanza (km) di ogni pixel dal mare e dalle lagune costiere.
    Si usa la mappa EUCROPMAP dell'area di ADDESTRAMENTO (piu' ampia, arriva fino all'Adriatico):
    li' la classe 0 ("nessun dato") e' il mare con le lagune esterne. Nell'area di analisi, piu'
    piccola, il mare quasi non c'e' e le distanze sarebbero sbagliate.
    Calcolo a 100 m (piu' che sufficiente per distanze in km), poi riportato sulla griglia a 20 m.
    """
    # masked=False: nella mappa lo 0 e' dichiarato "nodata", ma per noi e' proprio il mare
    with rioxarray.open_rasterio(PROJECT_ROOT / cfg["paths"]["aux"] / cfg["drought"]["cropmap_wide_file"],
                                 masked=False) as src:
        wide = src.squeeze("band", drop=True).load()
    wide = wide.rio.write_nodata(None)
    wide100 = wide.rio.reproject(grid.rio.crs, resolution=100, resampling=Resampling.mode)
    sea = wide100.values == 0
    dist = wide100.copy(data=distance_transform_edt(~sea) * 100 / 1000.0)
    return dist.rio.reproject_match(grid, resampling=Resampling.bilinear).values


def ta_correction(date):
    """Correzione dell'aria usata da TSEB per quella data (None se la stazione mancava)."""
    f = PROJECT_ROOT / "outputs" / date.replace("-", "") / "p2_05_ta_offset.yaml"
    if not f.exists():
        return None
    return float(yaml.safe_load(f.read_text(encoding="utf-8"))["ta_offset_K"])


def main():
    cfg = load_config(verbose=False)
    dc = cfg["drought"]
    sfx = cfg["tseb"]["run_suffix"]
    tab_dir = PROJECT_ROOT / cfg["paths"]["tables_common"]
    fig_dir = PROJECT_ROOT / cfg["paths"]["figures_common"]
    fig_dir.mkdir(parents=True, exist_ok=True)
    years = [dc["target_year"]] + dc["reference_years"]
    log(f"Confronto siccita': {dc['target_year']} contro {dc['reference_years']} "
        f"(altezza della chioma: {cfg['tseb']['height_source']})")

    # --- Griglia comune e maschere dalla mappa colturale ---
    grid = None
    stacks = {v: {} for v in ("ET", "EF", "CWSI")}
    rows = []
    for year in years:
        for date in dc["dates"][year]:
            proc = PROJECT_ROOT / "data" / "processed" / f"p1_{date.replace('-', '')}"
            et_file = proc / f"tseb_et_day{sfx}.tif"
            off = ta_correction(date)
            if not et_file.exists():
                log(f"  {date}: risultati TSEB mancanti, esclusa")
                continue
            if off is None:
                log(f"  {date}: nessuna correzione dell'aria dalla stazione, esclusa")
                continue
            et = read(et_file, grid)
            grid = et if grid is None else grid
            stacks["ET"][date] = et.values
            stacks["EF"][date] = read(proc / f"tseb_ef{sfx}.tif", grid).values
            stacks["CWSI"][date] = read(proc / f"cwsi_tseb{sfx}.tif", grid).values
            rows.append({"data": date, "anno": year, "correzione_Ta_K": off})

    crop = read(PROJECT_ROOT / cfg["paths"]["aux"] / cfg["tseb"]["cropmap_file"], grid).values
    land = ~np.isin(crop, cfg["tseb"]["exclude_classes"]) & np.isfinite(crop)
    px = abs(float(grid.x[1] - grid.x[0]))
    dist_km = distance_from_sea(cfg, grid)

    # --- 1. Statistiche per data (sulla terraferma) ---
    for r in rows:
        d = r["data"]
        et, ef, cw = stacks["ET"][d], stacks["EF"][d], stacks["CWSI"][d]
        r.update({"frazione_valida": np.isfinite(et[land]).mean(),
                  "ET_media_mm": np.nanmean(et[land]), "ET_mediana_mm": np.nanmedian(et[land]),
                  "EF_media": np.nanmean(ef[land]),
                  "CWSI_medio_vegetazione": np.nanmean(cw),
                  "quota_vegetazione_stressata": np.nanmean(cw[np.isfinite(cw)] > 0.3)})
    by_date = pd.DataFrame(rows).round(3)
    by_date.to_csv(tab_dir / f"p4_02_by_date{sfx}.csv", index=False)
    log("\nPer data (media sulla terraferma; stress solo sulla vegetazione):")
    log(by_date.to_string(index=False))

    # --- 2. Media per anno in ogni pixel, riferimento = mediana degli anni di riferimento ---
    def yearly_mean(var, year):
        arrs = [stacks[var][d] for d in dc["dates"][year] if d in stacks[var]]
        if not arrs:
            return np.full(grid.shape, np.nan), 0
        cube = np.stack(arrs)
        n_ok = np.isfinite(cube).sum(axis=0)
        with np.errstate(invalid="ignore"):
            mean = np.where(n_ok >= min(dc["min_dates"], len(arrs)), np.nanmean(cube, axis=0), np.nan)
        return mean, len(arrs)

    maps, year_rows = {}, []
    for var in ("ET", "EF", "CWSI"):
        per_year = {}
        for y in years:
            per_year[y], n = yearly_mean(var, y)
            mask = land if var != "CWSI" else np.isfinite(per_year[y])
            year_rows.append({"grandezza": var, "anno": y, "date_usate": n,
                              "media_area": np.nanmean(per_year[y][mask])})
        ref = np.nanmedian(np.stack([per_year[y] for y in dc["reference_years"]]), axis=0)
        maps[var] = {"target": per_year[dc["target_year"]], "ref": ref,
                     "anom": per_year[dc["target_year"]] - ref}
    by_year = pd.DataFrame(year_rows).round(3)
    by_year.to_csv(tab_dir / f"p4_02_by_year{sfx}.csv", index=False)
    log("\nMedia stagionale per anno (media delle date di ogni pixel, poi media dell'area):")
    log(by_year.pivot(index="anno", columns="grandezza", values="media_area").round(3).to_string())

    et_t, et_r = maps["ET"]["target"], maps["ET"]["ref"]
    ok = land & np.isfinite(et_t) & np.isfinite(et_r)
    anom_mm = np.nanmean(et_t[ok] - et_r[ok])
    anom_pct = 100 * (np.nanmean(et_t[ok]) / np.nanmean(et_r[ok]) - 1)
    ha = ok.sum() * px * px / 1e4
    # Volume d'acqua: mm/giorno x ettari x 10 m3/(mm ha) x giorni della stagione giugno-agosto
    vol = anom_mm * ha * 10 * dc["season_days"] / 1e6
    log(f"\nANOMALIA {dc['target_year']} (terraferma, {ha:,.0f} ha): ET {anom_mm:+.2f} mm/giorno "
        f"({anom_pct:+.0f}%), EF {np.nanmean(maps['EF']['anom'][ok]):+.3f}")
    log(f"Ordine di grandezza sull'estate ({dc['season_days']} giorni): {vol:+.1f} milioni di m3 "
        "(stima grezza: 5 fotografie per estate, non una serie continua)")

    # --- 3. Anomalia per distanza dal mare ---
    bins = dc["distance_bins_km"]
    df = pd.DataFrame({"dist_km": dist_km[ok], "anom_ET": (et_t - et_r)[ok],
                       "ET_ref": et_r[ok], "anom_EF": maps["EF"]["anom"][ok]})
    df["fascia_km"] = pd.cut(df["dist_km"], bins, include_lowest=True)
    by_dist = df.groupby("fascia_km", observed=True).agg(
        ettari=("anom_ET", lambda s: s.size * px * px / 1e4), ET_rif_mm=("ET_ref", "mean"),
        anomalia_ET_mm=("anom_ET", "mean"), anomalia_EF=("anom_EF", "mean"))
    by_dist["anomalia_ET_%"] = 100 * by_dist["anomalia_ET_mm"] / by_dist["ET_rif_mm"]
    by_dist = by_dist.round(3)
    by_dist.to_csv(tab_dir / f"p4_02_by_distance{sfx}.csv")
    log("\nAnomalia per distanza dal mare:")
    log(by_dist.to_string())

    # --- 3b. Il gradiente costiero e' reale o un artefatto? Due controlli ---
    # (1) Coerenza tra le date: per ogni data, rapporto tra l'ET dell'entroterra (4-8 km) e quella
    #     della costa (0-2 km). Se il gradiente del 2022 fosse un caso, comparirebbe in una data
    #     sola; se e' un effetto della siccita', il rapporto 2022 e' piu' basso in quasi tutte le date.
    coast = land & (dist_km <= 2)
    inland = land & (dist_km > 4) & (dist_km <= 8)
    grad_rows = []
    for r in rows:
        et = stacks["ET"][r["data"]]
        grad_rows.append({"data": r["data"], "anno": r["anno"],
                          "ET_costa_0_2km": np.nanmean(et[coast]),
                          "ET_entroterra_4_8km": np.nanmean(et[inland])})
    grad = pd.DataFrame(grad_rows)
    grad["rapporto_entroterra_costa"] = grad["ET_entroterra_4_8km"] / grad["ET_costa_0_2km"]
    grad = grad.round(3)
    grad.to_csv(tab_dir / f"p4_02_coast_gradient{sfx}.csv", index=False)
    log("\nGradiente costiero data per data (rapporto ET entroterra / costa):")
    log(grad.groupby("anno")["rapporto_entroterra_costa"].agg(["median", "min", "max", "count"])
        .round(3).to_string())

    # (2) Un indice che NON usa la temperatura dell'aria: il CWSI contestuale (Step 2d).
    #     Se anche lui mostra piu' stress nell'entroterra nel 2022 (rispetto agli anni di
    #     riferimento), il gradiente non puo' venire dalla correzione dell'aria.
    ctx = {}
    for r in rows:
        f = (PROJECT_ROOT / "data" / "processed" / f"p1_{r['data'].replace('-', '')}"
             / f"cwsi_context{sfx}.tif")
        if f.exists():
            ctx[r["data"]] = read(f, grid).values
    if ctx:
        stacks["CTX"] = ctx
        ctx_year = {y: yearly_mean("CTX", y)[0] for y in years}
        ctx_anom = ctx_year[dc["target_year"]] - np.nanmedian(
            np.stack([ctx_year[y] for y in dc["reference_years"]]), axis=0)
        sel = np.isfinite(ctx_anom)
        d2 = pd.DataFrame({"fascia_km": pd.cut(dist_km[sel], bins, include_lowest=True),
                           "anomalia_CWSI_contestuale": ctx_anom[sel]})
        ctx_tab = d2.groupby("fascia_km", observed=True)["anomalia_CWSI_contestuale"].agg(
            ["mean", "count"]).round(3)
        ctx_tab.to_csv(tab_dir / f"p4_02_coast_gradient_contextual{sfx}.csv")
        log("\nAnomalia del CWSI contestuale (non usa la temperatura dell'aria) per distanza dal mare:")
        log(ctx_tab.to_string())

    # --- 4. Per classe di copertura "stabile" (la mappa e' del 2022: i seminativi ruotano) ---
    groups = {"seminativi (tutte le colture annuali)": (crop >= 200) & (crop < 300),
              "foraggere (erba medica, pluriennale)": crop == 250,
              "prati": crop == 500, "bosco, frutteti, pioppeti": crop == 300}
    cls_rows = []
    for name, g in groups.items():
        sel = ok & g
        if sel.sum() < 100:
            continue
        cls_rows.append({"classe": name, "ettari": sel.sum() * px * px / 1e4,
                         "ET_2022": np.nanmean(et_t[sel]), "ET_rif": np.nanmean(et_r[sel]),
                         "anomalia_%": 100 * (np.nanmean(et_t[sel]) / np.nanmean(et_r[sel]) - 1)})
    by_class = pd.DataFrame(cls_rows).round(2)
    by_class.to_csv(tab_dir / f"p4_02_by_class{sfx}.csv", index=False)
    log("\nPer classe di copertura:")
    log(by_class.to_string(index=False))

    # --- Figure ---
    fig, axs = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)
    colors = {y: ("tab:red" if y == dc["target_year"] else c)
              for y, c in zip(years, ["tab:red", "tab:blue", "tab:green", "tab:purple"])}
    for ax, (col, label) in zip(axs, [("ET_media_mm", "ET media (mm/giorno)"),
                                      ("EF_media", "frazione evaporativa media"),
                                      ("CWSI_medio_vegetazione", "CWSI medio (vegetazione)")]):
        for y in years:
            sub = by_date[by_date["anno"] == y]
            doy = pd.to_datetime(sub["data"]).dt.dayofyear
            ax.plot(doy, sub[col], "o-", color=colors[y], label=str(y), lw=2.5 if y == dc["target_year"] else 1.2)
        ax.set_xticks([152, 166, 182, 196, 213, 227, 244], ["1 giu", "15 giu", "1 lug", "15 lug",
                                                             "1 ago", "15 ago", "1 set"])
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
    axs[0].legend()
    fig.suptitle("Delta del Po: le date di ogni estate (rosso = 2022)")
    fig.savefig(fig_dir / f"p4_02_timeseries{sfx}.png", dpi=150)

    fig, axs = plt.subplots(2, 3, figsize=(18, 11), constrained_layout=True)
    ext = [grid.x.min() / 1e3, grid.x.max() / 1e3, grid.y.min() / 1e3, grid.y.max() / 1e3]
    panels = [(maps["ET"]["target"], f"ET media {dc['target_year']} (mm/g)", "YlGnBu", (0, 6)),
              (maps["ET"]["ref"], "ET media, riferimento 2023-2025 (mm/g)", "YlGnBu", (0, 6)),
              (maps["ET"]["anom"], "anomalia ET (mm/g)", "RdBu", (-3, 3)),
              (maps["CWSI"]["target"], f"CWSI medio {dc['target_year']}", "RdYlBu_r", (0, 1)),
              (maps["CWSI"]["ref"], "CWSI medio, riferimento", "RdYlBu_r", (0, 1)),
              (maps["CWSI"]["anom"], "anomalia CWSI", "RdBu_r", (-0.6, 0.6))]
    for ax, (arr, title, cmap, lim) in zip(axs.flat, panels):
        a = np.where(land, arr, np.nan)
        ax.set_facecolor("0.85")
        im = ax.imshow(a, extent=ext, cmap=cmap, vmin=lim[0], vmax=lim[1])
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle(f"Siccita' {dc['target_year']} contro anni di riferimento (grigio = acqua, edificato, mare)")
    fig.savefig(fig_dir / f"p4_02_anomaly_maps{sfx}.png", dpi=130)
    log(f"Figure: outputs/figures/p4_02_timeseries{sfx}.png, p4_02_anomaly_maps{sfx}.png")

    (PROJECT_ROOT / "outputs" / f"p4_02_summary{sfx}.txt").write_text("\n".join(LINES) + "\n",
                                                                       encoding="utf-8")


if __name__ == "__main__":
    main()
