"""
p2_05_station_check.py - Confronto della temperatura dell'aria ERA5-Land con una stazione a terra.

Perche': TSEB dipende dalla differenza temperatura superficie - temperatura aria.
Nello Step 2c abbiamo visto che +1 K di aria vale circa +1 mm/giorno di ET sulle colture fitte.
Prima di correggere Ta serve una MISURA, non un'ipotesi.

Stazione: ARPAE-SIMC GIRALDA (Codigoro, FE), 44.814 N 12.248 E, -1 m s.l.m.,
nella stessa bonifica del Delta, ~7 km a sud dell'area di analisi. Dati scaricati da Dext3r
(temperatura media oraria a 2 m; licenza CC BY 4.0, citare ARPAE-SIMC).

Confronto corretto tra grandezze diverse:
  - la stazione da' la MEDIA dell'ora h -> h+1;
  - ERA5-Land da' valori ISTANTANEI alle ore piene;
  quindi la media oraria si confronta con (ERA5(h) + ERA5(h+1)) / 2.

Output:
  outputs/tables/p2_05_station_vs_era5.csv
  outputs/figures/p2_05_station_vs_era5.png
  outputs/p2_05_summary.txt  (con la correzione suggerita per tseb.ta_offset_K)

Come eseguirlo:
    python scripts/p2_05_station_check.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import yaml

from config_utils import load_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"
LINES = []


def log(text=""):
    print(text)
    LINES.append(text)


def read_dext3r(path):
    """
    Legge un CSV di Dext3r: intestazione con il nome della stazione, poi la tabella oraria,
    poi (dopo righe vuote) i metadati della stazione. Restituisce (serie, metadati).
    """
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("Inizio validit"))
    rows = []
    for l in lines[start + 1:]:
        if not l.strip():
            break                                    # fine della tabella
        t0, t1, val = l.split(",")
        rows.append((pd.Timestamp(t0).tz_convert(None), float(val) if val else np.nan))
    series = pd.Series(dict(rows), name="T_staz_C")   # indice = inizio dell'ora (UTC)

    meta_i = next(i for i, l in enumerate(lines) if l.startswith("Nome della stazione"))
    keys = lines[meta_i].split(",")
    vals = lines[meta_i + 1].split(",")
    meta = dict(zip(keys, vals))
    return series, meta


def main():
    cfg_all = load_config()
    st_cfg = cfg_all["stations"]
    proc = PROJECT_ROOT / cfg_all["paths"]["processed"]
    fig_dir = PROJECT_ROOT / cfg_all["paths"]["figures"]
    tab_dir = PROJECT_ROOT / cfg_all["paths"]["tables"]

    st_file = PROJECT_ROOT / st_cfg["file"]
    if not st_file.exists():
        raise SystemExit(f"Manca il file della stazione {st_cfg['file']}: scaricalo da Dext3r "
                         "(temperatura media oraria, giorno prima/dopo compresi).")
    series, meta = read_dext3r(st_file)
    # Il file puo' coprire un'estate intera: teniamo il giorno prima, quello del passaggio e il dopo
    day = pd.Timestamp(cfg_all["reference_date"])
    series = series[(series.index >= day - pd.Timedelta(days=1)) & (series.index < day + pd.Timedelta(days=2))]
    if series.empty:
        raise SystemExit(f"Il file {st_cfg['file']} non contiene il {cfg_all['reference_date']}.")
    # Un buco di UNA sola ora si riempie interpolando tra l'ora prima e l'ora dopo (la temperatura
    # dell'aria varia con continuita'); buchi piu' lunghi restano vuoti e fermano la correzione.
    filled = series.isna() & series.shift(1).notna() & series.shift(-1).notna()   # buchi isolati
    series = series.where(~filled, (series.shift(1) + series.shift(-1)) / 2)
    if filled.any():
        log("Ore mancanti riempite per interpolazione: " +
            ", ".join(f"{t:%d/%m %H}:00" for t in series.index[filled]))
    lon = float(meta["Longitudine (Gradi Centesimali)"])
    lat = float(meta["Latitudine (Gradi Centesimali)"])
    log(f"Stazione {meta['Nome della stazione']} ({meta['Comune']}), lon {lon}, lat {lat}")

    # ERA5 interpolato sulla stazione; poi media di due ore piene consecutive
    era = xr.open_dataset(proc / "era5land_filled.nc")
    t2m = (era["t2m"] - 273.15).interp(longitude=lon, latitude=lat).to_series()
    era_mid = (t2m + t2m.shift(-1)) / 2                 # media tra h e h+1, indicizzata su h
    df = pd.DataFrame({"T_staz_C": series, "T_era5_C": era_mid.reindex(series.index)})
    df["diff_staz_meno_era5"] = df["T_staz_C"] - df["T_era5_C"]
    df = df.round(2)
    df.to_csv(tab_dir / "p2_05_station_vs_era5.csv", index_label="inizio_ora_utc")

    log("\nOra UTC   stazione  ERA5   differenza")
    for t, r in df.iterrows():
        log(f"{t:%d/%m %H}:00   {r.T_staz_C:6.2f}  {r.T_era5_C:6.2f}  {r.diff_staz_meno_era5:+6.2f}")

    # Valore all'ora del passaggio: media delle due ore a cavallo (09-10 e 10-11 per le 09:58)
    t_pass = pd.Timestamp(f"{cfg_all['reference_date']} {cfg_all['tseb']['overpass_utc']}")
    centre = t_pass.round("h")                       # 09:58 -> 10:00
    win = df.loc[centre - pd.Timedelta(hours=1): centre]   # ore 09-10 e 10-11
    if win["T_staz_C"].isna().any():
        raise SystemExit(f"La stazione non ha misure attorno al passaggio del {cfg_all['reference_date']}: "
                         "nessuna correzione di Ta per questa data (TSEB usera' ERA5 cosi' com'e').")
    bias = win["diff_staz_meno_era5"].mean()
    day = df.between_time("07:00", "15:00")["diff_staz_meno_era5"]
    log(f"\nAl passaggio ({t_pass:%H:%M} UTC): stazione {win.T_staz_C.mean():.2f} C, "
        f"ERA5 {win.T_era5_C.mean():.2f} C -> ERA5 troppo fresca di {bias:.2f} K")
    log(f"Ore 07-15 UTC: differenza media {day.mean():+.2f} K (min {day.min():+.2f}, max {day.max():+.2f})")
    offset = round(bias * 2) / 2
    log(f"Correzione per TSEB: {offset:+.1f} K "
        "(arrotondata a 0,5 K; una sola stazione -> da trattare come ipotesi)")
    # Salvata per questa data: p2_04 la legge se nel config ta_offset_K e' "auto"
    (PROJECT_ROOT / cfg_all["paths"]["summaries"] / "p2_05_ta_offset.yaml").write_text(
        yaml.safe_dump({"ta_offset_K": float(offset), "bias_misurato_K": round(float(bias), 2),
                        "stazione": meta["Nome della stazione"],
                        "ora_passaggio_utc": cfg_all["tseb"]["overpass_utc"]},
                       allow_unicode=True), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(9, 4.5), constrained_layout=True)
    x = df.index + pd.Timedelta(minutes=30)             # centro dell'ora
    ax.plot(x, df["T_staz_C"], "o-", label=f"Stazione {meta['Nome della stazione']} (media oraria)")
    ax.plot(x, df["T_era5_C"], "s--", label="ERA5-Land interpolato sulla stazione")
    ax.axvline(t_pass, color="k", lw=1, ls=":", label="Passaggio Sentinel-3")
    ax.set_ylabel("Temperatura dell'aria a 2 m (C)")
    ax.set_title(f"ERA5-Land contro stazione a terra - {cfg_all['reference_date']}")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.savefig(fig_dir / "p2_05_station_vs_era5.png", dpi=150)
    log(f"Figura: {cfg_all['paths']['figures']}/p2_05_station_vs_era5.png")

    (PROJECT_ROOT / cfg_all["paths"]["summaries"] / "p2_05_summary.txt").write_text("\n".join(LINES) + "\n",
                                                                 encoding="utf-8")


if __name__ == "__main__":
    main()
