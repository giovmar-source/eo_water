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
    cfg_all = yaml.safe_load(CONFIG_FILE.read_text(encoding="utf-8"))
    st_cfg = cfg_all["stations"]
    proc = PROJECT_ROOT / cfg_all["paths"]["processed"]
    fig_dir = PROJECT_ROOT / cfg_all["paths"]["figures"]
    tab_dir = PROJECT_ROOT / cfg_all["paths"]["tables"]

    series, meta = read_dext3r(PROJECT_ROOT / st_cfg["file"])
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
    bias = win["diff_staz_meno_era5"].mean()
    day = df.between_time("07:00", "15:00")["diff_staz_meno_era5"]
    log(f"\nAl passaggio ({t_pass:%H:%M} UTC): stazione {win.T_staz_C.mean():.2f} C, "
        f"ERA5 {win.T_era5_C.mean():.2f} C -> ERA5 troppo fresca di {bias:.2f} K")
    log(f"Ore 07-15 UTC: differenza media {day.mean():+.2f} K (min {day.min():+.2f}, max {day.max():+.2f})")
    log(f"Correzione suggerita: tseb.ta_offset_K = {round(bias * 2) / 2:+.1f} "
        "(arrotondata a 0,5 K; una sola stazione -> da trattare come ipotesi)")

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
    log(f"Figura: outputs/figures/p2_05_station_vs_era5.png")

    (PROJECT_ROOT / "outputs" / "p2_05_summary.txt").write_text("\n".join(LINES) + "\n",
                                                                 encoding="utf-8")


if __name__ == "__main__":
    main()
