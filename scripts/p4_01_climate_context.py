"""
p4_01_climate_context.py - Quanto era anomalo il 2022? Precipitazioni e temperatura ERA5-Land
sul Delta del Po, 1991-2025, per scegliere gli anni di confronto con un criterio oggettivo.

Perche': per dire "la siccita' del 2022 ha ridotto l'ET del X%" serve un confronto con anni
NORMALI. Scegliere gli anni "a sensazione" e' un punto debole; qui li scegliamo con un dato:
  - precipitazione cumulata gennaio-luglio (la pioggia che ricarica suolo e falde prima
    dell'estate) in % della media 1991-2020 (il "trentennio climatologico" standard WMO);
  - anomalia della temperatura media giugno-agosto rispetto al 1991-2020.

Dato: ERA5-Land, medie mensili (reanalysis-era5-land-monthly-means), dal Climate Data Store.
  total_precipitation nelle medie mensili e' la pioggia MEDIA GIORNALIERA del mese (metri):
  il totale del mese = valore x giorni del mese x 1000 (mm).
Media spaziale sulle sole celle di terra (ERA5-Land non ha valori sul mare).

Verifica a terra: la pioggia di ERA5-Land a 9 km puo' essere molto lontana da quella locale.
Se c'e' il file giornaliero della stazione (Dext3r), si confrontano i totali gennaio-luglio e
gli anni di riferimento si scelgono sulla misura a terra.

Output:
  data/raw/climate/era5land_monthly_1991_2025.nc   (escluso da Git)
  outputs/tables/p4_01_climate_context.csv
  outputs/figures/p4_01_climate_context.png

Come eseguirlo:
    python scripts/p4_01_climate_context.py
"""

import calendar
import zipfile
from pathlib import Path

import cdsapi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from config_utils import load_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]
YEARS = list(range(1991, 2026))
NORMAL = (1991, 2020)


def download(cfg, out_file):
    """Scarica (una sola volta) le medie mensili di pioggia e temperatura."""
    if out_file.exists():
        print(f"{out_file.name} esiste gia': salto il download.")
        return
    a, m = cfg["training_aoi"], cfg["climate"]["margin_deg"]
    request = {
        "product_type": ["monthly_averaged_reanalysis"],
        "variable": ["total_precipitation", "2m_temperature"],
        "year": [str(y) for y in YEARS],
        "month": [f"{mo:02d}" for mo in range(1, 13)],
        "time": ["00:00"],
        "area": [a["north"] + m, a["west"] - m, a["south"] - m, a["east"] + m],   # N, W, S, E
        "data_format": "netcdf",
        "download_format": "unarchived",
    }
    print("Richiesta al Climate Data Store (medie mensili 1991-2025, puo' restare in coda) ...")
    cdsapi.Client().retrieve("reanalysis-era5-land-monthly-means", request).download(str(out_file))
    if zipfile.is_zipfile(out_file):                     # a volte il CDS consegna uno zip
        tmp = out_file.with_suffix("")
        with zipfile.ZipFile(out_file) as z:
            z.extractall(tmp)
        merged = xr.merge([xr.open_dataset(f).load() for f in sorted(tmp.glob("*.nc"))],
                          compat="override")
        out_file.unlink()
        merged.to_netcdf(out_file)


def monthly_series(path):
    """Serie mensili mediate sulle celle di terra: pioggia (mm/mese) e temperatura (°C)."""
    ds = xr.open_dataset(path)
    if "valid_time" in ds.dims:
        ds = ds.rename({"valid_time": "time"})
    tp = ds["tp"].mean(["latitude", "longitude"], skipna=True).to_series()
    t2 = (ds["t2m"].mean(["latitude", "longitude"], skipna=True) - 273.15).to_series()
    t = pd.to_datetime(tp.index)
    days = np.array([calendar.monthrange(d.year, d.month)[1] for d in t])
    df = pd.DataFrame({"P_mm": tp.values * days * 1000, "T_C": t2.values}, index=t)
    df["anno"], df["mese"] = df.index.year, df.index.month
    return df


def station_precip(path, max_missing_days=10):
    """
    Pioggia gennaio-luglio misurata alla stazione (CSV giornaliero di Dext3r).
    Un anno vale solo se mancano al massimo max_missing_days giorni: con troppi buchi
    il totale sarebbe sottostimato, e meglio non averlo che averlo sbagliato.
    """
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    rows = []
    for l in lines:
        parts = l.split(",")
        if len(parts) >= 3 and parts[0][:2] == "20":
            t = pd.Timestamp(parts[0]).tz_convert(None) + pd.Timedelta(hours=12)  # giorno di misura
            rows.append((t, float(parts[2]) if parts[2].strip() else np.nan))
    s = pd.Series(dict(rows))
    s = s[s.index.month <= 7]
    out = {}
    for year, g in s.groupby(s.index.year):
        missing = int(g.isna().sum())
        out[year] = g.sum() if missing <= max_missing_days else np.nan
    return pd.Series(out, name="P_stazione_gen_lug_mm")


def main():
    cfg = load_config(verbose=False)
    out_dir = PROJECT_ROOT / "data" / "raw" / "climate"
    out_dir.mkdir(parents=True, exist_ok=True)
    nc = out_dir / "era5land_monthly_1991_2025.nc"
    download(cfg, nc)

    m = monthly_series(nc)
    p_jj = m[m["mese"] <= 7].groupby("anno")["P_mm"].sum()               # gennaio-luglio
    t_ja = m[m["mese"].between(6, 8)].groupby("anno")["T_C"].mean()       # giugno-agosto
    normal = (p_jj.index >= NORMAL[0]) & (p_jj.index <= NORMAL[1])
    p_norm, t_norm = p_jj[normal].mean(), t_ja[normal].mean()

    df = pd.DataFrame({"P_gen_lug_mm": p_jj.round(0),
                       "P_%_normale": (100 * p_jj / p_norm).round(0),
                       "T_giu_ago_C": t_ja.round(2),
                       "anomalia_T_K": (t_ja - t_norm).round(2)})
    df["rango_secco_1991_2025"] = df["P_gen_lug_mm"].rank().astype(int)   # 1 = il piu' secco

    # Verifica con la pioggia misurata a terra (stazione ARPAE GIRALDA), dove disponibile
    st_file = PROJECT_ROOT / cfg["stations"]["precip_file"]
    if st_file.exists():
        df["P_stazione_gen_lug_mm"] = station_precip(st_file).round(0)
    df.to_csv(PROJECT_ROOT / cfg["paths"]["tables_common"] / "p4_01_climate_context.csv")

    print(f"\nNormale {NORMAL[0]}-{NORMAL[1]}: pioggia gennaio-luglio {p_norm:.0f} mm, "
          f"temperatura giugno-agosto {t_norm:.1f} C")
    recent = df.loc[cfg["date_search"]["years"]]
    print("\nAnni dell'era Sentinel-2 / Sentinel-3:")
    print(recent.to_string())

    # Anni di riferimento: se c'e' la stazione, gli anni con misura a terra e piu' pioggia del
    # 2022; altrimenti i piu' vicini alla norma ERA5 per pioggia e temperatura.
    if "P_stazione_gen_lug_mm" in recent and recent["P_stazione_gen_lug_mm"].notna().sum() > 1:
        st = recent["P_stazione_gen_lug_mm"].dropna()
        print("\nPioggia gennaio-luglio, stazione contro ERA5-Land (mm):")
        print(pd.DataFrame({"stazione": st, "ERA5-Land": recent.loc[st.index, "P_gen_lug_mm"]})
              .to_string())
        ref = [y for y in st.index if y != 2022 and st[y] > st.get(2022, np.inf)]
        print(f"\nAnni di riferimento (misurati a terra, piu' piovosi del 2022): {ref}")
    else:
        cand = recent.drop(index=2022, errors="ignore")
        score = (cand["P_%_normale"] - 100).abs() / 10 + cand["anomalia_T_K"].abs()
        ref = score.sort_values().index[:2].tolist()
        print(f"\nAnni di riferimento suggeriti (piu' vicini alla norma ERA5): {ref}")

    fig, axs = plt.subplots(1, 3, figsize=(20, 5), constrained_layout=True)
    colors = ["tab:red" if y == 2022 else ("tab:blue" if y in ref else "0.7") for y in df.index]
    axs[0].bar(df.index, df["P_%_normale"], color=colors)
    axs[0].axhline(100, color="k", lw=1)
    axs[0].set_ylabel(f"pioggia gennaio-luglio (% della media {NORMAL[0]}-{NORMAL[1]})")
    axs[0].set_title("Pioggia prima dell'estate")
    axs[1].bar(df.index, df["anomalia_T_K"], color=colors)
    axs[1].axhline(0, color="k", lw=1)
    axs[1].set_ylabel(f"anomalia temperatura giugno-agosto (K, rispetto a {NORMAL[0]}-{NORMAL[1]})")
    axs[1].set_title("Caldo estivo")
    if "P_stazione_gen_lug_mm" in df:
        both = df.dropna(subset=["P_stazione_gen_lug_mm"])
        x = np.arange(len(both))
        axs[2].bar(x - 0.2, both["P_gen_lug_mm"], 0.4, label="ERA5-Land (9 km)", color="0.6")
        axs[2].bar(x + 0.2, both["P_stazione_gen_lug_mm"], 0.4, label="stazione GIRALDA",
                   color="tab:green")
        axs[2].set_xticks(x, both.index)
        axs[2].set_ylabel("pioggia gennaio-luglio (mm)")
        axs[2].set_title("Verifica a terra: ERA5-Land contro stazione")
        axs[2].legend()
    fig.suptitle("Delta del Po (rosso = 2022, blu = anni di riferimento scelti)")
    out_fig = PROJECT_ROOT / cfg["paths"]["figures_common"] / "p4_01_climate_context.png"
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, dpi=150)
    print(f"Figura: {out_fig.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
