"""
step02_first_look.py - Primo sguardo: termico Sentinel-3 (1 km) vs ottico Sentinel-2 (20 m).

Cosa fa:
  1. Apre i due NetCDF scaricati da step02_download.py.
  2. Per ogni data calcola quanta parte dell'area e' libera da nuvole:
       - Sentinel-2: dalla classificazione della scena (SCL)
       - Sentinel-3: dal bit 14 del flag di qualita' confidence_in (= nuvola)
  3. Sceglie la data Sentinel-2 piu' limpida e la data Sentinel-3 limpida piu' vicina.
  4. Calcola l'NDVI (Sentinel-2) e la LST in gradi Celsius (Sentinel-3).
  5. Salva una figura affiancata e un riepilogo testuale.

Come eseguirlo:
    python scripts/step02_first_look.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "step02_first_look.yaml"

# Classi SCL di Sentinel-2 considerate non valide:
# 0 nessun dato, 1 saturato/difettoso, 3 ombra di nuvola, 8-9 nuvola, 10 cirro
SCL_INVALID = [0, 1, 3, 8, 9, 10]
# Bit 14 di confidence_in (Sentinel-3 SLSTR L2 LST) = pixel nuvoloso
CLOUD_BIT_S3 = 14
# Passaggi diurni di Sentinel-3 sull'Italia: circa 09-10 UTC (quelli notturni sono verso le 20-21 UTC).
# Per lo stress idrico serve la temperatura di giorno, quando le piante traspirano.
DAY_HOURS_UTC = (8, 12)

LINES = []


def log(text=""):
    print(text)
    LINES.append(text)


def open_cube(path):
    """Apre un NetCDF di openEO e rinomina la dimensione temporale in 'time'."""
    if not path.exists():
        raise FileNotFoundError(f"{path} non trovato: esegui prima step02_download.py")
    ds = xr.open_dataset(path).load()  # .load(): legge tutto e chiude il file (importante su Windows)
    ds.close()
    time_dim = next((d for d in ("t", "time") if d in ds.dims), None)
    if time_dim is None:
        raise ValueError(f"Nessuna dimensione temporale in {path.name}: dimensioni {dict(ds.sizes)}")
    if time_dim != "time":
        ds = ds.rename({time_dim: "time"})
    return ds


def s2_clear_fraction(ds):
    """Frazione di pixel validi (non nuvola/ombra/vuoto) per ogni data Sentinel-2."""
    scl = ds["SCL"]
    valid = scl.notnull() & ~scl.isin(SCL_INVALID)
    return valid.mean(dim=("y", "x")).to_series()


def s3_cloud_mask(ds):
    """True dove il bit 14 di confidence_in segnala nuvola."""
    conf = ds["confidence_in"]
    conf_int = conf.fillna(0).astype("int64")  # i flag vanno letti come bit di un intero
    return (conf_int & (1 << CLOUD_BIT_S3)) != 0


def s3_daytime(ds):
    """Tiene solo i passaggi diurni di Sentinel-3 (il prodotto LST contiene anche quelli notturni)."""
    hours = ds["time"].dt.hour
    return ds.sel(time=(hours >= DAY_HOURS_UTC[0]) & (hours < DAY_HOURS_UTC[1]))


def s3_clear_fraction(ds):
    """Frazione dei pixel di TERRA con LST valida e non nuvolosa, per ogni passaggio Sentinel-3.

    Il prodotto LST non ha valori sul mare: i pixel di mare non vanno contati come "nuvola".
    Consideriamo "terra" i pixel che hanno una LST valida in almeno un passaggio.
    """
    land = ds["LST"].notnull().any(dim="time")
    valid = ds["LST"].notnull() & ~s3_cloud_mask(ds) & land
    return (valid.sum(dim=("y", "x")) / land.sum()).to_series()


def extent(da):
    """Estensione [xmin, xmax, ymin, ymax] in metri, per affiancare le mappe correttamente."""
    x, y = da["x"].values, da["y"].values
    dx = abs(x[1] - x[0]) if x.size > 1 else 0
    dy = abs(y[1] - y[0]) if y.size > 1 else 0
    return [x.min() - dx / 2, x.max() + dx / 2, y.min() - dy / 2, y.max() + dy / 2]


def main():
    cfg = yaml.safe_load(open(CONFIG_FILE, encoding="utf-8"))
    raw = PROJECT_ROOT / cfg["paths"]["raw"]
    fig_dir = PROJECT_ROOT / cfg["paths"]["figures"]
    fig_dir.mkdir(parents=True, exist_ok=True)

    s2 = open_cube(raw / "s2_l2a.nc")
    s3 = open_cube(raw / "s3_lst.nc")

    log(f"Sentinel-2: {s2.sizes['time']} date, griglia {s2.sizes['y']} x {s2.sizes['x']} pixel")
    log(f"Sentinel-3: {s3.sizes['time']} passaggi, griglia {s3.sizes['y']} x {s3.sizes['x']} pixel")
    s3 = s3_daytime(s3)
    log(f"Sentinel-3: {s3.sizes['time']} passaggi diurni tenuti (tra le {DAY_HOURS_UTC[0]} e le {DAY_HOURS_UTC[1]} UTC)")

    # --- Quanto e' limpida ogni data? ---
    s2_clear = s2_clear_fraction(s2)
    s3_clear = s3_clear_fraction(s3)
    log("\nFrazione di area libera da nuvole (Sentinel-3: solo pixel di terra)")
    log("Sentinel-2:")
    for t, f in s2_clear.items():
        log(f"  {pd.Timestamp(t).date()}  {f:6.1%}")
    log("Sentinel-3 (ora UTC):")
    for t, f in s3_clear.items():
        log(f"  {pd.Timestamp(t):%Y-%m-%d %H:%M}  {f:6.1%}")

    # --- Scelta delle date ---
    t_s2 = s2_clear.idxmax()
    candidates = s3_clear[s3_clear >= 0.8]
    if candidates.empty:
        log("\nATTENZIONE: nessuna data Sentinel-3 con almeno 80% di area limpida; uso la migliore.")
        candidates = s3_clear
    # distanza in ore (non in giorni), cosi' tra due passaggi dello stesso giorno sceglie il piu' vicino
    gap_hours = pd.Series(np.abs((candidates.index - t_s2) / pd.Timedelta(hours=1)), index=candidates.index)
    t_s3 = gap_hours.idxmin()
    log(f"\nData Sentinel-2 scelta: {pd.Timestamp(t_s2).date()} ({s2_clear[t_s2]:.0%} limpida)")
    log(f"Passaggio Sentinel-3 scelto: {pd.Timestamp(t_s3):%Y-%m-%d %H:%M} UTC ({s3_clear[t_s3]:.0%} limpido)")

    # --- NDVI da Sentinel-2 ---
    scene = s2.sel(time=t_s2)
    # I valori Sentinel-2 L2A sono riflettanze moltiplicate per 10000.
    # Su acqua e zone molto scure possono essere negativi (correzione atmosferica + offset):
    # in quei casi l'NDVI esplode (es. -435), quindi teniamo solo riflettanze positive.
    red = scene["B04"].astype("float32") / 10000
    nir = scene["B08"].astype("float32") / 10000
    invalid = scene["SCL"].isin(SCL_INVALID) | scene["SCL"].isnull() | (red <= 0) | (nir <= 0)
    ndvi = ((nir - red) / (nir + red)).where(~invalid)

    # --- LST da Sentinel-3 (Kelvin -> gradi Celsius) ---
    s3_scene = s3.sel(time=t_s3)
    lst_c = (s3_scene["LST"] - 273.15).where(~s3_cloud_mask(s3_scene))

    log(f"\nNDVI: min {float(ndvi.min()):.2f}, mediana {float(ndvi.median()):.2f}, max {float(ndvi.max()):.2f}")
    log(f"LST:  min {float(lst_c.min()):.1f} C, mediana {float(lst_c.median()):.1f} C, max {float(lst_c.max()):.1f} C")
    log("Controllo: d'estate di giorno la LST dovrebbe stare circa tra 25 e 50 C;")
    log("valori intorno a 15-25 C potrebbero indicare un passaggio notturno.")

    # --- Figura affiancata ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    im0 = axes[0].imshow(ndvi.values, extent=extent(ndvi), cmap="YlGn", vmin=0, vmax=0.9)
    axes[0].set_title(f"Sentinel-2 NDVI, 20 m\n{pd.Timestamp(t_s2).date()}")
    fig.colorbar(im0, ax=axes[0], shrink=0.8, label="NDVI")

    im1 = axes[1].imshow(lst_c.values, extent=extent(lst_c), cmap="inferno")
    axes[1].set_title(f"Sentinel-3 LST diurna, ~1 km\n{pd.Timestamp(t_s3):%Y-%m-%d %H:%M} UTC")
    fig.colorbar(im1, ax=axes[1], shrink=0.8, label="Temperatura superficiale (°C)")

    # Stessa finestra geografica per entrambe le mappe (quella di Sentinel-2)
    xmin, xmax, ymin, ymax = extent(ndvi)
    for ax in axes:
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel("UTM 33N, x (m)")
        ax.set_ylabel("UTM 33N, y (m)")
        ax.ticklabel_format(style="plain")
        ax.tick_params(axis="x", labelrotation=30)
    fig.suptitle(cfg["aoi"]["description"])

    fig_path = fig_dir / "step02_ndvi_vs_lst.png"
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)
    log(f"\nFigura salvata: {fig_path.relative_to(PROJECT_ROOT)}")

    summary = PROJECT_ROOT / "outputs" / "step02_summary.txt"
    summary.write_text("\n".join(LINES), encoding="utf-8")
    print(f"Riepilogo salvato: {summary.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
