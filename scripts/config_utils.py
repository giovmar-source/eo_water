"""
config_utils.py - Lettura della configurazione con la data come parametro.

Tutti gli script della Fase 1 la usano cosi':

    from config_utils import load_config
    cfg = load_config()

e si lanciano scegliendo la data dalla riga di comando:

    python scripts/p1_02_download.py --date 2022-08-11

Senza --date si usa reference_date del file di configurazione.

Cosa fa load_config():
  1. legge config/phase1_sharpening.yaml;
  2. sostituisce {date} nei percorsi (es. data/raw/p1_{date} -> data/raw/p1_20220811),
     cosi' ogni data ha le sue cartelle e niente viene sovrascritto;
  3. applica le eccezioni della sezione "dates" (impostazioni valide solo per una data);
  4. se l'ora del passaggio Sentinel-3 e' "auto", la sceglie guardando i dati scaricati:
     il passaggio diurno (8-12 UTC) con piu' pixel senza nuvole;
  5. crea le cartelle di output.

Perche' un modulo comune: la stessa logica copiata in 12 script finirebbe per divergere.
"Una sola fonte di verita'" e' una regola base del software scientifico.
"""

import argparse
import copy
import json
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"
CLOUD_BIT_S3 = 14          # bit "nuvola" del campo confidence_in di Sentinel-3 SLSTR


def deep_update(base, extra):
    """Unisce due dizionari annidati: i valori di 'extra' sostituiscono quelli di 'base'."""
    for key, value in extra.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def pick_s3_pass(path, hours=(8, 12), target_hour=10.0):
    """
    Sceglie il passaggio Sentinel-3 diurno con piu' pixel limpidi.
    A parita' (entro il 5%), preferisce quello piu' vicino a target_hour (10 UTC, oppure
    l'ora di Landsat se quel giorno c'e' Landsat).
    Restituisce l'ora come testo "HH:MM".
    """
    import xarray as xr
    ds = xr.open_dataset(path)
    if "t" in ds.dims:
        ds = ds.rename({"t": "time"})
    times = pd.to_datetime(ds["time"].values)
    rows = []
    for i, t in enumerate(times):
        h = t.hour + t.minute / 60
        if not hours[0] <= h <= hours[1]:
            continue
        scene = ds.isel(time=i)
        flags = scene["confidence_in"].fillna(0).astype("int64")
        clear = scene["LST"].notnull() & ((flags & (1 << CLOUD_BIT_S3)) == 0)
        rows.append((t, int(clear.sum()), abs(h - target_hour)))
    ds.close()
    if not rows:
        raise SystemExit(f"Nessun passaggio Sentinel-3 tra le {hours[0]} e le {hours[1]} UTC in {path}")
    best_n = max(r[1] for r in rows)
    good = [r for r in rows if r[1] >= 0.95 * best_n]
    t, n, _ = min(good, key=lambda r: r[2])
    print("Sentinel-3, passaggi diurni: " +
          ", ".join(f"{r[0]:%H:%M} ({r[1]} pixel limpidi)" for r in rows) + f" -> scelto {t:%H:%M}")
    return f"{t:%H:%M}"


def load_config(verbose=True):
    """Configurazione per la data scelta con --date (default: reference_date)."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--date", help="data da elaborare, formato AAAA-MM-GG")
    args, _ = parser.parse_known_args()

    cfg = yaml.safe_load(CONFIG_FILE.read_text(encoding="utf-8"))
    date = args.date or str(cfg["reference_date"])
    pd.Timestamp(date)                          # errore chiaro se la data e' scritta male
    cfg["reference_date"] = date
    compact = date.replace("-", "")

    # Eccezioni valide solo per questa data (sezione "dates" del config)
    overrides = (cfg.get("dates") or {}).get(date, {})
    deep_update(cfg, copy.deepcopy(overrides))

    # Percorsi con {date}
    cfg["paths"] = {k: v.format(date=compact) for k, v in cfg["paths"].items()}
    if "stations" in cfg:
        cfg["stations"]["file"] = cfg["stations"]["file"].format(date=compact)

    # Ora del passaggio Sentinel-3 (serve a sharpening, TSEB e confronto con la stazione)
    s3_file = PROJECT_ROOT / cfg["paths"]["raw"] / "s3_lst.nc"
    if str(cfg["training"]["s3_time_utc"]).lower() == "auto":
        # Se quel giorno c'e' Landsat, si preferisce il passaggio S3 piu' vicino a Landsat:
        # la superficie si scalda in fretta la mattina, quindi conta la simultaneita'.
        target = 10.0
        meta = PROJECT_ROOT / cfg["paths"]["raw"] / "landsat_metadata.json"
        if meta.exists():
            t_ls = pd.Timestamp(json.loads(meta.read_text(encoding="utf-8"))["datetime_utc"])
            target = t_ls.hour + t_ls.minute / 60
        cfg["training"]["s3_time_utc"] = (pick_s3_pass(s3_file, target_hour=target)
                                          if s3_file.exists() else None)
    cfg["tseb"]["overpass_utc"] = cfg["training"]["s3_time_utc"]

    for key in ("tables", "figures", "summaries", "processed", "raw"):
        (PROJECT_ROOT / cfg["paths"][key]).mkdir(parents=True, exist_ok=True)

    if verbose:
        print(f"=== Data {date} | Sentinel-3 {cfg['training']['s3_time_utc'] or 'non ancora scaricato'} "
              f"| dati in {cfg['paths']['raw']} ===")
    return cfg


def ta_offset_for(cfg):
    """
    Correzione della temperatura dell'aria ERA5 per la data (K).
    - numero nel config  -> si usa quello;
    - "auto"             -> si legge il risultato del confronto con la stazione (p2_05);
                            se manca, 0 con un avviso (e va dichiarato).
    """
    value = cfg["tseb"].get("ta_offset_K", 0.0)
    if str(value).lower() != "auto":
        return float(value), "dal file di configurazione"
    f = PROJECT_ROOT / cfg["paths"]["summaries"] / "p2_05_ta_offset.yaml"
    if f.exists():
        info = yaml.safe_load(f.read_text(encoding="utf-8"))
        return float(info["ta_offset_K"]), f"dalla stazione {info['stazione']} (p2_05)"
    print("ATTENZIONE: nessun confronto con la stazione per questa data -> correzione Ta = 0 K")
    return 0.0, "nessuna stazione: NESSUNA correzione"


if __name__ == "__main__":
    c = load_config()
    print(yaml.safe_dump({"paths": c["paths"], "s3": c["training"]["s3_time_utc"]}, allow_unicode=True))
