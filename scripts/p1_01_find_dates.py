"""
p1_01_find_dates.py - Trova i giorni con passaggio comune di Sentinel-2, Sentinel-3 e Landsat 8/9.

Perche': per validare lo sharpening del termico serve un giorno in cui abbiamo
  - Sentinel-3 (termico a ~1 km, da affinare),
  - Sentinel-2 (ottico a 20 m, i predittori),
  - Landsat 8 o 9 (termico a 100 m, il riferimento indipendente),
tutti con poche nuvole. Questo script NON scarica immagini: interroga solo i cataloghi
STAC (gli "indici" degli archivi satellitari), quindi e' veloce e non consuma crediti.

Cosa produce:
  - una tabella a schermo con i giorni candidati, dal piu' limpido al piu' nuvoloso;
  - outputs/tables/p1_coincident_dates.csv con la stessa tabella.

Come eseguirlo:
    python scripts/p1_01_find_dates.py
"""

import time
from pathlib import Path

import pandas as pd
import yaml
from pystac_client import Client
from pystac_client.stac_api_io import StacApiIO
from urllib3.util.retry import Retry

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"

# Comportamento "educato" verso le API pubbliche:
#  - in caso di risposta 429 (troppe richieste) o errori temporanei del server,
#    riprova fino a 6 volte aspettando sempre di piu' (circa 2, 4, 8, 16... secondi);
#  - rispetta l'eventuale tempo di attesa indicato dal server (Retry-After).
RETRY = Retry(
    total=6,
    backoff_factor=2,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=None,          # riprova anche le richieste POST usate dalla ricerca STAC
    respect_retry_after_header=True,
)
PAGE_SIZE = 100                    # scene per pagina: meno pagine = meno richieste al server
PAUSE_BETWEEN_SEARCHES_S = 5       # pausa tra una ricerca e l'altra


def bbox_of(aoi):
    """Riquadro [ovest, sud, est, nord] come lo vogliono i cataloghi STAC."""
    return [aoi["west"], aoi["south"], aoi["east"], aoi["north"]]


def search(url, collection, bbox, period):
    """Cerca nel catalogo STAC tutte le scene che toccano l'area nel periodo indicato.

    Restituisce una tabella (DataFrame) con una riga per scena.
    """
    print(f"Cerco {collection} ...", end=" ", flush=True)
    try:
        catalog = Client.open(url, stac_io=StacApiIO(max_retries=RETRY))
        found = catalog.search(
            collections=[collection],
            bbox=bbox,
            datetime=f"{period['start']}/{period['end']}",
            limit=PAGE_SIZE,
            max_items=5000,
        )
        items = list(found.items())
    except Exception as error:  # rete assente, catalogo non raggiungibile, ecc.
        print("ERRORE")
        raise SystemExit(f"Impossibile interrogare {url}: {error}")

    rows = []
    for item in items:
        props = item.properties
        when = pd.Timestamp(item.datetime) if item.datetime else None
        if when is not None and when.tzinfo is not None:
            when = when.tz_convert(None)  # tutto in UTC "semplice", senza fuso orario allegato
        rows.append({
            "id": item.id,
            "datetime_utc": when,
            "date": when.normalize() if when is not None else None,
            "platform": props.get("platform"),
            # "eo:cloud_cover" e' il nome standard STAC; non tutti i prodotti lo hanno
            "cloud": props.get("eo:cloud_cover"),
        })
    print(f"{len(rows)} scene")
    return pd.DataFrame(rows)


def main():
    cfg = yaml.safe_load(open(CONFIG_FILE, encoding="utf-8"))
    bbox = bbox_of(cfg["training_aoi"])
    period = cfg["search_period"]
    thr = cfg["date_search"]
    cat = cfg["catalogs"]

    print(f"Area: {cfg['training_aoi']['name']}  {bbox}")
    print(f"Periodo: {period['start']} -> {period['end']}\n")

    s2 = search(cat["cdse_stac"], cat["s2_collection"], bbox, period)
    time.sleep(PAUSE_BETWEEN_SEARCHES_S)
    s3 = search(cat["cdse_stac"], cat["s3_collection"], bbox, period)
    ls = search(cat["landsat_stac"], cat["landsat_collection"], bbox, period)

    # --- Sentinel-2: piu' tasselli per giorno -> teniamo la copertura nuvolosa MASSIMA (caso peggiore) ---
    if s2["cloud"].isna().all():
        print("ATTENZIONE: il catalogo Sentinel-2 non riporta la copertura nuvolosa; filtro disattivato.")
        s2_day = s2.groupby("date").agg(s2_scenes=("id", "count"))
        s2_day["s2_cloud"] = float("nan")
    else:
        s2_day = s2.groupby("date").agg(s2_scenes=("id", "count"), s2_cloud=("cloud", "max"))
        s2_day = s2_day[s2_day["s2_cloud"] <= thr["max_cloud_s2"]]

    # --- Sentinel-3: solo passaggi diurni ---
    h0, h1 = thr["s3_day_hours_utc"]
    s3 = s3[(s3["datetime_utc"].dt.hour >= h0) & (s3["datetime_utc"].dt.hour < h1)]
    s3_day = s3.groupby("date").agg(
        s3_passes=("id", "count"),
        s3_times=("datetime_utc", lambda t: ", ".join(sorted(t.dt.strftime("%H:%M").unique()))),
    )

    # --- Landsat: solo 8 e 9 (il 7 ha un difetto di acquisizione dal 2003), poche nuvole ---
    ls = ls[ls["platform"].isin(cat["landsat_platforms"])]
    ls = ls[ls["cloud"] <= thr["max_cloud_landsat"]]
    ls_day = ls.groupby("date").agg(
        landsat=("platform", lambda p: "+".join(sorted(p.unique()))),
        landsat_cloud=("cloud", "min"),
        landsat_time=("datetime_utc", lambda t: t.min().strftime("%H:%M")),
    )

    # --- Giorni in cui ci sono tutti e tre (join = "unione per data") ---
    table = s2_day.join(s3_day, how="inner").join(ls_day, how="inner")
    table = table.sort_values(["landsat_cloud", "s2_cloud"])
    table.index = table.index.strftime("%Y-%m-%d")
    table.index.name = "data"

    print(f"\nGiorni con S2 <= {thr['max_cloud_s2']}% nuvole, S3 diurno e Landsat <= {thr['max_cloud_landsat']}% nuvole:")
    if table.empty:
        print("Nessuno. Prova ad alzare le soglie in config/phase1_sharpening.yaml.")
    else:
        with pd.option_context("display.width", 140, "display.max_columns", 20):
            print(table.round(1).to_string())

    out_dir = PROJECT_ROOT / cfg["paths"]["tables"]
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "p1_coincident_dates.csv"
    table.to_csv(out)
    print(f"\nTabella salvata: {out.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
