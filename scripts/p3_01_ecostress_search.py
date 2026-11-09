"""
p3_01_ecostress_search.py - Cerca e scarica i passaggi ECOSTRESS diurni sul Delta del Po.

Perche' ECOSTRESS: e' un radiometro termico della NASA montato sulla Stazione Spaziale (ISS),
con pixel da 70 m. E' un dato INDIPENDENTE dalla nostra catena (Sentinel-3 affinato con
Sentinel-2 + TSEB) e serve per validarla:
  - LST a 70 m       -> confronto con la nostra LST affinata
  - ET giornaliera   -> confronto con l'ET di TSEB (prodotto ECOSTRESS L3 JET)
  - ESI              -> indice di stress evaporativo (ET / ET potenziale), confronto col CWSI

Particolarita': la ISS non e' su un'orbita eliofila come Sentinel, quindi ECOSTRESS passa
a ORE DIVERSE ogni giorno (mattina, mezzogiorno, sera, notte). Per confrontarlo con
Sentinel-3 (ore 10 circa) servono passaggi a meta' mattina: li filtriamo per ora.

Prodotti (versione 002, "tiled" = ritagliati su tasselli MGRS da ~110 km, file GeoTIFF):
  ECO_L2T_LSTE : LST, cloud (maschera nuvole), water, QC
  ECO_L3T_JET  : ETdaily (ET giornaliera, mm)
  ECO_L4T_ESI  : ESI (indice di stress evaporativo)

Credenziali: NASA Earthdata, lette da _netrc nella cartella utente (mai nel codice).

Output:
  data/raw/ecostress/<data_ora>/<tassello>_<variabile>.tif  (ritagliati sull'area di analisi)
  outputs/tables/p3_01_ecostress_passes.csv  (passaggi, ora, copertura nuvolosa sull'area)

Come eseguirlo:
    python scripts/p3_01_ecostress_search.py
"""

from pathlib import Path

import earthaccess
import numpy as np
import pandas as pd
import rioxarray
import yaml

from config_utils import load_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"


def granule_info(g):
    """Ora di acquisizione, nome e link ai file di un granulo ECOSTRESS."""
    umm = g["umm"]
    t = pd.Timestamp(umm["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"]).tz_convert(None)
    name = umm["GranuleUR"]
    tile = name.split("_")[5]                       # es. 32TQR
    links = [u for u in g.data_links() if u.endswith(".tif")]
    return t, name, tile, links


def pick(links, suffix):
    """Il link del file che finisce con _<suffix>.tif (es. _LST.tif, _cloud.tif)."""
    hits = [u for u in links if u.endswith(f"_{suffix}.tif")]
    return hits[0] if hits else None


def clip_to_aoi(path, bbox):
    """Ritaglia un GeoTIFF (in UTM, 70 m) sull'area di analisi espressa in gradi."""
    da = rioxarray.open_rasterio(path, masked=True).squeeze("band", drop=True)
    return da.rio.clip_box(*bbox, crs="EPSG:4326").load()


def main():
    cfg_all = load_config()
    ec = cfg_all["ecostress"]
    a = cfg_all["analysis_aoi"]
    bbox = (a["west"], a["south"], a["east"], a["north"])
    raw = PROJECT_ROOT / "data" / "raw" / "ecostress"
    tmp = PROJECT_ROOT / "data" / "raw" / "ecostress" / "_tiles"   # tasselli interi (cancellabili)
    tmp.mkdir(parents=True, exist_ok=True)

    earthaccess.login(strategy="netrc")          # usa le credenziali salvate in _netrc

    # 1. Ricerca dei granuli LST nell'intervallo, poi filtro sulle ore diurne
    found = earthaccess.search_data(short_name="ECO_L2T_LSTE", version="002",
                                    bounding_box=bbox,
                                    temporal=(ec["start"], ec["end"]))
    print(f"Granuli LST trovati: {len(found)}")
    h0, h1 = ec["hours_utc"]
    def hour_of(g):
        t = granule_info(g)[0]
        return t.hour + t.minute / 60
    lst_granules = [g for g in found if h0 <= hour_of(g) <= h1]
    print(f"Di giorno tra le {h0} e le {h1} UTC: {len(lst_granules)}")

    rows = []
    for g in lst_granules:
        t, name, tile, links = granule_info(g)
        stamp = f"{t:%Y%m%dT%H%M}"
        files = {v: pick(links, v) for v in ("LST", "cloud", "water", "QC")}
        files = {k: v for k, v in files.items() if v}
        local = earthaccess.download(list(files.values()), local_path=str(tmp))
        local = {k: Path(p) for k, p in zip(files, local)}

        lst = clip_to_aoi(local["LST"], bbox)
        cloud = clip_to_aoi(local["cloud"], bbox)
        water = clip_to_aoi(local["water"], bbox) if "water" in local else lst * 0
        land = (water.values == 0)
        covered = land & np.isfinite(lst.values)
        cloudy = land & (cloud.values == 1)
        row = {"data_ora_utc": t, "tassello": tile, "granulo": name,
               "frazione_area_coperta": round(covered.mean(), 3),
               "nuvole_su_terra_%": round(100 * cloudy.sum() / max(land.sum(), 1), 1),
               "LST_mediana_C": round(float(np.nanmedian(lst.values[covered & ~cloudy])) - 273.15, 1)
               if (covered & ~cloudy).any() else np.nan}
        rows.append(row)
        print(f"  {t:%Y-%m-%d %H:%M} {tile}: area coperta {row['frazione_area_coperta']:.0%}, "
              f"nuvole {row['nuvole_su_terra_%']}%, LST mediana {row['LST_mediana_C']} C")

        # Salva i ritagli dei passaggi utili, poi scarica anche ET ed ESI di quel passaggio
        if row["frazione_area_coperta"] >= ec["min_coverage"] and \
                row["nuvole_su_terra_%"] <= ec["max_cloud_pct"]:
            out = raw / stamp
            out.mkdir(parents=True, exist_ok=True)
            for k, p in local.items():
                clip_to_aoi(p, bbox).rio.to_raster(out / f"{tile}_{k}.tif", compress="deflate")
            for short, var in (("ECO_L3T_JET", "ETdaily"), ("ECO_L4T_ESI", "ESI")):
                same = earthaccess.search_data(short_name=short, version="002",
                                               bounding_box=bbox,
                                               temporal=(t - pd.Timedelta(minutes=2),
                                                         t + pd.Timedelta(minutes=2)))
                same = [s for s in same if granule_info(s)[2] == tile]
                if not same:
                    print(f"    {short}: nessun granulo corrispondente")
                    continue
                url = pick(granule_info(same[0])[3], var)
                if url is None:
                    print(f"    {short}: file {var} non trovato")
                    continue
                p = Path(earthaccess.download([url], local_path=str(tmp))[0])
                clip_to_aoi(p, bbox).rio.to_raster(out / f"{tile}_{var}.tif", compress="deflate")
                print(f"    salvato {var}")

    df = pd.DataFrame(rows).sort_values(["data_ora_utc", "tassello"])
    out_csv = PROJECT_ROOT / cfg_all["paths"]["tables_common"] / "p3_01_ecostress_passes.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nTabella: {out_csv.relative_to(PROJECT_ROOT)}")
    print(df.drop(columns="granulo").to_string(index=False))
    print(f"\nI tasselli interi in {tmp.relative_to(PROJECT_ROOT)} si possono cancellare.")


if __name__ == "__main__":
    main()
