"""
p2_04_tseb.py - Bilancio energetico a due sorgenti (TSEB) e evapotraspirazione giornaliera a 20 m.

L'idea in una riga: l'energia che arriva dal Sole (radiazione netta Rn) deve andare da qualche
parte. Una parte scalda il suolo (G), una parte scalda l'aria (calore sensibile H), il resto
fa evaporare acqua (calore latente LE):

        Rn = G + H + LE        ->        LE = Rn - G - H   (LE e' l'ET in unita' di energia)

H si stima dalla differenza di temperatura tra superficie e aria: una superficie molto piu'
calda dell'aria cede tanto calore sensibile e quindi "evapora poco".

TSEB (Norman et al. 1995; Kustas & Norman 1999) divide il pixel in DUE sorgenti, suolo e
chioma, ognuna con la sua temperatura e i suoi flussi. La versione Priestley-Taylor (TSEB-PT):
  1. parte dall'ipotesi che la chioma traspiri al potenziale (coefficiente alpha = 1,26);
  2. ricava la temperatura della chioma, e quindi quella del suolo dalla temperatura
     radiometrica del pixel (la nostra LST a 20 m);
  3. se il suolo risulta "impossibile" (LE del suolo < 0, cioe' condensazione di giorno),
     riduce alpha per passi: la chioma e' sotto stress e traspira meno del potenziale.

Input (tutti sulla stessa griglia 20 m, EPSG:32633, ritagliati sull'AOI di analisi):
  - LST affinata (Step 1d, foresta casuale + correzione dei residui)  -> Tr_K
  - LAI di SNAP (Step 2b)                                            -> LAI, altezza chioma
  - meteo ERA5-Land (Step 2a) all'ora del passaggio, interpolato sui pixel
  - geometria solare calcolata per l'ora esatta del passaggio Sentinel-3

Output:
  data/processed/p1_20220719/tseb_*.tif   (Rn, H, LE, G, ET giornaliera, flag; esclusi da Git)
  outputs/figures/p2_04_tseb.png
  outputs/tables/p2_04_tseb_by_lai.csv
  outputs/p2_04_summary.txt

Come eseguirlo:
    python scripts/p2_04_tseb.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray  # noqa: F401  (attiva .rio)
import xarray as xr
import yaml
from pyproj import Transformer
from rasterio.enums import Resampling
from pyTSEB import TSEB
from pyTSEB import meteo_utils as met
from pyTSEB import net_radiation as rad
from pyTSEB import resistances as res

from p1_03_training_table import as_raster, open_nc

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"

LINES = []


def log(text=""):
    print(text)
    LINES.append(text)


# ----------------------------------------------------------------------------
# Lettura e preparazione degli input
# ----------------------------------------------------------------------------

def read_tif(path, bbox, crs):
    """Legge un GeoTIFF a una banda e lo ritaglia sull'AOI (bbox in gradi lon/lat)."""
    da = xr.open_dataarray(path, engine="rasterio").squeeze("band", drop=True)
    da = as_raster(da.rio.clip_box(*bbox, crs="EPSG:4326"), crs).load()
    return da.where(np.isfinite(da))


def pixel_lonlat(grid):
    """Longitudine e latitudine di ogni pixel (servono per meteo e Sole)."""
    xx, yy = np.meshgrid(grid.x.values, grid.y.values)
    to_wgs84 = Transformer.from_crs(grid.rio.crs, "EPSG:4326", always_xy=True)
    lon, lat = to_wgs84.transform(xx, yy)
    return lon, lat


def meteo_at_overpass(era5_file, lon, lat, t_utc):
    """
    Meteo ERA5-Land all'ora del passaggio, interpolato in modo bilineare su ogni pixel.

    Grandezze istantanee (t2m, d2m, vento, pressione): interpolate nel tempo all'istante esatto.
    Radiazione: ERA5-Land la CUMULA dalle 00 UTC (J/m2). Il flusso medio in un intervallo e'
    la differenza dei cumulati divisa per i secondi. Usiamo le due ore attorno al passaggio
    (dalle 09 alle 11 UTC per un passaggio alle 10) cosi' il valore e' centrato sull'istante.
    Il totale giornaliero e' il cumulato alle 00 UTC del giorno dopo (= 24 h di somma).
    """
    ds = xr.open_dataset(era5_file)
    t = pd.Timestamp(t_utc)
    h0 = t.floor("h") - pd.Timedelta(hours=1) if t.minute < 30 else t.floor("h")
    h1 = h0 + pd.Timedelta(hours=2)
    day_end = t.normalize() + pd.Timedelta(days=1)

    inst = ds[["t2m", "d2m", "u10", "v10", "sp"]].interp(time=t)
    sw = (ds["ssrd"].sel(time=h1) - ds["ssrd"].sel(time=h0)) / 7200.0
    lw = (ds["strd"].sel(time=h1) - ds["strd"].sel(time=h0)) / 7200.0
    sw_day = ds["ssrd"].sel(time=day_end)                     # J/m2 nelle 24 h

    fields = {
        "T_A_K": inst["t2m"],
        "Td_K": inst["d2m"],
        "u": np.hypot(inst["u10"], inst["v10"]),            # velocita' del vento (m/s)
        "p": inst["sp"] / 100.0,                             # Pa -> hPa (= mb)
        "S_dn": sw, "L_dn": lw, "S_day": sw_day,
    }
    # Ogni campo si porta dietro la sua coordinata "time" (istanti diversi): la togliamo
    fields = {k: v.drop_vars(["time", "number", "expver"], errors="ignore")
              for k, v in fields.items()}

    # Interpolazione bilineare su tutti i pixel in un colpo solo
    pts_lon = xr.DataArray(lon.ravel(), dims="pix")
    pts_lat = xr.DataArray(lat.ravel(), dims="pix")
    out = {}
    for name, field in fields.items():
        vals = field.interp(longitude=pts_lon, latitude=pts_lat).values
        out[name] = vals.reshape(lon.shape).astype("float64")

    # Pressione di vapore dal punto di rugiada (formula di Tetens/Bolton, hPa)
    td_c = out["Td_K"] - 273.15
    out["ea"] = 6.112 * np.exp(17.67 * td_c / (td_c + 243.5))
    return out


def canopy_height(lai, cfg, crop=None):
    """
    Altezza della chioma (non misurata: e' un'IPOTESI da dichiarare).
    Lineare tra h_min (suolo nudo) e un'altezza massima, raggiunta a LAI = lai_at_hmax.
    - senza mappa colturale: la stessa altezza massima per tutti (h_max_m);
    - con la mappa EUCROPMAP: altezza massima tipica di ogni coltura a meta' luglio
      (mais ~2,5 m, soia ~0,8 m, ...). Il LAI fa da indicatore dello stadio di crescita.
    """
    h_max = np.full_like(lai, cfg["h_max_m"])
    if crop is not None:
        table = {int(k): v for k, v in cfg["crop_h_max_m"].items()}
        h_max = pd.Series(crop).map(table).fillna(cfg["h_max_m"]).to_numpy()
    frac = np.clip(lai / cfg["lai_at_hmax"], 0, 1)
    return cfg["h_min_m"] + (h_max - cfg["h_min_m"]) * frac


def green_fraction(lai, fapar, sza, f_min=0.01):
    """
    Frazione verde della chioma f_g = FAPAR / FIPAR.
      FAPAR = luce assorbita dalle sole parti VERDI (fotosintetiche), dal processore SNAP;
      FIPAR = luce intercettata da TUTTA la chioma, verde o secca, calcolata dal LAI con la
              legge di Beer: FIPAR = 1 - exp(-k * LAI), k = 0,5 / cos(zenit) (foglie sferiche).
    Stoppie di frumento: intercettano un po' di luce ma non ne usano -> f_g ~ 0.
    Colture verdi: f_g ~ 1. In TSEB-PT la traspirazione potenziale e' moltiplicata per f_g.
    """
    k = 0.5 / np.cos(np.radians(sza))
    fipar = 1 - np.exp(-k * lai)
    fg = np.where(fipar > 0.02, fapar / np.maximum(fipar, 1e-6), 1.0)   # suolo nudo: irrilevante
    return np.clip(fg, f_min, 1.0)


# ----------------------------------------------------------------------------
# TSEB
# ----------------------------------------------------------------------------

def run_tseb(Tr_K, lai, mt, sza, cfg):
    """
    Esegue TSEB-PT su vettori 1D di pixel validi. Restituisce un dizionario di flussi (W/m2).
    """
    # 1. Radiazione solare diretta e diffusa, nel visibile (VIS) e nell'infrarosso vicino (NIR).
    #    La chioma assorbe diversamente luce diretta e diffusa, e VIS e NIR (le foglie
    #    assorbono il rosso, riflettono il NIR): per questo si divide in 4 parti.
    difvis, difnir, fvis, fnir = rad.calc_difuse_ratio(mt["S_dn"], sza, press=mt["p"])
    S_dir = mt["S_dn"] * (fvis * (1 - difvis) + fnir * (1 - difnir))
    S_dif = mt["S_dn"] * (fvis * difvis + fnir * difnir)

    # 2. Radiazione solare netta assorbita da chioma e suolo (modello di Campbell & Norman):
    #    piu' LAI -> piu' luce intercettata dalla chioma, meno al suolo.
    #    pyTSEB vuole le proprieta' ottiche come vettori lunghi quanto i pixel.
    optics = [np.full_like(lai, cfg[k]) for k in (
        "rho_leaf_vis", "tau_leaf_vis", "rho_leaf_nir", "tau_leaf_nir",
        "rho_soil_vis", "rho_soil_nir")]
    Sn_C, Sn_S = rad.calc_Sn_Campbell(lai, sza, S_dir, S_dif, fvis, fnir, *optics, x_LAD=1)

    # 3. Rugosita' aerodinamica: quanto la superficie "frena" il vento.
    #    z_0M = rugosita' per la quantita' di moto, d_0 = altezza di spostamento del piano zero.
    h_C = canopy_height(lai, cfg, mt.get("crop"))
    z_0M, d_0 = res.calc_roughness(
        lai, h_C, w_C=np.ones_like(lai), landcover=np.full_like(lai, res.CROP))
    z_0M = np.maximum(z_0M, cfg.get("z0_soil_m", 0.01))   # mai meno rugoso del suolo nudo

    # Quote del meteo. ERA5 da' Ta "a 2 m" e vento "a 10 m" sopra una superficie di riferimento:
    # sopra un mais di 2,5 m, 2 m dal suolo sarebbe DENTRO la chioma, dove la teoria dei profili
    # non vale. Ipotesi: le quote valgono sopra il piano di spostamento d_0 (la "superficie
    # efficace" che il vento vede), quindi z = quota + d_0.
    if cfg.get("heights_above_d0", False):
        z_u, z_T = cfg["z_u_m"] + d_0, cfg["z_T_m"] + d_0
    else:
        z_u, z_T = np.full_like(lai, cfg["z_u_m"]), np.full_like(lai, cfg["z_T_m"])

    # 4. TSEB-PT vero e proprio
    out = TSEB.TSEB_PT(
        Tr_K, np.full_like(lai, cfg["vza_deg"]), mt["T_A_K"], mt["u"], mt["ea"], mt["p"],
        Sn_C, Sn_S, mt["L_dn"], lai, h_C, cfg["emis_C"], cfg["emis_S"], z_0M, d_0,
        z_u, z_T, leaf_width=cfg["leaf_width_m"],
        z0_soil=cfg.get("z0_soil_m", 0.01), alpha_PT=1.26, x_LAD=1, f_c=1.0,
        f_g=mt.get("f_g", np.ones_like(lai)),   # frazione verde: solo il verde traspira
        calcG_params=[[1], 0.35])   # G = 35% della radiazione netta del suolo
    (flag, T_S, T_C, T_AC, Ln_S, Ln_C, LE_C, H_C, LE_S, H_S, G,
     R_S, R_x, R_A, u_star, L, n_iter, _alpha) = out
    # NB: l'alpha restituito da pyTSEB 2.5.4 NON e' affidabile come diagnostico: viene
    # reimpostato a 1,36 a ogni iterazione di stabilita', quindi i pixel che convergono
    # presto non mostrano la riduzione. Per lo stress useremo i flussi (Step 2d).

    Rn = Sn_C + Sn_S + Ln_C + Ln_S

    # Traspirazione POTENZIALE della chioma: la stessa formula di Priestley-Taylor che TSEB usa
    # come ipotesi iniziale (Norman 1995, eq. 14), con alpha = 1,26 e la stessa frazione verde.
    #   LE_C,pot = 1,26 * f_g * s/(s+gamma) * Rn_C
    # s = pendenza della curva di saturazione del vapore, gamma = costante psicrometrica.
    f_g = mt.get("f_g", np.ones_like(lai))
    s_slope = met.calc_delta_vapor_pressure(mt["T_A_K"]) * 10            # kPa/K -> mb/K
    c_p = met.calc_c_p(mt["p"], mt["ea"])
    gamma = met.calc_psicr(c_p, mt["p"], met.calc_lambda(mt["T_A_K"]))
    Rn_C = Sn_C + Ln_C
    LE_C_pot = 1.26 * f_g * s_slope / (s_slope + gamma) * Rn_C
    # Rapporto di traspirazione: 1 = al potenziale, 0 = chioma che non traspira
    T_ratio = np.clip(LE_C / np.where(LE_C_pot > 1, LE_C_pot, np.nan), 0, 1)
    return {
        "flag": flag, "T_S": T_S, "T_C": T_C, "Rn": Rn, "G": G,
        "H": H_C + H_S, "LE": LE_C + LE_S, "LE_C": LE_C, "LE_S": LE_S,
        "h_C": h_C, "Sn": Sn_C + Sn_S, "f_g": f_g,
        "Rn_C": Rn_C, "LE_C_pot": LE_C_pot, "T_ratio": T_ratio,
    }


def daily_et(LE, S_inst, S_day, lam):
    """
    Da LE istantaneo (W/m2) a ET giornaliera (mm/giorno).
    Ipotesi: il rapporto LE / radiazione solare resta costante durante il giorno
    (metodo della "frazione solare", Cammalleri et al. 2014). Allora:
        LE_giorno (J/m2) = LE / S_inst * S_giorno
    e 1 mm di acqua su 1 m2 = 1 kg, che richiede lambda = 2,45 MJ per evaporare.
    """
    return LE / S_inst * S_day / lam


def sensitivity(Tr, lai, mt, sza, cfg, n=30000):
    """
    Quanto cambia l'ET se la temperatura dell'aria o l'altezza della chioma sono sbagliate?
    TSEB dipende dalla DIFFERENZA superficie - aria: 1 K di errore su Ta si traduce in decine
    di W/m2 di H (e quindi di LE). Rifacciamo il calcolo su un campione con valori alternativi.
    """
    ok = np.flatnonzero((np.isfinite(Tr) & np.isfinite(lai) & np.isfinite(mt["T_A_K"])).ravel())
    idx = np.random.default_rng(0).choice(ok, min(n, ok.size), replace=False)
    sub = {k: v.ravel()[idx] for k, v in mt.items()}
    L = np.maximum(lai.ravel()[idx], 0.01)
    bins = [0, 0.5, 1, 2, 3, 8]
    rows = []
    variants = {"altezza da LAI, senza f_g": dict(cfg, heights_above_d0=False),
                "altezza da coltura, senza f_g": cfg,
                "altezza da coltura, con f_g": cfg}
    for vname, vcfg in variants.items():
        for dT in cfg["sens_ta_offsets_K"]:
            # dT e' in piu' rispetto alla correzione gia' applicata (ta_offset_K)
            sb = dict(sub, T_A_K=sub["T_A_K"] + dT)
            if vname.startswith("altezza da LAI"):
                sb.pop("crop", None)
            if vname.endswith("senza f_g"):
                sb.pop("f_g", None)
            fl = run_tseb(Tr.ravel()[idx], L, sb, sza.ravel()[idx], vcfg)
            et = daily_et(fl["LE"], sb["S_dn"], sb["S_day"], float(cfg["lambda_J_kg"]))
            med = pd.Series(et).groupby(pd.cut(L, bins, include_lowest=True), observed=True).median()
            row = {"variante": vname, "dTa_extra_K": dT}
            row.update({f"LAI {iv.left:g}-{iv.right:g}": round(v, 2) for iv, v in med.items()})
            row["LE=0 su LAI>3 (%)"] = round(100 * np.mean(et[L > 3] < 0.05), 1)
            rows.append(row)
    return pd.DataFrame(rows).set_index(["variante", "dTa_extra_K"])


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def run_case(name, Tr, lai, mt, sza, cfg):
    """Esegue TSEB su tutti i pixel validi e rimette i risultati in 2D."""
    ok = np.isfinite(Tr) & np.isfinite(lai)
    for v in ("T_A_K", "u", "ea", "p", "S_dn", "L_dn", "S_day"):
        ok &= np.isfinite(mt[v])
    log(f"  [{name}] pixel validi: {ok.sum():,} su {ok.size:,} ({100 * ok.mean():.1f}%)")

    sub = {k: v[ok] for k, v in mt.items()}
    # LAI minimo 0,01: con LAI = 0 alcune formule dividono per zero (suolo nudo di fatto)
    fl = run_tseb(Tr[ok], np.maximum(lai[ok], 0.01), sub, sza[ok], cfg)
    fl["ET_day"] = daily_et(fl["LE"], sub["S_dn"], sub["S_day"], float(cfg["lambda_J_kg"]))
    fl["EF"] = fl["LE"] / (fl["Rn"] - fl["G"])           # frazione evaporativa

    maps = {}
    for k, v in fl.items():
        full = np.full(Tr.shape, np.nan)
        full[ok] = v
        maps[k] = full
    bad = maps["flag"] >= TSEB.F_OUT_LIMS            # 252-255: nessuna soluzione fisica
    for k in maps:
        if k != "flag":
            maps[k][bad] = np.nan
    log(f"  [{name}] pixel senza soluzione (flag >= 252): {int(bad.sum()):,}")
    return maps


def main():
    cfg_all = yaml.safe_load(CONFIG_FILE.read_text(encoding="utf-8"))
    cfg = cfg_all["tseb"]
    a = cfg_all["analysis_aoi"]
    bbox = (a["west"], a["south"], a["east"], a["north"])
    crs = f"EPSG:{cfg_all['download']['output_crs']}"
    raw = PROJECT_ROOT / cfg_all["paths"]["raw"]
    proc = PROJECT_ROOT / cfg_all["paths"]["processed"]
    fig_dir = PROJECT_ROOT / cfg_all["paths"]["figures"]
    tab_dir = PROJECT_ROOT / cfg_all["paths"]["tables"]
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    date = cfg_all["reference_date"]
    t_utc = f"{date} {cfg['overpass_utc']}"
    log(f"TSEB-PT - {a['name']} - {t_utc} UTC")

    # --- Superficie ---
    Tr = read_tif(proc / cfg["lst_file"], bbox, crs)
    lai_da = read_tif(proc / cfg["lai_file"], bbox, crs)
    lai_da = lai_da.rio.reproject_match(Tr)          # stessa griglia (gia' lo e', per sicurezza)
    log(f"Griglia: {Tr.shape[0]} x {Tr.shape[1]} pixel a 20 m")

    # --- Sole all'ora esatta del passaggio S3 (non quella di S2) ---
    lon, lat = pixel_lonlat(Tr)
    t = pd.Timestamp(t_utc)
    ftime = t.hour + t.minute / 60.0
    sza, saa = met.calc_sun_angles(lat, lon, 0.0, t.dayofyear, ftime)   # stdlon 0 = ora UTC
    log(f"Zenit solare medio: {np.nanmean(sza):.1f} gradi")

    # --- Meteo ---
    mt = meteo_at_overpass(proc / "era5land_filled.nc", lon, lat, t_utc)
    log("Meteo medio sull'AOI all'ora del passaggio:")
    log(f"  Ta {np.nanmean(mt['T_A_K']) - 273.15:.1f} C | ea {np.nanmean(mt['ea']):.1f} hPa | "
        f"vento {np.nanmean(mt['u']):.1f} m/s | p {np.nanmean(mt['p']):.0f} hPa")
    log(f"  SW {np.nanmean(mt['S_dn']):.0f} W/m2 | LW {np.nanmean(mt['L_dn']):.0f} W/m2 | "
        f"SW giornaliera {np.nanmean(mt['S_day']) / 1e6:.1f} MJ/m2")
    # Correzione esplicita della temperatura dell'aria (default 0). Le celle ERA5-Land dell'AOI
    # sono mare e vengono riempite dalla costa: l'aria rischia di essere troppo fresca.
    # Si cambia SOLO con una misura di stazione a supporto, e va dichiarato.
    ta_off = float(cfg.get("ta_offset_K", 0.0))
    mt["T_A_K"] = mt["T_A_K"] + ta_off
    log(f"  correzione Ta applicata: {ta_off:+.1f} K")

    # --- Mappa colturale (EUCROPMAP 2022): altezza per coltura e esclusione di acqua/edificato ---
    crop_file = PROJECT_ROOT / cfg_all["paths"]["aux"] / cfg["cropmap_file"]
    crop = xr.open_dataarray(crop_file, engine="rasterio").squeeze("band", drop=True)
    crop = crop.rio.reproject_match(Tr, resampling=Resampling.nearest).values.astype("float64")
    mt["crop"] = crop

    # --- Frazione verde (stoppie e residui secchi non traspirano) ---
    fapar = read_tif(proc / cfg["fapar_file"], bbox, crs).rio.reproject_match(Tr).values
    lai_tmp = np.maximum(read_tif(proc / cfg["lai_file"], bbox, crs).values, 0.01)
    mt["f_g"] = green_fraction(lai_tmp, fapar, sza)
    mt["f_g"][~np.isfinite(mt["f_g"])] = 1.0
    log(f"Frazione verde f_g: mediana {np.nanmedian(mt['f_g']):.2f}, "
        f"pixel con f_g < 0,3: {100 * np.mean(mt['f_g'] < 0.3):.1f}%")
    excluded = np.isin(crop, cfg["exclude_classes"])
    log(f"Mappa colturale: esclusi {100 * excluded.mean():.1f}% dei pixel "
        f"(classi {cfg['exclude_classes']}: nessun dato, edificato, acqua, zone umide)")

    # --- Caso principale: LAI di SNAP ---
    log("\nCaso principale (LAI SNAP):")
    Tr_v, lai_v = Tr.values.astype("float64"), lai_da.values.astype("float64")
    lai_v[excluded] = np.nan          # TSEB e' un modello per suolo e vegetazione, non per l'acqua
    maps = run_case("SNAP", Tr_v, lai_v, mt, sza, cfg)

    # --- Sensibilita': LAI di BIOPAR (piu' alto sui campi fitti) ---
    log("\nSensibilita' (LAI BIOPAR):")
    bp = open_nc(raw / "biopar_lai.nc")["LAI"].isel(time=0, drop=True)
    bp = as_raster(bp, crs).rio.reproject_match(Tr)
    bp_v = bp.values.astype("float64")
    bp_v[~np.isfinite(lai_v)] = np.nan               # stessi pixel del caso principale
    maps_bp = run_case("BIOPAR", Tr_v, bp_v, mt, sza, cfg)

    # --- Sensibilita' a temperatura dell'aria e altezza della chioma (campione di pixel) ---
    log("\nSensibilita' (30.000 pixel a caso): ET mediana (mm/g) per classe di LAI")
    sens = sensitivity(Tr_v, lai_v, mt, sza, cfg)
    sens.to_csv(tab_dir / "p2_04_sensitivity.csv")
    log(sens.to_string())

    # --- Salvataggio GeoTIFF ---
    for k in ("Rn", "H", "LE", "G", "ET_day", "EF", "T_C", "T_S", "flag",
              "LE_C", "LE_S", "LE_C_pot", "T_ratio", "f_g", "h_C"):
        da = Tr.copy(data=maps[k].astype("float32"))
        da.attrs, da.encoding, da.name = {}, {}, k
        as_raster(da, crs).rio.to_raster(proc / f"tseb_{k.lower()}.tif", compress="deflate")
    log(f"\nGeoTIFF salvati in {proc.relative_to(PROJECT_ROOT)}/tseb_*.tif")

    # --- Controlli di plausibilita' ---
    valid = np.isfinite(maps["LE"])
    closure = maps["Rn"] - maps["G"] - maps["H"] - maps["LE"]
    log("\nControlli:")
    log(f"  chiusura Rn - G - H - LE: max |residuo| {np.nanmax(np.abs(closure)):.2e} W/m2 "
        "(TSEB chiude il bilancio per costruzione)")
    flags, counts = np.unique(maps["flag"][valid], return_counts=True)
    log("  flag di TSEB (0 = soluzione con chioma al potenziale; 3 = alpha ridotto e "
        "evaporazione del suolo azzerata; 5 = LE totale azzerato, H = Rn - G): " +
        ", ".join(f"{int(f)}: {100 * c / valid.sum():.1f}%" for f, c in zip(flags, counts)))

    # Statistiche per classe di LAI (attese: suolo nudo ~1 mm/g, mais irrigato 5-7 mm/g)
    bins = [0, 0.5, 1, 2, 3, 4, 8]
    df = pd.DataFrame({k: maps[k][valid] for k in ("Rn", "H", "LE", "G", "ET_day", "EF", "LE_C", "LE_S", "f_g")})
    df["Tr_C"] = Tr_v[valid] - 273.15
    df["ET_day_biopar"] = maps_bp["ET_day"][valid]
    df["LAI"] = lai_v[valid]
    df["classe_LAI"] = pd.cut(df["LAI"], bins, include_lowest=True)
    tab = df.groupby("classe_LAI", observed=True).agg(
        pixel=("LAI", "size"), Tr_C=("Tr_C", "median"), Rn=("Rn", "median"),
        H=("H", "median"), LE=("LE", "median"), EF=("EF", "median"),
        ET_mm=("ET_day", "median"), ET_p10=("ET_day", lambda s: s.quantile(0.1)),
        ET_p90=("ET_day", lambda s: s.quantile(0.9)), ET_biopar_mm=("ET_day_biopar", "median"))
    tab.round(2).to_csv(tab_dir / "p2_04_tseb_by_lai.csv")

    # Statistiche per coltura (le classi principali)
    names = {int(k): v for k, v in cfg_all["cropmap"]["legend"].items()}
    df["coltura"] = pd.Series(crop[valid]).map(names).fillna("altro").to_numpy()
    by_crop = df.groupby("coltura").agg(
        ettari=("LAI", lambda s: s.size * 400 / 1e4), LAI=("LAI", "median"),
        Tr_C=("Tr_C", "median"), H=("H", "median"), LE=("LE", "median"),
        LE_C=("LE_C", "median"), LE_S=("LE_S", "median"), f_g=("f_g", "median"),
        ET_mm=("ET_day", "median"), ET_p10=("ET_day", lambda s: s.quantile(0.1)),
        ET_p90=("ET_day", lambda s: s.quantile(0.9)),
        LE0_pct=("ET_day", lambda s: 100 * (s < 0.05).mean()))
    by_crop = by_crop[by_crop["ettari"] > 50].sort_values("ettari", ascending=False).round(2)
    by_crop.to_csv(tab_dir / "p2_04_tseb_by_crop.csv")
    log("\nMediane per coltura (EUCROPMAP 2022, classi > 50 ha):")
    log(by_crop.to_string())
    log("\nMediane per classe di LAI:")
    log(tab.round(2).to_string())

    d = df["ET_day_biopar"] - df["ET_day"]
    log(f"\nET giornaliera AOI: mediana {df['ET_day'].median():.2f} mm "
        f"(p5 {df['ET_day'].quantile(.05):.2f}, p95 {df['ET_day'].quantile(.95):.2f})")
    log(f"Sensibilita' al LAI (BIOPAR - SNAP): mediana {d.median():+.2f} mm/g, "
        f"p5 {d.quantile(.05):+.2f}, p95 {d.quantile(.95):+.2f}")

    # --- Figura ---
    fig, axs = plt.subplots(2, 3, figsize=(16, 10.5), constrained_layout=True)
    ext = [Tr.x.min() / 1e3, Tr.x.max() / 1e3, Tr.y.min() / 1e3, Tr.y.max() / 1e3]
    panels = [
        (Tr_v - 273.15, "LST affinata (C)", "inferno", None),
        (lai_v, "LAI SNAP (m2/m2)", "Greens", (0, 5)),
        (maps["Rn"], "Radiazione netta Rn (W/m2)", "magma", None),
        (maps["H"], "Calore sensibile H (W/m2)", "Reds", (0, 400)),
        (maps["ET_day"], "ET giornaliera (mm/giorno)", "YlGnBu", (0, 8)),
    ]
    for ax, (arr, title, cmap, lim) in zip(axs.flat, panels):
        vmin, vmax = lim if lim else np.nanpercentile(arr, [2, 98])
        im = ax.imshow(arr, extent=ext, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.set_xlabel("UTM 33N x (km)")
        ax.set_ylabel("y (km)")
        fig.colorbar(im, ax=ax, shrink=0.8)
    ax = axs.flat[5]
    sample = df.sample(min(len(df), 40000), random_state=0)
    ax.hexbin(sample["LAI"], sample["ET_day"], gridsize=60, bins="log", cmap="viridis",
              extent=(0, 6, 0, 9))
    ax.plot(tab.index.map(lambda iv: iv.mid), tab["ET_mm"], "o-", color="w", mec="k",
            label="mediana SNAP")
    ax.plot(tab.index.map(lambda iv: iv.mid), tab["ET_biopar_mm"], "s--", color="orange",
            mec="k", label="mediana con LAI BIOPAR")
    ax.set_xlabel("LAI SNAP (m2/m2)")
    ax.set_ylabel("ET giornaliera (mm)")
    ax.set_title("ET contro LAI")
    ax.legend(loc="upper left")
    fig.suptitle(f"TSEB-PT a 20 m - Delta del Po - {date}", fontsize=14)
    fig.savefig(fig_dir / "p2_04_tseb.png", dpi=150)
    log(f"Figura: {(fig_dir / 'p2_04_tseb.png').relative_to(PROJECT_ROOT)}")

    (PROJECT_ROOT / "outputs" / "p2_04_summary.txt").write_text(
        "\n".join(LINES) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
