"""
p1_04_models.py - Confronta i modelli che prevedono la LST (1 km) dalle bande Sentinel-2.

Domanda: quanto bene le bande Sentinel-2 spiegano la temperatura misurata da Sentinel-3?
Il modello migliore sara' quello usato per lo sharpening a 20 m.

Modelli, dal piu' semplice al piu' complesso:
  1. "Media"            : prevede sempre la LST media (riferimento minimo da battere)
  2. "Lineare NDVI"     : LST = a + b * NDVI
  3. "Lineare tutte"    : LST = a + b1*B02 + ... + b12*NDMI
  4. "Foresta casuale"  : centinaia di alberi di regressione, risposta media

Validazione incrociata SPAZIALE: l'area e' divisa in blocchi di 5 x 5 km; ogni volta
il modello si addestra su alcuni blocchi e viene verificato su blocchi DIVERSI.
Per confronto calcoliamo anche la validazione incrociata CASUALE (pixel mescolati),
che di solito da' risultati troppo ottimisti: e' una lezione importante.

Output:
  outputs/tables/p1_04_model_scores.csv
  outputs/figures/p1_04_models.png
  outputs/p1_04_summary.txt

Come eseguirlo:
    python scripts/p1_04_models.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PROJECT_ROOT / "config" / "phase1_sharpening.yaml"

LINES = []


def log(text=""):
    print(text)
    LINES.append(text)


def build_models(cfg_m, n_predictors):
    """Dizionario nome -> (modello, colonne da usare)."""
    rf = RandomForestRegressor(
        n_estimators=cfg_m["rf_trees"],
        min_samples_leaf=cfg_m["rf_min_samples_leaf"],
        random_state=cfg_m["random_seed"],
        n_jobs=-1,               # usa tutti i core del processore
    )
    all_cols = cfg_m["predictors"]
    return {
        "Media": (DummyRegressor(strategy="mean"), all_cols),
        "Lineare NDVI": (LinearRegression(), ["NDVI"]),
        # StandardScaler porta ogni predittore a media 0 e deviazione 1:
        # i coefficienti diventano confrontabili tra loro
        "Lineare tutte": (make_pipeline(StandardScaler(), LinearRegression()), all_cols),
        "Foresta casuale": (rf, all_cols),
    }


def cross_validate(model, X, y, splits):
    """Previsioni 'fuori campione': ogni pixel e' previsto da un modello che non l'ha mai visto."""
    pred = np.full(len(y), np.nan)
    for train_idx, test_idx in splits:
        model.fit(X.iloc[train_idx], y.iloc[train_idx])
        pred[test_idx] = model.predict(X.iloc[test_idx])
    return pred


def scores(y, pred):
    err = pred - y
    return {
        "RMSE_K": float(np.sqrt(np.mean(err ** 2))),  # errore tipico, penalizza gli errori grandi
        "MAE_K": float(np.mean(np.abs(err))),         # errore medio in valore assoluto
        "bias_K": float(np.mean(err)),                # tendenza a sovra/sottostimare
        "R2": float(1 - np.sum(err ** 2) / np.sum((y - y.mean()) ** 2)),  # quota di variabilita' spiegata
    }


def main():
    cfg = yaml.safe_load(open(CONFIG_FILE, encoding="utf-8"))
    m = cfg["models"]
    table = pd.read_csv(PROJECT_ROOT / cfg["paths"]["processed"] / "training_table.csv")
    X = table[m["predictors"]]
    y = table["LST_K"]
    log(f"Esempi: {len(table)}  |  Predittori: {len(m['predictors'])}  |  "
        f"LST: media {y.mean() - 273.15:.1f} °C, deviazione std {y.std():.2f} K")

    # --- Blocchi spaziali: ogni pixel riceve l'etichetta del blocco 5x5 km in cui cade ---
    block_m = m["cv_block_km"] * 1000
    blocks = (np.floor(table["x"] / block_m).astype(int).astype(str) + "_" +
              np.floor(table["y"] / block_m).astype(int).astype(str))
    log(f"Blocchi spaziali da {m['cv_block_km']} km: {blocks.nunique()}")

    spatial_splits = list(GroupKFold(n_splits=m["cv_folds"]).split(X, y, groups=blocks))
    random_splits = list(KFold(n_splits=m["cv_folds"], shuffle=True,
                               random_state=m["random_seed"]).split(X, y))

    rows, oof = [], {}
    for name, (model, cols) in build_models(m, len(m["predictors"])).items():
        for cv_name, splits in [("spaziale", spatial_splits), ("casuale", random_splits)]:
            pred = cross_validate(model, X[cols], y, splits)
            rows.append({"modello": name, "validazione": cv_name, **scores(y.values, pred)})
            if cv_name == "spaziale":
                oof[name] = pred

    results = pd.DataFrame(rows)
    log("\nRisultati (previsioni su dati mai visti dal modello):")
    log(results.round(3).to_string(index=False))

    # --- Quali predittori usa la foresta casuale? (importanza per permutazione) ---
    # Si "mescola" una colonna alla volta e si misura quanto peggiora l'errore:
    # piu' peggiora, piu' quel predittore e' importante. Calcolata sui blocchi esclusi.
    rf_model, rf_cols = build_models(m, len(m["predictors"]))["Foresta casuale"]
    imp = np.zeros(len(rf_cols))
    for train_idx, test_idx in spatial_splits:
        rf_model.fit(X.iloc[train_idx], y.iloc[train_idx])
        r = permutation_importance(rf_model, X.iloc[test_idx], y.iloc[test_idx],
                                   scoring="neg_root_mean_squared_error",
                                   n_repeats=10, random_state=m["random_seed"], n_jobs=-1)
        imp += r.importances_mean / len(spatial_splits)
    importance = pd.Series(imp, index=rf_cols).sort_values(ascending=False)
    log("\nImportanza dei predittori (aumento dell'RMSE in K se il predittore viene 'mescolato'):")
    for name, v in importance.items():
        log(f"  {name:5s} {v:6.3f}")

    # --- Salvataggi ---
    tab_dir = PROJECT_ROOT / cfg["paths"]["tables"]
    tab_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(tab_dir / "p1_04_model_scores.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), constrained_layout=True)
    lims = [y.min() - 273.15 - 0.5, y.max() - 273.15 + 0.5]
    for ax, name in zip(axes[:2], ["Lineare tutte", "Foresta casuale"]):
        s = results[(results.modello == name) & (results.validazione == "spaziale")].iloc[0]
        ax.scatter(y - 273.15, oof[name] - 273.15, s=10, alpha=0.6)
        ax.plot(lims, lims, "k--", lw=1, label="1:1")
        ax.set_xlim(lims); ax.set_ylim(lims)
        ax.set_xlabel("LST misurata da Sentinel-3 (°C)")
        ax.set_ylabel("LST prevista (°C)")
        ax.set_title(f"{name} - validazione spaziale\nRMSE {s.RMSE_K:.2f} K, R² {s.R2:.2f}")
        ax.grid(alpha=0.3); ax.legend()
    importance[::-1].plot.barh(ax=axes[2], color="tab:orange")
    axes[2].set_xlabel("Aumento dell'RMSE (K) se il predittore viene mescolato")
    axes[2].set_title("Importanza dei predittori (foresta casuale)")
    axes[2].grid(alpha=0.3, axis="x")
    fig.suptitle(f"Previsione della LST a 1 km dalle bande Sentinel-2 - {cfg['reference_date']}")
    fig_dir = PROJECT_ROOT / cfg["paths"]["figures"]
    fig.savefig(fig_dir / "p1_04_models.png", dpi=150)
    plt.close(fig)
    log(f"\nFigura: {cfg['paths']['figures']}/p1_04_models.png")
    log(f"Tabella: {cfg['paths']['tables']}/p1_04_model_scores.csv")
    (PROJECT_ROOT / "outputs" / "p1_04_summary.txt").write_text("\n".join(LINES), encoding="utf-8")


if __name__ == "__main__":
    main()
