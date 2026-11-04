# PROJECT LOG

Registro sintetico del progetto. Va allegato (o letto dalla cartella) all'inizio di ogni nuova conversazione.

## Progetto

**Nome provvisorio:** l'acqua dall'alto e dal basso — consumo idrico agricolo, stress idrico e deformazione del suolo durante le siccità europee (Sentinel-1/2/3).

- Fase 1: evapotraspirazione e stress idrico da fusione Sentinel-2 / Sentinel-3 (Doñana, Delta del Po), 2018–2025.
- Fase 2: deformazione del suolo con InSAR Sentinel-1 su area ridotta, elaborazione propria e confronto con EGMS.
- Fase 3: report, Zenodo, serie LinkedIn; integrazione leggera in EarthEmis tramite file standard.

Repository: https://github.com/giovmar-source/eo_water (privato per ora). Cartella locale: `C:\EO\eo_water`.

## Stato

| Fase | Step | Stato | Data |
| --- | --- | --- | --- |
| 0 | 1 – Account, Miniforge, Git, clone del repository | fatto | 2026-11-04 |
| 0 | 1 – Ambiente `eo_water`, `.cdsapirc`, verifica (`check_env.py`: tutto OK) | fatto | 2026-11-04 |
| 0 | 1 – Primo push su GitHub (commit 0169a39) | fatto | 2026-11-04 |
| 0 | 2 – Primo accesso ai dati: Sentinel-3 LST e Sentinel-2 sulla stessa area (openEO): figura NDVI 20 m vs LST 1 km | fatto | 2026-11-04 |
| 0 | 3 – Basi teoriche: dalla temperatura superficiale all'evapotraspirazione; scelta del metodo di sharpening | da iniziare | |

## Decisioni

| Data | Decisione | Motivazione |
| --- | --- | --- |
| 2026-11-04 | Progetto D (ET + InSAR) a fasi | Allineato alle competenze richieste (SAR, fusione, Python); Fase 1 autonoma se la Fase 2 non regge |
| 2026-11-04 | Ambiente Conda `eo_water`, Python 3.11, canale conda-forge | Riproducibilità, librerie geospaziali stabili su Windows |
| 2026-11-04 | Pipeline separata da EarthEmis; integrazione solo in Fase 3 via COG/GeoJSON/CSV | Dimostra due competenze distinte; non rischia di rompere il sito |
| 2026-11-04 | Editor: VS Code (estensioni Python e Jupyter) | Preferenza dell'utente, diffuso in azienda |
| 2026-11-04 | Aprire i raster con `with rioxarray.open_rasterio(...)` e `.load()` | Su Windows un file aperto non si può cancellare o sovrascrivere |

| 2026-11-04 | Sentinel-3 LST: usare solo passaggi diurni (08-12 UTC) | Il prodotto contiene anche passaggi notturni (~20-21 UTC, LST 12-25 °C) inutili per lo stress idrico |
| 2026-11-04 | Copertura nuvolosa S3 calcolata sui soli pixel di terra | La LST non esiste sul mare: contarlo come nuvola sottostima le date limpide |
| 2026-11-04 | Sentinel-2: riflettanza = DN/10000, scartare valori <= 0 | Valori negativi su acqua facevano esplodere l'NDVI (min -435) |

## File prodotti

- `environment.yml`, `check_env.py`, `.gitignore`, `PROJECT_LOG.md`, `.vscode/settings.json` (locale)
- `config/step02_first_look.yaml`, `scripts/step02_download.py`, `scripts/step02_first_look.py`

## Criticità aperte

- Nome definitivo del progetto e del repository.
- Area InSAR da decidere in Fase 2 dopo verifica in letteratura.
- Disponibilità di torri di flusso (eddy covariance) per la validazione dell'ET: da verificare.
- Codice di EarthEmis da analizzare prima della Fase 3.
- Verificare sul Sentinel-3 SLSTR Land Handbook il significato di tutti i bit di `confidence_in`: il bit 14 sembra segnalare come nuvola anche pixel su lagune e zone umide in una giornata limpida (17/07/2022).
