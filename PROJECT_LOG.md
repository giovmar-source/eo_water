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
| 0 | 3 – Basi teoriche: bilancio energetico, TSEB, metodi di sharpening, scelta del metodo | fatto | 2026-11-05 |
| 1 | 1a – Ricerca date con passaggio comune S2 + S3 + Landsat (STAC): trovato 19/07/2022 | fatto | 2026-11-05 |
| 1 | 1b – Download dati del 19/07/2022 (S2 11 bande, S3 LST, Landsat 9 LC09_L2SP_192029_20220719) | fatto | 2026-11-05 |
| 1 | 1c – Sharpening su una data con validazione Landsat: tabella, modelli, PSF, mappa a 20 m | fatto | 2026-11-05 |
| 1 | 2 – Dalla LST a 20 m all'evapotraspirazione (TSEB) e allo stress idrico sul 19/07/2022 | da iniziare | |

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
| 2026-11-05 | Modello di ET: bilancio energetico a due sorgenti (TSEB) | Separa suolo e vegetazione usando Sentinel-2; standard ESA Sen-ET |
| 2026-11-05 | Sharpening del termico con predittori Sentinel-2 come metodo principale; fusione temporale tipo STARFM con Landsat solo come confronto | Disponibile ogni giorno limpido; la LST varia troppo da un giorno all'altro per la fusione temporale |
| 2026-11-05 | Predittori: tutte le bande S2 a 20 m (incl. red-edge B05-B07 e SWIR B11-B12) + quota | Lo SWIR vede il contenuto d'acqua, riduce (non elimina) il ritardo dell'ottico sullo stress |
| 2026-11-05 | Validazione dello sharpening con Landsat 8/9 (termico 100 m) nei giorni con passaggio comune; test su campi con NDVI simile e LST diversa | Riferimento indipendente; verifica del rischio di circolarità |
| 2026-11-05 | Uso del suolo: Sentinel-2 come mappa dinamica per lo sharpening; mappe a livello di campo (LPIS/SIGPAC, EuroCrops, WorldCereal, da verificare) per maschere e stratificazione | CORINE troppo grossolano (25 ha) e poco aggiornato per i singoli campi |
| 2026-11-05 | Sentinel-3C (lanciato 15/09/2026, in collaudo) e LSTM (dal 2028) citati come evoluzione operativa | Stesso SLSTR a ~1 km: più passaggi, non più risoluzione; non coprono il 2018-2025 |
| 2026-11-05 | Area di addestramento ~40x40 km separata dall'area di analisi; un modello per data; validazione incrociata spaziale; dal modello lineare alla random forest | ~220 pixel da 1 km nell'area di analisi sono pochi per ~10 predittori |
| 2026-11-05 | Data di riferimento per il prototipo: 19/07/2022 (S2 ~0% nuvole; S3 09:19 e 09:58 UTC; Landsat 9 0,6% nuvole alle 09:58 UTC) | Unico giorno con i tre satelliti limpidi nell'estate 2022; S3 e Landsat quasi simultanei |
| 2026-11-05 | Ricerche STAC con pagine da 100, pausa e retry con backoff su 429 | Il catalogo CDSE ha un limite di richieste |
| 2026-11-05 | Validazione: separare distorsione tra sensori e capacità di ricostruire il disegno spaziale (bias, RMSE + correlazione, variabilità ricostruita/vera) | S3 vs Landsat 9 a 1 km il 19/07/2022: bias -1,2 °C, RMSE 2,3 °C, r 0,74; std interna ai pixel da 1 km 3,6 °C (mediana) |
| 2026-11-05 | Nodata sempre dichiarato (NaN) prima di riproiettare o aggregare | Landsat usa 0 = nessun dato: senza dichiararlo la media include 0 K |
| 2026-11-05 | PSF di Sentinel-3 modellata come gaussiana con sigma 500 m, applicata a S2 prima dell'aggregazione (`psf_sigma_m` in config) | Stimata in due modi indipendenti (correlazione S3-Landsat sfocato; R² dei modelli al variare di sigma) |
| 2026-11-05 | Tenere sia il modello lineare sia la foresta casuale fino alla validazione a scala fine con Landsat | A 1 km si equivalgono (R² 0,46 vs 0,43); conta la ricostruzione dei dettagli |
| 2026-11-05 | Modello principale per lo sharpening: foresta casuale + correzione dei residui (residuo interpolato bilineare, sfocato con la stessa PSF) | Validazione Landsat 9 a 100 m (68.316 celle): r 0,87, RMSE 2,65 K (2,25 K senza bias) contro r 0,50, RMSE 4,27 K di S3 interpolato; variabilità interna al km ricostruita all'89% (lineare 139%: estrapola fuori dall'intervallo di addestramento) |
| 2026-11-05 | Stress idrico valutato solo sui pixel vegetati (CWSI / TSEB); i campi caldi con NDMI basso a luglio sono suolo nudo (stoppie, arature), non stress | Caldo ≠ stress: il suolo nudo è caldo perché non traspira |

## File prodotti

- `environment.yml`, `check_env.py`, `.gitignore`, `PROJECT_LOG.md`, `.vscode/settings.json` (locale)
- `config/step02_first_look.yaml`, `scripts/step02_download.py`, `scripts/step02_first_look.py`
- `config/phase1_sharpening.yaml`, `scripts/p1_01_find_dates.py`, `outputs/tables/p1_coincident_dates.csv`, `scripts/p1_02_download.py`, `scripts/p1_03_training_table.py`, `scripts/p1_04_models.py`, `scripts/p1_05_sharpen.py`
- `GUIDA_PERSONALE.md` (locale, esclusa da Git): spiegazioni, codice commentato, domande da colloquio. Va aggiornata a ogni step.

## Criticità aperte

- **Risultato chiave 19/07/2022**: a 1 km le bande S2 spiegano solo ~35% della varianza della LST Sentinel-3 (R² spaziale 0,35; lineare NDVI = lineare tutte = foresta casuale), ma ~80% della LST Landsat 9 mediata a 1 km (R² 0,82). Il limite è nel dato Sentinel-3, non nella fisica. Nessuno spostamento di un pixel intero; la correlazione S3-Landsat cresce da 0,74 a 0,79 sfocando Landsat con una gaussiana di sigma ~500 m → l'impronta effettiva di Sentinel-3 è più ampia della cella da 1 km (PSF). Applicata la stessa sfocatura (sigma 500 m) a S2 prima dell'aggregazione: R² spaziale (lineare tutte) da 0,35 a 0,46; massimo a 500 m tra 0/250/500/750 m. Divario residuo con Landsat (0,82) attribuito a rumore della LST S3 e differenze tra sensori.

- Un solo giorno di validazione Landsat nell'estate 2022 sul Delta: ampliare (soglie, altri anni, Doñana) prima di trarre conclusioni generali.

- Rischio di circolarità dello sharpening ottico sullo stress precoce (stomi chiusi, foglie ancora verdi): da quantificare con Landsat.
- Ipotesi di frazione evaporativa costante per passare da ET istantanea a giornaliera: da dichiarare e, se possibile, verificare.

- Nome definitivo del progetto e del repository.
- Area InSAR da decidere in Fase 2 dopo verifica in letteratura.
- Disponibilità di torri di flusso (eddy covariance) per la validazione dell'ET: da verificare.
- Codice di EarthEmis da analizzare prima della Fase 3.
- Verificare sul Sentinel-3 SLSTR Land Handbook il significato di tutti i bit di `confidence_in`: il bit 14 sembra segnalare come nuvola anche pixel su lagune e zone umide in una giornata limpida (17/07/2022).
