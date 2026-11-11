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
| 1 | 2a – Meteo ERA5-Land del 19/07/2022 (celle mare riempite con la terra più vicina) | fatto | 2026-11-05 |
| 1 | 2b – LAI, FCOVER, FAPAR a 20 m con la rete neurale SNAP (8 bande + angoli), confronto con BIOPAR (3 bande) | fatto | 2026-11-05 |
| 1 | 2c – TSEB-PT (pyTSEB 2.5.4) a 20 m, ET giornaliera, sensibilità a Ta, altezza chioma e LAI | primo run, Ta da verificare | 2026-11-05 |
| 1 | 2c-bis – Verifica ERA5 con stazione ARPAE GIRALDA (Codigoro): ERA5 −2,1 K al passaggio → ta_offset_K +2 | fatto | 2026-11-06 |
| 1 | 2c-ter – Altezza della chioma per coltura (EUCROPMAP 2022), quote sopra d_0, esclusi acqua/edificato: LE = 0 dal 7% al 3%, ET mediana 2,5 mm/g | fatto | 2026-11-06 |
| 1 | 2c-quater – Frazione verde f_g = FAPAR/FIPAR (stoppie: traspirazione da 82 a 4 W/m², ET frumento 2,1 → 1,6 mm/g) | fatto | 2026-11-06 |
| 1 | 2d – Stress idrico: CWSI TSEB (1 − LE_C/LE_C,pot) e CWSI contestuale per classi di LAI; accordo sul rango 0,84 per LAI ≥ 2, 0,43–0,58 per vegetazione rada | fatto (1 data) | 2026-11-06 |

| 1 | 3a – Ricerca passaggi ECOSTRESS estate 2022 (8–12 UTC): 7 sereni sul tassello 32TQQ; scelto 11/08/2022 09:12 UTC (S2 stesso giorno, 8% nuvole) | fatto | 2026-11-09 |

| 1 | 3b – Pipeline parametrica (--date), config_utils, run_chain, output per data; test di non regressione sul 19/07 superato | fatto | 2026-11-09 |

| 1 | 3c – Catena sull'11/08/2022 + confronto ECOSTRESS a 500 m: LST r 0,69 (bias +3,1 K, ECOSTRESS 49 min prima), ET r 0,61 (bias −0,8 mm/g), stress nessun accordo (definizioni diverse) | fatto | 2026-11-09 |

| 1 | 3d – Nuova ricerca date 2019–2025 (apr–set) con nuvole misurate sull'area (QA Landsat, SCL S2): 33 giorni utilizzabili; selezionate 9 date per la validazione multi-data | fatto | 2026-11-09 |

| 1 | 3e – Validazione multi-data: 10 date Landsat 2019–2025, r mediano 0,51 → 0,86, RMSE senza bias 3,10 → 1,99 K, migliorato 10/10 | fatto | 2026-11-09 |

| 1 | 4a – Contesto climatico ERA5-Land 1991–2025 + verifica pioggia stazione GIRALDA: 2022 il più secco a terra (184 mm gen–lug); ERA5 sbaglia il 2023 | fatto | 2026-11-10 |
| 1 | 4b – Variante altezza chioma (--altezza lai/coltura) con suffisso dei file | fatto | 2026-11-10 |
| 1 | 4c – Serie 2022 contro 2023–2025 (20 date, --altezza lai): 18/20 complete; ET 2022 −42% (−1,6 mm/g), EF −0,30, CWSI 0,35 contro 0,05–0,14 | fatto (01/08/2022 sostituito dal 29/07: S3 senza dati) | 2026-11-11 |
| 1 | 4d – Sensibilità all'altezza della chioma: ET 2022 −39% (LAI) contro −36% (mappa colturale); gradiente costiero in entrambe; bosco e frutteti NON robusto (−31% contro −41%) | fatto | 2026-11-11 |

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
| 2026-11-05 | ERA5-Land: celle classificate mare (gran parte del Delta) riempite con la cella di terra più vicina | ERA5-Land è solo terra; le celle costiere (Ta ~28 °C alle 10 UTC) rappresentano il Delta meglio dell'entroterra (~32 °C) |
| 2026-11-05 | Parametri biofisici con il processore SNAP v2.1 (Weiss & Baret) reimplementato in Python, coefficienti scaricati dal repository ufficiale s2tbx; BIOPAR openEO (3 bande) solo come confronto | L'NDVI satura (NDVI 0,85-0,95 corrisponde a LAI 2-6+); SNAP usa red-edge e SWIR; BIOPAR openEO usa solo B03, B04, B08 |
| 2026-11-05 | FCOVER usa il dominio di validità degli input di LAI | Il file non esiste per FCOVER in SNAP 2.1; quelli di LAI, FAPAR e LAI_Cab sono identici |

| 2026-11-05 | TSEB-PT di pyTSEB 2.5.4 (pip, con radiative-transfer-models e Py6S) | Implementazione di riferimento (Nieto, Guzinski), la stessa del plugin Sen-ET di ESA |
| 2026-11-05 | Sole calcolato all'ora del passaggio S3 (09:58), non dagli angoli S2 | La LST è quella delle 09:58 |
| 2026-11-05 | Radiazione ERA5 media 09–11 UTC (differenza dei cumulati / 7200 s) | Valore centrato sull'istante del passaggio |
| 2026-11-05 | ET giornaliera con rapporto LE / radiazione solare costante (Cammalleri 2014) | Usa la radiazione giornaliera ERA5, più robusto della frazione evaporativa costante |
| 2026-11-05 | Nessuna correzione di Ta (ta_offset_K = 0) finché non c'è una stazione | Una correzione senza misura sarebbe calibrare il risultato che si vuole ottenere |
| 2026-11-05 | Non usare l'alpha restituito da pyTSEB come indicatore di stress | Viene reimpostato a 1,36 a ogni iterazione: non mostra la riduzione per molti pixel |
| 2026-11-06 | Ta ERA5 + 2,0 K su tutta l'area | Stazione ARPAE GIRALDA: ERA5 troppo fresca di 2,1 K alle 10 UTC; uno scarto (non il valore della stazione) conserva il gradiente costa-entroterra |
| 2026-11-06 | Mappa colturale JRC EUCROPMAP 2022 (10 m, CC BY 4.0), ricampionata a 20 m con la moda | Distingue mais / riso / soia / cereali trebbiati, che hanno altezze diverse a parità di LAI; dato categorico → moda, mai media |
| 2026-11-06 | Altezza massima per coltura a metà luglio, scalata dal LAI; quote meteo = quota ERA5 + d_0 | I pixel "secchi" a LAI > 4 erano per il 58% soia (bassa) a cui il LAI dava 1,8 m; con mais 2,5 m, 2 m dal suolo sarebbe dentro la chioma |
| 2026-11-06 | Cartella `data/aux` rinominata `data/ancillary` | AUX è un nome riservato in Windows |
| 2026-11-06 | Frazione verde f_g = FAPAR / FIPAR (FIPAR dal LAI con Beer, k = 0,5/cos zenit) | Separa stoppie (0,02) da colture verdi anche rade (mais giovane 1,09); una soglia sull'NDVI le confonderebbe |
| 2026-11-06 | Stress da due indici indipendenti (TSEB assoluto, contestuale relativo senza Ta) | La mappa è credibile dove concordano sul rango; il contestuale non dà livelli assoluti |
| 2026-11-06 | Stress calcolato solo per LAI ≥ 0,5 e f_g ≥ 0,5, esclusi acqua ed edificato | Sotto è suolo nudo o stoppie |
| 2026-11-09 | Validazione con ECOSTRESS su una nuova data (11/08/2022), non sul 19/07 | Il 19/07 ECOSTRESS passa solo di notte e al tramonto; l'11/08 passa alle 09:12 UTC con S2 sereno lo stesso giorno |
| 2026-11-09 | Data come parametro, percorsi {date}, S3 'auto', Ta 'auto' da p2_05, validazione Landsat/ECOSTRESS | Serve per validazione, serie temporale e articolo; una sola fonte di verità |
| 2026-11-09 | Output spostati in outputs/<data>/ | Più date senza sovrascrivere |
| 2026-11-09 | Confronto con ECOSTRESS a 500 m, non a 70–100 m | Scene ECOSTRESS 2022 sfocate e con strisce; le anomalie sotto il km concordano (r 0,56 a 300 m), ma la validazione a livello di campo resta Landsat |
| 2026-11-09 | Nuvole valutate sull'area (ritagli delle maschere COG), non sulla scena intera | La copertura di scena (110–185 km) scartava giorni limpidi sul Delta: da 1 a 33 giorni utilizzabili |
| 2026-11-09 | Passaggio S3 'auto' più vicino all'ora di Landsat quando c'è Landsat | Simultaneità: la superficie si scalda di 2–3 °C/ora la mattina |
| 2026-11-09 | Landsat scaricato intero con requests e ritagliato in locale (non più letto via HTTP con GDAL) | Su Windows la lettura GDAL via HTTP lasciava il processo appeso e bloccava run_chain |
| 2026-11-10 | Anni di riferimento 2023, 2024, 2025 (non 2020–2021) | Scelti sulla pioggia misurata a terra; GIRALDA non ha dati prima del 2022; ERA5-Land indicava il 2023 come il più secco, la stazione come il più piovoso |
| 2026-11-10 | Analisi principale con altezza dal LAI per tutti (--altezza lai), mappa colturale 2022 come sensibilità | Le colture ruotano: per un confronto tra anni la coerenza conta più della precisione sul singolo campo |
| 2026-11-10 | WaPOR L1 (FAO, 300 m, decadale) come confronto indipendente dell'anomalia 2022 | Serie continua, modello diverso, riferimento noto |
| 2026-11-11 | 01/08/2022 sostituito con 29/07/2022 | L'unico passaggio S3 del 01/08 non aveva pixel di temperatura sull'area |
| 2026-11-11 | Non affermare che gli alberi hanno sofferto meno | Il risultato cambia da −31% a −41% con la variante dell'altezza |

## File prodotti

- `environment.yml`, `check_env.py`, `.gitignore`, `PROJECT_LOG.md`, `.vscode/settings.json` (locale)
- `config/step02_first_look.yaml`, `scripts/step02_download.py`, `scripts/step02_first_look.py`
- `config/phase1_sharpening.yaml`, `scripts/p1_01_find_dates.py`, `outputs/tables/p1_coincident_dates.csv`, `outputs/tables/p3_03_validation_all_dates.csv`, `outputs/figures/p3_03_validation_all_dates.png`, `scripts/p1_02_download.py`, `scripts/p1_03_training_table.py`, `scripts/p1_04_models.py`, `scripts/p1_05_sharpen.py`, `scripts/p2_01_download_era5.py`, `scripts/p2_02_download_biophys_inputs.py`, `scripts/p2_03_biophysical.py`, `scripts/p2_04_tseb.py`, `scripts/p2_05_station_check.py`, `scripts/p2_06_download_cropmap.py`, `scripts/p2_07_water_stress.py`, `scripts/p3_01_ecostress_search.py`, `scripts/config_utils.py`, `scripts/run_chain.py`, `scripts/p3_02_compare_ecostress.py`, `scripts/p3_03_validation_summary.py`, `scripts/p4_01_climate_context.py`, `scripts/p4_02_drought_comparison.py`, `outputs/tables/p4_01_climate_context.csv`, `outputs/tables/p1_coincident_dates.csv`, `outputs/tables/p3_03_validation_all_dates.csv`, `outputs/figures/p3_03_validation_all_dates.png`, `outputs/tables/p3_01_ecostress_passes.csv`, `outputs/tables/p2_07_stress_by_crop.csv`, `outputs/tables/p2_07_stress_by_lai.csv`, `outputs/tables/p2_04_tseb_by_crop.csv`, `data/ancillary/stations/arpae_giralda_20220719.csv` (escluso da Git), `outputs/tables/p2_04_tseb_by_lai.csv`, `outputs/tables/p2_04_sensitivity.csv`
- `GUIDA_PERSONALE.md` (locale, esclusa da Git): spiegazioni, codice commentato, domande da colloquio. Va aggiornata a ogni step.

## Criticità aperte

- **Gradiente costiero dell'anomalia 2022 — confermato da 3 prove indipendenti**: rapporto ET entroterra/costa 0,76 nel 2022 (tutte le date < 0,85) contro circa 1,00 negli altri anni; presente con entrambe le varianti dell'altezza; presente nel CWSI contestuale, che non usa Ta (da −0,06 sulla costa a +0,08 a 6–8 km). Resta da escludere una differenza di Ta costa–entroterra solo nel 2022 (controllo ERA5). Interpretazioni da verificare: falda vicino alle lagune; gestione delle prese irrigue dei consorzi.
- La stazione GIRALDA ha buchi isolati di un'ora: ora riempiti per interpolazione (p2_05).

- ~~Validazione a livello di campo su un giorno solo~~: risolto, 10 date Landsat 2019–2025 (r mediano 0,86). ECOSTRESS 2022 sul Delta è sfocato: va bene solo a scala di qualche centinaio di metri. L'ESI di ECOSTRESS (PT-JPL) non è confrontabile con il CWSI termico di TSEB.

- **Stress 19/07/2022**: su LAI > 3 TSEB dà CWSI mediano 0,5–0,6 (77–89% della superficie > 0,3). Coerente con le temperature (campi fitti da 32,6 a 38,6 °C) ma da validare. Ipotesi da verificare: campi non irrigati per siccità o cuneo salino (dati consortili di bonifica, ECOSTRESS). TSEB è poco sensibile allo stress della vegetazione rada (il calore in più finisce nel suolo). L'indice contestuale mostra un gradiente ovest-est probabilmente dovuto alla temperatura dell'aria.

- **TSEB dopo le correzioni (06/11)**: le colture fitte (LAI > 3) hanno ancora un'ET più bassa di quelle con LAI 1–3 (circa 2,7–3 contro 4 mm/g). Può essere stress reale (luglio 2022: ondata di calore, siccità, cuneo salino) o un errore residuo: serve una validazione indipendente (ECOSTRESS, torri di flusso). Stoppie risolte con f_g, ma l'evaporazione del suolo trebbiato resta ~100 W/m² (residuo di TSEB-PT): probabilmente sovrastimata su suolo secco, dipende dalla resistenza del suolo R_S. Il riso è quasi assente nell'area (15 ha), contrariamente all'ipotesi iniziale.

- **TSEB 19/07/2022, problema principale aperto**: con Ta ERA5 (28,5 °C) i campi più fitti (LAI > 3) risultano con ET = 0 (H ≈ Rn), mentre dovrebbero essere i più umidi. Landsat conferma che sono davvero caldi (36–37 °C; l'1% più fresco di tutto il suolo è 33 °C): il problema è la differenza superficie − aria. Le celle ERA5 dell'AOI sono mare, riempite da celle costiere più fresche. Sensibilità: +1 K di Ta ≈ +1 mm/g sulle colture fitte; altezza massima della chioma 1,0 m invece di 1,8 m riduce i pixel con LE = 0 su LAI > 3 dal 75% al 30%. **Aggiornamento 06/11**: con Ta +2 K (stazione GIRALDA) i pixel con LE = 0 scendono dal 21% al 7%; LAI 1–3 a ~3,7–4 mm/g, ma LAI > 4 ancora 0,4 mm/g. Restano da fare: (1) ~~Ta di stazione~~ fatto, solo una stazione (cercarne altre, anche ARPAV, verso la costa); (2) mappa colturale per l'altezza (riso ~1 m, mais ~2,5 m, stoppie di frumento); (3) eventualmente meteo a 100 m di quota come in Sen-ET.

- **LAI SNAP vs BIOPAR (19/07/2022)**: accordo complessivo buono (r 0,93, scarto quadratico medio 0,61, inferiore all'errore dichiarato ~0,9 di ciascun metodo). Previsione NON confermata: BIOPAR non è più basso sulla vegetazione fitta, è più alto (+0,45 per NDVI > 0,85). Dentro la classe NDVI > 0,85 il LAI SNAP segue l'indice red-edge NDRE (Spearman 0,88 contro 0,65 di BIOPAR), BIOPAR segue di più l'NDVI (0,66 contro 0,58). Senza misure di LAI a terra non si può dire quale sia più vicino al vero: SNAP resta la scelta per motivi fisici (usa red-edge e SWIR); BIOPAR verrà usato come analisi di sensibilità in TSEB.

- Verificare ERA5-Land (Ta, umidità, vento) con stazioni meteo del Delta (ARPAV / ARPAE) il 19/07/2022 alle 10 UTC: gradiente costa-entroterra di ~4,7 °C e Delta classificato "mare" in ERA5-Land.

- **Risultato chiave 19/07/2022**: a 1 km le bande S2 spiegano solo ~35% della varianza della LST Sentinel-3 (R² spaziale 0,35; lineare NDVI = lineare tutte = foresta casuale), ma ~80% della LST Landsat 9 mediata a 1 km (R² 0,82). Il limite è nel dato Sentinel-3, non nella fisica. Nessuno spostamento di un pixel intero; la correlazione S3-Landsat cresce da 0,74 a 0,79 sfocando Landsat con una gaussiana di sigma ~500 m → l'impronta effettiva di Sentinel-3 è più ampia della cella da 1 km (PSF). Applicata la stessa sfocatura (sigma 500 m) a S2 prima dell'aggregazione: R² spaziale (lineare tutte) da 0,35 a 0,46; massimo a 500 m tra 0/250/500/750 m. Divario residuo con Landsat (0,82) attribuito a rumore della LST S3 e differenze tra sensori.

- Un solo giorno di validazione Landsat nell'estate 2022 sul Delta: ampliare (soglie, altri anni, Doñana) prima di trarre conclusioni generali.

- Rischio di circolarità dello sharpening ottico sullo stress precoce (stomi chiusi, foglie ancora verdi): da quantificare con Landsat.
- Ipotesi di frazione evaporativa costante per passare da ET istantanea a giornaliera: da dichiarare e, se possibile, verificare.

- Pubblicazione: per un articolo servono più date/anni, una domanda scientifica (es. siccità e risalita del cuneo salino 2022) e un confronto con dati indipendenti (ECOSTRESS, torri di flusso se disponibili).
- Nome definitivo del progetto e del repository.
- Area InSAR da decidere in Fase 2 dopo verifica in letteratura.
- Disponibilità di torri di flusso (eddy covariance) per la validazione dell'ET: da verificare.
- Codice di EarthEmis da analizzare prima della Fase 3.
- Verificare sul Sentinel-3 SLSTR Land Handbook il significato di tutti i bit di `confidence_in`: il bit 14 sembra segnalare come nuvola anche pixel su lagune e zone umide in una giornata limpida (17/07/2022).
