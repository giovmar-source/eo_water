"""
run_chain.py - Esegue tutta la catena della Fase 1 per una data, nell'ordine giusto.

    python scripts/run_chain.py --date 2022-08-11
    python scripts/run_chain.py --date 2022-08-11 --from p1_03   # riparte da uno step

Ordine:
  p1_02  download Sentinel-2, Sentinel-3 (openEO) e Landsat se passa quel giorno
  p2_01  meteo ERA5-Land
  p2_02  angoli Sentinel-2 e LAI BIOPAR (openEO)
  p1_03  tabella di addestramento a 1 km
  p1_04  confronto dei modelli
  p1_05  sharpening a 20 m e validazione (Landsat o ECOSTRESS)
  p2_03  LAI, FCOVER, FAPAR con la rete SNAP
  p2_05  confronto con la stazione (se c'e' il file della stazione per quella data)
  p2_04  TSEB ed ET giornaliera
  p2_07  stress idrico

Ogni step e' un processo Python separato: se uno fallisce, la catena si ferma e si vede
subito quale (e si puo' ripartire da li' con --from).
Gli step di download saltano i file gia' presenti, quindi rilanciare non riscarica nulla.
"""

import argparse
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
STEPS = ["p1_02_download", "p2_01_download_era5", "p2_02_download_biophys_inputs",
         "p1_03_training_table", "p1_04_models", "p1_05_sharpen", "p2_03_biophysical",
         "p2_05_station_check", "p2_04_tseb", "p2_07_water_stress"]
OPTIONAL = {"p2_05_station_check"}     # senza stazione si prosegue (TSEB usera' Ta senza correzione)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--date", required=True, help="AAAA-MM-GG")
    parser.add_argument("--from", dest="start", default=None,
                        help="prefisso dello step da cui ripartire, es. p1_03")
    args = parser.parse_args()

    steps = STEPS
    if args.start:
        idx = next((i for i, s in enumerate(STEPS) if s.startswith(args.start)), None)
        if idx is None:
            sys.exit(f"Step {args.start} sconosciuto. Possibili: {', '.join(STEPS)}")
        steps = STEPS[idx:]

    for step in steps:
        print(f"\n{'#' * 70}\n# {step}  ({args.date})\n{'#' * 70}", flush=True)
        result = subprocess.run([sys.executable, str(SCRIPTS / f"{step}.py"), "--date", args.date])
        if result.returncode != 0:
            if step in OPTIONAL:
                print(f"\n{step} non riuscito: proseguo senza (vedi messaggio sopra).")
                continue
            sys.exit(f"\nCatena interrotta a {step}. Correggi e riparti con: "
                     f"python scripts/run_chain.py --date {args.date} --from {step[:5]}")
    print(f"\nCatena completata per il {args.date}.")


if __name__ == "__main__":
    main()
