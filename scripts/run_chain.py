"""
run_chain.py - Esegue tutta la catena della Fase 1 per una data, nell'ordine giusto.

    python scripts/run_chain.py --date 2022-08-11
    python scripts/run_chain.py --date 2022-08-11 --from p1_03   # riparte da uno step
    python scripts/run_chain.py --date 2019-08-04 2020-06-19 --solo-sharpening   # piu' date

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
    parser.add_argument("--date", required=True, nargs="+", help="una o piu' date AAAA-MM-GG")
    parser.add_argument("--from", dest="start", default=None,
                        help="prefisso dello step da cui partire, es. p1_03")
    parser.add_argument("--to", dest="stop", default=None,
                        help="prefisso dell'ultimo step da eseguire, es. p1_05")
    parser.add_argument("--altezza", choices=["coltura", "lai"], default=None,
                        help="altezza della chioma in TSEB (default: quella del config)")
    parser.add_argument("--solo-sharpening", action="store_true",
                        help="solo download S2/S3/Landsat, sharpening e validazione (step p1_)")
    args = parser.parse_args()

    def index_of(prefix):
        idx = next((i for i, s in enumerate(STEPS) if s.startswith(prefix)), None)
        if idx is None:
            sys.exit(f"Step {prefix} sconosciuto. Possibili: {', '.join(STEPS)}")
        return idx

    first = index_of(args.start) if args.start else 0
    last = index_of(args.stop) if args.stop else len(STEPS) - 1
    steps = STEPS[first:last + 1]
    if args.solo_sharpening:
        steps = [s for s in steps if s.startswith("p1_")]

    failed = []
    for date in args.date:
        ok = True
        for step in steps:
            print(f"\n{'#' * 70}\n# {step}  ({date})\n{'#' * 70}", flush=True)
            cmd = [sys.executable, str(SCRIPTS / f"{step}.py"), "--date", date]
            if args.altezza:
                cmd += ["--altezza", args.altezza]
            result = subprocess.run(cmd)
            if result.returncode != 0:
                if step in OPTIONAL:
                    print(f"\n{step} non riuscito: proseguo senza (vedi messaggio sopra).")
                    continue
                print(f"\nCatena interrotta a {step} per il {date}. Per ripartire: "
                      f"python scripts/run_chain.py --date {date} --from {step[:5]}")
                failed.append((date, step))
                ok = False
                break                    # passa alla data successiva invece di fermare tutto
        if ok:
            print(f"\nCatena completata per il {date}.")
    if failed:
        print("\nDate non completate: " + ", ".join(f"{d} ({s})" for d, s in failed))
        sys.exit(1)


if __name__ == "__main__":
    main()
