# Public examples

These are transformed public measurements. Original authors and dataset DOIs are listed in `../DATA_SOURCES.json`; the data retain CC BY 4.0.

| File / field | Meaning |
|---|---|
| `levo_design.csv`: `sample_id` | Deposited design-run number, 15 runs |
| `dose_g_L` | Hydroxyapatite photocatalyst dose in g/L |
| `concentration_ppm` | Initial LEVO concentration in ppm |
| `pH` | Nominal solution pH |
| `response` | Deposited percentage removal used for the design fit |
| `banana_design.csv`: `sample_id` | Deposited run number, nine settings |
| `time_h` | Drying duration in hours |
| `temperature_C` | Drying temperature in °C |
| `moisture_percent` | Deposited moisture response, percent |
| `colour_deltaE` | Deposited colour difference response |

LEVO rows come from the `BBD Model` sheet of `LEVO_RSM.xlsx`. The article-cited deposit and the initially located deposit contain a workbook with the same verified SHA256. Neither resolves the separate confirmation discrepancy; that sheet is excluded. Banana rows come from `Data and Design`; coded time and temperature are converted by time = 36 + 12A and temperature = 50 + 5B. Deposited responses are retained; the ambiguous printed raw-mass formula is not recalculated.

Preparation is reproducible with `../scripts/reproduce_examples.py`; see `../REPRODUCIBILITY.md`.
