# Source

Notebooks in `notebooks/` are **generated**. Edit the cell files here, then rebuild:

```bash
python3 build_all.py          # from the repository root: rebuilds every notebook
```

or for one route:

```bash
cd src/route_e_ovary && python3 build.py
```

## Layout of a route

```
cells/c01.py … c0N.py   one file per notebook cell, in order
build.py                assembles cells/ into the .ipynb
simulate.py             writes mock data in the same shape as the real files
run.py                  executes every cell against a folder of data
```

Cell files are plain Python. The first line is a `#@title` comment, which Colab turns into the
cell heading, and `#@param` comments become the input boxes in Cell 1.

## Testing a route

Each route is checked against simulated data before use. The pattern:

```bash
cd src/route_e_ovary
python3 simulate.py /tmp/mock_effect effect    # mock data with a planted effect
python3 simulate.py /tmp/mock_null   null      # same, with no effect

# routes that download from GEO need a local stand-in for it
(cd /tmp/mock_effect && python3 -m http.server 8800 &)
ROUTE_E_GEO=http://127.0.0.1:8800/geo python3 run.py /tmp/out_effect http://127.0.0.1:8800/geo
```

Two things must hold: the planted effect is recovered at roughly the right size, and the null
run reports nothing. Several routes also have mocks for input that is wrong in a specific way
— shuffled sample order, facility codes instead of sample names, a SuperSeries with no files at
series level, several versions of one dataset in the same folder — where the notebook must stop
rather than produce a result.

`run.py` reads environment variables so tests do not touch Google Drive or the real servers:

| Route | Variables |
| --- | --- |
| route_a_depmap | `ROUTE_A_DATA` |
| route_b_fly_scatac | `ROUTE_B_DATA` |
| route_b_fly_bulk | `ROUTE_BB_DATA`, `ROUTE_BB_GTF` |
| route_c_rat_liver | `ROUTE_CR_DATA`, `ROUTE_CR_GEO` |
| route_d_salivary | `ROUTE_D_DATA`, `ROUTE_D_GEO`, `ROUTE_D_GTF` |
| route_e_ovary | `ROUTE_E_DATA`, `ROUTE_E_GEO`, `ROUTE_E_EBI`, `ROUTE_E_GTF` |

## Rebuilding the PDF guides

`src/route_a_depmap/pdf/` and `src/route_b_fly_scatac/pdf/` build the step-by-step guides from
HTML using Playwright and Chromium:

```bash
pip install playwright && playwright install chromium
cd src/route_b_fly_scatac/pdf && python3 build_pdf.py
```

`build_all.py` runs these too and copies the output to `docs/`.
