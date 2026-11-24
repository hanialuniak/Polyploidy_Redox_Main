#@title Thioredoxin arm across the switch - reads saved results, no re-run needed
import os, glob
import numpy as np
import pandas as pd
try:
    from IPython.display import display
except ImportError:
    display = print
pd.set_option("display.max_columns", None, "display.width", 250)

DATA_DIR = "/content/drive/MyDrive/Colab Notebooks/route_e_ovary"  #@param {type:"string"}
DATA_DIR = os.environ.get("ROUTE_E_DATA", DATA_DIR)
if DATA_DIR.startswith("/content/drive") and not os.path.isdir("/content/drive"):
    from google.colab import drive
    drive.mount("/content/drive")

# The thioredoxin/peroxiredoxin arm, split by where each enzyme works.
ARM = {
    "Jafrac1": "cytosolic peroxiredoxin (Prx2)",
    "Jafrac2": "ER peroxiredoxin",
    "Prx3": "mitochondrial peroxiredoxin",
    "Prx5": "mitochondrial/peroxisomal peroxiredoxin",
    "Prx2540-1": "peroxiredoxin",
    "Prx2540-2": "peroxiredoxin",
    "Trxr-1": "thioredoxin reductase (recycles Trx, NADPH-dependent)",
    "Trxr-2": "thioredoxin reductase",
    "Trx-2": "thioredoxin (reduces oxidised Prx)",
    "TrxT": "thioredoxin, testis-biased",
    "dhd": "thioredoxin, female germline (deadhead)",
    "Srx": "sulfiredoxin (rescues hyperoxidised Prx)",
}
REF = {
    "Gclc": "glutathione synthesis, rate-limiting",
    "Gclm": "glutathione synthesis, modifier",
    "Gr1": "glutathione reductase",
    "Cat": "catalase",
    "Sod1": "superoxide dismutase, cytosolic",
    "Ero1L": "ER oxidase, makes H2O2 per disulfide",
    "Hsc70-3": "BiP, ER chaperone (folding load marker)",
}

found = sorted(glob.glob(os.path.join(DATA_DIR, "*", "genes.csv")))
if not found:
    raise SystemExit(f"No genes.csv found under {DATA_DIR}/<GSE>/. "
                     "Run the ovary notebook once per dataset first, or fix DATA_DIR above.")
print("Reading:", ", ".join(os.path.basename(os.path.dirname(p)) for p in found))

tabs, libcols = {}, {}
for p in found:
    gse = os.path.basename(os.path.dirname(p))
    d = pd.read_csv(p, index_col=0)
    tabs[gse] = d
    libcols[gse] = [c for c in d.columns if c.startswith("lib_")]

rows = []
for g, role in list(ARM.items()) + list(REF.items()):
    r = {"gene": g, "arm": "thioredoxin/peroxiredoxin" if g in ARM else "reference", "role": role}
    vals = []
    for gse, d in tabs.items():
        if g in d.index:
            v = float(d.loc[g, "log2_change"])
            r[gse] = v
            r[f"{gse}_pct"] = float(d.loc[g, "pct_of_all_genes"]) * 100
            vals.append(v)
            for c in libcols[gse]:
                r[f"{gse}:{c[4:]}"] = float(d.loc[g, c])
        else:
            r[gse] = np.nan
    r["agree"] = ("down in both" if len(vals) == len(tabs) and all(v < 0 for v in vals) else
                  "up in both" if len(vals) == len(tabs) and all(v > 0 for v in vals) else
                  "mixed" if len(vals) == len(tabs) else "not measured in all")
    rows.append(r)
tab = pd.DataFrame(rows).set_index("gene")

main = ["arm", "role"] + list(tabs) + [f"{g}_pct" for g in tabs] + ["agree"]
print("\nlog2 change, endocycling vs mitotic follicle cells. "
      "Negative = lower after the switch. pct = percentile among all genes.")
display(tab[main].round(2))

print("\nPer library (Rust atlas has three):")
percol = [c for c in tab.columns if ":" in c]
if percol:
    display(tab.loc[list(ARM), percol].round(2))

# Does the arm move as a group, and does it move more than the genome as a whole?
print("\n--- Is the thioredoxin arm falling as a group? ---")
for gse, d in tabs.items():
    arm = [g for g in ARM if g in d.index]
    v = d.loc[arm, "log2_change"]
    n_down = int((v < 0).sum())
    print(f"{gse}: {n_down} of {len(arm)} measured enzymes down, mean {v.mean():+.2f}, median {v.median():+.2f}")
    pp = os.path.join(DATA_DIR, gse, "programmes.csv")
    if os.path.exists(pp):
        pr = pd.read_csv(pp)
        row = pr[pr.programme.str.contains("thioredoxin", case=False, na=False)]
        if len(row):
            r0 = row.iloc[0]
            extra = ""
            if "random_median" in pr.columns and r0.random_median == r0.random_median:
                extra = f", random-gene median {r0.random_median:+.2f}"
            q = r0.get("q_fdr", np.nan)
            print(f"    as a programme: {r0.log2_change:+.2f}{extra}"
                  f"{'' if q != q else f', q = {q:.3f}'}")

print("""
How to read this. The capacity explanation needs the arm to fall as a group, not just Jafrac1.
  - Most enzymes down in both atlases, Trxr-1 and Trx-2 among them: reducing capacity is falling,
    so a brighter HyPer7 signal may reflect slower probe reduction rather than more H2O2.
  - Only Jafrac1 down, the thioredoxin reductases flat: capacity is largely intact and the
    Jafrac1 fall is a narrower change than it looked.
  - Jafrac2 is the one to watch separately: it is the ER peroxiredoxin, so it sits in the same
    compartment as Ero1L. Flat Jafrac2 against a rising translocon is a local imbalance;
    rising Jafrac2 says the cells are keeping up with the folding load.
Percentiles matter as much as the log2 values: a -0.5 that sits at the 3rd percentile of all
genes is a bigger statement than a -0.5 in a library where everything drifted down.
""")
