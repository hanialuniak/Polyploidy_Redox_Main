"""Run the notebook's code cells on mock DepMap/CCLE files with a planted answer."""
import json, os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
scenario, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
os.environ["ROUTE_A_DATA"] = outdir
rng = np.random.default_rng(3 if scenario == "effect" else 4)
nb = json.load(open("RouteA_ploidy_redox.ipynb"))
codes = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
ns = {}
exec(codes[0], ns)  # settings
gene_cell = [c for c in codes if "PROGRAMMES = {" in c][0]
tmp = dict(ns); tmp["re"] = __import__("re"); exec(gene_cell, tmp)
prog = tmp["PROGRAMMES"]; prolif = tmp["PROLIFERATION"]
genes = list(dict.fromkeys([g for v in prog.values() for g in v] + prolif + ["GMNN", "CDT1"] + [f"BG{i}" for i in range(3000)]))
n, tissues = 1200, [f"T{i}" for i in range(29)] + ["CNS/Brain"]
ids = [f"ACH-{i:06d}" for i in range(n)]
tissue = rng.choice(tissues, n)
t_effect = {t: rng.normal(0, 0.6) for t in tissues}
tp = np.array([t_effect[t] for t in tissue])
ploidy = np.clip(2.6 + 0.8 * tp + rng.normal(0, 0.7, n), 1.5, 6)
wgd = (ploidy > 2.9).astype(float)
keap1 = rng.random(n) < 0.04
tp53 = rng.random(n) < 0.5
if scenario == "tp53":
    ploidy = ploidy + 0.6 * tp53
prolif_z = 0.5 * (ploidy - ploidy.mean()) + rng.normal(0, 1, n)
base = rng.normal(3, 1.5, len(genes))
X = base + rng.normal(0, 0.6, (n, len(genes)))
col = {g: i for i, g in enumerate(genes)}
for g in prolif: X[:, col[g]] += 0.8 * prolif_z
bg = [g for g in genes if g.startswith("BG")]
for g in rng.choice(bg, 600, replace=False):  # many ordinary genes also follow proliferation
    X[:, col[g]] += rng.uniform(0.2, 0.8) * prolif_z
for name, members in prog.items():
    act = rng.normal(0, 0.3, n) + 1.0 * tp
    if name.startswith("REDOX"):
        act += 0.4 * prolif_z
    if name in ("REDOX: NRF2 targets", "REDOX: thioredoxin/peroxiredoxin"):
        act += 2.0 * keap1
    if scenario == "effect" and name in ("REDOX: NRF2 targets", "REDOX: glutathione"):
        act += 0.35 * (ploidy - ploidy.mean())
    if scenario == "tp53" and name in ("REDOX: NRF2 targets", "REDOX: glutathione"):
        act += 0.35 * (ploidy - ploidy.mean()) * (~tp53) - 1.0 * tp53
    if name == "CONTROL: DNA damage":
        act -= 1.0 * tp53
    for g in members: X[:, col[g]] += act * 0.7
X[:, col["GMNN"]] += 0.6 * prolif_z
X[:, col["CDT1"]] += 0.6 * prolif_z
if scenario == "effect":
    X[:, col["GMNN"]] -= 0.5 * (ploidy - ploidy.mean())
# keap1 mutants also polyploid: confound
ploidy[keap1] += 0.8
pd.DataFrame({"ModelID": ids, "CCLEName": [f"LINE{i}_{t.replace('/','')}" for i, t in enumerate(tissue)],
              "OncotreeLineage": tissue}).to_csv(f"{outdir}/Model.csv", index=False)
sig = pd.DataFrame({"Ploidy": ploidy, "WGD": wgd, "Aneuploidy": 1.0}, index=ids); sig.loc[ids[:20], "Ploidy"] = np.nan
sig.to_csv(f"{outdir}/OmicsSignatures.csv")
expr = pd.DataFrame(np.maximum(X, 0), index=ids, columns=[f"{g} ({i+1})" for i, g in enumerate(genes)])
expr.to_csv(f"{outdir}/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
dam = pd.DataFrame(0.0, index=ids, columns=["KEAP1 (9817)", "CUL3 (8452)", "TP53 (7157)"]); dam.loc[np.array(ids)[keap1], "KEAP1 (9817)"] = 1
dam.loc[np.array(ids)[tp53], "TP53 (7157)"] = 1
dam.to_csv(f"{outdir}/OmicsSomaticMutationsMatrixDamaging.csv")
hot = pd.DataFrame(0.0, index=ids, columns=["NFE2L2 (4780)", "KRAS (3845)"]); hot.to_csv(f"{outdir}/OmicsSomaticMutationsMatrixHotspot.csv")
ko = pd.DataFrame(rng.normal(-0.3, 0.3, (n, 3)), index=ids, columns=["GMNN (51053)", "CDT1 (81620)", "TP53 (7157)"])
if scenario == "effect":
    ko["GMNN (51053)"] -= 0.15 * (ploidy - ploidy.mean())
ko.iloc[::9] = np.nan
ko.to_csv(f"{outdir}/CRISPRGeneEffect.csv")
sub = rng.choice(n, 420, replace=False)
prot = (X[sub] - X[sub].mean(0)) * 0.8 + rng.normal(0, 0.3, (420, len(genes)))
names = [f"LINE{i}_{tissue[i].replace('/','')}_TenPx{rng.integers(1,43):02d}" for i in sub]
names[1] = names[0].rsplit("_TenPx", 1)[0] + "_TenPx99"
prot_df = pd.DataFrame(prot.T, columns=names)
prot_df.insert(0, "Gene_Symbol", genes); prot_df.insert(0, "Protein_Id", range(len(genes)))
prot_df.iloc[::7, 5:40] = np.nan
prot_df.to_csv(f"{outdir}/protein_quant_current_normalized.csv.gz", index=False)
os.environ["ROUTE_A_DATA"] = outdir
ns = {"__name__": "__main__"}
import time
for i, c in enumerate(codes):
    t = time.time(); exec(c, ns); print(f"--- cell {i+1} ok in {time.time()-t:.1f}s")
