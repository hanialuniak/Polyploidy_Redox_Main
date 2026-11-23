import os, re, sys, numpy as np, pandas as pd
f = sys.argv[1]; out = sys.argv[2]
d = pd.read_excel(f, sheet_name="Liver Proteomics")
d = d[d["Gene name"].notna()].copy()
d["gene"] = d["Gene name"].astype(str).str.split().str[0]
d = d.drop_duplicates("gene")
L2, L4, S, L = "large (2N) / small (2N)", "large (4N) / small (4N)", "4N (small) / 2N (small)", "4N (large) / 2N (large)"
d["ploidy_small"] = d[S]; d["ploidy_large"] = d[L]; d["dna_conc"] = d[S] - d[L2]; d["size_2n"] = d[L2]; d["size_4n"] = d[L4]
SETS = {
 "Glutathione S-transferases": "Gsta1 Gsta2 Gsta3 Gsta4 Gstm1 Gstm2 Gstm3 Gstm5 Gstm6 Gstm7 Gstp1 Gstp2 Gstk1 Gsto1 Gstt1 Gstt2 Gstt3 Gstz1 Mgst1 Mgst3",
 "Glutathione synthesis & recycling": "Gclc Gclm Gss Gsr Gpx1 Gpx4 Glrx Glrx2 Glrx3 Glrx5",
 "Thioredoxin / peroxiredoxin": "Txn1 Txnrd1 Txn2 Txnrd2 Prdx1 Prdx2 Prdx3 Prdx4 Prdx5 Prdx6 Srxn1 Txnl1",
 "H2O2 / ROS-producing enzymes": "Cyp2e1 Xdh Maoa Maob Ero1a Ero1lb Acox1 Acox2 Acox3 Hao1 Hao2 Dao Aox1 Nox4",
 "NADPH supply": "G6pdx Pgd Tkt Taldo1 Me1 Me2 Idh1 Idh2 Nnt Mthfd1 Aldh1l1",
 "NRF2 targets (other)": "Nqo1 Hmox1 Ephx1 Cbr1 Cbr3 Akr7a5 Akr1b8 Ftl1 Fth1 Blvrb Ugdh",
 "Superoxide dismutase & catalase": "Sod1 Sod2 Sod3 Cat",
}
rng = np.random.default_rng(7)
def comp_label(r):
    if r.Mitochondria: return "mito"
    if r.ER: return "ER"
    if r.Nucleus and not r.Cytoplasm: return "nucleus"
    return "cyto"
d["comp"] = d.apply(comp_label, axis=1)

def set_test(dd, genes, col, n=5000, by_comp=False, stat=np.median):
    dd = dd.reset_index(drop=True)
    pos = {g: i for i, g in enumerate(dd.gene)}; ii = np.array([pos[g] for g in genes if g in pos])
    if len(ii) < 3: return np.nan, np.nan, len(ii)
    key = pd.qcut(dd.UniqID.rank(method="first"), 5, labels=False).astype(str)
    if by_comp: key = key + "|" + dd.comp
    pools = {k: g.index.to_numpy() for k, g in dd.groupby(key)}
    ks = key.to_numpy()[ii]; v = dd[col].to_numpy()
    null = np.array([stat(v[[rng.choice(pools[k]) for k in ks]]) for _ in range(n)])
    obs = stat(v[ii]); c = np.median(null)
    return obs, (1 + np.sum(np.abs(null - c) >= abs(obs - c))) / (n + 1), len(ii)

lines = []
def say(s=""): print(s); lines.append(s)
base = d[d.UniqID >= 2]
cols = ["ploidy_small", "ploidy_large", "dna_conc"]
say("TEST 1 - Robustness: stricter protein filter (>=5 peptides), mean instead of median, compartment-matched random sets")
rows = []
for name, gs in SETS.items():
    gl = gs.split()
    for col in cols:
        r = {"set": name, "comparison": col}
        r["base_med"], r["base_p"], r["n"] = set_test(base, gl, col)
        r["pep5_med"], r["pep5_p"], r["n5"] = set_test(d[d.UniqID >= 5], gl, col)
        r["mean"], r["mean_p"], _ = set_test(base, gl, col, stat=np.mean)
        r["comp_med"], r["comp_p"], _ = set_test(base, gl, col, by_comp=True)
        rows.append(r)
t1 = pd.DataFrame(rows); pd.set_option("display.width", 260)
for col in cols:
    say(f"\n  {col}"); say(t1[t1.comparison == col].drop(columns="comparison").set_index("set").round(3).to_string())
t1.to_csv(os.path.join(out, "test1_robustness.csv"), index=False)

say("\nTEST 2 - Leave-one-out: is a set driven by one protein? (ploidy_small)")
for name in ["Glutathione S-transferases", "Glutathione synthesis & recycling", "Thioredoxin / peroxiredoxin", "H2O2 / ROS-producing enzymes"]:
    gl = [g for g in SETS[name].split() if g in set(base.gene)]
    meds = []
    for g in gl:
        m, p, _ = set_test(base, [x for x in gl if x != g], "ploidy_small", n=2000); meds.append((g, m, p))
    worst = max(meds, key=lambda x: x[2])
    say(f"  {name}: medians without each protein {min(m for _, m, _ in meds):+.2f} to {max(m for _, m, _ in meds):+.2f}; "
        f"largest p {worst[2]:.3f} (when {worst[0]} is left out)")
    frac = (base.set_index("gene").loc[gl, "ploidy_small"] < 0).mean()
    say(f"     {int(round(frac * len(gl)))} of {len(gl)} proteins lower in 4n (all proteins: {(base.ploidy_small < 0).mean():.0%})")

say("\nTEST 3 - Remove what cell size alone predicts. Small 4n cells are 2x bigger than small 2n, so take out each protein's")
say("size scaling (mean of the 2n and 4n size ratios) and test what remains.")
size = (base.size_2n + base.size_4n) / 2
b = np.polyfit(size, base.ploidy_small, 1)[0]
base = base.assign(ploidy_small_resid=base.ploidy_small - b * size)
say(f"  Across all proteins, ploidy_small = {b:+.3f} x size scaling + rest")
for name, gs in SETS.items():
    m, p, n = set_test(base, gs.split(), "ploidy_small_resid")
    say(f"  {name:36s} n={n:2d}  after removing size: {m:+.3f}  p={p:.4f}")

say("\nTEST 4 - Reproducibility: do small-cell and large-cell ploidy comparisons agree, protein by protein?")
say(f"  All proteins: r = {base.ploidy_small.corr(base.ploidy_large):.2f}; proteins with >=10 peptides: "
    f"r = {base[base.UniqID >= 10].ploidy_small.corr(base[base.UniqID >= 10].ploidy_large):.2f}")
allg = sorted({g for v in SETS.values() for g in v.split()} & set(base.gene))
x = base.set_index("gene").loc[allg]
say(f"  Redox enzymes ({len(allg)}): r = {x.ploidy_small.corr(x.ploidy_large):.2f}; lower in 4n in both: "
    f"{int(((x.ploidy_small < 0) & (x.ploidy_large < 0)).sum())}, higher in both: {int(((x.ploidy_small > 0) & (x.ploidy_large > 0)).sum())}")

say("\nTEST 5 - Unbiased screen: every GO biological process with >=10 proteins, ranked by ploidy_small (4n vs 2n, same DNA per volume)")
go = {}
for g, s in zip(base.gene, base["Gene ontology (biological process)"].fillna("")):
    for term in re.findall(r"([^;]+?) \[GO:\d+\]", s):
        go.setdefault(term.strip(), []).append(g)
go = {k: v for k, v in go.items() if 10 <= len(v) <= 300}
z = (base.ploidy_small - base.ploidy_small.mean()) / base.ploidy_small.std()
zs = pd.Series(z.to_numpy(), index=base.gene)
from scipy.stats import mannwhitneyu
rows = []
for term, gs in go.items():
    inn = zs.loc[gs]; outn = zs.drop(gs)
    p = mannwhitneyu(inn, outn).pvalue
    rows.append({"GO term": term, "n": len(gs), "median_log2": base.set_index("gene").loc[gs, "ploidy_small"].median(), "p": p})
g5 = pd.DataFrame(rows).sort_values("p")
p = g5.p.to_numpy(); o = np.argsort(p); q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
g5["q"] = np.minimum(q, 1)[np.argsort(o)] if False else np.minimum(q, 1)
g5 = g5.reset_index(drop=True); g5["rank"] = np.arange(1, len(g5) + 1)
say(f"  {len(g5)} GO terms tested; {int((g5.q < 0.05).sum())} with q < 0.05")
say(g5.head(25)[["rank", "GO term", "n", "median_log2", "q"]].round(4).to_string(index=False))
red = g5[g5["GO term"].str.contains("glutathione|oxidative|oxidant|redox|hydrogen peroxide|reactive oxygen|xenobiotic|NADP", case=False)]
say("\n  Redox-related GO terms and their rank:")
say(red[["rank", "GO term", "n", "median_log2", "q"]].round(4).to_string(index=False))
g5.to_csv(os.path.join(out, "test5_go_screen_ploidy_small.csv"), index=False)
open(os.path.join(out, "follow_up_results.txt"), "w").write("\n".join(lines))
