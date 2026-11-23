import sys, numpy as np, pandas as pd
f = sys.argv[1]
d = pd.read_excel(f, sheet_name="Liver Proteomics")
d = d[d["Gene name"].notna()].copy()
d["gene"] = d["Gene name"].astype(str).str.split().str[0]
d = d.drop_duplicates("gene")
L2, L4, S4S2, L4L2 = "large (2N) / small (2N)", "large (4N) / small (4N)", "4N (small) / 2N (small)", "4N (large) / 2N (large)"
d["ploidy_same_size"] = d[S4S2] - d[L2]          # small 4n vs large 2n: same size, 2x DNA per volume
d["ploidy_same_dna_conc"] = d[S4S2]               # small 4n vs small 2n: same DNA per volume, 4n twice as big
d["size_2n"] = d[L2]                               # large vs small 2n: genome diluted
d["size_4n"] = d[L4]
COMP = {"ploidy_same_dna_conc": "4n vs 2n, same DNA per volume (small 4n / small 2n)",
        "ploidy_same_dna_conc_large": "4n vs 2n, large cells (large 4n / large 2n)",
        "ploidy_same_size": "same size, 2x DNA per volume (small 4n / large 2n)",
        "size_2n": "2x bigger, same ploidy 2n (large / small)",
        "size_4n": "2x bigger, same ploidy 4n (large / small)"}
d["ploidy_same_dna_conc_large"] = d[L4L2]
d = d[d.UniqID >= 2].reset_index(drop=True)
SETS = {
 "Glutathione synthesis & recycling": "Gclc Gclm Gss Gsr Gpx1 Gpx4 Glrx Glrx2 Glrx3 Glrx5",
 "Glutathione S-transferases": "Gsta1 Gsta2 Gsta3 Gsta4 Gstm1 Gstm2 Gstm3 Gstm5 Gstm6 Gstm7 Gstp1 Gstp2 Gstk1 Gsto1 Gstt1 Gstt2 Gstt3 Gstz1 Mgst1 Mgst3",
 "Thioredoxin / peroxiredoxin": "Txn1 Txnrd1 Txn2 Txnrd2 Prdx1 Prdx2 Prdx3 Prdx4 Prdx5 Prdx6 Srxn1 Txnl1",
 "Superoxide dismutase & catalase": "Sod1 Sod2 Sod3 Cat",
 "NADPH supply": "G6pdx Pgd Tkt Taldo1 Me1 Me2 Idh1 Idh2 Nnt Mthfd1 Aldh1l1",
 "NRF2 targets (other)": "Nqo1 Hmox1 Ephx1 Cbr1 Cbr3 Akr7a5 Akr1b8 Ftl1 Fth1 Blvrb Ugdh",
 "H2O2 / ROS-producing enzymes": "Cyp2e1 Xdh Maoa Maob Ero1a Ero1lb Acox1 Acox2 Acox3 Hao1 Hao2 Dao Aox1 Nox4",
}
POS = {"CONTROL: histones (track DNA)": d.loc[d.Histone == True, "gene"].tolist(),
       "CONTROL: mitochondrial proteins": d.loc[d.Mitochondria == True, "gene"].tolist()}
present = set(d.gene)
sets = {k: [g for g in v.split() if g in present] for k, v in SETS.items()}
sets.update({k: v for k, v in POS.items()})
rng = np.random.default_rng(1)
# random sets matched on number of peptides (measurement precision)
d["pbin"] = pd.qcut(d.UniqID.rank(method="first"), 10, labels=False)
pool = {b: g.index.to_numpy() for b, g in d.groupby("pbin")}
idx = {g: i for i, g in enumerate(d.gene)}
rows = []
for name, genes in sets.items():
    ii = np.array([idx[g] for g in genes]); bins = d.pbin.to_numpy()[ii]
    null_ix = np.stack([np.array([rng.choice(pool[b]) for b in bins]) for _ in range(5000)])
    for col, lab in COMP.items():
        v = d[col].to_numpy(); obs = np.median(v[ii]); null = np.median(v[null_ix], axis=1)
        p = (1 + np.sum(np.abs(null - np.median(null)) >= abs(obs - np.median(null)))) / (len(null) + 1)
        rows.append({"set": name, "n": len(ii), "comparison": col, "median_log2": obs, "random_median": np.median(null), "p_vs_random": p})
res = pd.DataFrame(rows)
redox = ~res.set.str.startswith("CONTROL")
for col in COMP:
    m = redox & (res.comparison == col); p = res.loc[m, "p_vs_random"].to_numpy(); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]; qq = np.empty_like(p); qq[o] = np.minimum(q, 1)
    res.loc[m, "q_fdr"] = qq
pd.set_option("display.width", 250)
print("Proteins analysed:", len(d))
for k, v in sets.items():
    print(f"  {k}: {len(v)} found" + ("" if k.startswith("CONTROL") else f" ({', '.join(v)})"))
piv = res.pivot(index="set", columns="comparison", values="median_log2")[list(COMP)]
print("\nMedian log2 ratio per set:"); print(piv.round(3).to_string())
pv = res.pivot(index="set", columns="comparison", values="p_vs_random")[list(COMP)]
print("\np vs random protein sets:"); print(pv.round(4).to_string())
qv = res.pivot(index="set", columns="comparison", values="q_fdr")[list(COMP)]
print("\nq (FDR across redox sets):"); print(qv.round(3).to_string())
genes = sorted({g for k, v in sets.items() if not k.startswith("CONTROL") for g in v})
gt = d.set_index("gene").loc[genes, ["UniqID"] + list(COMP)]
# per-protein z vs all proteins
for col in COMP:
    gt[col + "_pct"] = d[col].rank(pct=True).to_numpy()[[idx[g] for g in genes]]
gt.to_csv("redox_enzymes_per_protein.csv")
res.to_csv("redox_sets.csv", index=False)
print("\nExtreme single enzymes (top/bottom 5% of all proteins) in ploidy_same_size:")
x = gt[(gt.ploidy_same_size_pct > .95) | (gt.ploidy_same_size_pct < .05)]
print(x[["UniqID", "ploidy_same_size", "ploidy_same_dna_conc", "size_2n"]].round(2).sort_values("ploidy_same_size").to_string())
print("\nExtreme in ploidy_same_dna_conc:")
x = gt[(gt.ploidy_same_dna_conc_pct > .95) | (gt.ploidy_same_dna_conc_pct < .05)]
print(x[["UniqID", "ploidy_same_dna_conc", "ploidy_same_dna_conc_large", "ploidy_same_size"]].round(2).sort_values("ploidy_same_dna_conc").to_string())
print("\nAll-protein spread (sd):", d[list(COMP)].std().round(3).to_dict())
print("Correlation size_2n vs ploidy_same_size across all proteins:", round(d.size_2n.corr(d.ploidy_same_size), 3))
