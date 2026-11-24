#@title Cell 7 - Compartment-resolved screen, gene table and saving
COMPART = {
    ("nucleus", "H2O2 source (chromatin demethylation)"): "Su(var)3-3",
    ("nucleus", "thioredoxins and oxidative DNA repair"): "dhd TrxT Rrp1 Ogg1",
    ("ER", "H2O2 source"): "Ero1L",
    ("ER", "H2O2 removal"): "Jafrac2",
    ("ER", "folding load"): "Hsc70-3 Pdi ERp60 Calr Sec61alpha Sec61beta Xbp1 Ire1 PEK crc Srp54",
    ("cytosol", "H2O2 removal"): "Jafrac1 Prx2540-1 Prx2540-2 PHGPx",
    ("cytosol", "thioredoxin system"): "Trxr-1 Trx-2 Srx",
    ("cytosol", "glutathione supply"): "Gclc Gclm Gss1 Gss2 Gr1 Grx1 Grx1t",
    ("cytosol", "superoxide dismutase"): "Sod1 Sod3",
    ("mitochondria", "H2O2 removal"): "Prx3 Prx5 Trxr-2",
    ("mitochondria", "superoxide dismutase"): "Sod2",
    ("peroxisome", "H2O2 removal"): "Cat",
    ("plasma membrane", "H2O2 source"): "Duox Nox",
}
resolved, missing = {}, []
for (comp, role), v in COMPART.items():
    keep = []
    for g in v.split():
        r = resolve(g, PRESENT)
        keep.append(r) if r else missing.append(g)
    resolved[(comp, role)] = keep
if missing:
    print("Not detected in this dataset:", ", ".join(missing))

pct = fc.rank(pct=True)
rows, crows = [], []
for (comp, role), gs in resolved.items():
    for g in gs:
        rows.append({"gene": g, "compartment": comp, "role": role,
                     "gland_specific": float(fc.loc[g]), "pct_of_all_genes": float(pct.loc[g]) * 100,
                     "gland_late_vs_early": float(fc_sg.loc[g]), "other_late_vs_early": float(fc_oth.loc[g])})
    if not gs:
        continue
    obs = float(fc.loc[gs].median())
    null = np.array([fc.loc[[rng.choice(pools[bins[g]]) for g in gs]].median() for _ in range(N_RANDOM_SETS)])
    c = float(np.median(null))
    crows.append({"compartment": comp, "role": role, "genes": len(gs), "gland_specific": obs,
                  "random_median": c, "vs_random": obs - c,
                  "p_vs_random": ((1 + np.sum(np.abs(null - c) >= abs(obs - c))) / (N_RANDOM_SETS + 1)
                                  if len(gs) >= 3 else np.nan),
                  "median_pct": float(pct.loc[gs].median()) * 100})
gene_tab = pd.DataFrame(rows).set_index("gene")
comp_tab = pd.DataFrame(crows).sort_values(["compartment", "role"])
print("\nEvery redox gene by where it works. Positive = rises in the gland beyond the embryo-wide change.")
display(gene_tab.sort_values(["compartment", "role", "gland_specific"]).round(2))
display(comp_tab.round(3))

word = lambda p: "up" if p > 65 else "down" if p < 35 else "flat"
print("\n--- Source against sink, within each compartment (gland-specific) ---")
for comp in ["ER", "nucleus", "cytosol", "mitochondria"]:
    src = [g for (c, ro), gs in resolved.items() if c == comp and "source" in ro for g in gs]
    snk = [g for (c, ro), gs in resolved.items() if c == comp and ("removal" in ro or "thioredoxin" in ro) for g in gs]
    if not src or not snk:
        print(f"  {comp}: not enough detected genes on both sides ({len(src)} source, {len(snk)} removal)")
        continue
    s, k = float(fc.loc[src].median()), float(fc.loc[snk].median())
    sp, kp = float(pct.loc[src].median()) * 100, float(pct.loc[snk].median()) * 100
    print(f"  {comp}: sources {s:+.2f} ({sp:.0f}th pct, {word(sp)}) [{', '.join(src)}]; "
          f"removal {k:+.2f} ({kp:.0f}th pct, {word(kp)}) [{', '.join(snk)}]")
    if word(sp) == "up" and word(kp) != "up":
        print(f"      Production rises and removal does not follow: an imbalance local to the {comp}.")
    elif word(sp) == "up" and word(kp) == "up":
        print("      Both rise together, so capacity is tracking production.")
    elif word(sp) != "up" and word(kp) == "down":
        print(f"      Removal falls without a rise in production: less reducing capacity in the {comp}.")
    else:
        print("      Neither side moves clearly against the genome.")

fig, ax = plt.subplots(figsize=(8, 0.34 * len(comp_tab) + 1.4))
lab = [f"{c} - {r}" for c, r in zip(comp_tab.compartment, comp_tab.role)]
ax.barh(lab, comp_tab.vs_random, color=["#b2432f" if v < 0 else "#2f6db2" for v in comp_tab.vs_random])
ax.axvline(0, color="black", lw=0.8); ax.invert_yaxis()
ax.set_xlabel("gland-specific log2 change, minus random genes")
ax.set_title(f"{GSE}: salivary gland, stages 13-16 vs 10-12", loc="left", fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(GDIR, "chart_compartments.png"), dpi=150); plt.show()

res.to_csv(os.path.join(GDIR, "programmes.csv"), index=False)
gene_tab.to_csv(os.path.join(GDIR, "compartments_genes.csv"))
comp_tab.to_csv(os.path.join(GDIR, "compartments_summary.csv"), index=False)
L.to_csv(os.path.join(GDIR, "pseudobulk_log2cpm.csv"))
print("\nSaved to", GDIR)
print("""
Reading this. The gland-specific contrast already removes anything happening across the whole embryo between
these windows, which is the main confound. Two things it cannot remove: the gland is differentiating into a
secretory tissue AND entering endocycles over the same hours, so this comparison cannot separate those, exactly
as the larval salivary gland could not. What it can do is date the change - whether the ER folding load and any
redox shift are already underway in the early window or only appear in the late one.
""")
