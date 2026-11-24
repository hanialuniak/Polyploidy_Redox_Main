#@title Cell 8 - Compartment-resolved redox screen, including the nucleus
# The sets above pool enzymes that work in different places. Jafrac1 (cytosol) and Jafrac2 (ER)
# sit in one "thioredoxin / peroxiredoxin" set and can move in opposite directions, which hides
# the compartment where H2O2 is actually being handled. This cell splits every redox gene by
# where it works and adds the nucleus, which no earlier set covered.
#
# Nuclear H2O2 in flies: Su(var)3-3 is the fly LSD1 / KDM1A, an FAD-dependent H3K4 demethylase
# that releases one H2O2 per demethylation, on chromatin. Flies have no NOX4, so the other
# nuclear source described in mammals is not available here. Rrp1 is the fly APEX1, whose repair
# role is specific to oxidative DNA damage; dhd and TrxT are nuclear-acting thioredoxins.
# Note also that the perinuclear space is continuous with the ER lumen and H2O2 crosses the
# nuclear pore freely, so ER-generated H2O2 is topologically next door to the nucleoplasm.
pd.set_option("display.max_columns", None, "display.width", 250)
ALIASES.update({
    "Su(var)3-3": ["dLsd1", "Lsd1", "Su-var(3)3", "Suvar(3)3", "SuVar3-3", "CG17149"],
    "Rrp1": ["APEX1", "Rrp-1", "CG3178"],
    "Ogg1": ["dOgg1", "CG1795"],
})

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

present = set(L.index)
resolved, missing = {}, []
for (comp, role), v in COMPART.items():
    keep = []
    for g in v.split():
        r = resolve(g, present)
        if r:
            keep.append(r)
        else:
            missing.append(g)
    resolved[(comp, role)] = keep
if missing:
    print("Not detected in this dataset:", ", ".join(missing))

word = lambda p: "up" if p > 65 else "down" if p < 35 else "flat"
gene_rows, comp_rows = [], []
for t, r in CONTRASTS:
    name = f"{t} vs {r}"
    fc = gmean[t] - gmean[r]
    pct = fc.rank(pct=True)
    for (comp, role), gs in resolved.items():
        for g in gs:
            gene_rows.append({"contrast": name, "gene": g, "compartment": comp, "role": role,
                              "log2_change": float(fc.loc[g]),
                              "pct_of_all_genes": float(pct.loc[g]) * 100})
        if not gs:
            continue
        obs = float(fc.loc[gs].median())
        null = np.array([fc.loc[[rng.choice(pools[bins[g]]) for g in gs]].median()
                         for _ in range(N_RANDOM_SETS)])
        c = float(np.median(null))
        p = ((1 + np.sum(np.abs(null - c) >= abs(obs - c))) / (N_RANDOM_SETS + 1)
             if len(gs) >= 3 else np.nan)
        comp_rows.append({"contrast": name, "compartment": comp, "role": role, "genes": len(gs),
                          "log2_change": obs, "random_median": c, "vs_random": obs - c,
                          "p_vs_random": p, "median_pct": float(pct.loc[gs].median()) * 100})

gene_tab = pd.DataFrame(gene_rows)
comp_tab = pd.DataFrame(comp_rows)
print("\nEvery redox gene by where it works. Negative = lower in the test group.")
display(gene_tab.pivot_table(index=["compartment", "role", "gene"], columns="contrast",
                             values="log2_change").round(2))
print("\nBy compartment. p is only computed where a group has at least three detected genes; "
      "for smaller groups read the percentile.")
display(comp_tab.sort_values(["contrast", "compartment", "role"]).round(3))

print("\n--- Source against sink, within each compartment ---")
for t, r in CONTRASTS:
    fc = gmean[t] - gmean[r]
    pct = fc.rank(pct=True)
    print(f"\n{t} vs {r}:")
    for comp in ["ER", "nucleus", "cytosol", "mitochondria"]:
        src = [g for (c, ro), gs in resolved.items() if c == comp and "source" in ro for g in gs]
        snk = [g for (c, ro), gs in resolved.items()
               if c == comp and ("removal" in ro or "thioredoxin" in ro) for g in gs]
        if not src or not snk:
            print(f"  {comp}: not enough detected genes on both sides "
                  f"({len(src)} source, {len(snk)} removal) - cannot compare")
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

print("""
The salivary gland comparison is between different tissues, so a compartment that moves here may
reflect secretory identity rather than ploidy. The S2 contrast is the control: endocycles are
induced without changing cell identity, so a compartment that moves in the gland and NOT in S2
cells is workload, not ploidy. Read the two contrasts together, never the gland alone.
""")
gene_tab.to_csv(os.path.join(GDIR, "compartments_genes.csv"), index=False)
comp_tab.to_csv(os.path.join(GDIR, "compartments_summary.csv"), index=False)
print("Saved compartments_genes.csv and compartments_summary.csv to", GDIR)
