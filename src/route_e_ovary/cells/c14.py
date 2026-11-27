#@title Cell 14 - Is any category enriched among the calibrated shared genes, or does it only look that way?
# Cell 13 produced a few hundred calibrated shared genes, and reading down that list by eye suggests two
# themes: chromatin, satellite-DNA and lamina proteins falling, and a few metal- and NADPH-handling genes
# rising. Spotting a theme in a ranked list by eye is exactly how this project has gone wrong before, so the
# themes are written out as explicit sets here and tested against expression-matched random sets of the same
# size, drawn from the genes that were actually eligible to be called shared. If a set is not enriched beyond
# those, the theme was an artefact of reading a list.
#
# Stated plainly for the record: the chromatin and metal sets below were assembled AFTER seeing the Cell 13
# output, so this is a test of a hypothesis generated from the same data. It can disqualify a theme but it
# cannot confirm one. Only a different dataset can do that.
N_RANDOM_SETS_14 = 2000      #@param {type:"integer"}
Q_CUT = 0.10                 #@param {type:"number"}

THEMES = {
    "Chromatin architecture & satellite DNA": "HmgD HmgZ D1 Dp1 Df31 His2Av HP1b HP1c HP4 HP5 HP6 Su(var)205 "
                                              "Su(var)3-9 HipHop ver mod(mdg4) CtBP jigr1 wac ball",
    "Nuclear lamina & envelope": "Lam LamC Mtor Nup154 Nup98-96 otefin bocks Ndc1 Koi klar",
    "Chromosome segregation": "polo aurB pav Incenp Det Sak Cks30A Bub3 Nnf1a Mad1 BubR1 cid Spc105R dgt2",
    "Replication machinery & repair": "RPA2 RfC4 RfC38 PCNA dUTPase mre11 rad50 nbs spn-A Rrp1 Ogg1 mus209 Fen1",
    "Metal handling & metallothionein": "MtnA MtnB MtnC MtnD Atox1 Ctr1A Ctr1B ATP7 Mvl Fer1HCH Fer2LCH Tsf1",
    "NADPH supply (re-test)": "Zw Pgd Men Men-b Idh Idh3a Nmdmc Tkt Taldo",
    "Peroxide removal (the original hypothesis)": "Jafrac1 Jafrac2 Prx3 Prx5 Prx2540-1 Prx2540-2 Trxr-1 Trx-2 "
                                                  "dhd Cat Sod1 Sod2 PHGPx",
}
THEMES = {k: [g for g in resolve_set(v) if g in de3.index] if "resolve_set" in dir()
          else [resolve(g, set(de3.index)) for g in v.split() if resolve(g, set(de3.index))]
          for k, v in THEMES.items()}

# The background is the set of genes that could have been called shared at all: tested on both sides, not
# flagged as ambient, and with a germline q-value. Anything else is not a fair comparison.
elig = de3[de3.germline_q.notna() & (de3.bleed_through_risk == "ok")
           & de3.lineage.astype(str).ne("not in germline data")].copy()
hit = set(elig.index[elig.shared_calibrated.fillna(False)])
print(f"Background: {len(elig):,} genes eligible to be called shared. Of those, {len(hit):,} were "
      f"({len(hit) / len(elig):.1%}). A set with no enrichment should land near that fraction.")

lvl14 = Lp.mean(axis=1).reindex(elig.index)
bins14 = pd.qcut(lvl14.rank(method="first"), 10, labels=False)
pool14 = {b: list(bins14.index[bins14 == b]) for b in range(10)}
rng14 = np.random.default_rng(11)
rows14 = []
for name, gs in THEMES.items():
    gs = [g for g in gs if g in elig.index]
    if len(gs) < 4:
        rows14.append({"theme": name, "genes_tested": len(gs), "note": "too few genes in this dataset"})
        continue
    obs_hit = sum(g in hit for g in gs)
    null = np.array([sum(rng14.choice(pool14[bins14[g]]) in hit for g in gs)
                     for _ in range(N_RANDOM_SETS_14)])
    rows14.append({"theme": name, "genes_tested": len(gs), "shared": obs_hit,
                   "fraction": obs_hit / len(gs),
                   "expected_matched": null.mean(),
                   "p_vs_matched": (1 + int((null >= obs_hit).sum())) / (N_RANDOM_SETS_14 + 1),
                   "members": ", ".join(g for g in gs if g in hit)[:70]})
en = pd.DataFrame(rows14)
if "p_vs_matched" in en:
    o14 = en.dropna(subset=["p_vs_matched"]).sort_values("p_vs_matched")
    en.loc[o14.index, "q_BH"] = np.minimum.accumulate(
        (o14.p_vs_matched.to_numpy() * len(o14) / np.arange(1, len(o14) + 1))[::-1])[::-1].clip(max=1.0)
display(en.round(3))

# A set whose expected hit count is near zero is "significantly enriched" on the strength of one gene, which
# is not a theme. Three members is the minimum for calling something a pattern; a single gene is reported as
# the single gene it is, and judged on its own q-value from Cell 13.
sig = en[(en.get("q_BH", pd.Series(dtype=float)).fillna(1) < 0.10) & (en.get("shared", 0).fillna(0) >= 3)]
thin = en[(en.get("q_BH", pd.Series(dtype=float)).fillna(1) < 0.10) & (en.get("shared", 0).fillna(0) < 3)]
for _, r in thin.iterrows():
    print(f"NOT a theme: {r.theme} has only {int(r.shared)} shared gene(s) ({r.members}). Against an expected "
          f"{r.expected_matched:.2f} that is formally significant, but one or two genes is not a programme - "
          f"read it as that gene and nothing more.")
if len(sig):
    print("Enriched beyond expression-matched random sets of the same size:")
    for _, r in sig.iterrows():
        print(f"   {r.theme}: {int(r.shared)}/{int(r.genes_tested)} shared against "
              f"{r.expected_matched:.1f} expected, q = {r.q_BH:.3f}")
        print(f"      {r.members}")
else:
    print("NOTHING is enriched beyond matched random sets. The themes read out of the Cell 13 list were an "
          "artefact of reading a ranked list, and should be dropped.")

if "Peroxide removal (the original hypothesis)" in set(en.theme):
    r = en[en.theme == "Peroxide removal (the original hypothesis)"].iloc[0]
    if pd.notna(r.get("p_vs_matched")):
        print(f"\nThe original hypothesis, for the record: {int(r.shared)}/{int(r.genes_tested)} peroxide-removal "
              f"genes shared, {r.expected_matched:.1f} expected by chance, p = {r.p_vs_matched:.3f}. "
              f"{'Still nothing.' if r.p_vs_matched > 0.1 else 'Worth a second look.'}")
en.to_csv(os.path.join(GDIR, "route_e_theme_enrichment.csv"), index=False)
print("\nSaved to", os.path.join(GDIR, "route_e_theme_enrichment.csv"))
