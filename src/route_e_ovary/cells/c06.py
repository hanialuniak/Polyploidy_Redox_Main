#@title Cell 6 - Redox programmes: endocycling vs mitotic follicle cells
PROGRAMMES = {
    "REDOX: glutathione S-transferases": "GstD1 GstD2 GstD3 GstD5 GstD9 GstE1 GstE3 GstE6 GstE7 GstE9 GstE12 GstE14 GstS1 GstO1 GstO2 GstO3 GstT1 GstZ1 GstZ2 Mgstl",
    "REDOX: glutathione synthesis & recycling": "Gclc Gclm Gss1 Gss2 Gr1 Grx1 Grx1t PHGPx GILT1",
    "REDOX: thioredoxin / peroxiredoxin": "Trxr-1 Trxr-2 Trx-2 TrxT dhd Jafrac1 Jafrac2 Prx3 Prx5 Prx2540-1 Prx2540-2 Srx",
    "REDOX: SOD and catalase": "Sod1 Sod2 Sod3 Cat",
    "REDOX: NADPH supply": "Zw Pgd Men Men-b Idh Idh3a Nmdmc Tkt Taldo",
    "REDOX: cnc (NRF2) targets": "cnc Keap1 GstD1 Gclc Gclm Cyp6a2 Cyp6a8 Ugt35a Jafrac1 Prx3",
    "REDOX: H2O2 / ROS-producing enzymes": "Ero1L Duox Nox Cyp4e3 Mav Uro Aox1 Pxd Pxt",
    "REDOX: mitochondrial respiration": "ND-75 ND-49 ND-30 ND-42 SdhA SdhB UQCR-C1 UQCR-C2 Cyt-c-p COX4 COX5A ATPsynbeta blw",
    "EXPLORE: ER folding load": "Xbp1 Hsc70-3 Ire1 PEK crc Pdi ERp60 Calr Ero1L Sec61alpha Sec61beta Srp54",
    "EXPLORE: ribosome / protein synthesis": "RpL3 RpL4 RpL5 RpL7 RpL8 RpL10 RpL11 RpL13 RpS2 RpS3 RpS4 RpS6 RpS8 RpS9 RpS11 RpS15",
    "EXPLORE: licensing and re-replication": "geminin dup Cdc6 Orc1 Mcm2 Mcm3 dpa Mcm5 Mcm6 Mcm7 CycE E2f1 Myc",
    "EXPLORE: p53 / DNA damage": "p53 Ku80 Irbp lok mei-41 tefu rad50 Xrp1 spn-A mus304",
    "CONTROL: heat shock": "Hsp70Aa Hsp83 Hsp68 Hsp27 Hsp26 Hsp23 DnaJ-1 Hsf",
    "CONTROL: immediate-early genes": "Hr38 sr kay Jra puc CrebA",
    "CONTROL: housekeeping": "Act5C alphaTub84B RpL32 Gapdh1 Tbp eIF4A Pgk Ef1alpha48D",
}
# Some genes are listed under a different name in the annotation, so each one is looked up under its synonyms too.
# FlyBase calls Geminin "geminin" (CG3183); "Gem" is only a synonym, which is why earlier runs could not find it.
ALIASES = {"geminin": ["Gem", "CG3183"], "dup": ["Cdt1"], "dpa": ["Mcm4"], "Zw": ["G6PD", "G6pd"],
           "Ero1L": ["Ero1", "CG1333"], "PHGPx": ["CG12013"], "Gss1": ["Gss"], "Gr1": ["Gr", "GR"],
           "Trxr-1": ["Trxr1"], "Jafrac1": ["Prx2"], "Men-b": ["Men-b", "CG5889"]}


def resolve(g, present):
    if g in present:
        return g
    for a in ALIASES.get(g, []):
        if a in present:
            return a
    return None


_renamed = {}
_resolved = {}
for n, v in PROGRAMMES.items():
    out = []
    for g in v.split():
        r = resolve(g, gpos)
        if r:
            out.append(r)
            if r != g:
                _renamed[g] = r
    _resolved[n] = out
PROGRAMMES = _resolved
if _renamed:
    print("Found under a synonym:", ", ".join(f"{k} -> {v}" for k, v in _renamed.items()))
_lic = [n for n in PROGRAMMES if "licensing" in n]
if _lic:
    print("Licensing set as tested:", ", ".join(PROGRAMMES[_lic[0]]))
for n, v in PROGRAMMES.items():
    print(f"{n}: {len(v)} genes found")
PROGRAMMES = {n: v for n, v in PROGRAMMES.items() if len(v) >= 3}

# Pseudobulk: add up the counts of each group, per library, then normalise by median of ratios.
# This avoids yolk protein genes (very high in endocycling cells) dominating the scale.
pb, lab = [], []
for lib in sorted(cells.library.unique()):
    for r in (MIT, ENDO):
        m = ((cells.library == lib) & (cells.role == r) & cells.low_stress).to_numpy()
        if m.sum() >= 50:
            pb.append(np.asarray(X[np.flatnonzero(m)].sum(axis=0)).ravel()); lab.append((lib, r))
pbd = pd.DataFrame(np.array(pb).T, index=genes, columns=pd.MultiIndex.from_tuples(lab))
pbd = pbd[(pbd / pbd.sum() * 1e6 >= 5).all(axis=1)]
lg = np.log(pbd.replace(0, np.nan)); ref = lg.mean(axis=1)
size = np.exp(lg[ref.notna()].sub(ref[ref.notna()], axis=0).median())
Lp = np.log2(pbd / size / pbd.sum().mean() * 1e6 + 1)
n_lib = len({l for l, _ in lab if (l, MIT) in lab and (l, ENDO) in lab})
print(f"\nStress-matched cells only. Pseudobulk groups: {len(lab)} ({n_lib} librar{'y' if n_lib == 1 else 'ies'} have both cell groups), "
      f"{len(Lp):,} genes expressed in all of them")
fc = Lp.xs(ENDO, level=1, axis=1).mean(axis=1) - Lp.xs(MIT, level=1, axis=1).mean(axis=1)
per_lib = {l: Lp[(l, ENDO)] - Lp[(l, MIT)] for l in sorted({a for a, b in lab}) if (l, ENDO) in lab and (l, MIT) in lab}

level = Lp.mean(axis=1)
bins = pd.qcut(level.rank(method="first"), 10, labels=False)
avoid = {g for v in PROGRAMMES.values() for g in v}
pools = {b: [g for g in bins.index[bins == b] if g not in avoid] for b in range(10)}
controls = [n for n in PROGRAMMES if n.startswith("CONTROL")]
rows = []
for n, gs in PROGRAMMES.items():
    gs = [g for g in gs if g in fc.index]
    if len(gs) < 3:
        continue
    obs = fc.loc[gs].median()
    null = np.array([fc.loc[[rng.choice(pools[bins[g]]) for g in gs]].median() for _ in range(N_RANDOM_SETS)])
    c = np.median(null)
    r = {"programme": n, "genes": len(gs), "log2_change": obs, "random_median": c,
         "p_vs_random": (1 + np.sum(np.abs(null - c) >= abs(obs - c))) / (N_RANDOM_SETS + 1),
         "beyond_all_controls": all(np.sign(obs - fc.loc[[x for x in PROGRAMMES[k] if x in fc.index]].median())
                                    == np.sign(obs - c) for k in controls)}
    for l, d in per_lib.items():
        r[f"lib_{l[:14]}"] = d.loc[gs].median()
    rows.append(r)
res = pd.DataFrame(rows)
m = ~res.programme.str.startswith("CONTROL")
p = res.loc[m, "p_vs_random"].to_numpy(); o = np.argsort(p)
q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
qq = np.empty_like(p); qq[o] = np.minimum(q, 1); res.loc[m, "q_fdr"] = qq
display(res.round(3))
print("log2_change: endocycling minus mitotic follicle cells. Positive means higher in endocycling cells.")
libcols = [c for c in res.columns if c.startswith("lib_")]
print("\nVerdict:")
for _, r in res[m].iterrows():
    agree = len(libcols) < 2 or len({np.sign(r[c] - res.loc[m, c].median()) for c in libcols}) == 1
    v = ("HIGHER" if r.log2_change > r.random_median else "LOWER") + " in endocycling follicle cells"
    if r.q_fdr >= 0.05:
        v = "no difference from random genes"
    elif not r.beyond_all_controls:
        v += " (but not beyond every stress/housekeeping control)"
    elif not agree:
        v += " (but the libraries disagree)"
    print(f"   {r.programme}: {v} ({r.log2_change:+.2f})")
if len(libcols) < 2:
    print("\nOnly one library has both cell groups, so there is no biological replication. The comparison with random "
          "gene sets is the only safeguard here: read the result as a lead, not as a measured effect size.")
print("\nNote: in the follicle cell lineage, ploidy and developmental stage rise together. A difference here means "
      "endocycling cells differ from mitotic ones, not that ploidy alone caused it.")
