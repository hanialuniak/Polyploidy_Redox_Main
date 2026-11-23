#@title Cell 5 - Redox programmes: endocycling vs dividing, against random genes and stress controls
PROGRAMMES = {
    "REDOX: glutathione S-transferases": "GstD1 GstD2 GstD3 GstD5 GstD9 GstE1 GstE3 GstE6 GstE7 GstE9 GstE12 GstE14 GstS1 GstO1 GstO2 GstO3 GstT1 GstZ1 GstZ2 Mgstl",
    "REDOX: glutathione synthesis & recycling": "Gclc Gclm Gss1 Gss2 Gr1 Grx1 Grx1t PHGPx GILT1",
    "REDOX: thioredoxin / peroxiredoxin": "Trxr-1 Trxr-2 Trx-2 TrxT dhd Jafrac1 Jafrac2 Prx3 Prx5 Prx2540-1 Prx2540-2 Srx",
    "REDOX: SOD and catalase": "Sod1 Sod2 Sod3 Cat",
    "REDOX: NADPH supply": "Zw Pgd Men Men-b Idh Idh3a Nmdmc Tkt Taldo",
    "REDOX: cnc (NRF2) targets": "cnc Keap1 GstD1 Gclc Gclm Cyp6a2 Cyp6a8 Ugt35a Jafrac1 Prx3",
    "REDOX: H2O2 / ROS-producing enzymes": "Ero1L Duox Nox Cyp4e3 Mav Uro Aox1 Pxd Pxt",
    "REDOX: mitochondrial respiration": "ND-75 ND-49 ND-30 ND-42 SdhA SdhB UQCR-C1 UQCR-C2 Cyt-c-p COX4 COX5A ATPsynbeta blw",
    "EXPLORE: ER folding load (Zhou et al. 2026 model)": "Xbp1 Hsc70-3 Ire1 PEK crc Pdi ERp60 Calr Ero1L Sec61alpha Sec61beta Srp54 Hsp83",
    "EXPLORE: ribosome / protein synthesis": "RpL3 RpL4 RpL5 RpL7 RpL8 RpL10 RpL11 RpL13 RpS2 RpS3 RpS4 RpS6 RpS8 RpS9 RpS11 RpS15",
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


_present = set(L.index)
_renamed, _resolved = {}, {}
for n, v in PROGRAMMES.items():
    out = []
    for g in v.split():
        r = resolve(g, _present)
        if r:
            out.append(r)
            if r != g:
                _renamed[g] = r
    _resolved[n] = out
PROGRAMMES = _resolved
if _renamed:
    print("Found under a synonym:", ", ".join(f"{k} -> {v}" for k, v in _renamed.items()))
for n, v in PROGRAMMES.items():
    print(f"{n}: {len(v)} genes")
PROGRAMMES = {n: v for n, v in PROGRAMMES.items() if len(v) >= 3}
level = L.mean(axis=1)
bins = pd.qcut(level.rank(method="first"), 10, labels=False)
avoid = {g for v in PROGRAMMES.values() for g in v}
pools = {b: [g for g in bins.index[bins == b] if g not in avoid] for b in range(10)}
controls = [n for n in PROGRAMMES if n.startswith("CONTROL")]
rows = []
for t, r in CONTRASTS:
    fc = gmean[t] - gmean[r]
    for n, gs in PROGRAMMES.items():
        obs = fc.loc[gs].median()
        null = np.array([fc.loc[[rng.choice(pools[bins[g]]) for g in gs]].median() for _ in range(N_RANDOM_SETS)])
        c = np.median(null)
        p = (1 + np.sum(np.abs(null - c) >= abs(obs - c))) / (N_RANDOM_SETS + 1)
        beyond = [fc.loc[gs].median() - fc.loc[PROGRAMMES[k]].median() for k in controls]
        rows.append({"contrast": f"{t} vs {r}", "programme": n, "genes": len(gs), "log2_change": obs,
                     "random_median": c, "p_vs_random": p,
                     "beyond_all_controls": all(np.sign(b) == np.sign(obs - c) for b in beyond)})
res = pd.DataFrame(rows)
for ct in res.contrast.unique():
    m = (res.contrast == ct) & ~res.programme.str.startswith("CONTROL")
    p = res.loc[m, "p_vs_random"].to_numpy(); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    qq = np.empty_like(p); qq[o] = np.minimum(q, 1); res.loc[m, "q_fdr"] = qq
display(res.round(3))

print("\nVerdict (a programme counts only if it differs from random genes in EVERY contrast, in the same direction):")
verdict = {}
for n in PROGRAMMES:
    if n.startswith("CONTROL"): continue
    s = res[res.programme == n]
    rel = np.sign(s.log2_change - s.random_median)
    if (s.q_fdr < 0.05).all() and abs(rel.sum()) == len(s):
        v = ("HIGHER" if rel.iloc[0] > 0 else "LOWER") + " in endocycling cells, consistently"
        if not s.beyond_all_controls.all(): v += " (but not beyond every stress/housekeeping control)"
    elif (s.q_fdr < 0.05).any():
        v = "differs in some contrasts only"
    else:
        v = "no difference from random genes"
    verdict[n] = v
    print(f"   {n}: {v}")
if GSE == "GSE121956":
    print("\nReminder: salivary glands, brain and discs are different tissues. A difference here can come from tissue "
          "identity (for example secretion), not ploidy. GSE121955 (same cells, endocycles induced) is the cleaner test.")
