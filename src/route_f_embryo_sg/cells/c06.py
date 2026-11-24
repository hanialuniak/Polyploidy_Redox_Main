#@title Cell 6 - Redox programmes: what changes in the gland between the two windows
PROGRAMMES = {
    "REDOX: glutathione S-transferases": "GstD1 GstD2 GstD3 GstD5 GstD9 GstE1 GstE3 GstE6 GstE7 GstE9 GstE12 GstE14 GstS1 GstO1 GstO2 GstO3 GstT1 GstZ1 GstZ2 Mgstl",
    "REDOX: glutathione synthesis & recycling": "Gclc Gclm Gss1 Gss2 Gr1 Grx1 Grx1t PHGPx GILT1",
    "REDOX: thioredoxin / peroxiredoxin": "Trxr-1 Trxr-2 Trx-2 TrxT dhd Jafrac1 Jafrac2 Prx3 Prx5 Prx2540-1 Prx2540-2 Srx",
    "REDOX: SOD and catalase": "Sod1 Sod2 Sod3 Cat",
    "REDOX: NADPH supply": "Zw Pgd Men Men-b Idh Idh3a Nmdmc Tkt Taldo",
    "REDOX: cnc (NRF2) targets": "cnc Keap1 GstD1 Gclc Gclm Cyp6a2 Cyp6a8 Ugt35a Jafrac1 Prx3",
    "REDOX: H2O2-producing enzymes": "Ero1L Duox Nox Aox1 Uro Cyp4e3",
    "REDOX: mitochondrial respiration": "ND-75 ND-49 ND-30 ND-42 SdhA SdhB UQCR-C1 UQCR-C2 Cyt-c-p COX4 COX5A ATPsynbeta blw",
    "EXPLORE: ER folding load": "Xbp1 Hsc70-3 Ire1 PEK crc Pdi ERp60 Calr Ero1L Sec61alpha Sec61beta Srp54",
    "EXPLORE: ribosome / protein synthesis": "RpL3 RpL4 RpL5 RpL7 RpL8 RpL10 RpL11 RpL13 RpS2 RpS3 RpS4 RpS6 RpS8 RpS9 RpS11 RpS15",
    "EXPLORE: licensing and re-replication": "geminin dup Cdc6 Orc1 Mcm2 Mcm3 dpa Mcm5 Mcm6 Mcm7 CycE E2f1 Myc",
    "EXPLORE: p53 / DNA damage": "p53 Ku80 Irbp lok mei-41 tefu rad50 Xrp1 spn-A mus304",
    "CONTROL: heat shock": "Hsp70Aa Hsp83 Hsp68 Hsp27 Hsp26 Hsp23 DnaJ-1 Hsf",
    "CONTROL: immediate-early genes": "Hr38 sr kay Jra puc",
    "CONTROL: housekeeping": "Act5C alphaTub84B RpL32 Gapdh1 Tbp eIF4A Pgk Ef1alpha48D",
}
_res, _renamed = {}, {}
for n, v in PROGRAMMES.items():
    out = []
    for g in v.split():
        r = resolve(g, PRESENT)
        if r:
            out.append(r)
            if r != g:
                _renamed[g] = r
    _res[n] = out
PROGRAMMES = {n: v for n, v in _res.items() if len(v) >= 3}
if _renamed:
    print("Found under a synonym:", ", ".join(f"{k} -> {v}" for k, v in _renamed.items()))
for n, v in PROGRAMMES.items():
    print(f"{n}: {len(v)} genes")

level = L.mean(axis=1)
bins = pd.qcut(level.rank(method="first"), 10, labels=False)
avoid = {g for v in PROGRAMMES.values() for g in v}
pools = {b: [g for g in bins.index[bins == b] if g not in avoid] for b in range(10)}
controls = [n for n in PROGRAMMES if n.startswith("CONTROL")]
rows = []
for label, series in [("gland-specific (gland minus other cells)", fc),
                      ("salivary gland, late vs early", fc_sg),
                      ("other cells, late vs early", fc_oth)]:
    for n, gs in PROGRAMMES.items():
        obs = series.loc[gs].median()
        null = np.array([series.loc[[rng.choice(pools[bins[g]]) for g in gs]].median() for _ in range(N_RANDOM_SETS)])
        c = np.median(null)
        p = (1 + np.sum(np.abs(null - c) >= abs(obs - c))) / (N_RANDOM_SETS + 1)
        beyond = all(np.sign(obs - series.loc[PROGRAMMES[k]].median()) == np.sign(obs - c) for k in controls)
        rows.append({"contrast": label, "programme": n, "genes": len(gs), "log2_change": obs,
                     "random_median": c, "p_vs_random": p, "beyond_all_controls": beyond,
                     "pct_of_all_genes": series.rank(pct=True).loc[gs].median() * 100})
res = pd.DataFrame(rows)
for ct in res.contrast.unique():
    m = (res.contrast == ct) & ~res.programme.str.startswith("CONTROL")
    p = res.loc[m, "p_vs_random"].to_numpy(); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    qq = np.empty_like(p); qq[o] = np.minimum(q, 1); res.loc[m, "q_fdr"] = qq
display(res[res.contrast == res.contrast.iloc[0]].round(3))
print("\nThe gland-specific contrast is the one that matters: it removes changes happening everywhere in the "
      "embryo between these windows. The other two are shown in the saved table for reference.")
