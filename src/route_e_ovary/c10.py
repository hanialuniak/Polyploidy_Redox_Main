#@title Cell 10 - The same switch in the germline: nurse cells as an independent endocycling lineage
# Follicle cells endocycle AND differentiate into an eggshell-secreting epithelium at the same time, so this
# notebook cannot tell those apart. The germline can. Cystocytes divide mitotically in the germarium, then
# become endocycling nurse cells - the same switch, in the same egg chambers, in the same libraries, with a
# completely different differentiation programme (maternal mRNA and protein, no eggshell).
#   Change seen in BOTH lineages  -> tracks endocycling.
#   Change seen only in follicle cells -> tracks secretory differentiation.
fc_follicle = fc.copy()
res_follicle = res.copy()

gl_cells = cells.cluster.map(zt[gl_col] >= 0.5).fillna(False).to_numpy()
if gl_cells.sum() < 100:
    raise StopHere(f"Only {int(gl_cells.sum())} germline cells found, too few for this comparison. "
                   "The germline marker column is " + repr(gl_col) + "; check the z-score table in Cell 4.")
g_clusters = [cl for cl in zt.index if zt.loc[cl, gl_col] >= 0.5]
print(f"Germline clusters: {g_clusters}  ({int(gl_cells.sum()):,} cells)")

# Split the germline on M PHASE alone. Nurse cells endocycle, so they run S phase repeatedly and express
# PCNA and the Mcm genes; what they never do is divide. So S-phase genes cannot tell the two apart and
# M-phase and cytokinesis genes can. (The follicle-cell marker set mixes the two, which is why it was wrong
# here.) bam marks germarial cysts, but it is on in only a few cells, so it is a check and not a requirement:
# vas, nos and Sxl are pan-germline and are higher in nurse cells simply because those cells are enormous.
MPHASE = "polo aurB pav Incenp CycB Det Sak Cks30A".split()
EARLY_GL = "bam bgcn".split()
#@markdown Leave empty to split automatically. If the automatic split looks wrong in the table below, put the
#@markdown germarium cluster ids here, comma separated, and run this cell again.
GERMARIUM_CLUSTERS = ""  #@param {type:"string"}

gm = pd.Series(score(MPHASE)[0]).groupby(cells.cluster.to_numpy()).mean().reindex(g_clusters)
ge = pd.Series(score(EARLY_GL)[0]).groupby(cells.cluster.to_numpy()).mean().reindex(g_clusters)
n_counts = cells.groupby("cluster").n_counts.median().reindex(g_clusters)
zm = (gm - gm.mean()) / (gm.std() or 1)
ze = (ge - ge.mean()) / (ge.std() or 1)

forced = [c.strip() for c in GERMARIUM_CLUSTERS.split(",") if c.strip()]
if forced:
    grole = {cl: ("mitotic germline (germarium)" if cl in forced else "endocycling germline (nurse cells)")
             for cl in g_clusters}
    print("Using the clusters you named as germarium:", forced)
else:
    # Rank on M phase and cut at the largest gap, which puts the dividing clusters on one side however many
    # clusters there are. With two clusters this is simply the split between them.
    srt = gm.sort_values()
    gaps = srt.diff().iloc[1:]
    cutm = srt.loc[gaps.idxmax()] if len(srt) > 1 and gaps.notna().any() else gm.median()
    grole = {cl: ("mitotic germline (germarium)" if gm[cl] >= cutm else "endocycling germline (nurse cells)")
             for cl in g_clusters}
    if len(g_clusters) > 2 and (gm.max() - gm.min()) < 1e-9:
        raise StopHere("The germline clusters do not differ in M-phase gene expression at all, so they cannot "
                       "be split. Raise N_CLUSTERS in Cell 1.")
display(pd.DataFrame({"role": pd.Series(grole), "Mphase_z": zm.round(2), "bam_z": ze.round(2),
                      "median_counts": n_counts.round(0),
                      "cells": cells.cluster.value_counts().reindex(g_clusters)}))
_g = [cl for cl in g_clusters if grole[cl].startswith("mitotic")]
_n = [cl for cl in g_clusters if grole[cl].startswith("endocycling")]
if _g and _n and ze.notna().any() and ze[_g].mean() < ze[_n].mean():
    print("NOTE: bam is not higher in the group called germarium. bam is on in only a handful of cells, so this "
          "often happens and is not by itself a problem; but if the M-phase separation is also weak, set "
          "GERMARIUM_CLUSTERS by hand from the table above.")
if _g and _n and n_counts[_n].mean() < n_counts[_g].mean():
    print("NOTE: the group called nurse cells does not have more counts per cell than the germarium group. "
          "Nurse cells are far larger, so check the split before trusting the result.")
cells["gl_role"] = cells.cluster.map(grole)
GMIT, GENDO = "mitotic germline (germarium)", "endocycling germline (nurse cells)"
n_ok = [r for r in (GMIT, GENDO) if (cells.gl_role == r).sum() >= 30]
if len(n_ok) < 2:
    raise StopHere("Could not separate dividing germarial germline from endocycling nurse cells "
                   f"(found {cells.gl_role.value_counts().to_dict()}). Raising N_CLUSTERS in Cell 1 often "
                   "splits them; if the dataset has no germarium, this comparison cannot be made here.")

# Same stress matching as the follicle comparison
cut_g = min(np.quantile(cells.stress[cells.gl_role == GMIT], 0.75),
            np.quantile(cells.stress[cells.gl_role == GENDO], 0.75))
low_g = cells.stress <= cut_g
pb, lab = [], []
for lib in sorted(cells.library.unique()):
    for r in (GMIT, GENDO):
        m = ((cells.library == lib) & (cells.gl_role == r) & low_g).to_numpy()
        if m.sum() >= 30:
            pb.append(np.asarray(X[np.flatnonzero(m)].sum(axis=0)).ravel()); lab.append((lib, r))
if len({l for l, r in lab if (l, GMIT) in lab and (l, GENDO) in lab}) == 0:
    raise StopHere("No library has at least 30 cells of both germline groups after stress matching.")
pbg = pd.DataFrame(np.array(pb).T, index=genes, columns=pd.MultiIndex.from_tuples(lab))
pbg = pbg[(pbg / pbg.sum() * 1e6 >= 5).all(axis=1)]
lg = np.log(pbg.replace(0, np.nan)); refg = lg.mean(axis=1)
sizeg = np.exp(lg[refg.notna()].sub(refg[refg.notna()], axis=0).median())
Lg = np.log2(pbg / sizeg / pbg.sum().mean() * 1e6 + 1)
fc_gl = Lg.xs(GENDO, level=1, axis=1).mean(axis=1) - Lg.xs(GMIT, level=1, axis=1).mean(axis=1)
print(f"\nGermline pseudobulk: {len(lab)} groups, {len(Lg):,} genes expressed in all of them")

# Positive control: nurse cells must look less mitotic and more endoreplicated than germarial cysts
# The control has to be strict. Germarial cysts divide; nurse cells never do, so M-phase genes should fall
# steeply between them. A shallow fall means both groups are endocycling - that is, early against late nurse
# cells - and that comparison is NOT the mitotic-to-endocycle switch this cell exists to test.
mitg = [g for g in MPHASE if g in fc_gl.index]
drop = float(fc_gl.loc[mitg].median()) if mitg else float("nan")
print(f"\nPositive control, M-phase genes in nurse cells vs germarium: {drop:+.2f} ({len(mitg)} genes)")
_bam = resolve("bam", set(genes))
if _bam and _bam in gpos:
    bi = gpos[_bam]
    for r in (GMIT, GENDO):
        col = X[np.flatnonzero((cells.gl_role == r).to_numpy())][:, bi]
        print(f"   cells expressing {_bam}: {r} {100 * float((col.toarray().ravel() > 0).mean()):.1f}%")
GERMLINE_USABLE = bool(mitg) and drop <= -0.8
if not GERMLINE_USABLE:
    print("""
WARNING - the germline comparison did not pass its control, so DO NOT read the verdicts below as
follicle-versus-germline. M-phase genes barely fall, which means both germline groups are almost certainly
endocycling already: this is early against late nurse cells, not germarium against nurse cells. A true
germarium is a small structure and many ovary libraries barely sample it.
What to try, in order:
  1. Raise N_CLUSTERS in Cell 1 and run again - a germarial cluster is often hidden inside a larger one.
  2. Run the other atlas. GSE136162 has three libraries and samples the germarium better.
  3. If neither separates a dividing germline population, this dataset cannot answer the question and the
     follicle-cell result has to stand on its own.
""")

levg = Lg.mean(axis=1)
binsg = pd.qcut(levg.rank(method="first"), 10, labels=False)
avoidg = {g for v in PROGRAMMES.values() for g in v}
poolsg = {b: [g for g in binsg.index[binsg == b] if g not in avoidg] for b in range(10)}
rows = []
for n, gs in PROGRAMMES.items():
    gs = [g for g in gs if g in fc_gl.index]
    if len(gs) < 3:
        continue
    obs = float(fc_gl.loc[gs].median())
    null = np.array([fc_gl.loc[[rng.choice(poolsg[binsg[g]]) for g in gs]].median() for _ in range(N_RANDOM_SETS)])
    c = float(np.median(null))
    rows.append({"programme": n, "genes": len(gs), "germline": obs, "random_median": c,
                 "p_vs_random": (1 + np.sum(np.abs(null - c) >= abs(obs - c))) / (N_RANDOM_SETS + 1),
                 "pct_of_all_genes": float(fc_gl.rank(pct=True).loc[gs].median()) * 100})
resg = pd.DataFrame(rows)
foll = res_follicle.set_index("programme")["log2_change"]
resg["follicle"] = resg.programme.map(foll)
resg["both_same_direction"] = np.sign(resg.germline) == np.sign(resg.follicle)
print("\nProgrammes, germline against follicle cells:")
display(resg[["programme", "genes", "follicle", "germline", "pct_of_all_genes",
              "p_vs_random", "both_same_direction"]].round(3))

KEY = ["Jafrac1", "Jafrac2", "Ero1L", "Gclc", "Trx-2", "Trxr-1", "Rrp1", "Hsc70-3", "Sec61beta", "Cat", "geminin"]
key = [resolve(g, set(fc_gl.index)) for g in KEY]
key = [(k, g) for k, g in zip(KEY, key) if g and g in fc_follicle.index]
tab_k = pd.DataFrame({"follicle": [float(fc_follicle.loc[g]) for _, g in key],
                      "germline": [float(fc_gl.loc[g]) for _, g in key],
                      "germline_pct": [float(fc_gl.rank(pct=True).loc[g]) * 100 for _, g in key]},
                     index=[k for k, _ in key])
fpct = fc_follicle.rank(pct=True)
tab_k["follicle_pct"] = [float(fpct.loc[g]) * 100 for _, g in key]


def verdict(fp, gp):
    """Judge each lineage against its own genome rather than against zero."""
    w = lambda p: "up" if p > 65 else "down" if p < 35 else "flat"
    a, b = w(fp), w(gp)
    if a == "flat" and b == "flat":
        return "neither"
    if a == b:
        return "BOTH lineages - tracks endocycling"
    if b == "flat":
        return "follicle only - likely differentiation"
    if a == "flat":
        return "germline only"
    return "opposite directions"


tab_k["verdict"] = [verdict(f, g) for f, g in zip(tab_k.follicle_pct, tab_k.germline_pct)]
tab_k = tab_k[["follicle", "follicle_pct", "germline", "germline_pct", "verdict"]]
hdr = ("\nGene by gene. Each lineage is judged against its own genome, so a value near zero reads as flat "
       "rather than as agreement:")
if not GERMLINE_USABLE:
    hdr += ("\n*** The control above did not pass, so the 'germline' column is early vs late nurse cells, not "
            "the switch. Every verdict in this table is unreliable until that is fixed. ***")
print(hdr)
display(tab_k.round(2))
print("""
How to read this. A change that appears in BOTH lineages tracks the mitotic-to-endocycle switch itself, because
the only thing nurse cells and follicle cells share here is that switch. A change in follicle cells ALONE is
most likely secretory differentiation, which is what the salivary gland comparison could never exclude.
Two cautions: nurse cells are germline, so germline-restricted genes such as dhd are not comparable across the
two lineages; and they are far larger than follicle cells, which is why the pseudobulk uses median-of-ratios.
""")
fc_gl.to_csv(os.path.join(GDIR, "germline_fold_changes.csv"))
resg.to_csv(os.path.join(GDIR, "germline_programmes.csv"), index=False)
tab_k.to_csv(os.path.join(GDIR, "germline_vs_follicle_genes.csv"))
print("Saved germline tables to", GDIR)
