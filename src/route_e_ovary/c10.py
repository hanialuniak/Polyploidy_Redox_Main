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

# Rank WITHIN the germline, exactly as the follicle cells were ranked within their own lineage.
EARLY_GL = "bam nos vas Sxl stg".split()          # germarium: dividing cystocytes
gmit, _ = score(MARKERS["mitotic (cycling)"].split())
gear, _ = score(EARLY_GL)
gm = pd.Series(gmit).groupby(cells.cluster.to_numpy()).mean().reindex(g_clusters)
ge = pd.Series(gear).groupby(cells.cluster.to_numpy()).mean().reindex(g_clusters)
gm = (gm - gm.mean()) / (gm.std() or 1)
ge = (ge - ge.mean()) / (ge.std() or 1)
grole = {}
for cl in g_clusters:
    if gm[cl] > 0 and ge[cl] > 0:
        grole[cl] = "mitotic germline (germarium)"
    elif gm[cl] < 0 and ge[cl] < 0:
        grole[cl] = "endocycling germline (nurse cells)"
    else:
        grole[cl] = "germline, unclear stage"
display(pd.DataFrame({"role": pd.Series(grole), "mitotic_z": gm.round(2), "early_z": ge.round(2),
                      "cells": cells.cluster.value_counts().reindex(g_clusters)}))
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
mitg = [g for g in MARKERS["mitotic (cycling)"].split() if g in fc_gl.index]
print(f"Positive control, mitotic genes in nurse cells vs germarium: {fc_gl.loc[mitg].median():+.2f} "
      f"({len(mitg)} genes)")
if fc_gl.loc[mitg].median() > -0.2:
    print("WARNING: mitotic genes do not fall in the nurse cell group. The germline split is probably wrong, "
          "so treat everything below as uninterpretable.")

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
print("\nGene by gene. Each lineage is judged against its own genome, so a value near zero reads as flat "
      "rather than as agreement:")
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
