#@title Cell 11 - Unbiased: what actually changes across the switch, every gene
# Every cell so far asked about genes chosen in advance. This one asks the opposite question: take the whole
# transcriptome, rank it by how much it changes from mitotic to endocycling, and look at what comes out. That
# is the only way to find out what the switch is mostly doing, rather than what we went looking for.
#
# Three things are reported for every gene, and the third is the one that matters:
#   1. log2 fold change, endocycling vs mitotic follicle cells (pseudobulk, median-of-ratios).
#   2. Whether every library agrees on the direction. With 3 libraries this is real replication; with 1 it is
#      not, and the column says so.
#   3. The same change in the GERMLINE. Shared with the germline -> the endocycle. Follicle only -> secretory
#      differentiation. This is what stopped the redox leads being believed, and it applies to every gene here.
#
# No p-values at gene level. The groups are cells, not animals, so per-cell tests would report thousands of
# significant genes off a handful of libraries. Effect size plus direction agreement across libraries plus the
# germline comparison is the honest version of the same question.

MIN_DETECT = 0.10          #@param {type:"number"}
TOP_N = 40                 #@param {type:"integer"}

_mitx = np.flatnonzero(((cells.role == MIT) & cells.low_stress).to_numpy())
_endx = np.flatnonzero(((cells.role == ENDO) & cells.low_stress).to_numpy())
order = [g for g in fc_follicle.index if g in gpos]
cols = [gpos[g] for g in order]
det_m = pd.Series(np.asarray((X[_mitx][:, cols] > 0).mean(axis=0)).ravel(), index=order)
det_e = pd.Series(np.asarray((X[_endx][:, cols] > 0).mean(axis=0)).ravel(), index=order)

de = pd.DataFrame({"log2FC": fc_follicle.reindex(order),
                   "detect_mitotic": det_m, "detect_endo": det_e})
de["expr_level"] = [float(Lp.loc[g].mean()) for g in order]
de["expr_rank_pct"] = de.expr_level.rank(pct=True) * 100

libcols2 = sorted(per_lib)
if len(libcols2) > 1:
    sgn = pd.DataFrame({f"lib_{str(l)[:14]}": per_lib[l].reindex(order) for l in libcols2})
    de["libraries_agree"] = (np.sign(sgn).nunique(axis=1) == 1)
    de = pd.concat([de, sgn.round(2)], axis=1)
else:
    de["libraries_agree"] = np.nan

# The germline column, and the verdict that goes with it.
if "fc_gl" in dir():
    gl_pct = fc_gl.rank(pct=True) * 100
    f_pct = fc_follicle.rank(pct=True) * 100
    de["germline_log2FC"] = fc_gl.reindex(order)
    de["germline_pct"] = gl_pct.reindex(order)
    de["follicle_pct"] = f_pct.reindex(order)
    # Cell 10's verdict judges each lineage against its own genome, which is right for a handful of chosen
    # genes but far too permissive across 10,000: a gene moving 4.4 in follicle cells and 0.3 in the germline
    # would read as "both". Across the whole transcriptome the germline change has to be a real change too,
    # so a minimum effect size is required on that side as well.
    MIN_GL_EFFECT = 0.3

    def lineage_call(f_lfc, g_lfc, fp, gp):
        if pd.isna(g_lfc) or pd.isna(f_lfc):
            return "not in germline data"
        f_moves = abs(f_lfc) >= MIN_GL_EFFECT and (fp > 65 or fp < 35)
        g_moves = abs(g_lfc) >= MIN_GL_EFFECT and (gp > 65 or gp < 35)
        if not f_moves and not g_moves:
            return "neither"
        if f_moves and g_moves:
            return ("BOTH lineages - tracks endocycling" if np.sign(f_lfc) == np.sign(g_lfc)
                    else "opposite directions")
        return "follicle only - likely differentiation" if f_moves else "germline only"

    de["lineage"] = [lineage_call(f, g, fp, gp) for f, g, fp, gp in
                     zip(de.log2FC, de.germline_log2FC, de.follicle_pct, de.germline_pct)]
else:
    de["lineage"] = "Cell 10 not run"

# Annotate anything that was in a curated set, so a familiar gene is recognisable in the ranked list.
member = {}
for nm, gs in (PROGRAMMES.items() if "PROGRAMMES" in dir() else []):
    for g in gs:
        member.setdefault(g, []).append(nm)
de["in_set"] = [", ".join(member.get(g, [])) for g in de.index]

keep = de[(de[["detect_mitotic", "detect_endo"]].max(axis=1) >= MIN_DETECT)].copy()
print(f"{len(de):,} genes in the comparison; {len(keep):,} detected in at least {MIN_DETECT:.0%} of cells in "
      f"one group or the other. The rest are too sparse to rank.")
if len(libcols2) > 1:
    agree_only = keep[keep.libraries_agree]
    print(f"Of those, {len(agree_only):,} change in the same direction in all {len(libcols2)} libraries. "
          f"Rankings below use that subset.")
    rank_on = agree_only
else:
    print("Only one library has both cell groups, so there is no direction check and nothing below is "
          "replicated. Treat the ranking as a list of candidates, not results.")
    rank_on = keep

show = ["log2FC", "detect_mitotic", "detect_endo", "expr_rank_pct", "germline_log2FC", "lineage", "in_set"]
show = [c for c in show if c in rank_on.columns]
print(f"\n=== TOP {TOP_N} UP in endocycling follicle cells ===")
display(rank_on.nlargest(TOP_N, "log2FC")[show].round(2))
print(f"\n=== TOP {TOP_N} DOWN in endocycling follicle cells ===")
display(rank_on.nsmallest(TOP_N, "log2FC")[show].round(2))

# The subset that matters most: large changes that the germline makes too.
if "germline_log2FC" in rank_on.columns:
    both = rank_on[rank_on.lineage.astype(str).str.startswith("BOTH")]
    both = both.reindex(both.log2FC.abs().sort_values(ascending=False).index)
    print(f"\n=== Changes shared with the GERMLINE: {len(both):,} genes ===")
    print("These are the candidates for tracking the endocycle itself rather than secretory differentiation.")
    display(both.head(TOP_N)[show].round(2))
    n_fol = int((rank_on.lineage.astype(str).str.startswith("follicle only")).sum())
    print(f"\nFor scale: {len(both):,} genes shared with the germline against {n_fol:,} follicle-only. "
          "A gene appearing here is necessary but not sufficient - the germline split is a two-group "
          "comparison with its own noise, so treat this as the shortlist to check, not a result.")

out = os.path.join(GDIR, "route_e_unbiased_de.csv")
de.sort_values("log2FC").to_csv(out)
print("\nSaved the full ranked table (every gene, unfiltered) to", out)
print("Columns: log2FC, detection in each group, expression rank, per-library fold changes, the germline "
      "value and the lineage verdict. Open it in a spreadsheet and sort on whatever you like.")
