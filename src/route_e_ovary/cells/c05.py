#@title Cell 5 - Pick mitotic and endocycling follicle cells, and check them
fc_col = [c for c in zt.columns if c.startswith("follicle cell")][0]
gl_col = [c for c in zt.columns if c.startswith("germline")][0]
mu_col = [c for c in zt.columns if c.startswith("muscle")][0]
mit_col = [c for c in zt.columns if c.startswith("mitotic")][0]
early_col = [c for c in zt.columns if c.startswith("early follicle")][0]
mid_col = [c for c in zt.columns if c.startswith("mid follicle")][0]
late_col = [c for c in zt.columns if c.startswith("late follicle")][0]
amp_col = [c for c in zt.columns if c.startswith("vitelline")][0]

esc_col = [c for c in zt.columns if c.startswith("escort")][0]
# Follicle cells are defined by EXCLUDING the other lineages. The general epithelial markers (CadN, Fas2, br)
# fade as follicle cells mature, so requiring them would throw away exactly the late, most polyploid cells.
is_fc = (zt[gl_col] < 0.5) & (zt[mu_col] < 0.5) & (zt[esc_col] < 0.5)
stage_raw = zt[[mid_col, late_col, amp_col]].mean(axis=1) - zt[early_col]
# Rank the follicle clusters among themselves, so the split is early vs late WITHIN the lineage
f = zt.index[is_fc]
stage = (stage_raw - stage_raw[f].mean()) / (stage_raw[f].std() or 1)
mitr = (zt[mit_col] - zt.loc[f, mit_col].mean()) / (zt.loc[f, mit_col].std() or 1)
role = {}
for cl in zt.index:
    if not is_fc[cl]:
        role[cl] = "other cell type"
    elif mitr[cl] > 0 and stage[cl] < 0:
        role[cl] = "mitotic follicle cell"
    elif mitr[cl] < 0 and stage[cl] > 0:
        role[cl] = "endocycling follicle cell"
    else:
        role[cl] = "follicle cell, unclear stage"
cells["role"] = cells.cluster.map(role)
summary = pd.DataFrame({"role": pd.Series(role), "mitotic_z": mitr.round(2),
                        "stage_score": stage.round(2), "cells": cells.cluster.value_counts()})
display(summary.sort_values(["role", "stage_score"]))
print(cells.role.value_counts().to_string())
MIT, ENDO = "mitotic follicle cell", "endocycling follicle cell"
if not {MIT, ENDO} <= set(cells.role.unique()):
    raise StopHere("Could not find both mitotic and endocycling follicle cell clusters. Send me the tables above; "
                   "raising N_CLUSTERS in Cell 1 often separates them.")

# Checks: the two groups must differ in the expected way
for name, gs in [("mitotic genes (expected much lower in endocycling cells)", MARKERS["mitotic (cycling)"]),
                 ("early markers CadN, Fas2, ct (expected lower in endocycling cells)", MARKERS["early follicle (stages 1-6)"]),
                 ("yolk and mid markers (expected higher in endocycling cells)", MARKERS["mid follicle (stages 7-10, endocycling)"])]:
    s, n = score(gs.split())
    a, b = s[(cells.role == MIT).to_numpy()].mean(), s[(cells.role == ENDO).to_numpy()].mean()
    print(f"{name}: mitotic {a:.2f}, endocycling {b:.2f}, difference {b - a:+.2f}")
print(f"\nMedian counts per cell: mitotic {cells.n_counts[cells.role == MIT].median():,.0f}, "
      f"endocycling {cells.n_counts[cells.role == ENDO].median():,.0f}")

# Dissociation stress: pulling a tissue apart switches on heat shock and immediate-early genes, and it hits
# large fragile cells hardest. If one group is more stressed, that alone shifts many genes.
DISSOC = "Hsp70Aa Hsp70Ab Hsp68 Hsp26 Hsp23 Hsp27 Hsp83 DnaJ-1 Hsf Jra kay Hr38 sr puc".split()
cells["stress"], _ = score(DISSOC)[0], None
cells["stress"] = score(DISSOC)[0]
for r in (MIT, ENDO):
    v = cells.stress[cells.role == r]
    print(f"Dissociation stress score, {r}: median {v.median():.2f}")
cut = min(np.quantile(cells.stress[cells.role == MIT], 0.75), np.quantile(cells.stress[cells.role == ENDO], 0.75))
cells["low_stress"] = cells.stress <= cut
print(f"Cells kept for the stress-matched comparison (score <= {cut:.2f}): "
      f"{int(((cells.role == MIT) & cells.low_stress).sum()):,} mitotic, "
      f"{int(((cells.role == ENDO) & cells.low_stress).sum()):,} endocycling")
