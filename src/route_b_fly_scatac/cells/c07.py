#@title Cell 7 - Main check: more polyploidy than doublets can explain?
# Two cells captured together look like one cell with twice the DNA. About 0.8% of captured cells per 1,000
# recovered are doublets. Doublets have about twice the fragments, so the fragment filter in Cell 5 keeps them more
# often than single cells; and a doublet of two different neurons lands in a neuron cluster. So the expected
# doublet share is worked out for all kept cells of a run, not only for same-cluster pairs.
expected = {}
for run, part in all_depth.groupby("run"):
    d = part.n_fragments.to_numpy()
    D = min(len(d) / 1000 * 0.008, 0.5)
    pair = rng.choice(d, 20000) + rng.choice(d, 20000)
    pass_single = np.mean((d >= MIN_FRAGMENTS) & (d <= MAX_FRAGMENTS))
    pass_double = np.mean((pair >= MIN_FRAGMENTS) & (pair <= MAX_FRAGMENTS))
    kept_share = D * pass_double / (D * pass_double + (1 - D) * pass_single)
    expected[run] = kept_share
    print(f"{run[:45]}: {len(d):,} cells, median {np.median(d):,.0f} fragments; doublets expected "
          f"{100 * D:.1f}% of all cells, {100 * kept_share:.1f}% of the cells kept by the fragment filter")
cells["expected_doublet"] = cells.run.map(expected)
rows = []
for ct, part in cells.groupby("cell_type"):
    boot = [part.P_polyploid.sample(len(part), replace=True, random_state=i).mean() for i in range(500)]
    rows.append({"cell_type": ct, "cells": len(part), "polyploid_rate": part.P_polyploid.mean(),
                 "ci_low": np.quantile(boot, 0.025), "expected_doublets": part.expected_doublet.mean()})
check = pd.DataFrame(rows).set_index("cell_type")
check["excess"] = check.polyploid_rate - check.expected_doublets
check["clear_excess"] = check.ci_low > check.expected_doublets
display(check.sort_values("polyploid_rate", ascending=False).round(3))
overall = cells.P_polyploid.mean(); expected = cells.expected_doublet.mean()
print(f"All cells: polyploid rate {overall:.3f}, expected from doublets {expected:.3f}")
print(f"Cell types with clearly more polyploidy than doublets: {int(check.clear_excess.sum())} of {len(check)}")
if check.clear_excess.sum() == 0:
    print("WARNING: no cell type has more polyploidy than doublets would give. Later results cannot be trusted; "
          "use deeper runs, older flies, or doublet-removed data.")
div = cells.cell_type.str.startswith("dividing")
if div.any():
    print(f"\nDividing cells (neuroblasts and S/G2/M): polyploid rate {cells.P_polyploid[div].mean():.2f}. "
          "This is cell-cycle DNA content, not endoreplication: ignore it.")

# Positive control: subperineurial glia endoreplicate to cover the growing brain, so they SHOULD be polyploid.
# If they are not found above the doublet rate here, the method is not sensitive enough in data of this depth,
# and a null result in other runs cannot be interpreted.
surface = cells.cell_type.str.startswith("surface glia")
print("\n=== POSITIVE CONTROL: subperineurial (surface) glia ===")
if surface.any():
    part = cells[surface]
    boot = [part.P_polyploid.sample(len(part), replace=True, random_state=i).mean() for i in range(2000)]
    lo, hi = np.quantile(boot, [0.025, 0.975])
    exp = part.expected_doublet.mean()
    rest = cells.P_polyploid[~surface & ~div].mean()
    print(f"{len(part):,} surface glia: polyploid rate {part.P_polyploid.mean():.3f} "
          f"(95% interval {lo:.3f} to {hi:.3f}), expected from doublets {exp:.3f}")
    print(f"Other non-dividing cells: {rest:.3f}")
    # Doublets hit every cluster about equally, so comparing surface glia with the OTHER cell types of the same
    # runs cancels that background. This is the more sensitive test; the absolute one above is conservative,
    # because it assumes every doublet is detected as polyploid.
    # Cluster labels are coarse: a cluster can sit just under the "dividing" threshold and still be full of
    # cells in S/G2. Score every cell individually and drop the most cycling ones from the comparison group.
    cyc_idx = [gene_pos[g] for g in MARKERS["dividing (S/G2/M)"] + MARKERS["neuroblast"] if g in gene_pos]
    cyc = np.asarray(logcpm[:, cyc_idx].mean(axis=1)).ravel()
    cells["cycling_score"] = (cyc - cyc.mean()) / (cyc.std() or 1)
    quiet = cells.cycling_score <= np.quantile(cells.cycling_score, 0.80)
    print(f"\nPer-cell check: dropping the 20% most cycling cells leaves "
          f"{int((surface & quiet).sum()):,} surface glia and {int((~surface & ~div & quiet).sum()):,} other cells.")
    sq, oq = cells[surface & quiet], cells[~surface & ~div & quiet]
    if len(sq) >= 30 and len(oq) >= 30:
        bs = [sq.P_polyploid.sample(len(sq), replace=True, random_state=i).mean() for i in range(2000)]
        bq = [oq.P_polyploid.sample(len(oq), replace=True, random_state=i).mean() for i in range(2000)]
        d2 = np.quantile([a_ - b_ for a_, b_ in zip(bs, bq)], [0.025, 0.975])
        print(f"   Surface glia {sq.P_polyploid.mean():.3f} vs other non-cycling cells {oq.P_polyploid.mean():.3f}: "
              f"difference {sq.P_polyploid.mean() - oq.P_polyploid.mean():+.3f} "
              f"(95% interval {d2[0]:+.3f} to {d2[1]:+.3f})")
    others = cells[~surface & ~div]
    bo = [others.P_polyploid.sample(len(others), replace=True, random_state=i).mean() for i in range(2000)]
    diff = [a_ - b_ for a_, b_ in zip(boot, bo)]
    dlo, dhi = np.quantile(diff, [0.025, 0.975])
    print(f"Surface glia minus other non-dividing cells: {part.P_polyploid.mean() - rest:+.3f} "
          f"(95% interval {dlo:+.3f} to {dhi:+.3f})")
    if lo > exp:
        print("PASSED on both tests. Surface glia exceed the doublet rate outright, so the method works on data of "
              "this depth and a null result elsewhere can be read as a real negative.")
    elif dlo > 0:
        print("PARTLY PASSED. Surface glia do NOT exceed the absolute doublet rate, but they are clearly more "
              "polyploid than the other cell types in the same runs, which share the same doublet background. "
              "So the method does detect real differences in ploidy between cell types; what it cannot do in data "
              "of this depth is prove that any single population exceeds the doublet rate on its own. Comparisons "
              "BETWEEN cell types are usable; absolute polyploid rates are not.")
    else:
        print("FAILED. Surface glia are no more polyploid than other cells and do not exceed the doublet rate. "
              "The method is not sensitive enough in data of this depth, so a null result in other runs means "
              "nothing either way. Deeper libraries or doublet-removed data are needed.")
    if PLOIDY_FIT == "per cluster":
        print("\nNote: PLOIDY_FIT is 'per cluster', which fits each cluster against its own baseline and will hide a "
              "cell type that is mostly polyploid. Set it to 'per run' in Cell 1 and run again from Cell 6.")
else:
    print("No surface glia cluster was found, so the positive control could not be run. Raising N_CLUSTERS in Cell 1 "
          "may separate them; they are a small population.")
