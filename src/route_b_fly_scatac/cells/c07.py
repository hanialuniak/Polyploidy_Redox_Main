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
surface = cells.cell_type.str.startswith("surface glia")
if surface.any():
    print(f"Bonus check: surface glia polyploid rate {cells.P_polyploid[surface].mean():.2f} "
          f"(expected high), other cells {cells.P_polyploid[~surface].mean():.2f}")
else:
    print("No surface glia cluster (expected: they are very few). The doublet check above is the main test.")
