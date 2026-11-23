#@title Cell 12 - Per cell type, charts and saving
per_type = []
for ct, part in cells.groupby("cell_type"):
    m = (cells.cell_type == ct).to_numpy()
    if m.sum() < 100 or cells.P_polyploid[m].std() < 0.02:
        continue
    for n in scores.columns:
        per_type.append({"cell_type": ct, "programme": n, "cells": int(m.sum()),
                         "partial_rho": float(corr(x[m, None], residual(scores[n].to_numpy())[m, None])[0])})
per_type = pd.DataFrame(per_type)
if not per_type.empty:
    display(per_type.pivot(index="programme", columns="cell_type", values="partial_rho").round(2))

fig, ax = plt.subplots(figsize=(7, 4.5))
order = cells.groupby("cell_type").P_polyploid.mean().sort_values()
ax.barh(order.index, order.values, color=["#c26a2e" if s.startswith("surface") else "#3b6ea8" for s in order.index])
ax.set_xlabel("Mean probability of being polyploid"); ax.set_title("Polyploidy by cell type", loc="left")
plt.tight_layout(); plt.savefig(os.path.join(DATA_DIR, "chart_ploidy_by_type.png"), dpi=150); plt.show()

for name, frame in {"route_b_cells.csv": cells, "route_b_programmes.csv": result, "route_b_genes.csv": genes_result,
                    "route_b_per_type.csv": per_type}.items():
    frame.to_csv(os.path.join(DATA_DIR, name), index=False)
print("Saved results to", DATA_DIR)
