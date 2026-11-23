#@title Cell 14 - Charts
plot = results[results.data.isin(["RNA", "Protein"])]
names = list(dict.fromkeys(plot.programme))
fig, ax = plt.subplots(figsize=(8, 0.45 * len(names) + 1.5))
for label, offset, colour in (("RNA", -0.18, "#3b6ea8"), ("Protein", 0.18, "#c26a2e")):
    part = plot[plot.data == label]
    if part.empty:
        continue
    y = np.array([names.index(n) for n in part.programme]) + offset
    ax.scatter(part.partial_rho, y, color=colour, label=label, zorder=3)
    has = part.random_low.notna().to_numpy()
    ax.hlines(y[has], part.random_low[has], part.random_high[has], color="#999999", lw=4, alpha=0.5)
ax.axvline(0, color="#444444", lw=0.8)
ax.set_yticks(range(len(names)))
ax.set_yticklabels([n.replace("REDOX: ", "").replace("CONTROL: ", "control: ") for n in names])
ax.invert_yaxis()
ax.set_xlabel("Partial correlation with ploidy (within tissue, adjusted for proliferation)")
ax.set_title("Programmes vs ploidy. Grey bars: range reached by random gene sets", loc="left", fontsize=11)
ax.legend(frameon=False, loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(DATA_DIR, "chart_programmes.png"), dpi=150)
plt.show()

fig, ax = plt.subplots(figsize=(7, 5))
part = gene_results[gene_results.data == "RNA"]
hit = part.q_fdr < FDR_MAX
ax.scatter(part.partial_rho[~hit], -np.log10(part.q_fdr[~hit]), s=12, color="#bbbbbb")
ax.scatter(part.partial_rho[hit], -np.log10(part.q_fdr[hit]), s=16, color="#3b6ea8")
for _, row in part.head(12).iterrows():
    ax.annotate(row["gene"], (row["partial_rho"], -np.log10(row["q_fdr"])), fontsize=8, xytext=(3, 2),
                textcoords="offset points")
ax.axhline(-np.log10(FDR_MAX), color="#999999", lw=0.8, ls="--")
ax.set_xlabel("Partial correlation with ploidy (RNA)")
ax.set_ylabel("-log10 false discovery rate")
ax.set_title("Gene by gene (Screen 6). Blue: false discovery rate below 0.05", loc="left", fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(DATA_DIR, "chart_genes.png"), dpi=150)
plt.show()
