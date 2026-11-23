#@title Cell 8 - Charts and saving
show = [n for n in scores.columns if not n.startswith("CONTROL")][:7]
fig, axes = plt.subplots(1, len(show), figsize=(3 * len(show), 3), sharey=True)
for ax, n in zip(np.atleast_1d(axes), show):
    for reg, col in (("OL", "#c26a2e"), ("CB", "#3b6ea8")):
        m = (d.region == reg).to_numpy()
        means = pd.Series(scores[n].to_numpy()[m]).groupby(d.age_days[m].to_numpy()).mean()
        ax.plot(means.index, means.values, "o-", color=col, label=reg)
    ax.set_xscale("symlog", base=2, linthresh=1); ax.set_title(n.split(": ")[1], fontsize=9); ax.set_xlabel("age (days)")
np.atleast_1d(axes)[0].set_ylabel("programme score"); np.atleast_1d(axes)[0].legend(frameon=False)
plt.tight_layout(); plt.savefig(os.path.join(DATA_DIR, "chart_programmes_by_age.png"), dpi=150); plt.show()
for name, frame in {"bulk_programmes.csv": result, "bulk_genes.csv": gene_result, "bulk_scores.csv": scores.assign(**d[["region", "age_days", "sex"]])}.items():
    frame.to_csv(os.path.join(DATA_DIR, name), index=False)
print("Saved to", DATA_DIR)
