#@title Cell 6 - Gene by gene, chart and saving
genes = list(dict.fromkeys(g for n, v in PROGRAMMES.items() if not n.startswith("CONTROL") for g in v))
tab = pd.DataFrame({f"{t} vs {r}": (gmean[t] - gmean[r]).loc[genes] for t, r in CONTRASTS})
tab["mean_log2_change"] = tab.mean(axis=1)
tab["same_direction_all"] = (np.sign(tab.iloc[:, :len(CONTRASTS)]).nunique(axis=1) == 1)
tab["programme"] = [", ".join(n.split(": ")[1] for n, v in PROGRAMMES.items() if g in v) for g in genes]
for g in info.group.unique():
    tab[f"mean_{g}"] = gmean[g].loc[genes]
display(tab.sort_values("mean_log2_change").round(2))

fig, ax = plt.subplots(figsize=(8, 0.35 * len(PROGRAMMES) + 1.5))
names = list(PROGRAMMES)
for k, ct in enumerate(res.contrast.unique()):
    s = res[res.contrast == ct].set_index("programme").loc[names]
    y = np.arange(len(names)) + (k - (len(CONTRASTS) - 1) / 2) * 0.25
    ax.scatter(s.log2_change - s.random_median, y, label=ct, s=40)
ax.axvline(0, color="black", lw=0.8); ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8)
ax.invert_yaxis(); ax.set_xlabel("log2 change, endocycling vs dividing, minus random genes")
ax.set_title(f"{GSE}: programmes in endocycling vs dividing cells", loc="left", fontsize=10); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(GDIR, "chart_programmes.png"), dpi=150); plt.show()
res.to_csv(os.path.join(GDIR, "programmes.csv"), index=False); tab.to_csv(os.path.join(GDIR, "genes.csv"))
print("Saved to", GDIR)
