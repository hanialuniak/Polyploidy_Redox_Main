#@title Cell 6 - Gene by gene, compared with the mouse hepatocyte proteins, chart and saving
genes = [g for k, v in SETS.items() if not k.startswith("CONTROL") for g in v]
tab = pd.DataFrame({"change_per_doubling": gene_b.loc[genes], "t": gene_t.loc[genes],
                    "pct_of_all_genes": gene_b.rank(pct=True).loc[genes]})
grp = samples.groupby("ploidy").gsm.apply(list)
for p_, cols in grp.items():
    tab[f"mean_{int(p_)}c"] = e.loc[genes, cols].mean(axis=1)
tab["set"] = [next(k for k, v in SETS.items() if g in v) for g in genes]
pd.set_option("display.max_rows", 200)
display(tab.sort_values("t").round(2))

mouse = os.path.join(DATA_DIR, "redox_enzymes_per_protein.csv")
if os.path.exists(mouse):
    mo = pd.read_csv(mouse).set_index("gene")
    both = tab.join(mo[["ploidy_same_dna_conc", "ploidy_same_dna_conc_large"]], how="inner")
    if len(both) >= 5:
        r = both.change_per_doubling.corr(both.ploidy_same_dna_conc)
        agree = int(((both.change_per_doubling < 0) & (both.ploidy_same_dna_conc < 0)).sum())
        print(f"\nRat mRNA vs mouse protein (Lanz, 4n vs 2n small cells), {len(both)} shared enzymes: r = {r:+.2f}; "
              f"lower with ploidy in both: {agree}")
else:
    print("\n(To compare with the mouse hepatocyte proteins, copy redox_enzymes_per_protein.csv from route_c_liver into this folder.)")

fig, ax = plt.subplots(figsize=(7, 4))
main = res[~res.model.str.contains("zonation")].set_index("set")
ax.barh(main.index, main.change_per_doubling, color=["#b2432f" if v < 0 else "#2f6db2" for v in main.change_per_doubling])
for i, (s, r_) in enumerate(main.iterrows()):
    ax.text(r_.change_per_doubling, i, f"  p={r_.p_vs_random:.3f}", va="center", fontsize=8)
ax.axvline(0, color="black", lw=0.8); ax.invert_yaxis()
ax.set_xlabel("Median log2 change per doubling of ploidy"); ax.set_title(f"{GSE}: redox enzyme genes vs hepatocyte ploidy", loc="left", fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(DATA_DIR, "chart_rat_redox_sets.png"), dpi=150); plt.show()
res.to_csv(os.path.join(DATA_DIR, "rat_redox_sets.csv"), index=False)
tab.to_csv(os.path.join(DATA_DIR, "rat_redox_genes.csv"))
print("Saved to", DATA_DIR)
