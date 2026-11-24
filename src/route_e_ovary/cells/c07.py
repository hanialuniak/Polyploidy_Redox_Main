#@title Cell 7 - Gene by gene, the graded stage trend, chart and saving
gl = list(dict.fromkeys(g for n, v in PROGRAMMES.items() if not n.startswith("CONTROL") for g in v if g in fc.index))
tab2 = pd.DataFrame({"log2_change": fc.loc[gl], "pct_of_all_genes": fc.rank(pct=True).loc[gl],
                     "mean_mitotic": Lp.xs(MIT, level=1, axis=1).mean(axis=1).loc[gl],
                     "mean_endocycling": Lp.xs(ENDO, level=1, axis=1).mean(axis=1).loc[gl]})
for l, d in per_lib.items():
    tab2[f"lib_{l[:14]}"] = d.loc[gl]
tab2["programme"] = [", ".join(n.split(": ")[1] for n, v in PROGRAMMES.items() if g in v) for g in gl]
display(tab2.sort_values("log2_change").round(2))

# Graded check: do the programmes change step by step across the follicle cell clusters, ordered by stage?
order = [cl for cl in summary.index if summary.loc[cl, "role"].startswith("follicle cell") or
         summary.loc[cl, "role"].endswith("follicle cell")]
order = sorted(order, key=lambda cl: stage[cl])
grade = {}
for n, gs in PROGRAMMES.items():
    s, _ = score(gs)
    grade[n.split(": ")[1]] = pd.Series(s).groupby(cells.cluster.to_numpy()).mean().reindex(order)
grade = pd.DataFrame(grade).T
grade.columns = [f"c{cl} ({summary.loc[cl, 'role'].replace(' follicle cell', '')}, stage {stage[cl]:+.1f})" for cl in order]
display(grade.round(2))
print("Follicle cell clusters left to right by stage score (earliest first). Values are mean log expression per cell.")

# Licensing genes one by one, in stage order: does Geminin change before or after the switch?
lic_name = next((n for n in PROGRAMMES if "licensing" in n), None)
if lic_name:
    lic = PROGRAMMES[lic_name]
    rows_l = {}
    for g in lic:
        v = np.asarray(Ln[:, gpos[g]].todense()).ravel()
        rows_l[g] = pd.Series(v).groupby(cells.cluster.to_numpy()).mean().reindex(order)
    lic_tab = pd.DataFrame(rows_l).T
    lic_tab.columns = grade.columns
    switch = [i for i, cl in enumerate(order) if summary.loc[cl, "role"] == ENDO]
    first_endo = min(switch) if switch else None
    print(f"\nLicensing genes across the follicle clusters, earliest stage first "
          f"({lic_name.split(': ')[1]}):")
    display(lic_tab.round(2))
    if first_endo is not None:
        print(f"The switch sits between column {first_endo - 1} and {first_endo}: columns 0-{first_endo - 1} are "
              f"mitotic clusters, {first_endo} onwards are endocycling.")

    print("\nEndocycling minus mitotic, per gene (same numbers as the table above, summarised):")
    lic_sum = pd.DataFrame({"log2_change": fc.reindex(lic),
                            "pct_of_all_genes": fc.rank(pct=True).reindex(lic)})
    for l, d in per_lib.items():
        lic_sum[f"lib_{l[:14]}"] = d.reindex(lic)
    libc = [c for c in lic_sum.columns if c.startswith("lib_")]
    lic_sum["same_direction"] = (np.sign(lic_sum[libc]).nunique(axis=1) == 1) if len(libc) > 1 else True
    display(lic_sum.sort_values("log2_change").round(2))
    for g in ("geminin", "Gem", "dup", "Cdt1"):
        if g in fc.index:
            d = fc[g]
            print(f"   {g}: {d:+.2f} log2 ({'lower' if d < 0 else 'higher'} in endocycling cells, "
                  f"percentile {fc.rank(pct=True)[g]:.2f} of all genes)")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for g in lic_tab.index:
        key = g in ("geminin", "Gem")
        ax.plot(range(lic_tab.shape[1]), lic_tab.loc[g], marker="o", ms=4,
                lw=2.5 if key else 1, color="#b2432f" if key else "#9aa7b4",
                label=g if key or lic_tab.loc[g].std() > 0.3 else None, zorder=3 if key else 1)
    if first_endo is not None:
        ax.axvline(first_endo - 0.5, color="black", ls="--", lw=1)
        ax.text(first_endo - 0.45, ax.get_ylim()[1], " switch", va="top", fontsize=8)
    ax.set_xticks(range(lic_tab.shape[1])); ax.set_xticklabels([c.split(" (")[0] for c in lic_tab.columns], fontsize=8)
    ax.set_xlabel("follicle cell clusters, earliest stage first"); ax.set_ylabel("mean log expression per cell")
    ax.set_title("Licensing genes across the mitotic-to-endocycle switch", loc="left", fontsize=10)
    ax.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(os.path.join(GDIR, "chart_licensing.png"), dpi=150); plt.show()
    lic_tab.to_csv(os.path.join(GDIR, "licensing_by_stage.csv"))

fig, ax = plt.subplots(figsize=(8, 0.34 * len(res) + 1.4))
col = ["#b2432f" if v < 0 else "#2f6db2" for v in (res.log2_change - res.random_median)]
ax.barh(res.programme, res.log2_change - res.random_median, color=col)
for i, r in res.reset_index().iterrows():
    ax.text(r.log2_change - r.random_median, i, f"  q={r.q_fdr:.3f}" if r.q_fdr == r.q_fdr else "  (control)",
            va="center", fontsize=8)
ax.axvline(0, color="black", lw=0.8); ax.invert_yaxis()
ax.set_xlabel("log2 change, endocycling vs mitotic follicle cells, minus random genes")
ax.set_title(f"{GSE}: programmes across the mitotic-to-endocycle switch", loc="left", fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(GDIR, "chart_programmes.png"), dpi=150); plt.show()
res.to_csv(os.path.join(GDIR, "programmes.csv"), index=False)
tab2.to_csv(os.path.join(GDIR, "genes.csv")); grade.to_csv(os.path.join(GDIR, "stage_trend.csv"))
cells[["library", "barcode", "cluster", "role", "n_counts", "n_genes"]].to_csv(os.path.join(GDIR, "cells.csv"), index=False)
print("Saved to", GDIR)
