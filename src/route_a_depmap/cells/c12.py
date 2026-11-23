#@title Cell 12 - Screen 7: Geminin and replication licensing
def partial_assoc(table, a, scores, extra=None, n_perm=N_PERMUTATIONS):
    """Partial rank correlation of vector a with each column of scores, within tissue, adjusted for
    proliferation (and optionally ploidy); p-value by shuffling a within tissue."""
    cols = [adjustment_basis(table)]
    if extra is not None:
        cols.append(np.asarray(extra, dtype=float)[:, None])
    basis, _ = np.linalg.qr(np.column_stack(cols))
    a_rank = pd.Series(a).rank().to_numpy()
    xa = residual(a_rank, basis)
    xsa = residual(shuffled_within_tissue(a_rank, table["tissue"].to_numpy(), n_perm), basis)
    rows = []
    for name in scores.columns:
        s = scores[name].to_numpy()
        r = residual(np.where(np.isnan(s), np.nanmean(s), s), basis)
        rho = float(corr(xa[:, None], r[:, None])[0])
        rows.append({"against": name, "partial_rho": rho,
                     "p_perm": perm_p(rho, corr(xsa, np.repeat(r[:, None], xsa.shape[1], axis=1)))})
    out = pd.DataFrame(rows)
    out["q_fdr"] = fdr(out["p_perm"])
    return out


def rho_only(table, A, scores, extra=None):
    """Partial correlations of many vectors (columns of A) with each score, no shuffling (fast)."""
    cols = [adjustment_basis(table)]
    if extra is not None:
        cols.append(np.asarray(extra, dtype=float)[:, None])
    basis, _ = np.linalg.qr(np.column_stack(cols))
    RA = residual(pd.DataFrame(A).rank().to_numpy(), basis)
    S = scores.to_numpy(float)
    RS = residual(np.where(np.isnan(S), np.nanmean(S, axis=0), S), basis)
    RA = (RA - RA.mean(0)) / RA.std(0)
    RS = (RS - RS.mean(0)) / RS.std(0)
    return RA.T @ RS / len(RA)


def matched_genes(gene, table, n=200):
    """Genes like `gene` in average expression and in how closely they follow proliferation."""
    data = expression.loc[table.index]
    avoid = set(g for v in PROGRAMMES.values() for g in v) | set(PROLIFERATION) | set(GOI)
    pool = [g for g in data.columns if g not in avoid and data[g].std() > 0]
    prolif = programme_scores(data, {"p": PROLIFERATION})[0]["p"].rank()
    level = data[pool + [gene]].mean()
    follow = data[pool + [gene]].rank().corrwith(prolif)
    lb = pd.qcut(level.rank(method="first"), 5, labels=False)
    fb = pd.qcut(follow.rank(method="first"), 5, labels=False)
    same = [g for g in pool if lb[g] == lb[gene] and fb[g] == fb[gene]]
    return list(rng.choice(same, min(n, len(same)), replace=False))


table = rna["table"]
measures = {}
for g in GOI:
    if g in expression.columns:
        measures[f"{g} RNA"] = expression.loc[table.index, g].to_numpy(float)
    if crispr is not None and g in crispr.columns:
        measures[f"{g} knockout effect"] = crispr[g].reindex(table.index).to_numpy(float)
if {"GMNN", "CDT1"} <= set(expression.columns):
    measures["GMNN minus CDT1 (RNA)"] = (expression.loc[table.index, "GMNN"] - expression.loc[table.index, "CDT1"]).to_numpy(float)

redox_scores = rna["scores"][[c for c in rna["scores"].columns if c.startswith("REDOX")]]
geminin = []
for label, values in measures.items():
    ok = ~np.isnan(values)
    if ok.sum() < 25:
        continue
    sub = table[ok]
    ploidy_scores = pd.DataFrame({"ploidy": sub["ploidy"].to_numpy(), "whole-genome doubling": sub["wgd"].fillna(0).to_numpy()})
    a = partial_assoc(sub, values[ok], ploidy_scores).assign(measure=label, adjusted_for_ploidy="-")
    b = partial_assoc(sub, values[ok], redox_scores.loc[sub.index].reset_index(drop=True)).assign(measure=label, adjusted_for_ploidy="no")
    c = partial_assoc(sub, values[ok], redox_scores.loc[sub.index].reset_index(drop=True),
                      extra=sub["ploidy"].rank().to_numpy()).assign(measure=label, adjusted_for_ploidy="yes")
    part = pd.concat([a, b, c], ignore_index=True).assign(lines=int(ok.sum()))
    gene = label.split()[0]
    if label.endswith(" RNA") and gene in expression.columns:
        others = matched_genes(gene, sub)
        A = expression.loc[sub.index, others].to_numpy(float)
        targets = pd.concat([ploidy_scores, redox_scores.loc[sub.index].reset_index(drop=True)], axis=1)
        null_no = rho_only(sub, A, targets)
        null_yes = rho_only(sub, A, redox_scores.loc[sub.index].reset_index(drop=True), sub["ploidy"].rank().to_numpy())
        for i, row in part.reset_index(drop=True).iterrows():
            null = (null_yes[:, list(redox_scores.columns).index(row["against"])] if row["adjusted_for_ploidy"] == "yes"
                    else null_no[:, list(targets.columns).index(row["against"])])
            part.loc[part.index[i], "matched_low"] = np.quantile(null, 0.025)
            part.loc[part.index[i], "matched_high"] = np.quantile(null, 0.975)
            part.loc[part.index[i], "p_vs_matched_genes"] = perm_p(row["partial_rho"], null)
    geminin.append(part)
geminin = pd.concat(geminin, ignore_index=True) if geminin else pd.DataFrame()
for c in ("matched_low", "matched_high", "p_vs_matched_genes"):
    if not geminin.empty and c not in geminin.columns:
        geminin[c] = np.nan
if geminin.empty:
    print("None of the genes of interest were found.")
else:
    for label in measures:
        part = geminin[geminin.measure == label]
        if part.empty:
            continue
        print(f"\n{label} ({part.lines.iloc[0]} lines)")
        cols = ["against", "adjusted_for_ploidy", "partial_rho", "p_perm", "q_fdr"]
        cols += [c for c in ("matched_low", "matched_high", "p_vs_matched_genes") if part[c].notna().any()]
        display(part[cols].round(3).reset_index(drop=True))
