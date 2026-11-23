#@title Cell 5 - Scoring and adjustment tools (no results yet)
def zscore(frame):
    return (frame - frame.mean()) / frame.std(ddof=0).replace(0, np.nan)


def programme_scores(data, programmes, min_genes=5):
    """Each programme's score = average z-score of its genes measured in at least 80% of lines."""
    z = zscore(data[data.columns[data.notna().mean() >= 0.8]])
    scores, used = {}, {}
    for name, genes in programmes.items():
        found = [g for g in genes if g in z.columns]
        if len(found) >= min_genes:
            scores[name], used[name] = z[found].mean(axis=1, skipna=True), found
    return pd.DataFrame(scores), used


def adjustment_basis(table):
    """What every test removes first: tissue of origin, and (per Cell 1) proliferation and TP53 status."""
    columns = [pd.get_dummies(table["tissue"]).astype(float).to_numpy()]
    if ADJUST_PROLIFERATION:
        prolif = programme_scores(expression.loc[table.index], {"p": PROLIFERATION})[0]["p"].rank().to_numpy()
        columns.append(((prolif - prolif.mean()) / prolif.std())[:, None])
    if TP53_MODE == "Adjust for TP53 status" and table["tp53_mutant"].nunique() > 1:
        columns.append(table["tp53_mutant"].astype(float).to_numpy()[:, None])
    q, _ = np.linalg.qr(np.column_stack(columns))
    return q


def residual(values, basis):
    values = np.asarray(values, dtype=float)
    return values - basis @ (basis.T @ values)


def corr(a, b):
    a, b = a - a.mean(axis=0), b - b.mean(axis=0)
    return (a * b).sum(axis=0) / np.sqrt((a**2).sum(axis=0) * (b**2).sum(axis=0))


def shuffled_within_tissue(values, groups, n):
    out = np.tile(np.asarray(values, dtype=float)[:, None], (1, n))
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        for j in range(n):
            out[idx, j] = out[rng.permutation(idx), j]
    return out


def perm_p(observed, null):
    return (1 + np.sum(np.abs(null) >= abs(observed))) / (len(null) + 1)


def fdr(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
    out = np.empty_like(p)
    out[order] = np.minimum(np.minimum.accumulate(adjusted[::-1])[::-1], 1)
    return out
print("Tools ready.")
