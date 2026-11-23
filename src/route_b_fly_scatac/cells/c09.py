#@title Cell 9 - Scores and tools
lognorm = logcpm.tocsc()


def gene_z(idx):
    """Per-gene accessibility (log CPM), z-scored within cluster."""
    m = np.asarray(lognorm[:, idx].todense(), float)
    out = np.empty_like(m)
    for cl in cells.cluster.unique():
        r = (cells.cluster == cl).to_numpy()
        sd = m[r].std(0); out[r] = (m[r] - m[r].mean(0)) / np.where(sd > 0, sd, 1)
    return out


def programme_score(gs):
    return gene_z([gene_pos[g] for g in gs]).mean(1)


scores = pd.DataFrame({n: programme_score(v) for n, v in PROGRAMMES.items()})
lfpc = np.log(cells.fragments_per_copy.to_numpy()); lfpc = (lfpc - lfpc.mean()) / lfpc.std()
dummies = pd.get_dummies(cells.cluster + "|" + cells.run).astype(float).to_numpy()
basis, _ = np.linalg.qr(np.column_stack([dummies, lfpc, lfpc**2]))
strata = (cells.cluster + "|" + cells.run).to_numpy()


def residual(v):
    v = np.asarray(v, float); return v - basis @ (basis.T @ v)


def corr(a, b):
    a, b = a - a.mean(0), b - b.mean(0)
    return (a * b).sum(0) / np.sqrt((a**2).sum(0) * (b**2).sum(0))


def perm_p(obs, null):
    return (1 + np.sum(np.abs(null) >= abs(obs))) / (len(null) + 1)


def fdr(p):
    p = np.asarray(p, float); o = np.argsort(p); a = p[o] * len(p) / np.arange(1, len(p) + 1)
    out = np.empty_like(p); out[o] = np.minimum(np.minimum.accumulate(a[::-1])[::-1], 1); return out


x = residual(cells.P_polyploid.to_numpy())
xs = np.tile(cells.P_polyploid.to_numpy()[:, None], (1, N_PERMUTATIONS))
for s in np.unique(strata):
    idx = np.flatnonzero(strata == s)
    for j in range(N_PERMUTATIONS):
        xs[idx, j] = xs[rng.permutation(idx), j]
xs = residual(xs)
print("Tools ready.")
