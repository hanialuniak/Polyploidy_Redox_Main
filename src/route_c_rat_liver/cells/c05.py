#@title Cell 5 - Redox enzyme sets vs ploidy, against random genes, with and without zonation adjustment
SETS = {
 "Glutathione S-transferases": "Gsta1 Gsta2 Gsta3 Gsta4 Gsta5 Gstm1 Gstm2 Gstm3 Gstm4 Gstm5 Gstm6 Gstm7 Gstp1 Gstk1 Gsto1 Gstt1 Gstt2 Gstt3 Gstz1 Mgst1 Mgst2 Mgst3",
 "Glutathione synthesis & recycling": "Gclc Gclm Gss Gsr Gpx1 Gpx4 Glrx Glrx2 Glrx3 Glrx5",
 "Thioredoxin / peroxiredoxin": "Txn1 Txnrd1 Txn2 Txnrd2 Prdx1 Prdx2 Prdx3 Prdx4 Prdx5 Prdx6 Srxn1 Txnl1",
 "H2O2 / ROS-producing enzymes": "Cyp2e1 Xdh Maoa Maob Ero1a Ero1b Acox1 Acox2 Acox3 Hao1 Hao2 Dao Aox1 Nox4",
 "NADPH supply": "G6pd Pgd Tkt Taldo1 Me1 Me2 Idh1 Idh2 Nnt Mthfd1 Aldh1l1",
 "NRF2 targets (other)": "Nqo1 Hmox1 Ephx1 Cbr1 Cbr3 Akr7a3 Akr1b8 Ftl1 Fth1 Blvrb Ugdh",
 "Superoxide dismutase & catalase": "Sod1 Sod2 Sod3 Cat",
 "CONTROL: cell cycle": " ".join(CHECKS["cell cycle (positive control: up in 8c, reported by the authors)"]),
}
present = set(e.index)
SETS = {k: [g for g in v.split() if g in present] for k, v in SETS.items()}
for k, v in SETS.items():
    print(f"{k}: {len(v)} genes found")


def slopes(X):
    """Per-gene coefficient of log2(ploidy): change per doubling of ploidy."""
    Y = e.to_numpy().T
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    dof = max(len(Y) - X.shape[1], 1)
    se = np.sqrt((resid**2).sum(0) / dof * np.linalg.inv(X.T @ X)[0, 0])
    return pd.Series(beta[0], index=e.index), pd.Series(beta[0] / np.where(se > 0, se, np.nan), index=e.index)


block = pd.get_dummies(samples.animal).astype(float).to_numpy() if PAIRED else np.ones((len(lp), 1))
X0 = np.column_stack([lp, block])                     # ploidy first, then one baseline per animal
zc = (zone - zone.mean()) / (zone.std() or 1)
X1 = np.column_stack([X0, zc])
models = {"ploidy" + (", within animal" if PAIRED else ""): X0}
if np.linalg.matrix_rank(X1) == X1.shape[1] and len(lp) > X1.shape[1] + 2:
    models["ploidy, zonation-adjusted" + (", within animal" if PAIRED else "")] = X1
level = e.mean(axis=1)
bins = pd.qcut(level.rank(method="first"), 10, labels=False)
avoid = {g for v in SETS.values() for g in v}
pools = {b: [g for g in bins.index[bins == b] if g not in avoid] for b in range(10)}
rows = []
for mname, X in models.items():
    b, t = slopes(X)
    for name, gs in SETS.items():
        if len(gs) < 3:
            continue
        obs = b.loc[gs].median()
        null = np.array([b.loc[[rng.choice(pools[bins[g]]) for g in gs]].median() for _ in range(N_RANDOM_SETS)])
        c = np.median(null)
        p = (1 + np.sum(np.abs(null - c) >= abs(obs - c))) / (N_RANDOM_SETS + 1)
        rows.append({"model": mname, "set": name, "genes": len(gs), "change_per_doubling": obs,
                     "random_median": c, "p_vs_random": p, "share_down": float((b.loc[gs] < 0).mean())})
    if mname.startswith("ploidy") and "zonation" not in mname:
        gene_b, gene_t = b, t
res = pd.DataFrame(rows)
for mname in res.model.unique():
    m = (res.model == mname) & ~res.set.str.startswith("CONTROL")
    p = res.loc[m, "p_vs_random"].to_numpy(); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    qq = np.empty_like(p); qq[o] = np.minimum(q, 1); res.loc[m, "q_fdr"] = qq
pd.set_option("display.width", 200); pd.set_option("display.max_columns", 12)
display(res.round(3))
if any("zonation" in m for m in models) and abs(r_zone) > 0.9:
    print(f"NOTE: zonation and ploidy are almost perfectly linked here (r = {r_zone:+.2f}). The adjusted model then cannot tell "
          "them apart, so a real ploidy effect can vanish after adjustment. Read the adjusted rows as 'not separable', not as 'no effect'.")
print("change_per_doubling: median log2 change of the set's genes per doubling of ploidy (2c -> 4c -> 8c).")
print(f"All genes: {(gene_b < 0).mean():.0%} go down with ploidy.")
