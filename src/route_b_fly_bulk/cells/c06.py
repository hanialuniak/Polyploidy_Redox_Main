#@title Cell 6 - Screens: does the optic lobe change differently with age from the central brain?
use = info.region.isin(["OL", "CB"]).to_numpy()
d = info[use].reset_index(drop=True)
age = np.log2(d.age_days.to_numpy() + 1)   # +1 because the first time point is day 0
ol = (d.region == "OL").to_numpy(float)
cols = [np.ones(len(d)), ol, age - age.mean(), (age - age.mean()) * ol]
if d.sex.nunique() > 1:
    cols.append((d.sex == "M").to_numpy(float))
D = np.column_stack(cols)
dof = len(d) - np.linalg.matrix_rank(D)
pinv = np.linalg.pinv(D)
c_var = np.linalg.pinv(D.T @ D)[3, 3]
from scipy import stats


def interaction_t(Y):
    """t statistic for (age slope in optic lobe) minus (age slope in central brain), per column of Y."""
    Y = np.asarray(Y, float)
    Y = Y[:, None] if Y.ndim == 1 else Y
    B = pinv @ Y
    resid = Y - D @ B
    s2 = (resid**2).sum(0) / dof
    return B[3] / np.sqrt(s2 * c_var), B[3]


Z = expressed.loc[:, np.flatnonzero(use)].T
Z = (Z - Z.mean()) / Z.std().replace(0, np.nan)
scores = pd.DataFrame({n: Z[v].mean(axis=1).to_numpy() for n, v in PROGRAMMES.items()})
t, b = interaction_t(scores.to_numpy())
result = pd.DataFrame({"programme": scores.columns, "genes": [len(PROGRAMMES[n]) for n in scores.columns],
                       "interaction": b, "t": t, "p": 2 * stats.t.sf(np.abs(t), dof)})
for fam in ("REDOX", "EXPLORE"):
    m = result.programme.str.startswith(fam)
    p = result.loc[m, "p"].to_numpy(); o = np.argsort(p)
    q = np.minimum(np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1], 1)
    qq = np.empty_like(q); qq[o] = q; result.loc[m, "q_fdr"] = qq

level = expressed.mean(axis=1)
bins = pd.qcut(level.rank(method="first"), 20, labels=False)
avoid = {g for v in PROGRAMMES.values() for g in v}
by_bin = {b_: [g for g in level.index[bins == b_] if g not in avoid] for b_ in range(20)}
for n in [n for n in scores.columns if not n.startswith("CONTROL")]:
    draws = np.column_stack([Z[[rng.choice(by_bin[bins[g]]) for g in PROGRAMMES[n]]].mean(axis=1) for _ in range(N_RANDOM_SETS)])
    null, _ = interaction_t(draws)
    i = result.programme == n; obs = result.loc[i, "t"].iloc[0]
    result.loc[i, "random_low"], result.loc[i, "random_high"] = np.quantile(null, [0.025, 0.975])
    result.loc[i, "p_vs_random"] = (1 + np.sum(np.abs(null) >= abs(obs))) / (len(null) + 1)
    worst, which = 0.0, ""
    for c in [c for c in scores.columns if c.startswith("CONTROL")]:
        tc, _ = interaction_t(scores[n].to_numpy() - scores[c].to_numpy())
        pc = 2 * stats.t.sf(abs(tc[0]), dof) if np.sign(tc[0]) == np.sign(obs) else 1.0
        if pc >= worst:
            worst, which = pc, c.replace("CONTROL: ", "")
    result.loc[i, "p_vs_controls"], result.loc[i, "closest_control"] = worst, which
ok = (result.q_fdr < FDR_MAX) & (result.p_vs_random < RANDOM_P_MAX) & (result.p_vs_controls < CONTROL_P_MAX)
result["verdict"] = np.where(result.programme.str.startswith("CONTROL"), "(control)",
                    np.where(ok, np.where(result.t > 0, "passes: rises more with age in optic lobe",
                                          "passes: falls more with age in optic lobe"),
                    np.where(result.q_fdr < FDR_MAX, "not specific", "no difference")))
print(f"{len(d)} samples (optic lobe and central brain), {dof} degrees of freedom")
display(result[["programme", "interaction", "t", "q_fdr", "p_vs_random", "p_vs_controls", "verdict"]].round(3))
print("\nFull table (random-set range and closest control):")
display(result.round(3))
