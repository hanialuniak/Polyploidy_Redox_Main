#@title Cell 7 - Screen 2: is it stronger than random genes?
def test_random(ctx):
    """For each redox programme, draw N_RANDOM_SETS random gene sets of the same size, each gene matched
    on average expression, and see how strongly they correlate with ploidy."""
    data, used, x, basis, result = ctx["data"], ctx["used"], ctx["x"], ctx["basis"], ctx["result"]
    usable = data.columns[data.notna().mean() >= 0.8]
    level = data[usable].mean()
    bins = pd.qcut(level.rank(method="first"), 20, labels=False)
    avoid = set().union(*used.values()) | set(PROLIFERATION)
    pool = [g for g in usable if g not in avoid]
    by_bin = {b: [g for g in level.index[bins == b] if g not in avoid] or pool for b in range(20)}
    z = zscore(data[usable])
    for name, genes in used.items():
        if not name.startswith("REDOX"):
            continue
        null = []
        for _ in range(N_RANDOM_SETS):
            draw = [rng.choice(by_bin[int(bins[g])]) for g in genes]
            null.append(corr(x[:, None], residual(z[draw].mean(axis=1).fillna(0).to_numpy(), basis)[:, None])[0])
        i = result.programme == name
        rho = float(result.loc[i, "partial_rho"].iloc[0])
        result.loc[i, "random_low"], result.loc[i, "random_high"] = np.quantile(null, [0.025, 0.975])
        result.loc[i, "p_vs_random"] = perm_p(rho, np.array(null))


test_random(rna)
display(rna["result"].loc[rna["result"].programme.str.startswith("REDOX"),
                          ["programme", "partial_rho", "random_low", "random_high", "p_vs_random"]].round(3))
