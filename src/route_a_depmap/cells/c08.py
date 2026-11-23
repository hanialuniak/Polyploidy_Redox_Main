#@title Cell 8 - Screen 3: is it specific to redox, or general stress?
def test_controls(ctx):
    """For each redox programme and each stress control, test whether (redox minus control) still goes with
    ploidy. The worst (largest) p-value across the four controls is kept."""
    scores, x, xs, basis, result = ctx["scores"], ctx["x"], ctx["xs"], ctx["basis"], ctx["result"]
    controls = [c for c in scores.columns if c.startswith("CONTROL")]
    for name in [n for n in scores.columns if n.startswith("REDOX")]:
        i = result.programme == name
        rho = float(result.loc[i, "partial_rho"].iloc[0])
        worst_p, worst = 0.0, ""
        for control in controls:
            diff = zscore(scores[[name]])[name] - zscore(scores[[control]])[control]
            r = residual(diff.fillna(0).to_numpy(), basis)
            observed = float(corr(x[:, None], r[:, None])[0])
            p = perm_p(observed, corr(xs, np.repeat(r[:, None], xs.shape[1], axis=1)))
            if np.sign(observed) != np.sign(rho):
                p = 1.0
            if p >= worst_p:
                worst_p, worst = p, control.replace("CONTROL: ", "")
        result.loc[i, "p_vs_controls"], result.loc[i, "closest_control"] = worst_p, worst


test_controls(rna)
display(rna["result"][["programme", "partial_rho"]].round(3).join(
    rna["result"][["p_vs_controls", "closest_control"]]))
