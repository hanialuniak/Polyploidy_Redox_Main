#@title Cell 10 - Screens 1-3: programmes vs ploidy, random genes, stress controls
rows = []
for n in scores.columns:
    r = residual(scores[n].to_numpy())
    rho = float(corr(x[:, None], r[:, None])[0])
    rows.append({"programme": n, "genes": len(PROGRAMMES[n]), "partial_rho": rho,
                 "p_perm": perm_p(rho, corr(xs, np.repeat(r[:, None], xs.shape[1], 1)))})
result = pd.DataFrame(rows)
for family in ("REDOX", "EXPLORE"):
    m = result.programme.str.startswith(family)
    if m.any():
        result.loc[m, "q_fdr"] = fdr(result.loc[m, "p_perm"])

level = np.asarray(lognorm.mean(axis=0)).ravel()
ok = np.flatnonzero(level > 0)
bins = pd.qcut(pd.Series(level[ok]).rank(method="first"), 20, labels=False).to_numpy()
bin_of = dict(zip(ok, bins))
avoid = {gene_pos[g] for v in PROGRAMMES.values() for g in v}
by_bin = {b: [g for g, bb in zip(ok, bins) if bb == b and g not in avoid] for b in range(20)}
for n in [n for n in scores.columns if n.startswith(("REDOX", "EXPLORE"))]:
    idx = [gene_pos[g] for g in PROGRAMMES[n]]
    null = [corr(x[:, None], residual(gene_z([rng.choice(by_bin[bin_of.get(i, 0)]) for i in idx]).mean(1))[:, None])[0]
            for _ in range(N_RANDOM_SETS)]
    i = result.programme == n
    result.loc[i, "random_low"], result.loc[i, "random_high"] = np.quantile(null, [0.025, 0.975])
    result.loc[i, "p_vs_random"] = perm_p(result.loc[i, "partial_rho"].iloc[0], np.array(null))

controls = [n for n in scores.columns if n.startswith("CONTROL")]
for n in [n for n in scores.columns if n.startswith(("REDOX", "EXPLORE"))]:
    i = result.programme == n; rho = result.loc[i, "partial_rho"].iloc[0]; worst, which = 0.0, ""
    for c in controls:
        r = residual(scores[n].to_numpy() - scores[c].to_numpy())
        obs = float(corr(x[:, None], r[:, None])[0])
        p = 1.0 if np.sign(obs) != np.sign(rho) else perm_p(obs, corr(xs, np.repeat(r[:, None], xs.shape[1], 1)))
        if p >= worst:
            worst, which = p, c.replace("CONTROL: ", "")
    result.loc[i, "p_vs_controls"], result.loc[i, "closest_control"] = worst, which
display(result.round(3))
passes = (result.q_fdr < FDR_MAX) & (result.p_vs_random < RANDOM_P_MAX) & (result.p_vs_controls < CONTROL_P_MAX)
result["verdict"] = np.where(result.programme.str.startswith("CONTROL"), "(control)",
                    np.where(passes, np.where(result.partial_rho > 0, "passes: more open in polyploid cells",
                                              "passes: less open in polyploid cells"),
                    np.where(result.q_fdr < FDR_MAX, "not specific", "no link")))
print()
for _, r in result[~result.programme.str.startswith("CONTROL")].iterrows():
    print(f"{r.programme}: {r.verdict} (partial rho {r.partial_rho:+.3f})")
