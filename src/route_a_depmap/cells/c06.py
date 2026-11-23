#@title Cell 6 - Screen 1: does each programme go with ploidy, within tissue? (RNA)
def test_link(data, label, exposure="ploidy"):
    """Partial rank correlation of programme score with ploidy, after removing tissue and proliferation.
    The p-value comes from shuffling ploidy within each tissue N_PERMUTATIONS times."""
    table = lines[lines.index.isin(data.index)].dropna(subset=[exposure])
    data = data.loc[table.index]
    basis = adjustment_basis(table)
    x = residual(table[exposure].rank().to_numpy(), basis)
    xs = residual(shuffled_within_tissue(table[exposure].rank().to_numpy(), table["tissue"].to_numpy(),
                                         N_PERMUTATIONS), basis)
    scores, used = programme_scores(data, PROGRAMMES)
    rows = []
    for name in scores.columns:
        s = scores[name].to_numpy()
        ok = ~np.isnan(s)
        r = residual(np.where(ok, s, np.nanmean(s)), basis)
        rho = float(corr(x[ok, None], r[ok, None])[0])
        null = corr(xs[ok], np.repeat(r[ok, None], xs.shape[1], axis=1))
        rows.append({"data": label, "programme": name, "lines": int(ok.sum()), "genes_used": len(used[name]),
                     "partial_rho": rho, "p_perm": perm_p(rho, null)})
    result = pd.DataFrame(rows)
    redox = result.programme.str.startswith("REDOX")
    result.loc[redox, "q_fdr"] = fdr(result.loc[redox, "p_perm"])
    return {"result": result, "scores": scores, "used": used, "x": x, "xs": xs, "basis": basis,
            "table": table, "data": data}


rna = test_link(expression, "RNA")
display(rna["result"].round(3))
