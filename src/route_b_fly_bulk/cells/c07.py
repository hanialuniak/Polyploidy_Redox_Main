#@title Cell 7 - Gene by gene, including Geminin, dup, Cdc6 and Orc1
genes = list(dict.fromkeys([g for v in PROGRAMMES.values() for g in v] + [g for g in ("geminin", "dup", "Cdc6", "Orc1") if g in expressed.index]))
tg, bg = interaction_t(Z[genes].to_numpy())
gene_result = pd.DataFrame({"gene": genes, "interaction": bg, "t": tg, "p": 2 * stats.t.sf(np.abs(tg), dof)})
p = gene_result.p.to_numpy(); o = np.argsort(p)
q = np.minimum(np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1], 1)
gene_result.loc[o, "q_fdr"] = q
gene_result["group"] = [", ".join(n.split(": ")[1] for n, v in PROGRAMMES.items() if g in v) for g in genes]
gene_result = gene_result.sort_values("p")
display(gene_result.head(15).round(3))
display(gene_result[gene_result.gene.isin(["geminin", "dup", "Cdc6", "Orc1"])].round(3))
