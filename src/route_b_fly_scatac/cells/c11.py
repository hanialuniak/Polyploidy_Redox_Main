#@title Cell 11 - Screen 4: gene by gene, including Geminin and dup
tested = list(dict.fromkeys([g for v in PROGRAMMES.values() for g in v] + GOI))
Z = np.column_stack([residual(col) for col in gene_z([gene_pos[g] for g in tested]).T])
Zs = (Z - Z.mean(0)) / np.where(Z.std(0) > 0, Z.std(0), 1)
xz = (x - x.mean()) / x.std(); xsz = (xs - xs.mean(0)) / xs.std(0)
rho = xz @ Zs / len(xz); null = xsz.T @ Zs / len(xz)
genes_result = pd.DataFrame({"gene": tested, "partial_rho": rho,
                             "p_perm": (1 + (np.abs(null) >= np.abs(rho)).sum(0)) / (null.shape[0] + 1)})
genes_result["q_fdr"] = fdr(genes_result.p_perm)
genes_result["group"] = [", ".join(n.split(": ")[1] for n, v in PROGRAMMES.items() if g in v) or "gene of interest"
                         for g in tested]
genes_result = genes_result.assign(_a=genes_result.partial_rho.abs()).sort_values(["p_perm", "_a"], ascending=[True, False]).drop(columns="_a")
display(genes_result.head(15).round(3))
display(genes_result[genes_result.gene.isin(GOI)].round(3))
