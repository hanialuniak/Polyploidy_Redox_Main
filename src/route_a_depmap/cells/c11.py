#@title Cell 11 - Screen 6: gene by gene
def gene_screen(ctx, label):
    """Every redox and control gene, plus GENES_OF_INTEREST, tested on its own against ploidy (same adjustment
    and shuffles as Screen 1). False discovery rate is computed across all genes tested."""
    data, x, xs, basis = ctx["data"], ctx["x"], ctx["xs"], ctx["basis"]
    genes = list(dict.fromkeys([g for v in PROGRAMMES.values() for g in v] + GOI))
    genes = [g for g in genes if g in data.columns and data[g].notna().mean() >= 0.8]
    z = zscore(data[genes])
    R = residual(z.fillna(z.mean()).fillna(0).to_numpy(), basis)
    R = (R - R.mean(0)) / np.where(R.std(0) > 0, R.std(0), 1)
    xz = (x - x.mean()) / x.std()
    xsz = (xs - xs.mean(0)) / xs.std(0)
    rho = xz @ R / len(xz)
    null = xsz.T @ R / len(xz)
    p = (1 + (np.abs(null) >= np.abs(rho)).sum(0)) / (null.shape[0] + 1)
    member = {g: ", ".join(n.split(": ")[1] for n, v in PROGRAMMES.items() if g in v) for g in genes}
    out = pd.DataFrame({"data": label, "gene": genes, "partial_rho": rho, "p_perm": p,
                        "programmes": [member[g] or "gene of interest" for g in genes]})
    out["q_fdr"] = fdr(out["p_perm"])
    return out.assign(_abs=out.partial_rho.abs()).sort_values(["p_perm", "_abs"], ascending=[True, False]).drop(columns="_abs")


gene_results = [gene_screen(rna, "RNA")]
if prot is not None:
    gene_results.append(gene_screen(prot, "Protein"))
gene_results = pd.concat(gene_results, ignore_index=True)
both = gene_results.pivot_table(index="gene", columns="data", values=["partial_rho", "q_fdr"])
print("Top 15 genes in RNA:")
display(gene_results[gene_results.data == "RNA"].head(15).round(3))
if prot is not None:
    agree = both[(both[("q_fdr", "RNA")] < FDR_MAX) & (both[("q_fdr", "Protein")] < FDR_MAX) &
                 (np.sign(both[("partial_rho", "RNA")]) == np.sign(both[("partial_rho", "Protein")]))]
    print(f"\nGenes significant in both RNA and protein, same direction: {len(agree)}")
    display(agree.round(3))
