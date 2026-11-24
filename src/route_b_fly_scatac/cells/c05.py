#@title Cell 5 - Group cells into clusters and name them from marker genes
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans

keep = (cells.n_fragments >= MIN_FRAGMENTS) & (cells.n_fragments <= MAX_FRAGMENTS)
cells = cells[keep.to_numpy()].reset_index(drop=True)
overlaps = overlaps.loc[cells.cell_id]
gene_counts = gene_counts[np.flatnonzero(keep.to_numpy())]
print(f"Cells with {MIN_FRAGMENTS:,}-{MAX_FRAGMENTS:,} autosomal fragments: {len(cells):,}")

detected = np.asarray((gene_counts > 0).mean(axis=0)).ravel()
use = np.flatnonzero(detected > 0.01)
tf = sp.diags(1 / np.maximum(np.asarray(gene_counts[:, use].sum(axis=1)).ravel(), 1)) @ gene_counts[:, use]
tfidf = tf.multiply(np.log1p(len(cells) / np.maximum(detected[use] * len(cells), 1))).tocsr()
tfidf.data = np.log1p(tfidf.data * 1e4)
svd = TruncatedSVD(n_components=30, random_state=SEED).fit_transform(tfidf)[:, 1:]
svd = svd / np.linalg.norm(svd, axis=1, keepdims=True)
cells["cluster"] = KMeans(n_clusters=N_CLUSTERS, n_init=5, random_state=SEED).fit_predict(svd).astype(str)

MARKERS = {"neuron": ["elav", "nSyb", "brp", "fne"], "glia": ["repo"],
           "surface glia": ["moody", "Mdr65", "Gli", "Nrx-IV"],          # subperineurial glia: known polyploid
           "astrocyte": ["alrm", "Eaat1"], "ensheathing glia": ["Eaat2"], "cortex glia": ["wrapper"],
           "neuroblast": ["dpn", "ase", "mira", "wor", "insc", "pros"],  # larval and pupal brains only
           "dividing (S/G2/M)": ["PCNA", "RnrS", "RnrL", "Mcm2", "Mcm5", "dpa", "CycB", "polo", "stg", "aurB"]}
gene_pos = {g: i for i, g in enumerate(genes.gene)}

# Polyploid cells have more fragments. Thin every cell to the same number of cut sites in genes
# so that accessibility scores do not rise with depth alone.
in_genes = np.asarray(gene_counts.sum(axis=1)).ravel()
TARGET = int(np.quantile(in_genes, 0.1))
thin = gene_counts.tocsr(copy=True).astype(np.int64)
prob = np.repeat(TARGET / np.maximum(in_genes, 1), np.diff(thin.indptr))
thin.data = rng.binomial(thin.data, np.minimum(prob, 1.0)); thin.eliminate_zeros()
enough = in_genes >= TARGET
print(f"Thinned every cell to {TARGET:,} cut sites in genes; {int((~enough).sum()):,} cells below that are left out")
cells = cells[enough].reset_index(drop=True)
overlaps = overlaps.loc[cells.cell_id]
gene_counts, thin, svd = gene_counts[enough], thin[enough], svd[enough]
logcpm = thin.astype(float); logcpm.data = np.log1p(logcpm.data)
marker_table = {}
for label, ms in MARKERS.items():
    idx = [gene_pos[m] for m in ms if m in gene_pos]
    if idx:
        v = np.asarray(logcpm[:, idx].mean(axis=1)).ravel()
        marker_table[label] = pd.Series(v).groupby(cells.cluster).mean()
marker_table = pd.DataFrame(marker_table)
z = (marker_table - marker_table.mean()) / marker_table.std()
names = {}
for cl, row in z.iterrows():
    if row.get("neuroblast", -9) > 1.2 or row.get("dividing (S/G2/M)", -9) > 1.5:
        names[cl] = "dividing cells"          # neuroblasts and their daughters: 2-4C because they are in cycle
    elif row.get("surface glia", -9) > 1.5:
        names[cl] = "surface glia"
    elif row.get("glia", -9) > 0.5:
        sub = max(("astrocyte", "ensheathing glia", "cortex glia"), key=lambda k: row.get(k, -9))
        names[cl] = "glia" + (f": {sub.replace(' glia', '')}" if row.get(sub, -9) > 1 else "")
    else:
        names[cl] = "neuron"
cells["cell_type"] = cells.cluster.map(names) + " (c" + cells.cluster + ")"
summary = z.round(2).assign(name=pd.Series(names), cells=cells.cluster.value_counts())
display(summary.sort_values("name"))
if (pd.Series(names) == "dividing cells").any():
    n_div = int((cells.cluster.map(names) == "dividing cells").sum())
    print(f"\n{n_div:,} cells are in clusters that look like neuroblasts or cells in S/G2/M. A cell in S or G2 has "
          "more than 2C of DNA and will be called polyploid by Cell 6, so these are reported separately in Cell 7 "
          "and must not be read as endoreplication.")
