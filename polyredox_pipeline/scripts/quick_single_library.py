"""Quick single-library check: ploidy from ATAC fragments vs redox scores from RNA.

For one 10x Multiome library (for example the public 10x human brain data). With one
donor this is a feasibility check, not a test that generalises. Edit the four paths below.
"""
from importlib import resources

import numpy as np, pandas as pd
from scipy.stats import spearmanr
from polyredox import doublets, fragments, pipeline, ploidy, scoring

FRAGMENTS = "atac_fragments.tsv.gz"   # from Cell Ranger ARC
RNA_DIR = "filtered_feature_bc_matrix"
N_RECOVERED = 3233                    # nuclei recovered in the channel (3,233 for brain 3k)
EXCLUDE_BED = None                    # e.g. a UCSC simple-repeats BED for GRCh38
GMT = resources.files("polyredox") / "data" / "starter_genesets.gmt"

# cells.tsv: cell_id (= barcode), cell_class (e.g. Exc, Inh, Oligo), cell_type
cells = pd.read_csv("cells.tsv", sep="\t")[["cell_id", "cell_class", "cell_type"]]

# 1-2. Fragment overlaps and ploidy posteriors
counts = fragments.library_overlaps(FRAGMENTS, cells["cell_id"].to_numpy(dtype=object),
                                    exclude=fragments.load_exclude_regions(EXCLUDE_BED))
overlaps = pd.DataFrame(counts, columns=[f"f{i}" for i in range(1, 7)], index=cells["cell_id"])
calls, _ = ploidy.call_ploidy(overlaps, cells["cell_class"].to_numpy(dtype=object), min_fragments=3000)
cells = cells.join(calls, on="cell_id").assign(library="lib1", donor="donor1", n_recovered=N_RECOVERED)

# 3. Is there more polyploidy than same-class doublets would give?
cells = doublets.expected_rates(cells)
print(cells.dropna(subset=["P_polyploid"]).groupby("cell_class")[["P_polyploid", "expected_doublet"]].mean().round(3))

# 4. Redox and control scores on depth-matched RNA
frame = cells.dropna(subset=["P_polyploid"]).reset_index(drop=True)
X, genes, kept = pipeline.load_rna({"format": "mtx", "path": RNA_DIR}, frame["cell_id"].to_numpy(dtype=object))
frame = frame.set_index("cell_id").loc[kept].reset_index()
target = int(np.quantile(np.asarray(X.sum(axis=1)).ravel(), 0.1))
X, keep = pipeline.thin_counts(X, target, seed=0)
frame = frame[keep].reset_index(drop=True)
sets, _, _ = scoring.resolve_sets(scoring.read_gmt(str(GMT)), genes)
names = list(sets)
scores = scoring.score_sets(scoring.rank_matrix(X), [sets[n] for n in names])
scores = pd.DataFrame(pipeline.zscore_by(scores, frame["cell_type"].to_numpy(dtype=object)), columns=names)

# 5. Correlation within each cell type, redox next to the stress controls
rows = []
for cell_type, index in frame.groupby("cell_type").groups.items():
    if len(index) < 100:
        continue
    for name in names:
        rho, p = spearmanr(frame.loc[index, "P_polyploid"], scores.loc[index, name])
        rows.append({"cell_type": cell_type, "programme": name, "rho": rho, "p": p, "n": len(index)})
print(pd.DataFrame(rows).pivot(index="programme", columns="cell_type", values="rho").round(3))
