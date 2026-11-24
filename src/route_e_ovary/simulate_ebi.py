"""Mock the EBI Single Cell Expression Atlas layout: one .mtx with FlyBase IDs, plus _rows and _cols."""
import gzip, io, os, shutil, sys, numpy as np, scipy.sparse as sp
from scipy.io import mmread, mmwrite
src, out = sys.argv[1], sys.argv[2]
old = os.path.join(src, "geo", "series", "GSE146nnn", "GSE146040", "suppl")
# empty GEO for GSE136162 (series, matrix and samples all absent)
os.makedirs(os.path.join(out, "geo", "series", "GSE136nnn", "GSE136162", "suppl"), exist_ok=True)
d = os.path.join(out, "ebi", "E-GEOD-136162")
os.makedirs(d, exist_ok=True)
Ms, cols = [], []
rows = None
for n in sorted(os.listdir(old)):
    if n.endswith("matrix.mtx.gz"):
        with gzip.open(os.path.join(old, n), "rb") as h:
            Ms.append(sp.csr_matrix(mmread(h)))
    if n.endswith("barcodes.tsv.gz"):
        lib = n.split("_barcodes")[0]
        cols += [f"{lib}-{l.strip()}" for l in gzip.open(os.path.join(old, n), "rt")]
    if n.endswith("features.tsv.gz") and rows is None:
        rows = [l.split("\t")[0] for l in gzip.open(os.path.join(old, n), "rt")]
M = sp.hstack(Ms).tocoo()                       # genes x cells, as EBI stores it
b = io.BytesIO(); mmwrite(b, M.astype(int), field="integer")
open(f"{d}/E-GEOD-136162.aggregated_filtered_counts.mtx", "wb").write(b.getvalue())
open(f"{d}/E-GEOD-136162.aggregated_filtered_counts.mtx_rows", "w").write("".join(f"{g}\t{g}\n" for g in rows))
open(f"{d}/E-GEOD-136162.aggregated_filtered_counts.mtx_cols", "w").write("".join(c + "\n" for c in cols))
open(f"{d}/ExpDesign-E-GEOD-136162.tsv", "w").write("Assay\tFactor\n" + "".join(f"{c}\tovary\n" for c in cols))
# a matching annotation so FlyBase IDs can be turned back into gene names
feat = [l.rstrip("\n").split("\t") for l in gzip.open(os.path.join(old, sorted(os.listdir(old))[1]), "rt")]
with gzip.open(os.path.join(out, "dm.gtf.gz"), "wt") as h:
    for fb, sym, *_ in feat:
        h.write(f'2L\tFlyBase\tgene\t1\t100\t.\t+\t.\tgene_id "{fb}"; gene_name "{sym}";\n')
print("built EBI mock:", M.shape)
