#@title Cell 4 - Stream the matrix into pseudobulk (a few minutes; the full matrix is never held in memory)
rows = [l.split("\t")[0].strip() for l in open(ROWS) if l.strip()]          # genes
cols = [l.split("\t")[0].strip() for l in open(COLS) if l.strip()]          # cells
print(f"Matrix index: {len(rows):,} rows x {len(cols):,} columns")

with open(MTX) as h:
    header = h.readline()
    line = h.readline()
    while line.startswith("%"):
        line = h.readline()
    n1, n2, nnz = (int(x) for x in line.split())
print(f"Matrix header says {n1:,} x {n2:,} with {nnz:,} non-zero entries")
# EBI writes genes as rows, but check rather than assume
if (n1, n2) == (len(rows), len(cols)):
    gene_axis, cell_axis = 0, 1
elif (n1, n2) == (len(cols), len(rows)):
    gene_axis, cell_axis = 1, 0
    print("   matrix is transposed relative to the index files; handled")
else:
    raise StopHere(f"Matrix is {n1} x {n2} but the index files have {len(rows)} genes and {len(cols)} cells. "
                   "These files do not belong together.")

groups = sorted(set(keep_ids.values()))
gidx = {g: i for i, g in enumerate(groups)}
col_group = np.full(len(cols), -1, dtype=np.int16)
for j, c in enumerate(cols):
    g = keep_ids.get(c)
    if g is not None:
        col_group[j] = gidx[g]
print(f"Columns matched to the annotation table: {int((col_group >= 0).sum()):,} of {len(cols):,}")
if (col_group >= 0).sum() < 0.2 * len(cols):
    raise StopHere("Fewer than a fifth of the matrix columns matched the annotation table, so the cell "
                   "identifiers do not line up. Send me the first few of each and I will fix the matching.")

sums = np.zeros((len(rows), len(groups)), dtype=np.float64)
per_cell = np.zeros(len(cols), dtype=np.float64)
read = 0
with open(MTX) as h:
    h.readline()
    line = h.readline()
    while line.startswith("%"):
        line = h.readline()
    for line in h:                                   # one non-zero entry per line: row col value
        a, b, v = line.split()
        i, j = int(a) - 1, int(b) - 1
        if gene_axis == 1:
            i, j = j, i
        g = col_group[j]
        val = float(v)
        per_cell[j] += val
        if g >= 0:
            sums[i, g] += val
        read += 1
        if read % 20_000_000 == 0:
            print(f"   {read:,} entries read")
print(f"Finished: {read:,} entries")

pb = pd.DataFrame(sums, index=rows, columns=groups)
pb = pb.groupby(level=0).sum()
n_cells = pd.Series(col_group[col_group >= 0]).map({i: g for g, i in gidx.items()}).value_counts().reindex(groups)
print("\nCells and total counts per group:")
display(pd.DataFrame({"cells": n_cells, "counts": pb.sum().round(0)}))
low = n_cells[n_cells < 25]
if len(low):
    raise StopHere(f"These groups have too few cells to pseudobulk: {list(low.index)}.")

# median-of-ratios, so that a few very high genes in one group cannot set the scale
pb = pb[(pb / pb.sum() * 1e6 >= 5).all(axis=1)]
lg = np.log(pb.replace(0, np.nan)); ref = lg.mean(axis=1)
size = np.exp(lg[ref.notna()].sub(ref[ref.notna()], axis=0).median())
L = np.log2(pb / size / pb.sum().mean() * 1e6 + 1)
print(f"\nGenes expressed in every group: {len(L):,}")
print("Size factors:", dict(size.round(2)))
