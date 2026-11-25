#@title Cell 4 - Stream the matrix into pseudobulk (a few minutes; the full matrix is never held in memory)
if DEC:
    # column 2 of the decorated file is the gene symbol; fall back to the id where it is blank
    rows = []
    for l in open(DEC):
        if not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        rows.append((f[1].strip() if len(f) > 1 and f[1].strip() else f[0].strip()))
else:
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

# The salivary gland groups are small (a few hundred cells), so a CPM floor applied to every group throws
# away genes that are perfectly well measured in the large ones. Require a few actual counts in each group
# instead, which keeps abundant genes like Sec61beta and the mitotic set that a CPM floor removed.
before = len(pb)
ALL_GENES = set(pb.index)   # kept so Cell 5 can say whether a missing gene was absent or filtered out
# Do NOT require expression in every group: a gene switched on across the window is near zero in the early
# gland group by definition, so an all-groups rule deletes exactly what this comparison exists to find.
# Ask instead that the gene is really expressed SOMEWHERE in the gland and SOMEWHERE in the background,
# with enough total signal to give a stable ratio.
sg_cols = [c for c in pb.columns if c.startswith("salivary gland")]
oth_cols = [c for c in pb.columns if c not in sg_cols]
keep = ((pb[sg_cols].max(axis=1) >= MIN_COUNTS_PER_GROUP) &
        (pb[oth_cols].max(axis=1) >= MIN_COUNTS_PER_GROUP) &
        (pb.sum(axis=1) >= MIN_TOTAL_COUNTS))
pb = pb[keep]
print(f"\nGene filter: {before:,} -> {len(pb):,} genes (at least {MIN_COUNTS_PER_GROUP} counts in one gland "
      f"group AND one background group, and {MIN_TOTAL_COUNTS} counts in total)")
zero_sg = int((pb[sg_cols] == 0).any(axis=1).sum())
print(f"   {zero_sg:,} of them are absent from one gland window - usually the point, but their fold changes "
      f"rest on few counts, so read the gene tables with that in mind")
# Diagnostic: what are these cells actually expressing? If the annotation is right, the top genes in the
# late gland window should look like a secretory gland. If it is wrong, nothing below can be trusted.
sg_late = [c for c in sg_cols if "late" in c][0]
raw = pd.DataFrame(sums, index=rows, columns=groups).groupby(level=0).sum()
print(f"\nTop 20 genes by counts in '{sg_late}':")
display(raw[sg_late].sort_values(ascending=False).head(20).astype(int).to_frame("counts"))
DIAG = "Sec61beta Sec61alpha Srp54 Sec31 zetaCOP Sar1 Sec23 CrebA fkh sage Pasilla hth CycB polo stg Hsc70-3 Pdi".split()
have = [g for g in DIAG if g in raw.index]
print("Raw counts for the control genes (before any filtering):")
display(raw.loc[have].astype(int))
gone = [g for g in DIAG if g not in raw.index]
if gone:
    print("Not present in the matrix under these names:", ", ".join(gone))

lg = np.log(pb.replace(0, np.nan)); ref = lg.mean(axis=1)
size = np.exp(lg[ref.notna()].sub(ref[ref.notna()], axis=0).median())
L = np.log2(pb / size / pb.sum().mean() * 1e6 + 1)
print(f"\nGenes expressed in every group: {len(L):,}")
print("Size factors:", dict(size.round(2)))
