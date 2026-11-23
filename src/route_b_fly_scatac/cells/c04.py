#@title Cell 4 - Read fragments: overlap counts and gene accessibility (10-30 minutes per run)
K = 6
SHIFT = np.uint64(32)


def count_overlaps(codes, starts, ends, n, k=K):
    """Per cell, f_x = number of fragments whose left end is covered by x same-cell fragments."""
    codes, starts, ends = (np.asarray(a, dtype=np.uint64) for a in (codes, starts, ends))
    if codes.size == 0:
        return np.zeros((n, k), dtype=np.int64)
    key = (codes << SHIFT) | starts
    order = np.lexsort((ends, key))
    key, ends, codes = key[order], ends[order], codes[order]
    keep = np.r_[True, (key[1:] != key[:-1]) | (ends[1:] != ends[:-1])]
    key, ends, codes = key[keep], ends[keep], codes[keep]
    end_key = np.sort((codes << SHIFT) | ends)
    depth = np.searchsorted(key, key, side="right") - np.searchsorted(end_key, key, side="right")
    depth = np.clip(depth, 1, k).astype(np.int64)
    return np.bincount(codes.astype(np.int64) * k + depth - 1, minlength=n * k).reshape(n, k)


chrom_bins = {}
offset = 0
for c in sorted(set(genes.chrom)):
    size = int(genes.loc[genes.chrom == c, "w_end"].max()) + 10 * BIN
    chrom_bins[c] = (offset, size // BIN + 1)
    offset += size // BIN + 1
n_bins = offset

overlaps, bin_blocks = [], []
for run in RUN_LIST:
    run_cells = cells[cells.run == run]
    index = pd.Index(run_cells["barcode"])
    f_counts = np.zeros((len(index), K), dtype=np.int64)
    acc = sp.csr_matrix((len(index), n_bins), dtype=np.float32)
    reader = pd.read_csv(paths[(run, "fragments.tsv.gz")], sep="\t", comment="#", header=None, usecols=[0, 1, 2, 3],
                         names=["chrom", "start", "end", "barcode"], dtype={"chrom": str, "barcode": str},
                         chunksize=5_000_000)
    buffer, current = [], None

    def flush(chrom, parts):
        global acc
        if not parts:
            return
        c = np.concatenate([p[0] for p in parts]); s = np.concatenate([p[1] for p in parts]); e = np.concatenate([p[2] for p in parts])
        if chrom in AUTOSOMES:
            f_counts[:] += count_overlaps(c, s, e, len(index))
        if chrom in chrom_bins:
            base, nb = chrom_bins[chrom]
            b = np.r_[np.minimum(s // BIN, nb - 1), np.minimum((e - 1) // BIN, nb - 1)].astype(np.int32) + base
            cc = np.r_[c, c].astype(np.int32)
            acc = acc + sp.csr_matrix((np.ones(b.size, dtype=np.float32), (cc, b)), shape=(len(index), n_bins))

    for chunk in reader:
        chunk["chrom"] = chunk["chrom"].map(norm_chrom)
        code = index.get_indexer(chunk["barcode"])
        chunk = chunk[code >= 0].assign(code=code[code >= 0])
        for chrom, part in chunk.groupby("chrom", sort=False):
            if chrom != current:
                flush(current, buffer)
                buffer, current = [], chrom
            buffer.append((part["code"].to_numpy(np.int32), part["start"].to_numpy(np.int64), part["end"].to_numpy(np.int64)))
    flush(current, buffer)
    overlaps.append(pd.DataFrame(f_counts, columns=[f"f{i}" for i in range(1, K + 1)], index=run_cells["cell_id"]))
    bin_blocks.append(acc)
    print(f"   {run}: {len(index):,} cells, {int(f_counts.sum()):,} autosomal fragments")

overlaps = pd.concat(overlaps).reindex(cells["cell_id"])
cut_bins = sp.vstack(bin_blocks).tocsr()

g_rows, g_cols = [], []
for gi, g in genes.iterrows():
    base, nb = chrom_bins[g.chrom]
    b0, b1 = base + g.w_start // BIN, base + min(g.w_end // BIN, nb - 1)
    g_rows.extend(range(b0, b1 + 1)); g_cols.extend([gi] * (b1 - b0 + 1))
membership = sp.csr_matrix((np.ones(len(g_rows), dtype=np.float32), (g_rows, g_cols)), shape=(n_bins, len(genes)))
gene_counts = (cut_bins @ membership).tocsr()
cells["n_fragments"] = overlaps.sum(axis=1).to_numpy()
cells["n_cuts"] = np.asarray(cut_bins.sum(axis=1)).ravel()
all_depth = cells[["run", "n_fragments"]].copy()   # every called cell, kept for the doublet check in Cell 7
print(f"Gene accessibility matrix: {gene_counts.shape[0]:,} cells x {gene_counts.shape[1]:,} genes")
