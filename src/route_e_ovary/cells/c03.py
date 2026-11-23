#@title Cell 3 - Load the count matrices (cells x genes)
def read_mtx_set(mtx, bar, feat):
    """A 10x triplet: matrix.mtx(.gz), barcodes.tsv(.gz), features/genes.tsv(.gz)."""
    from scipy.io import mmread
    op = gzip.open if mtx.endswith(".gz") else open
    with op(mtx, "rb") as h:
        M = sp.csr_matrix(mmread(h))
    rd = lambda p: [l.rstrip("\n").split("\t") for l in (gzip.open(p, "rt") if p.endswith(".gz") else open(p))]
    cells_ = [r[0] for r in rd(bar)]
    rows = rd(feat)
    names = [r[1] if len(r) > 1 else r[0] for r in rows]           # column 2 is the gene symbol
    if M.shape[0] == len(names) and M.shape[1] == len(cells_):
        M = M.T.tocsr()
    return M, cells_, names


def read_dense(path):
    op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    first = op.readline(); op.close()
    sep = "," if first.count(",") > first.count("\t") else "\t"
    d = pd.read_csv(path, sep=sep, index_col=0, low_memory=False)
    d = d.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    if d.shape[0] < d.shape[1]:                                     # genes should be the long side
        d = d.T
    return sp.csr_matrix(d.to_numpy(dtype=np.float32).T), list(d.columns), list(d.index)


def to_symbols(gs):
    """EBI files name genes by FlyBase ID; turn them into gene names using the annotation."""
    if np.mean([str(g).startswith("FBgn") for g in gs[:500]]) <= 0.5:
        return list(gs), False
    if not os.path.exists(GTF_PATH):
        raise StopHere(f"The matrix uses FlyBase IDs and the annotation is not at {GTF_PATH}. It is the same file "
                       "used in Route B; fix GTF_PATH in Cell 1.")
    global _ID2SYM
    if "_ID2SYM" not in globals():
        _ID2SYM = {}
        with gzip.open(GTF_PATH, "rt") as h:
            for line in h:
                if "\tgene\t" not in line:
                    continue
                a = re.search(r'gene_id "([^"]+)"', line); b = re.search(r'gene_name "([^"]+)"', line)
                if a and b:
                    _ID2SYM[a.group(1)] = b.group(1)
    return [_ID2SYM.get(re.sub(r"\.\d+$", "", str(g)), str(g)) for g in gs], True


def rank_of(name):
    """Lower is better. GEO per-sample files keep the library structure, so they come first.
    The EBI entry ships the same data three times over; take the filtered raw counts."""
    n = name.lower()
    if "aggregated" in n or n.startswith("e-geod") or n.startswith("e-mtab"):
        if "filtered_normalised" in n:
            return 3
        return 1 if "filtered" in n else 2
    return 0


found_blocks = []
by_dir = {}
for p in files:
    by_dir.setdefault(os.path.dirname(p), []).append(p)
for d, ps in sorted(by_dir.items()):
    triplets = {}
    for p in ps:
        n = os.path.basename(p).lower()
        key = re.sub(r"((barcodes|features|genes|matrix)\.(tsv|mtx)|\.mtx(_rows|_cols)?)(\.gz)?$", "", n)
        kind = ("mtx" if n.endswith((".mtx", ".mtx.gz")) or ("matrix" in n and ".mtx" in n)
                else "bar" if "barcode" in n or n.endswith((".mtx_cols", ".mtx_cols.gz"))
                else "feat" if "feature" in n or "genes" in n or n.endswith((".mtx_rows", ".mtx_rows.gz"))
                else None)
        if kind:
            triplets.setdefault(key, {})[kind] = p
    for key, t in sorted(triplets.items()):
        if {"mtx", "bar", "feat"} <= set(t):
            M, cs, gs = read_mtx_set(t["mtx"], t["bar"], t["feat"])
            lib = (key.strip("_-. ") or os.path.basename(d)) or "library1"
            found_blocks.append({"rank": rank_of(lib), "lib": lib, "M": M, "genes": gs, "cells": cs})
if not found_blocks:                                                # no matrix triplets: try a dense table
    for p in sorted(files, key=os.path.getsize, reverse=True):
        if re.search(r"\.(csv|tsv|txt|tab)(\.gz)?$", p, re.I):
            try:
                M, cs, gs = read_dense(p)
            except Exception as ex:
                print(f"   {os.path.basename(p)}: could not read ({ex})"); continue
            if M.shape[0] > 200 and M.shape[1] > 2000:
                lib = re.sub(r"\.(csv|tsv|txt|tab)(\.gz)?$", "", os.path.basename(p), flags=re.I)
                found_blocks.append({"rank": rank_of(lib), "lib": lib, "M": M, "genes": gs, "cells": cs})
if not found_blocks:
    print("Files found:", [os.path.basename(f) for f in files][:40])
    raise StopHere("None of the files could be read as a single-cell matrix. Send me this output.")

print("Matrices found in your folder:")
for b in sorted(found_blocks, key=lambda b: b["rank"]):
    print(f"   {b['lib']}: {b['M'].shape[0]:,} cells x {b['M'].shape[1]:,} genes")
best = min(b["rank"] for b in found_blocks)
chosen = [b for b in found_blocks if b["rank"] == best]
if best > 0 and len(chosen) > 1:                                    # EBI ships one dataset in several forms
    chosen = chosen[:1]
if len(chosen) < len(found_blocks):
    print(f"\nSeveral versions of this dataset are present. Using "
          f"{'the per-sample files' if best == 0 else chosen[0]['lib']} and ignoring the rest, so that one dataset "
          "is not counted twice.")
blocks, groups = [], []
for b in chosen:
    gs, converted = to_symbols(b["genes"])
    blocks.append((b["M"], gs))
    groups.append(pd.DataFrame({"library": b["lib"], "barcode": b["cells"]}))
    print(f"Using {b['lib']}: {b['M'].shape[0]:,} cells x {b['M'].shape[1]:,} genes"
          + (" (FlyBase IDs converted to gene names)" if converted else ""))

shared = sorted(set.intersection(*[set(g) for _, g in blocks]))
if len(shared) < 3000:
    raise StopHere(f"Only {len(shared)} gene names are shared across the libraries. Send me this message.")
X = sp.vstack([M[:, [{g: i for i, g in enumerate(gs)}[s] for s in shared]] for M, gs in blocks]).tocsr()
cells = pd.concat(groups, ignore_index=True)
genes = pd.Index(shared)
X.data = np.nan_to_num(X.data)
cells["n_counts"] = np.asarray(X.sum(axis=1)).ravel()
cells["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()
mito = [i for i, g in enumerate(genes) if g.startswith(("mt:", "mt-"))]
cells["pct_mito"] = 100 * np.asarray(X[:, mito].sum(axis=1)).ravel() / np.maximum(cells.n_counts, 1) if mito else 0.0
keep = (cells.n_genes >= MIN_GENES) & (cells.n_counts >= MIN_COUNTS) & (cells.pct_mito < 30)
print(f"\nCells before filtering: {len(cells):,}  |  kept: {int(keep.sum()):,}")
cells, X = cells[keep.to_numpy()].reset_index(drop=True), X[np.flatnonzero(keep.to_numpy())]
detected = np.asarray((X > 0).sum(axis=0)).ravel()
gkeep = detected >= 10
X, genes = X[:, np.flatnonzero(gkeep)].tocsr(), genes[gkeep]
print(f"Genes kept (detected in >=10 cells): {len(genes):,}  |  libraries: {cells.library.nunique()}")
display(cells.groupby("library").agg(cells=("barcode", "size"), median_counts=("n_counts", "median")))
