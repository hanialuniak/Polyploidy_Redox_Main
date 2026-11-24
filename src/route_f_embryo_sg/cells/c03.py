#@title Cell 3 - Find the salivary gland cells and the two stage windows
meta = pd.read_csv(META, sep="\t", dtype=str, low_memory=False)
meta.columns = [c.strip() for c in meta.columns]
print(f"Annotation table: {len(meta):,} rows x {meta.shape[1]} columns")


def pick_column(patterns, exclude=()):
    """Find the column whose name matches, preferring annotated values over ontology URLs."""
    for pat in patterns:
        hits = [c for c in meta.columns
                if re.search(pat, c, re.I) and "ontology" not in c.lower()
                and not any(re.search(x, c, re.I) for x in exclude)]
        hits = [c for c in hits if meta[c].notna().any() and meta[c].nunique() > 1]
        if hits:
            return sorted(hits, key=lambda c: meta[c].nunique())[-1]
    return None


cell_col = pick_column([r"inferred cell type.*author", r"inferred cell type", r"\bcell type\b", r"annotat"])
stage_col = pick_column([r"developmental stage", r"\bstage\b", r"\bage\b", r"time.?point"])
id_col = next((c for c in meta.columns if re.search(r"^assay$|cell.?id|barcode|run", c, re.I)), meta.columns[0])
print(f"Cell-type column: {cell_col!r}\nStage column:     {stage_col!r}\nCell id column:   {id_col!r}")
if cell_col is None or stage_col is None:
    print("\nColumns available:"); print("\n".join("   " + c for c in meta.columns))
    raise StopHere("Could not identify the cell-type or stage column. Send me the list above and I will set them by name.")

print(f"\nCell types annotated ({meta[cell_col].nunique()}):")
display(meta[cell_col].value_counts().head(40))
print(f"\nStage values ({meta[stage_col].nunique()}):")
display(meta[stage_col].value_counts())


def window_of(v):
    """Stages 10-12 are the early window, 13-16 the late one. Written to accept 'stage 13-16' or 'stage 14'."""
    nums = [int(n) for n in re.findall(r"\d+", str(v))]
    nums = [n for n in nums if 1 <= n <= 17]
    if not nums:
        return None
    lo, hi = min(nums), max(nums)
    if hi <= 12:
        return "early (stages 10-12)"
    if lo >= 13:
        return "late (stages 13-16)"
    return None                      # spans both windows: cannot be assigned


meta["window"] = meta[stage_col].map(window_of)
meta["is_sg"] = meta[cell_col].astype(str).str.contains(TISSUE_PATTERN, case=False, na=False)
meta["tissue"] = np.where(meta.is_sg, "salivary gland", "other embryonic cells")
print("\nCells per group:")
display(pd.crosstab(meta.tissue, meta.window))
if not meta.is_sg.any():
    raise StopHere(f"No cell type matched {TISSUE_PATTERN!r}. Look at the cell-type list above and set "
                   f"TISSUE_PATTERN in Cell 1 to something that appears there.")
counts = meta.groupby(["tissue", "window"]).size()
missing = [k for k in [("salivary gland", "early (stages 10-12)"), ("salivary gland", "late (stages 13-16)")]
           if counts.get(k, 0) < 25]
if missing:
    raise StopHere(f"Too few salivary gland cells in {missing}. The comparison needs both windows. "
                   f"Check the stage values printed above.")
meta = meta[meta.window.notna()].copy()
keep_ids = dict(zip(meta[id_col].astype(str), meta.tissue + " | " + meta.window))
print(f"\n{len(keep_ids):,} cells assigned to a tissue and a window.")
