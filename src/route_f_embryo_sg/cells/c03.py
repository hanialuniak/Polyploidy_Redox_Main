#@title Cell 3 - Find the salivary gland cells and the two stage windows
meta = pd.read_csv(META, sep="\t", dtype=str, low_memory=False)
meta.columns = [c.strip() for c in meta.columns]
print(f"Annotation table: {len(meta):,} rows x {meta.shape[1]} columns")


#@markdown Leave these empty to find the columns automatically. If that fails, the cell prints every column
#@markdown name; copy the right ones in here and run the cell again.
CELL_TYPE_COLUMN = ""  #@param {type:"string"}
STAGE_COLUMN = ""  #@param {type:"string"}

# Column names arrive differently depending on where the file came from: "Sample Characteristic[cell type]",
# "inferred_cell_type_-_authors_labels", and so on. Match against a flattened version of each name so that
# spaces, underscores, hyphens and brackets all read the same.
flat = {c: re.sub(r"[^a-z0-9]+", " ", c.lower()).strip() for c in meta.columns}


def pick_column(patterns, exclude=()):
    """Find the column whose name matches, preferring annotated values over ontology URLs."""
    for pat in patterns:
        hits = [c for c in meta.columns
                if re.search(pat, flat[c]) and "ontology" not in flat[c]
                and not any(re.search(x, flat[c]) for x in exclude)]
        hits = [c for c in hits if meta[c].notna().any() and meta[c].nunique() > 1]
        if hits:
            return sorted(hits, key=lambda c: meta[c].nunique())[-1]
    return None


cell_col = CELL_TYPE_COLUMN.strip() or pick_column(
    [r"inferred cell type.*author", r"inferred cell type", r"\bcell type\b", r"annotat", r"cluster"])
stage_col = STAGE_COLUMN.strip() or pick_column(
    [r"developmental stage", r"\bstage\b", r"\bage\b", r"time ?point", r"\bhours?\b"])
id_col = next((c for c in meta.columns if re.search(r"^assay$|^id$|cell id|barcode|^run$", flat[c])), meta.columns[0])
for nm, c in (("CELL_TYPE_COLUMN", cell_col), ("STAGE_COLUMN", stage_col)):
    if c and c not in meta.columns:
        raise StopHere(f"{nm} is set to {c!r}, which is not a column in this table.")
print(f"Cell-type column: {cell_col!r}\nStage column:     {stage_col!r}\nCell id column:   {id_col!r}")
if cell_col is None or stage_col is None:
    print("\nColumns available:"); print("\n".join("   " + c for c in meta.columns))
    raise StopHere("Could not identify the cell-type or stage column. Send me the list above and I will set them by name.")

print(f"\nCell types annotated ({meta[cell_col].nunique()}):")
display(meta[cell_col].value_counts().head(40))
print(f"\nStage values ({meta[stage_col].nunique()}):")
display(meta[stage_col].value_counts())


# Some tables give hours after egg laying rather than stages. At 25 C the boundary between stage 12 and
# stage 13 falls at about 9h20, so that is where the two windows are cut.
IN_HOURS = bool(re.search(r"\bhours?\b|\bh ael\b|\btime\b", flat.get(stage_col, "")) or
                meta[stage_col].astype(str).str.contains(r"\bh\b|hour", case=False, na=False).mean() > 0.5)
if IN_HOURS:
    print("\nThe stage column looks like hours after egg laying, so the windows are cut at 9h20 "
          "(the stage 12 to stage 13 boundary at 25 C).")


def window_of(v):
    """Stages 10-12 are the early window, 13-16 the late one. Accepts 'stage 13-16', 'stage 14' or hours."""
    nums = [float(n) for n in re.findall(r"\d+(?:\.\d+)?", str(v))]
    if IN_HOURS:
        if not nums:
            return None
        mid = sum(nums) / len(nums)
        return "early (stages 10-12)" if mid < 9.33 else "late (stages 13-16)"
    nums = [int(n) for n in nums if 1 <= n <= 17]
    if not nums:
        return None
    lo, hi = min(nums), max(nums)
    if hi <= 12:
        return "early (stages 10-12)"
    if lo >= 13:
        return "late (stages 13-16)"
    return None                      # spans both windows: cannot be assigned


meta["window"] = meta[stage_col].map(window_of)
if meta.window.notna().sum() == 0:
    print("\nStage values seen:", list(meta[stage_col].dropna().unique())[:25])
    raise StopHere("No cell could be assigned to an early or late window from the stage column. The values are "
                   "printed above; send them to me and I will set the rule for this table.")
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
