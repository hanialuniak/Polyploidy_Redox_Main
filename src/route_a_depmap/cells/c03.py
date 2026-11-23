#@title Cell 3 - Load and match the data
def read_model_matrix(path):
    """Read a DepMap cell-line x gene file (works with old and new layouts)."""
    table = pd.read_csv(path, low_memory=False)
    gene_columns = [c for c in table.columns if re.search(r"\(\d+\)\s*$", str(c))]
    if "ModelID" in table.columns:
        for flag in ("IsDefaultEntryForModel", "is_default_entry"):
            if flag in table.columns:
                table = table[table[flag].astype(str).str.lower().isin(["yes", "true", "1"])]
        table = table.drop_duplicates("ModelID").set_index("ModelID")
    else:
        table = table.set_index(table.columns[0])
    if not gene_columns:
        gene_columns = [c for c in table.columns if pd.api.types.is_numeric_dtype(table[c])]
    table = table[gene_columns].astype("float32")
    table.columns = [re.sub(r"\s*\(\d+\)\s*$", "", str(c)) for c in table.columns]
    table = table.loc[:, ~table.columns.duplicated()]
    table.index = table.index.astype(str)
    return table


def find_column(table, *names):
    lower = {str(c).lower(): c for c in table.columns}
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    raise StopHere(f"Expected a column like {names[0]}; found {list(table.columns)[:12]}")


models = pd.read_csv(paths["models"], low_memory=False)
models = models.set_index(find_column(models, "ModelID"))
models.index = models.index.astype(str)
lineage_col = find_column(models, "OncotreeLineage", "Lineage")
ccle_col = find_column(models, "CCLEName", "CCLE_Name")

signatures = pd.read_csv(paths["signatures"])
signatures = signatures.set_index("ModelID" if "ModelID" in signatures.columns else signatures.columns[0])
signatures.index = signatures.index.astype(str)

print("Reading expression (about 1 minute) ...")
expression = read_model_matrix(paths["expression"])
damaging = read_model_matrix(paths["damaging"])
hotspot = read_model_matrix(paths["hotspot"])


def mutated(matrix, genes):
    present = [g for g in genes if g in matrix.columns]
    return (matrix[present].fillna(0) > 0).any(axis=1) if present else pd.Series(False, index=matrix.index)


nrf2_mutant = pd.concat([mutated(damaging, ["KEAP1", "CUL3"]), mutated(hotspot, ["NFE2L2", "KEAP1"])], axis=1).any(axis=1)
lines = pd.DataFrame({
    "ploidy": pd.to_numeric(signatures[find_column(signatures, "Ploidy")], errors="coerce"),
    "wgd": pd.to_numeric(signatures[find_column(signatures, "WGD")], errors="coerce"),
}).dropna(subset=["ploidy"])
lines = lines[lines.index.isin(expression.index) & lines.index.isin(models.index)]
lines["tissue"] = models.loc[lines.index, lineage_col].astype(str)
lines["nrf2_mutant"] = nrf2_mutant.reindex(lines.index).fillna(False).astype(bool)
tp53_mutant = pd.concat([mutated(damaging, ["TP53"]), mutated(hotspot, ["TP53"])], axis=1).any(axis=1)
lines["tp53_mutant"] = tp53_mutant.reindex(lines.index).fillna(False).astype(bool)
print(f"Cell lines with ploidy and expression: {len(lines)}")
if EXCLUDE_NRF2_MUTANTS:
    print(f"Leaving out {int(lines.nrf2_mutant.sum())} lines with KEAP1, CUL3 or NFE2L2 mutations")
    lines = lines[~lines.nrf2_mutant]
print(f"TP53-mutant lines: {int(lines.tp53_mutant.sum())} of {len(lines)}")
if TP53_MODE == "TP53-normal lines only":
    print(f"Keeping only the {int((~lines.tp53_mutant).sum())} TP53-normal lines")
    lines = lines[~lines.tp53_mutant]
if LINEAGE == "CNS/Brain only":
    lines = lines[lines.tissue == "CNS/Brain"]
sizes = lines.tissue.value_counts()
lines = lines[lines.tissue.isin(sizes[sizes >= MIN_LINES_PER_TISSUE].index)]
print(f"Cell lines analysed: {len(lines)} from {lines.tissue.nunique()} tissues")
if len(lines) < 25:
    raise StopHere(f"Only {len(lines)} cell lines remain. Choose 'All lineages' in Cell 1 and rerun.")

proteins = None
if USE_PROTEINS:
    raw = pd.read_csv(paths["proteins"], low_memory=False)
    gene_col = find_column(raw, "Gene_Symbol", "GeneSymbol")
    sample_cols = [c for c in raw.columns if re.match(r"^.+_TenPx\d+$", str(c))]
    if not sample_cols:
        print("WARNING: protein file layout not recognised; proteins will be skipped.")
    by_gene = raw[[gene_col] + sample_cols].dropna(subset=[gene_col]).groupby(gene_col).mean()
    by_gene.columns = [re.sub(r"_TenPx\d+$", "", c) for c in by_gene.columns]
    by_gene = by_gene.T.groupby(level=0).mean()
    name_to_id = pd.Series(models.index, index=models[ccle_col].astype(str)).groupby(level=0).first()
    by_gene.index = by_gene.index.map(name_to_id)
    proteins = by_gene[by_gene.index.notna() & by_gene.index.isin(lines.index)]
    print(f"Cell lines that also have protein data: {len(proteins)}")

crispr = None
if USE_CRISPR:
    crispr = read_model_matrix(paths["crispr"])
    crispr = crispr[crispr.index.isin(lines.index)]
    print(f"Cell lines that also have CRISPR knockout data: {len(crispr)}")
