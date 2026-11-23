#@title Cell 3 - Cells and genes
def norm_chrom(c):
    return re.sub(r"^chr", "", str(c))


cells = []
for run in RUN_LIST:
    sc = pd.read_csv(paths[(run, "singlecell.csv")])
    sc = sc[sc["barcode"] != "NO_BARCODE"]
    flag = "is__cell_barcode" if "is__cell_barcode" in sc.columns else None
    called = sc[sc[flag] == 1] if flag else sc
    cells.append(pd.DataFrame({"run": run, "barcode": called["barcode"].astype(str), "n_recovered": len(called)}))
cells = pd.concat(cells, ignore_index=True)
cells["cell_id"] = cells["run"] + ":" + cells["barcode"]
print(f"Cells called by Cell Ranger: {len(cells):,} in {len(RUN_LIST)} run(s)")

rows = []
with gzip.open(gtf, "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if len(f) < 9 or f[2] != "gene":
            continue
        name = re.search(r'gene_name "([^"]+)"', f[8]) or re.search(r'gene_id "([^"]+)"', f[8])
        rows.append((norm_chrom(f[0]), int(f[3]) - 1, int(f[4]), f[6], name.group(1)))
genes = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gene"]).drop_duplicates("gene")
genes = genes[genes.chrom.isin(AUTOSOMES | {"X", "4"})].reset_index(drop=True)
genes["w_start"] = np.where(genes.strand == "+", genes.start - PROMOTER, genes.start).clip(0)
genes["w_end"] = np.where(genes.strand == "-", genes.end + PROMOTER, genes.end)
print(f"Genes read from the annotation: {len(genes):,}")
