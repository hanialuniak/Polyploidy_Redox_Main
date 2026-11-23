#@title Cell 3 - Samples, ploidy labels and gene-level expression
def read_matrix(path):
    meta, rows, inside = {}, [], False
    with gzip.open(path, "rt", errors="replace") as h:
        for line in h:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                inside = True; continue
            if line.startswith("!series_matrix_table_end"):
                break
            if inside:
                rows.append(line.split("\t"))
            elif line.startswith("!Sample_"):
                k, *v = line.split("\t")
                meta.setdefault(k, []).append([x.strip('"') for x in v])
    table = pd.DataFrame(rows[1:], columns=[c.strip('"') for c in rows[0]]) if len(rows) > 1 else None
    return meta, table


samples, tables = [], []
for m in matrices:
    meta, tab = read_matrix(m)
    gsm = meta["!Sample_geo_accession"][0]
    text = [" | ".join(v[i] for vs in meta.values() for v in vs if i < len(v)) for i in range(len(gsm))]
    samples.append(pd.DataFrame({"gsm": gsm, "title": meta["!Sample_title"][0], "text": text}))
    if tab is not None:
        tables.append(tab.set_index(tab.columns[0]))
samples = pd.concat(samples, ignore_index=True)
if not tables or tables[0].shape[0] == 0:
    raise StopHere("The series matrix holds no expression values (the data are only in supplementary files). Send me this message.")
expr = pd.concat(tables, axis=1)
expr = expr.apply(pd.to_numeric, errors="coerce")
expr.index = expr.index.astype(str).str.strip('"')


def ploidy_of(s):
    m = re.search(r"(?<![\d.])(2|4|8|16)\s*-?\s*[cCnN](?![a-zA-Z])", s) or re.search(r"(?i)(di|tetra|octa)ploid", s)
    if not m:
        return np.nan
    g = m.group(1).lower()
    return {"di": 2, "tetra": 4, "octa": 8}.get(g, float(g) if g.isdigit() else np.nan)


samples["ploidy"] = [ploidy_of(t) if not np.isnan(ploidy_of(t)) else ploidy_of(x) for t, x in zip(samples.title, samples.text)]
sheet = os.path.join(DATA_DIR, "sample_labels.csv")
if samples.ploidy.isna().any() and os.path.exists(sheet):
    samples = samples.drop(columns="ploidy").merge(pd.read_csv(sheet)[["gsm", "ploidy"]], on="gsm")
    print("Using your edited sample_labels.csv")
samples[["gsm", "title", "ploidy"]].to_csv(sheet, index=False)
display(samples[["gsm", "title", "ploidy"]])
if samples.ploidy.isna().any():
    raise StopHere("Some samples have no ploidy label. Open sample_labels.csv in your Drive folder, type 2, 4 or 8 in "
                   "the ploidy column, save, and run this cell again.")
samples = samples[samples.gsm.isin(expr.columns)].reset_index(drop=True)
expr = expr[samples.gsm]
print(samples.groupby("ploidy").size().rename("samples").to_string())

# Probe -> gene symbol
if os.path.exists(annot):
    with gzip.open(annot, "rt", errors="replace") as h:
        txt = h.read()
else:
    txt = open(full, errors="replace").read()
if "!platform_table_begin" in txt:
    txt = txt.split("!platform_table_begin", 1)[1].split("!platform_table_end", 1)[0]
lines = [l for l in txt.splitlines() if l and not l.startswith(("!", "#", "^"))]
ann = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", dtype=str, on_bad_lines="skip")
idcol = ann.columns[0]
symcol = next((c for c in ann.columns if c.lower().replace("_", " ") in ("gene symbol", "symbol", "genesymbol")), None)
if symcol is None:
    symcol = next((c for c in ann.columns if "assignment" in c.lower() or "symbol" in c.lower()), None)
if symcol is None:
    raise StopHere(f"No gene symbol column in the annotation. Columns: {list(ann.columns)[:15]}. Send me this message.")
sym = ann.set_index(ann[idcol].astype(str))[symcol].astype(str)
if "assignment" in symcol.lower():
    sym = sym.str.split("//").str[1]
sym = sym.str.split("///").str[0].str.strip()
sym = sym[(sym != "") & (sym != "nan") & (sym != "---")]
if np.nanmax(expr.to_numpy()) > 100:
    expr = np.log2(expr.clip(lower=0) + 1)
    print("Values were not log-scaled: log2 applied.")
e = expr.loc[expr.index.intersection(sym.index)].copy()
e["gene"] = sym.loc[e.index].to_numpy()
e["mean"] = e[samples.gsm].mean(axis=1)
e = e.sort_values("mean", ascending=False).drop_duplicates("gene").set_index("gene")[samples.gsm]
e = e[e.notna().all(axis=1)]
# Make samples comparable: quantile normalisation (every sample gets the same distribution of values)
ranks = e.rank(axis=0, method="first").astype(int) - 1
ref = np.sort(e.to_numpy(), axis=0).mean(axis=1)
e = pd.DataFrame(ref[ranks.to_numpy()], index=e.index, columns=e.columns)
print("Samples quantile-normalised.")
# Animal / replicate number from titles such as 2c_3 (the same number = cells sorted from the same liver)
samples["animal"] = samples.title.str.extract(r"[_\s-](\d+)\s*$", expand=False)
e = e[e.mean(axis=1) > np.quantile(e.mean(axis=1), 0.25)]        # drop the 25% least expressed genes
print(f"Genes kept: {len(e):,} (most expressed probe per gene; bottom quarter by expression removed)")
