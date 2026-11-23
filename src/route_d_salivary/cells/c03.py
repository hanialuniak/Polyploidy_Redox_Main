#@title Cell 3 - Read samples and the expression table
def read_matrix(path):
    meta = {}
    with gzip.open(path, "rt", errors="replace") as h:
        for line in h:
            if line.startswith("!series_matrix_table_begin"): break
            if line.startswith("!Sample_"):
                k, *v = line.rstrip("\n").split("\t")
                meta.setdefault(k, []).append([x.strip('"') for x in v])
    gsm = meta["!Sample_geo_accession"][0]
    text = [" | ".join(v[i] for vs in meta.values() for v in vs if i < len(v)) for i in range(len(gsm))]
    return pd.DataFrame({"gsm": gsm, "title": meta["!Sample_title"][0], "text": text})


samples = pd.concat([read_matrix(m) for m in matrices], ignore_index=True)


def read_any(name, handle):
    low = name.lower()
    raw = handle.read()
    if low.endswith(".gz"):
        raw = gzip.decompress(raw); low = low[:-3]
    if low.endswith((".xlsx", ".xls")):
        sheets = pd.read_excel(io.BytesIO(raw), sheet_name=None, header=None,
                               engine="openpyxl" if low.endswith(".xlsx") else None)
        df = max(sheets.values(), key=lambda d: d.shape[0] * d.shape[1])
        # the header is the first row whose later cells are mostly text
        hdr = next((i for i in range(min(10, len(df))) if df.iloc[i, 1:].map(lambda v: isinstance(v, str)).mean() > 0.5), None)
        if hdr is not None:
            df.columns = df.iloc[hdr].astype(str); df = df.iloc[hdr + 1:].reset_index(drop=True)
        return df
    txt = raw.decode("utf-8", errors="replace")
    lines = [l for l in txt.splitlines() if l.strip() and not l.startswith(("#", "__"))]
    sep = "," if lines and lines[0].count(",") > lines[0].count("\t") else "\t"
    return pd.read_csv(io.StringIO("\n".join(lines)), sep=sep, header=None if not re.search(r"[A-Za-z]", lines[0].split(sep)[-1]) else 0)


def as_numeric_table(df):
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    skip = r"(?i)^(length|len|start|end|chr|chrom|strand|gene_?length|effective)|fold|log2|pval|p\.?adj|padj|fdr|stat|lfcse|basemean"
    num = [c for c in df.columns[1:] if pd.to_numeric(df[c], errors="coerce").notna().mean() > 0.9 and not re.search(skip, c)]
    idc = next((c for c in df.columns if c not in num), df.columns[0])
    out = df.set_index(df[idc].astype(str))[num].apply(pd.to_numeric, errors="coerce")
    return out[~out.index.str.startswith(("__", "N_", "no_feature", "ambiguous", "not_aligned", "too_low", "alignment_not"))]


tables = []
for path in got:
    name = os.path.basename(path)
    if name.lower().endswith(".tar"):
        per = {}
        with tarfile.open(path) as tf:
            for m in tf.getmembers():
                if not m.isfile(): continue
                g = re.search(r"GSM\d+", m.name)
                try:
                    t = as_numeric_table(read_any(m.name, tf.extractfile(m)))
                except Exception:
                    continue
                if g and t.shape[1] >= 1:
                    col = t.columns[0] if t.shape[1] == 1 else next((c for c in t.columns if re.search(r"count|read", c, re.I)), t.columns[-1])
                    per[g.group(0)] = t[col]
        if per:
            tables.append(pd.DataFrame(per)); print(f"   {name}: {len(per)} per-sample files")
    elif re.search(r"\.(txt|tsv|csv|tab|xlsx?)(\.gz)?$", name, re.I):
        try:
            t = as_numeric_table(read_any(name, open(path, "rb")))
            if t.shape[1] >= 2:
                tables.append(t); print(f"   {name}: {t.shape[0]:,} genes x {t.shape[1]} sample columns")
        except Exception as ex:
            print(f"   {name}: could not read ({ex})")
if not tables:
    raise StopHere("None of the downloaded files could be read as an expression table. Send me this output.")
# prefer the raw count table; differential-expression result tables have few sample columns
expr = max(tables, key=lambda t: t.shape[1])
print(f"Using the table with {expr.shape[1]} sample columns: {list(expr.columns)}")


def match(col):
    c = str(col)
    for i, r in samples.iterrows():
        if r.gsm in c or c == r.title or (len(r.title) > 3 and (r.title in c or c in r.title)):
            return i
    return None


def group_of(text):
    for label, pat, role in GROUP_RULES.get(GSE, []):
        if re.search(pat, text):
            return label, role
    return None, None


hits = [match(c) for c in expr.columns]
by_name = [group_of(str(c)) for c in expr.columns]
if all(h is not None for h in hits):
    info = samples.loc[hits].reset_index(drop=True); info["column"] = list(expr.columns)
elif all(g[0] for g in by_name):
    info = pd.DataFrame({"column": list(expr.columns), "gsm": "", "title": [str(c) for c in expr.columns], "text": ""})
    print("Samples labelled from the column names of the expression table.")
else:
    print("Column names:", list(expr.columns)[:12]); print("GEO titles:", list(samples.title)[:12])
    if expr.shape[1] == len(samples):
        info = samples.copy(); info["column"] = list(expr.columns)
        print("WARNING: columns matched to GEO samples by ORDER. Check the table below; edit sample_labels.csv if wrong.")
    else:
        raise StopHere(f"{expr.shape[1]} expression columns vs {len(samples)} GEO samples and the names do not match. Send me this output.")

lab = [group_of(f"{c} | {t} | {x}") for c, t, x in zip(info.column.astype(str), info.title, info.text)]
info["group"], info["role"] = [l[0] for l in lab], [l[1] for l in lab]
sheet = os.path.join(GDIR, "sample_labels.csv")
if os.path.exists(sheet) and info.group.isna().any():
    info = info.drop(columns=["group", "role"]).merge(pd.read_csv(sheet)[["column", "group", "role"]], on="column")
    print("Using your edited sample_labels.csv")
info[["column", "gsm", "title", "group", "role"]].to_csv(sheet, index=False)
display(info[["column", "gsm", "title", "group", "role"]])
if info.group.isna().any():
    raise StopHere("Some samples have no group. Open sample_labels.csv in the dataset folder, fill in group and role "
                   "(test = endocycling, reference = dividing), save, and run this cell again.")
info = info[info.role.isin(["test", "reference"])].reset_index(drop=True)
expr = expr[info.column]

# Gene names from the fly annotation (FBgn -> symbol)
id2sym = {}
with gzip.open(GTF_PATH, "rt") as h:
    for line in h:
        if "\tgene\t" not in line: continue
        gid = re.search(r'gene_id "([^"]+)"', line); gnm = re.search(r'gene_name "([^"]+)"', line)
        if gid and gnm: id2sym[gid.group(1)] = gnm.group(1)
ids = expr.index.str.replace(r"\.\d+$", "", regex=True)
sym = [id2sym.get(i, i) for i in ids]
expr = expr.groupby(pd.Index(sym)).sum()
print(f"\nGenes: {expr.shape[0]:,}  |  samples per group:\n" + info.groupby(["role", "group"]).size().to_string())
