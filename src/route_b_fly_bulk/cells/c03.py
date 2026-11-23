#@title Cell 3 - Build the expression table and read sample information
# Gene IDs to symbols, from the fly annotation
id2sym = {}
with gzip.open(GTF_PATH, "rt") as h:
    for line in h:
        if "\tgene\t" in line:
            gid = re.search(r'gene_id "([^"]+)"', line); sym = re.search(r'gene_name "([^"]+)"', line)
            if gid:
                id2sym[gid.group(1)] = sym.group(1) if sym else gid.group(1)

# Sample information from the series matrix
meta = {}
with gzip.open(os.path.join(DATA_DIR, f"{GSE}_series_matrix.txt.gz"), "rt", errors="replace") as h:
    for line in h:
        if line.startswith("!Sample_"):
            key, *vals = line.rstrip("\n").split("\t")
            meta.setdefault(key, []).append([v.strip('"') for v in vals])
samples = pd.DataFrame({"gsm": meta["!Sample_geo_accession"][0], "title": meta["!Sample_title"][0]})
samples["text"] = [" | ".join(str(v[i]) for k, vs in meta.items() for v in vs if i < len(v)) for i in range(len(samples))]


def read_table(path):
    opener = gzip.open if path.endswith(".gz") else open
    if re.search(r"\.xlsx?$", path):
        return pd.read_excel(path)
    with opener(path, "rt", errors="replace") as h:
        head = h.read(5000)
    sep = "," if head.count(",") > head.count("\t") else "\t"
    return pd.read_csv(path, sep=sep, comment=None, low_memory=False)


def as_matrix(table):
    """Gene x sample numbers from a table whose first text column holds gene IDs or symbols."""
    text_cols = [c for c in table.columns if not pd.api.types.is_numeric_dtype(table[c])]
    num_cols = [c for c in table.columns if pd.api.types.is_numeric_dtype(table[c])]
    if not text_cols or not num_cols:
        return None
    key = max(text_cols, key=lambda c: table[c].astype(str).str.match(r"^FBgn\d+").mean())
    genes = table[key].astype(str).map(lambda g: id2sym.get(g.split(".")[0], g))
    m = table[num_cols].set_index(genes).groupby(level=0).sum()
    return m.loc[:, [c for c in m.columns if not re.search(r"length|^start|^end|chr", str(c), re.I)]]


candidates = [p for p in glob.glob(os.path.join(DATA_DIR, f"{GSE}*")) if not re.search(r"series_matrix|\.tar$", p)]
expr = None
for p in candidates:
    try:
        m = as_matrix(read_table(p))
    except Exception as error:
        print(f"   could not read {os.path.basename(p)}: {type(error).__name__}")
        continue
    if m is not None and m.shape[1] >= 6 and m.shape[0] > 5000:
        expr = m
        print(f"   using {os.path.basename(p)}: {m.shape[0]:,} genes x {m.shape[1]} samples")
        break
if expr is None and os.path.isdir(os.path.join(DATA_DIR, "raw")):
    cols = {}
    for p in sorted(glob.glob(os.path.join(DATA_DIR, "raw", "*"))):
        try:
            m = as_matrix(read_table(p))
        except Exception:
            continue
        if m is not None and m.shape[0] > 5000:
            cols[os.path.basename(p).split("_")[0]] = m.iloc[:, -1]
    if len(cols) >= 6:
        expr = pd.DataFrame(cols).fillna(0)
        print(f"   built from {len(cols)} per-sample files: {expr.shape[0]:,} genes")
if expr is None:
    raise StopHere("No expression table found among the GEO files. Send me the file list printed in Cell 2.")

# Match each expression column to a GEO sample: first by any shared ID, otherwise by order
def match(col):
    c = str(col)
    code = re.sub(r"^\D+", "", c)
    for i, r in samples.iterrows():
        if c in (r.gsm, r.title) or r.gsm in c or (len(r.title) > 3 and r.title in c):
            return i
        if len(code) >= 4 and re.search(rf"\b{code}\b", r.text):
            return i
    return None


hits = [match(c) for c in expr.columns]
if sum(h is not None for h in hits) == expr.shape[1]:
    info = samples.loc[hits].reset_index(drop=True)
    info["column"] = list(expr.columns)
    print("Matched every expression column to a GEO sample by name.")
    by_order = False
elif expr.shape[1] == len(samples):
    info = samples.copy().reset_index(drop=True)
    info["column"] = list(expr.columns)
    print(f"The expression columns ({expr.columns[0]}, ...) do not name the GEO samples.")
    print(f"There are {expr.shape[1]} columns and {len(samples)} GEO samples, so they are matched IN ORDER.")
    print("This is checked below: replicates of the same group should resemble each other most.")
    by_order = True
else:
    raise StopHere(f"The table has {expr.shape[1]} sample columns but GEO lists {len(samples)} samples, and "
                   "the column names do not match. Send me this message.")

# Labels from titles such as 0dFemaleOL1 (age, sex, region, replicate) or from the descriptions
t = info.title.astype(str)
low = (info.title + " | " + info.text).str.lower()
info["region"] = t.str.extract(r"(OL|CB|VNC)\d*$", expand=False)
alt = np.select([low.str.contains("optic"), low.str.contains("central brain"), low.str.contains("ventral")],
                ["OL", "CB", "VNC"], "unknown")
info["region"] = info["region"].fillna(pd.Series(alt))
info["age_days"] = t.str.extract(r"^(\d+)\s*d", expand=False).astype(float)
info["age_days"] = info["age_days"].fillna(low.str.extract(r"(\d+)\s*(?:d\b|day)", expand=False).astype(float))
tl = t.str.lower()       # the title is the most reliable source (e.g. 21dMaleCB2); descriptions can mention both sexes
info["sex"] = np.select([tl.str.contains("female"), tl.str.contains("male"),
                         low.str.contains("female") & ~low.str.contains(r"\bmale"), low.str.contains(r"\bmale")],
                        ["F", "M", "F", "M"], "unknown")
info["group"] = info.title.str.replace(r"\d+$", "", regex=True)

sheet = os.path.join(DATA_DIR, "sample_labels.csv")
auto_ok = not (info.region.isin(["unknown"]) | info.region.isna() | info.age_days.isna() | (info.sex == "unknown")).any()
if os.path.exists(sheet) and not auto_ok:
    info = pd.read_csv(sheet)
    print("Using your edited sample_labels.csv")
else:
    info[["column", "gsm", "title", "group", "region", "age_days", "sex"]].to_csv(sheet, index=False)
display(info[["column", "gsm", "title", "region", "age_days", "sex"]])
bad = info.region.isin(["unknown"]) | info.region.isna() | info.age_days.isna()
if bad.any():
    raise StopHere(f"{int(bad.sum())} samples could not be labelled. Open sample_labels.csv in your Drive folder, fill in "
                   "region (OL, CB or VNC), age_days and sex, save, and run this cell again.")
print(info.groupby(["region", "age_days"]).size().unstack(fill_value=0))
print(info.groupby(["region", "sex"]).size().unstack(fill_value=0))

# Check: is each sample most similar to a replicate of its own group?
L = np.log2(expr[info.column].astype(float) + 1)
L = L[(L > 3).mean(axis=1) >= 0.8]
L = L.loc[L.var(axis=1).sort_values().index[-2000:]]      # the 2,000 most variable expressed genes
C = np.corrcoef((L.sub(L.mean(axis=1), axis=0)).T.to_numpy())
np.fill_diagonal(C, -1)
nearest = C.argmax(axis=1)
same_group = (info.group.to_numpy()[nearest] == info.group.to_numpy()).mean()
same_region = (info.region.to_numpy()[nearest] == info.region.to_numpy()).mean()
print(f"\nCheck: nearest neighbour is a replicate of the same group for {100 * same_group:.0f}% of samples, "
      f"and from the same region for {100 * same_region:.0f}%")
if same_region < 0.8:
    msg = "Samples do not group by region as expected."
    if by_order:
        msg += " The column order probably does not follow GEO's sample order."
    raise StopHere(msg + " Do not go on; send me this output.")
print("Labels look right: samples group by region.")
