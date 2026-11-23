#@title Cell 2 - Download the single-cell files from GEO
BASE = os.environ.get("ROUTE_E_GEO", "https://ftp.ncbi.nlm.nih.gov/geo")


def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"   already here: {os.path.basename(dest)}"); return True
    r = subprocess.run(["curl", "-fsSL", "-m", "3600", "-o", dest + ".part", url], capture_output=True, text=True)
    if r.returncode != 0:
        if os.path.exists(dest + ".part"): os.remove(dest + ".part")
        return False
    os.replace(dest + ".part", dest)
    print(f"   downloaded {os.path.basename(dest)} ({os.path.getsize(dest) / 1e6:.1f} MB)"); return True


def listing(url, tries=3):
    for _ in range(tries):
        r = subprocess.run(["curl", "-fsSL", "-m", "180", url], capture_output=True, text=True)
        if r.returncode == 0:
            return [h for h in re.findall(r'href="([^"?/][^"]*)"', r.stdout) if not h.startswith("http")]
    return []


def series_dir(acc):
    return f"{BASE}/series/{acc[:-3]}nnn/{acc}/"


def data_files(acc):
    """Supplementary files of a series: at series level, under its sub-series, or per sample."""
    sdir = series_dir(acc) + "suppl/"
    names = [h for h in listing(sdir) if not h.endswith("/") and h != "filelist.txt"]
    if names:
        return [(sdir + n, n) for n in names]
    out = []
    mdir = series_dir(acc) + "matrix/"
    for mn in [h for h in listing(mdir) if h.endswith("_series_matrix.txt.gz")]:
        mp = os.path.join(GDIR, mn)
        if not fetch(mdir + mn, mp):
            continue
        with gzip.open(mp, "rt", errors="replace") as h:
            head = h.read(400000)
        subs = re.findall(r"SuperSeries of: (GSE\d+)", head)
        for sub in dict.fromkeys(subs):                       # a SuperSeries keeps its files in the sub-series
            print(f"   {acc} is a SuperSeries; looking inside {sub}")
            out += data_files(sub)
        if not subs:
            gsms = dict.fromkeys(re.findall(r"GSM\d+", head))
            print(f"   no files at series level; checking {len(gsms)} samples one by one")
            for g in gsms:
                gd = f"{BASE}/samples/{g[:-3]}nnn/{g}/suppl/"
                out += [(gd + n, n) for n in listing(gd, tries=1) if not n.endswith("/") and n != "filelist.txt"]
    return out


EBI = os.environ.get("ROUTE_E_EBI", "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments")


def ebi_files(acc):
    """Some studies put only raw sequencing on GEO. EBI's Single Cell Expression Atlas holds the ready-made matrix."""
    e = "E-GEOD-" + re.sub(r"\D", "", acc)
    d = f"{EBI}/{e}/"
    names = [n for n in listing(d) if re.search(r"\.(mtx|mtx_rows|mtx_cols|tsv)$", n)]
    keep = [n for n in names if ".mtx" in n or "cell_metadata" in n or n.startswith("ExpDesign")]
    if keep:
        print(f"   GEO has no processed files; using EBI Single Cell Expression Atlas entry {e}")
    return [(d + n, n) for n in keep]


found = data_files(GSE) or ebi_files(GSE)
print(f"Files to download: {len(found)}")
for _, n in found[:12]:
    print("   ", n)
got = [os.path.join(GDIR, n) for u, n in found if fetch(u, os.path.join(GDIR, n))]
if not got:
    got = [os.path.join(GDIR, n) for n in sorted(os.listdir(GDIR))
           if not n.endswith((".part", ".csv", ".png")) and "series_matrix" not in n]
    if got:
        print("GEO did not answer; using files already in your folder:", [os.path.basename(g) for g in got])
if not got:
    raise StopHere(f"No data files found for {GSE} at the series, sub-series or sample level. "
                   "GEO may be busy: try running this cell again. If it keeps failing, send me this message.")
for p in [g for g in got if g.lower().endswith(".tar")]:        # unpack per-sample bundles
    with tarfile.open(p) as tf:
        tf.extractall(os.path.join(GDIR, "unpacked"))
    print(f"   unpacked {os.path.basename(p)}")
files = []
for root, _, names in os.walk(GDIR):
    files += [os.path.join(root, n) for n in names
              if not n.endswith((".part", ".tar")) and "series_matrix" not in n]
print(f"\n{len(files)} data files available.")
