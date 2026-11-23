#@title Cell 2 - Download from GEO (sample descriptions and the expression files)
BASE = os.environ.get("ROUTE_D_GEO", "https://ftp.ncbi.nlm.nih.gov/geo")


def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"   already here: {os.path.basename(dest)}"); return True
    r = subprocess.run(["curl", "-fsSL", "-m", "3600", "-o", dest + ".part", url], capture_output=True, text=True)
    if r.returncode != 0:
        if os.path.exists(dest + ".part"): os.remove(dest + ".part")
        return False
    os.replace(dest + ".part", dest)
    print(f"   downloaded {os.path.basename(dest)} ({os.path.getsize(dest) / 1e6:.1f} MB)"); return True


def listing(url):
    r = subprocess.run(["curl", "-fsSL", "-m", "120", url], capture_output=True, text=True)
    return [h for h in re.findall(r'href="([^"?/][^"]*)"', r.stdout) if not h.startswith("http")] if r.returncode == 0 else []


stub = GSE[:-3] + "nnn"
mdir, sdir = f"{BASE}/series/{stub}/{GSE}/matrix/", f"{BASE}/series/{stub}/{GSE}/suppl/"
matrices = [os.path.join(GDIR, n) for n in ([h for h in listing(mdir) if h.endswith("_series_matrix.txt.gz")]
                                           or [f"{GSE}_series_matrix.txt.gz"]) if fetch(mdir + n, os.path.join(GDIR, n))]
if not matrices:
    raise StopHere(f"Could not download the sample descriptions from {mdir}. Check GEO opens in your browser, then run again.")
suppl = [h for h in listing(sdir) if not h.endswith("/")]
print("Supplementary files on GEO:", suppl or "none listed")
got = [os.path.join(GDIR, n) for n in suppl if fetch(sdir + n, os.path.join(GDIR, n))]
if not got:   # GEO listing failed this time: use files downloaded on an earlier run
    got = [os.path.join(GDIR, n) for n in sorted(os.listdir(GDIR))
           if n.startswith(GSE) and "series_matrix" not in n and not n.endswith(".part")
           and re.search(r"\.(tar|txt|tsv|csv|tab|xlsx?)(\.gz)?$", n, re.I)]
    if got:
        print("Using files already in your folder:", [os.path.basename(g) for g in got])
if not got:
    raise StopHere("No expression files could be downloaded. Send me this message and the list above.")
if not os.path.exists(GTF_PATH):
    raise StopHere(f"Gene annotation not found at {GTF_PATH}. It is the same file used in Route B; fix GTF_PATH in Cell 1.")
