#@title Cell 2 - Download the GEO series matrix and the platform annotation
BASE = os.environ.get("ROUTE_CR_GEO", "https://ftp.ncbi.nlm.nih.gov/geo")


def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"   already here: {os.path.basename(dest)}")
        return True
    r = subprocess.run(["curl", "-fsSL", "-m", "1800", "-o", dest + ".part", url], capture_output=True, text=True)
    if r.returncode != 0:
        if os.path.exists(dest + ".part"):
            os.remove(dest + ".part")
        return False
    os.replace(dest + ".part", dest)
    print(f"   downloaded {os.path.basename(dest)} ({os.path.getsize(dest) / 1e6:.1f} MB)")
    return True


def listing(url):
    r = subprocess.run(["curl", "-fsSL", "-m", "120", url], capture_output=True, text=True)
    return re.findall(r'href="([^"]+)"', r.stdout) if r.returncode == 0 else []


stub = GSE[:-3] + "nnn"
mdir = f"{BASE}/series/{stub}/{GSE}/matrix/"
names = [h for h in listing(mdir) if h.endswith("_series_matrix.txt.gz")] or [f"{GSE}_series_matrix.txt.gz"]
matrices = []
for n in names:
    dest = os.path.join(DATA_DIR, n)
    if fetch(mdir + n, dest):
        matrices.append(dest)
if not matrices:
    raise StopHere(f"Could not download the series matrix from {mdir}. Check that GEO opens in your browser and run again.")

with gzip.open(matrices[0], "rt", errors="replace") as h:
    head = h.read(200000)
gpl = re.search(r'!Series_platform_id\s+"?(GPL\d+)', head).group(1)
gstub = gpl[:-3] + "nnn" if len(gpl) > 6 else "GPLnnn"
annot = os.path.join(DATA_DIR, f"{gpl}.annot.gz")
full = os.path.join(DATA_DIR, f"{gpl}_full.txt")
got = fetch(f"{BASE}/platforms/{gstub}/{gpl}/annot/{gpl}.annot.gz", annot)
if not got:
    got = fetch(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=data", full)
if not got:
    raise StopHere(f"Could not download the annotation for platform {gpl}. Send me this message.")
print("Platform:", gpl)
