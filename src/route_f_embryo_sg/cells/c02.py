#@title Cell 2 - Download the matrix and the cell annotations
EBI = os.environ.get("ROUTE_F_EBI", "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments")
BASE = os.environ.get("ROUTE_F_GEO", "https://ftp.ncbi.nlm.nih.gov/geo")


def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"   already here: {os.path.basename(dest)} ({os.path.getsize(dest) / 1e6:.0f} MB)"); return True
    r = subprocess.run(["curl", "-fsSL", "-m", "7200", "-o", dest + ".part", url], capture_output=True, text=True)
    if r.returncode != 0:
        if os.path.exists(dest + ".part"):
            os.remove(dest + ".part")
        return False
    os.replace(dest + ".part", dest)
    print(f"   downloaded {os.path.basename(dest)} ({os.path.getsize(dest) / 1e6:.0f} MB)"); return True


def listing(url, tries=3):
    for _ in range(tries):
        r = subprocess.run(["curl", "-fsSL", "-m", "180", url], capture_output=True, text=True)
        if r.returncode == 0:
            return [h for h in re.findall(r'href="([^"?/][^"]*)"', r.stdout) if not h.startswith("http")]
    return []


# EBI's Single Cell Expression Atlas holds the ready-made matrix and the per-cell annotations.
# GEO carries only raw sequencing for this study, which is far too large for Colab.
acc = "E-GEOD-" + re.sub(r"\D", "", GSE)
d = f"{EBI}/{acc}/"
names = listing(d)
want = [n for n in names if re.search(r"\.mtx$|\.mtx_rows$|\.mtx_cols$", n)]
# prefer the filtered counts over raw, and take one matrix only
mtx = [n for n in want if n.endswith(".mtx")]
pick = next((m for m in mtx if "filtered" in m), mtx[0] if mtx else None)
want = [n for n in want if pick and n.startswith(pick.rsplit(".mtx", 1)[0])]
design = [n for n in names if n.startswith("ExpDesign") or "cell_metadata" in n or n.endswith("condensed-sdrf.tsv")]
found = [(d + n, n) for n in want + design]
if not found:
    raise StopHere(f"Nothing downloadable found at {d}. Open that address in your browser: if it is empty, the "
                   f"study may be listed under a different accession. Send me what you see.")
print(f"Files to download ({len(found)}):")
for _, n in found:
    print("   ", n)
got = {n: os.path.join(GDIR, n) for u, n in found if fetch(u, os.path.join(GDIR, n))}
if not got:
    got = {n: os.path.join(GDIR, n) for n in sorted(os.listdir(GDIR)) if not n.endswith(".part")}
    if got:
        print("Download failed; using files already in your folder:", list(got))
need = [k for k in ("mtx", "mtx_rows", "mtx_cols") if not any(n.endswith("." + k) for n in got)]
if need or not any(n.startswith("ExpDesign") or "cell_metadata" in n or "sdrf" in n for n in got):
    raise StopHere(f"Missing {need or 'the cell annotation table'}. Check your folder {GDIR} and run again.")
MTX = next(p for n, p in got.items() if n.endswith(".mtx"))
ROWS = next(p for n, p in got.items() if n.endswith(".mtx_rows"))
COLS = next(p for n, p in got.items() if n.endswith(".mtx_cols"))
META = next(p for n, p in got.items() if n.startswith("ExpDesign") or "cell_metadata" in n or "sdrf" in n)
print("\nUsing matrix:", os.path.basename(MTX), "| annotations:", os.path.basename(META))
