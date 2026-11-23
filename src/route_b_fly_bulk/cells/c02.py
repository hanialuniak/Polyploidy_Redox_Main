#@title Cell 2 - Download the GEO files (a few minutes the first time)
GEO = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{GSE[:-3]}nnn/{GSE}"


def fetch(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        print(f"   {os.path.basename(path)}: already downloaded")
        return True
    r = subprocess.run(["curl", "-fsSL", "-m", "3600", "-o", path + ".part", url], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(path + ".part") or os.path.getsize(path + ".part") == 0:
        print(f"   {os.path.basename(path)}: FAILED ({r.stderr.strip()[:120]})")
        if os.path.exists(path + ".part"):
            os.remove(path + ".part")
        return False
    os.replace(path + ".part", path)
    print(f"   {os.path.basename(path)}: {os.path.getsize(path) / 1e6:.1f} MB")
    return True


fetch(f"{GEO}/matrix/{GSE}_series_matrix.txt.gz", os.path.join(DATA_DIR, f"{GSE}_series_matrix.txt.gz"))
page = subprocess.run(["curl", "-fsSL", "-m", "120", f"{GEO}/suppl/"], capture_output=True, text=True).stdout
names = sorted(set(re.findall(r'href="(GSE[^"/]+)"', page)))
print("\nSupplementary files on GEO:", names or "(could not list; files already in your folder will be used)")
tables = [n for n in names if not n.endswith(".tar") and re.search(r"count|fpkm|tpm|rpkm|express|matrix|txt|csv|tsv|xls", n, re.I)]
for n in tables:
    fetch(f"{GEO}/suppl/{n}", os.path.join(DATA_DIR, n))
if not tables and f"{GSE}_RAW.tar" in names:
    print("No combined table; downloading the per-sample archive instead")
    fetch(f"{GEO}/suppl/{GSE}_RAW.tar", os.path.join(DATA_DIR, f"{GSE}_RAW.tar"))
for tar in glob.glob(os.path.join(DATA_DIR, "*.tar")):
    out = os.path.join(DATA_DIR, "raw")
    if not os.path.isdir(out):
        with tarfile.open(tar) as t:
            t.extractall(out)
    print(f"   unpacked {os.path.basename(tar)}: {len(os.listdir(out))} files")
if not os.path.exists(os.path.join(DATA_DIR, f"{GSE}_series_matrix.txt.gz")):
    raise StopHere("The sample information file could not be downloaded. Run this cell again; if it still fails, send me the output.")
if not os.path.exists(GTF_PATH):
    raise StopHere(f"Gene annotation not found at {GTF_PATH}. Fix GTF_PATH in Cell 1.")
print("Done.")
