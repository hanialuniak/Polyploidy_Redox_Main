#@title Cell 2 - List the runs on the Aerts lab site (about a minute)
def listing(url):
    r = subprocess.run(["curl", "-fsSL", "-m", "60", url], capture_output=True, text=True)
    if r.returncode != 0:
        return None, r.stderr.strip()[:200]
    links = re.findall(r'href="([^"?#]+)"', r.stdout)
    out = []
    for h in links:
        full = urllib.parse.urljoin(url, h)
        if full.startswith(url) and full != url:
            out.append(full)
    return sorted(set(out)), ""


def size_of(url):
    r = subprocess.run(["curl", "-fsSLI", "-m", "60", url], capture_output=True, text=True)
    m = re.findall(r"(?i)content-length:\s*(\d+)", r.stdout)
    return int(m[-1]) / 1e6 if m else float("nan")


top, err = listing(BASE_URL)
if top is None:
    raise StopHere(f"Could not open {BASE_URL} ({err}). Check the address opens in your browser, then run again.")
found, queue, seen = [], [(u, 0) for u in top], set()
while queue:
    u, depth = queue.pop(0)
    if u in seen:
        continue
    seen.add(u)
    name = u.rstrip("/").split("/")[-1]
    if not u.endswith("/"):
        if name in WANTED or name.endswith(("fragments.tsv.gz", "singlecell.csv")):
            found.append(u)
        continue
    if depth < 4:
        sub, _ = listing(u)
        queue += [(s, depth + 1) for s in (sub or [])]

rows = []
for u in found:
    rel = u[len(BASE_URL):]
    folder = re.sub(r"(^|/)outs$", "", os.path.dirname(rel)).strip("/")   # one sequencing run per folder
    run = folder.replace("/", "_") or "top"
    rows.append({"run": run, "file": rel, "MB": round(size_of(u))})
files = pd.DataFrame(rows)
if files.empty:
    print("Top-level entries:", [t[len(BASE_URL):] for t in top][:60])
    raise StopHere("No fragments.tsv.gz or singlecell.csv files found under this address. Send me the list above.")
pd.set_option("display.max_rows", 500)
display(files)
summary = files.groupby("run").agg(files=("file", "size"), MB=("MB", "sum"))
display(summary)
print("\nRuns that have both files:", [r for r, g in files.groupby("run")
                                       if any(f.endswith("fragments.tsv.gz") for f in g.file)
                                       and any(f.endswith("singlecell.csv") for f in g.file)])
print("Type the runs you want into RUNS in Cell 1 (comma-separated), run Cell 1 again, then Cell 3.")
