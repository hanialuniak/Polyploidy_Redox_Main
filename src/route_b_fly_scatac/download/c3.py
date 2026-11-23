#@title Cell 3 - Download the chosen runs into your Drive folder
chosen = [r.strip().strip("/") for r in RUNS.split(",") if r.strip()]
if not chosen:
    raise StopHere("RUNS in Cell 1 is empty. Type run names from the Cell 2 list, run Cell 1 again, then this cell.")
unknown = [r for r in chosen if r not in set(files.run)]
if unknown:
    raise StopHere(f"Not in the Cell 2 list: {unknown}. Copy the names exactly.")
need = files[files.run.isin(chosen)]
print(f"About {need.MB.sum() / 1000:.1f} GB to download. Check your Google Drive has that much free space.")
for run, g in need.groupby("run"):
    os.makedirs(os.path.join(DATA_DIR, run), exist_ok=True)
    for name in WANTED:
        hits = [f for f in g.file if f.endswith(name)]
        if not hits:
            print(f"   {run}: no {name} on the site")
            continue
        url = BASE_URL + hits[0]
        dest = os.path.join(DATA_DIR, run, name)
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"   {run}/{name}: already downloaded")
            continue
        r = subprocess.run(["curl", "-fsSL", "-C", "-", "-m", "7200", "-o", dest + ".part", url],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"   {run}/{name}: FAILED ({r.stderr.strip()[:150]}). Run this cell again to resume.")
            continue
        os.replace(dest + ".part", dest)
        print(f"   {run}/{name}: {os.path.getsize(dest) / 1e6:.0f} MB")
print("\nDone. In the main notebook (RouteB_fly_ploidy_redox), set RUNS in Cell 1 to:", ", ".join(chosen))
