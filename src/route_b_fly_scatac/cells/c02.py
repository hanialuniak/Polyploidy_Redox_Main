#@title Cell 2 - Check the data files in your Drive folder
paths, missing = {}, []
for run in RUN_LIST:
    for name in ("fragments.tsv.gz", "singlecell.csv"):
        p = os.path.join(DATA_DIR, run, name)
        if os.path.exists(p):
            print(f"   found   {run}/{name}  ({os.path.getsize(p) / 1e6:.0f} MB)")
            paths[(run, name)] = p
        else:
            print(f"   MISSING {run}/{name}")
            missing.append(f"{run}/{name}")
gtf = os.path.join(DATA_DIR, GTF_FILE)
if os.path.exists(gtf):
    print(f"   found   {GTF_FILE}  ({os.path.getsize(gtf) / 1e6:.0f} MB)")
else:
    print(f"   MISSING {GTF_FILE}")
    missing.append(GTF_FILE)
if missing:
    raise StopHere(f"{len(missing)} file(s) missing. Check folder and file names (see the guide), then run this cell again.")
print("All files ready.")
