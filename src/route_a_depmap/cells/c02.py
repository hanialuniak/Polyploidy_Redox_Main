#@title Cell 2 - Check the data files in your Drive folder
NEEDED = {
    "models": ["Model.csv"],
    "signatures": ["OmicsSignatures.csv"],
    "expression": ["OmicsExpressionProteinCodingGenesTPMLogp1.csv"],
    "damaging": ["OmicsSomaticMutationsMatrixDamaging.csv"],
    "hotspot": ["OmicsSomaticMutationsMatrixHotspot.csv"],
}
if USE_CRISPR:
    NEEDED["crispr"] = ["CRISPRGeneEffect.csv"]
if USE_PROTEINS:
    NEEDED["proteins"] = ["protein_quant_current_normalized.csv.gz", "protein_quant_current_normalized.csv"]

present = sorted(os.listdir(DATA_DIR))
paths, missing = {}, []
for key, names in NEEDED.items():
    found = next((n for n in names if os.path.exists(os.path.join(DATA_DIR, n))), None)
    if found:
        paths[key] = os.path.join(DATA_DIR, found)
        print(f"   found   {found}  ({os.path.getsize(paths[key]) / 1e6:.0f} MB)")
    else:
        stem = names[0].split(".")[0]
        similar = [p for p in present if p.startswith(stem)]
        missing.append(names[0])
        print(f"   MISSING {names[0]}" + (f"  (similar name in folder: {similar[0]} - rename it)" if similar else ""))
if missing:
    raise StopHere(f"{len(missing)} file(s) missing from {DATA_DIR}. Add or rename them, then run this cell again. "
                   "(Untick USE_CRISPR or USE_PROTEINS in Cell 1 to run without those files.)")
print("All files ready.")
