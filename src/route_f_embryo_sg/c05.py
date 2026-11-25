#@title Cell 5 - Gene symbols, the two contrasts, and the positive control
if DEC:
    print("Symbols came from the decorated rows file, so no GTF conversion is needed.")
elif not os.path.exists(GTF_PATH):
    raise StopHere(f"Gene annotation not found at {GTF_PATH}. It is the same file Route B and Route D use; "
                   "fix GTF_PATH in Cell 1.")
id2sym = {}
if not DEC:
    with gzip.open(GTF_PATH, "rt") as h:
        for line in h:
            if "\tgene\t" not in line:
                continue
            gid = re.search(r'gene_id "([^"]+)"', line); gnm = re.search(r'gene_name "([^"]+)"', line)
            if gid and gnm:
                id2sym[gid.group(1)] = gnm.group(1)
L.index = [id2sym.get(re.sub(r"\.\d+$", "", i), i) for i in L.index]
L = L.groupby(level=0).sum()
print(f"Genes after converting to symbols: {len(L):,}")

SG, OTH = "salivary gland", "other embryonic cells"
EARLY, LATE = "early (stages 10-12)", "late (stages 13-16)"
need = [f"{t} | {w}" for t in (SG, OTH) for w in (EARLY, LATE)]
miss = [c for c in need if c not in L.columns]
if miss:
    raise StopHere(f"Missing pseudobulk groups: {miss}. Available: {list(L.columns)}")

fc_sg = L[f"{SG} | {LATE}"] - L[f"{SG} | {EARLY}"]          # change across the window, in salivary gland
fc_oth = L[f"{OTH} | {LATE}"] - L[f"{OTH} | {EARLY}"]       # the same change in every other cell type
fc = fc_sg - fc_oth                                          # what is specific to the salivary gland
print("\nThree contrasts are available:")
print("   fc_sg   late vs early, salivary gland only")
print("   fc_oth  late vs early, all other cells - the embryo-wide maturation background")
print("   fc      fc_sg minus fc_oth - the SALIVARY-GLAND-SPECIFIC change, used for the screens below")

ALIASES = {"geminin": ["Gem", "CG3183"], "dup": ["Cdt1"], "Zw": ["G6PD", "G6pd"], "Ero1L": ["Ero1", "CG1333"],
           "PHGPx": ["CG12013"], "Gss1": ["Gss"], "Gr1": ["Gr", "GR"], "Trxr-1": ["Trxr1"], "Jafrac1": ["Prx2"],
           "Su(var)3-3": ["dLsd1", "Lsd1", "Su-var(3)3", "Suvar(3)3", "CG17149"], "Rrp1": ["APEX1", "CG3178"],
           "Sec61alpha": ["Sec61"], "SrpRbeta": ["SrpRb"]}


def resolve(g, present):
    if g in present:
        return g
    for a in ALIASES.get(g, []):
        if a in present:
            return a
    return None


PRESENT = set(L.index)
CONTROLS = {
    "POSITIVE - CrebA secretory pathway targets (must RISE in the gland)":
        "Sec61beta Sec61alpha SrpRalpha SrpRbeta Spase25 p24-1 zetaCOP Sec23 Sec31 Srp54 Sar1",
    "POSITIVE - salivary gland identity (should stay on)": "fkh sage CrebA Pasilla hth",
    "NEGATIVE - mitotic genes (gland leaves the mitotic cycle early)": "CycB polo aurB stg Cdk1 CycA pav tum",
}
print("\n=== Positive control: does the gland build its secretory machinery across these windows? ===")
for name, gs in CONTROLS.items():
    g = [r for x in gs.split() if (r := resolve(x, PRESENT))]
    lost = [x for x in gs.split() if not resolve(x, PRESENT)]
    if lost:
        # say WHY each one is missing, so a failed control can be diagnosed
        filtered = [x for x in lost if resolve(x, ALL_GENES)]
        absent = [x for x in lost if x not in filtered]
        if filtered:
            print(f"   {name.split(' - ')[0]}: removed by the gene filter -> {', '.join(filtered)} "
                  f"(lower MIN_COUNTS_PER_GROUP in Cell 1 to keep them)")
        if absent:
            print(f"   {name.split(' - ')[0]}: not in this matrix at all -> {', '.join(absent)}")
    if len(g) < 3:
        print(f"{name}: only {len(g)} genes found, skipped"); continue
    print(f"\n{name}  ({len(g)} genes: {', '.join(g)})")
    print(f"   salivary gland, late vs early: {fc_sg.loc[g].median():+.2f}")
    print(f"   other cells, late vs early:    {fc_oth.loc[g].median():+.2f}")
    print(f"   gland-specific:                {fc.loc[g].median():+.2f}")
crebA = [r for x in CONTROLS[list(CONTROLS)[0]].split() if (r := resolve(x, PRESENT))]
if len(crebA) >= 3 and fc_sg.loc[crebA].median() <= 0:
    print("\nWARNING: the CrebA secretory targets do not rise in the gland across these windows. That is the one "
          "thing this comparison must show, so the staging or the cell-type annotation is probably wrong. "
          "Treat everything below as uninterpretable until this is resolved.")
