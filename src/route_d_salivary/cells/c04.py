#@title Cell 4 - Normalise, and check samples and the positive controls
vals = expr.to_numpy(float)
is_counts = np.nanmax(vals) > 1000 and np.allclose(vals[np.isfinite(vals)][:5000] % 1, 0)
if is_counts:
    keep = (expr / expr.sum() * 1e6 >= 1).sum(axis=1) >= info.groupby("group").size().min()
    x = expr[keep]
    lg = np.log(x.replace(0, np.nan)); ref = lg.mean(axis=1)
    ok = ref.notna()
    size = np.exp((lg[ok].sub(ref[ok], axis=0)).median())               # median-of-ratios size factors
    L = np.log2(x / size / x.sum().mean() * 1e6 + 1)
    print("Raw counts detected: median-of-ratios normalisation, log2 CPM.")
else:
    L = np.log2(expr.clip(lower=0) + 1)
    L = L[(L > 1).sum(axis=1) >= info.groupby("group").size().min()]
    print("Already-normalised values detected: log2(x + 1).")
print(f"Expressed genes kept: {len(L):,}")

C = np.corrcoef(L.sub(L.mean(axis=1), axis=0).T.to_numpy()); np.fill_diagonal(C, -1)
same = np.mean(info.group.to_numpy()[C.argmax(axis=1)] == info.group.to_numpy())
print(f"Nearest-neighbour sample is from the same group for {same:.0%} of samples")

CONTRASTS = [(t, r) for t in info.loc[info.role == "test", "group"].unique() for r in info.loc[info.role == "reference", "group"].unique()]
gmean = {g: L[info.column[info.group == g]].mean(axis=1) for g in info.group.unique()}
POS = {"mitotic genes (Myb-MuvB targets): expected DOWN in endocycling cells":
       "CycB polo aurB stg Cdk1 CycA Det pav tum feo sti asp Klp61F Incenp BubR1 Ncd80 Nuf2 CycB3 mad2",
       "E2F1 S-phase genes (expected down in salivary glands, not in S2 iECs)":
       "PCNA RnrS RnrL Mcm2 Mcm3 Mcm5 Mcm6 Mcm7 dpa DNApol-alpha50 Orc1 dup Cdc6 E2f1"}
for name, gs in POS.items():
    g = [x for x in gs.split() if x in L.index]
    print(f"\n{name} ({len(g)} genes)")
    for t, r in CONTRASTS:
        print(f"   {t} vs {r}: median log2 change {(gmean[t] - gmean[r]).loc[g].median():+.2f}")
