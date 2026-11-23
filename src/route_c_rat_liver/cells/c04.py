#@title Cell 4 - Checks: do samples group by ploidy, the cell-cycle control, and liver zonation
lp = np.log2(samples.ploidy.to_numpy())
C = np.corrcoef((e.sub(e.mean(axis=1), axis=0)).T.to_numpy()); np.fill_diagonal(C, -1)
nn = C.argmax(axis=1)
print(f"Nearest-neighbour sample has the same ploidy for {np.mean(samples.ploidy.to_numpy()[nn] == samples.ploidy.to_numpy()):.0%} of samples")
PAIRED = samples.animal.notna().all() and samples.animal.nunique() < len(samples)
if PAIRED:
    same_animal = np.mean(samples.animal.to_numpy()[nn] == samples.animal.to_numpy())
    print(f"Nearest-neighbour sample has the same replicate number for {same_animal:.0%} of samples "
          f"(chance about {1 / samples.animal.nunique():.0%}).")
    print("Samples with the same number are treated as coming from the same animal: every comparison is made within animal.")

Z = e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1).replace(0, 1), axis=0)


def sample_score(genes):
    g = [x for x in genes if x in Z.index]
    return (Z.loc[g].mean(), len(g)) if len(g) >= 3 else (pd.Series(np.nan, index=Z.columns), len(g))


CHECKS = {"cell cycle (positive control: up in 8c, reported by the authors)":
          "Mki67 Top2a Ccnb1 Cdk1 Ccna2 Bub1 Aurkb Plk1 Cdc20 Ccnb2 Kif11 Cenpf Nusap1 Prc1".split(),
          "pericentral (central vein) genes": "Glul Cyp2e1 Cyp1a2 Oat Lect2 Axin2 Rhbg Cyp7a1 Slc1a2 Cyp27a1 Lgr5".split(),
          "periportal genes": "Cps1 Ass1 Arg1 Hal Sds Pck1 Cyp2f4 Gls2 Hsd17b13 Aldob Igfbp2 Hamp Cyp8b1".split()}
chk = {}
for name, gs in CHECKS.items():
    s, n = sample_score(gs); chk[name] = s
    by = s.groupby(samples.ploidy.to_numpy()).mean()
    print(f"{name} ({n} genes): " + ", ".join(f"{int(k)}c {v:+.2f}" for k, v in by.items()))
zone = (chk["pericentral (central vein) genes"] - chk["periportal genes"]).to_numpy()
r_zone = np.corrcoef(zone, lp)[0, 1]
print(f"\nZonation score (pericentral minus periportal) vs ploidy: r = {r_zone:+.2f}")
if abs(r_zone) > 0.5:
    print("Ploidy groups differ in zonation. Cell 5 therefore also reports results adjusted for the zonation score.")
