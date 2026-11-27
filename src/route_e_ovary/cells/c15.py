#@title Cell 15 - The held-out test: was the chromatin theme real, or did I build the set from the answers?
# Cell 14 says the chromatin and satellite-DNA set is enriched at q = 0.001, 10 of 17 genes shared against 2
# expected. That number cannot be taken at face value, and the reason is a bookkeeping problem rather than a
# statistical one: the set was written AFTER reading the Cell 13 output, and most of the genes in it were put
# there BECAUSE they appeared as hits. Testing whether genes chosen from the hit list are enriched in the hit
# list is circular, and would give q = 0.001 even if the theme were meaningless.
#
# The fix is to split each set in two. The genes below were visible in the printed Cell 11-13 output and so
# could have influenced how the set was written; they are the discovery half and carry no evidential weight.
# Everything else in the set is held out: it was added because it belongs to the same biology, not because it
# was seen to be a hit. If the held-out half is enriched too, the theme generalises beyond the genes it was
# built from. If it is not, the theme is an artefact of my reading and should be dropped.
N_RANDOM_SETS_15 = 4000      #@param {type:"integer"}

# Every gene name that appeared in output I had already read when the sets in Cell 14 were written.
SEEN = set("""Vm34Ca Yp1 Yp2 Yp3 Vm26Aa Vm26Ab sosie CG8303 veil ImpL2 Rtnl1 CrebA CG9338 wbl sm CG4670
Tsp74F CG6983 CG7194 luna spri CG11200 CG8602 CG42268 fzr Reg-5 CG5367 CG5112 Gdh CG1637 Picot nrv1 Inos
CG5885 SsRbeta CG11899 CG4325 l(1)G0320 Gale Sdc stai CG9932 stg Uvrag CG10512 ct hdc Fas2 CG1943 Nnf1a
SPARC msi fax eya CG12560 Rrp1 Cbs ball fru HmgD dgt2 E23 Df31 CG31808 CG42566 piwi apt Obp99a jigr1
CG31705 mre11 CG4404 CG1265 Aldh Bub3 CG31915 CG9135 Det sle CG13252 AsnS mirr corn CG34348 CG43313 Dip-B
Lam CG17271 spi MYPT-75D CG8417 Atox1 Lsd-2 CG9953 Hnf4 CG8177 CG9132 AnxB10 CG5567 Spp Aprt wac CG8320
D1 CG4572 CG31098 Lmpt ogre Dp1 Npc2a HipHop Mal-A5 Tudor-SN Aos1 RpS5b Cdk1 RPA2 HP4 Sptr CG8490 Men
Rbfox1 CG7872 dUTPase Vha55 CtBP MtnA RfC4 Vha68-2 peb Srp9 CG31937 GstE1 mod(mdg4) Karybeta3 geminin
Jafrac1 Jafrac2 Ero1L Gclc Trx-2 Trxr-1 Cat Hsc70-3 Sec61beta""".split())

rows15 = []
for name, gs in THEMES.items():
    gs = [g for g in gs if g in elig.index]
    disc = [g for g in gs if g in SEEN]
    held = [g for g in gs if g not in SEEN]
    for half, members in (("discovery (no weight)", disc), ("HELD OUT", held)):
        if len(members) < 3:
            rows15.append({"theme": name, "half": half, "genes": len(members),
                           "note": "too few genes to test"})
            continue
        obs_h = sum(g in hit for g in members)
        null = np.array([sum(rng14.choice(pool14[bins14[g]]) in hit for g in members)
                         for _ in range(N_RANDOM_SETS_15)])
        rows15.append({"theme": name, "half": half, "genes": len(members), "shared": obs_h,
                       "expected": null.mean(),
                       "p": (1 + int((null >= obs_h).sum())) / (N_RANDOM_SETS_15 + 1),
                       "members": ", ".join(g for g in members if g in hit)[:60]})
ho = pd.DataFrame(rows15)
display(ho.round(3))

print("\nVerdict on each theme, reading ONLY the held-out half:")
for name in THEMES:
    r = ho[(ho.theme == name) & (ho.half == "HELD OUT")]
    if not len(r) or pd.isna(r.iloc[0].get("p")):
        n = len(r) and int(r.iloc[0].genes)
        print(f"   {name}: cannot be tested - only {n} held-out gene(s) in this dataset. The Cell 14 number "
              f"for it is circular and should not be quoted.")
        continue
    r = r.iloc[0]
    if r.p < 0.05 and int(r.shared) >= 3:
        print(f"   {name}: SURVIVES. {int(r.shared)}/{int(r.genes)} held-out genes shared against "
              f"{r.expected:.1f} expected, p = {r.p:.4f}. Members: {r.members}")
    elif r.p < 0.05:
        print(f"   {name}: too thin to call. Only {int(r.shared)} held-out gene(s) shared "
              f"({r.members}) against {r.expected:.1f} expected, p = {r.p:.4f}. Formally significant "
              f"because the expectation is near zero, but one or two genes is not a theme - read it as "
              f"those genes alone.")
    else:
        print(f"   {name}: DOES NOT survive. {int(r.shared)}/{int(r.genes)} held-out genes shared against "
              f"{r.expected:.1f} expected, p = {r.p:.3f}. The Cell 14 result for this set came from the genes "
              f"it was built out of, so drop the theme.")

print("""
What a surviving theme is and is not. It means: genes belonging to that biology, which were NOT the ones that
suggested the theme, also change in both lineages more often than expression-matched random genes do. That is
a real internal replication and it is the strongest thing these data can say. It is still one dataset, one
weak germline split of 197 cells, and a hypothesis generated and tested on the same cells - so it is grounds
for designing an experiment, not for a claim in a paper.
""")
ho.to_csv(os.path.join(GDIR, "route_e_theme_heldout.csv"), index=False)
print("Saved to", os.path.join(GDIR, "route_e_theme_heldout.csv"))
