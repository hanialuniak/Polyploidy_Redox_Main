#@title Cell 10 - Screen 5: same tests on proteins
prot = None
if proteins is not None and len(proteins) >= 25:
    prot = test_link(proteins, "Protein")
    test_random(prot)
    test_controls(prot)
    display(prot["result"].loc[prot["result"].programme.str.startswith("REDOX"),
                               ["programme", "lines", "genes_used", "partial_rho", "q_fdr", "p_vs_random",
                                "p_vs_controls"]].round(3))
else:
    print("Proteins skipped: switched off in Cell 1 or fewer than 25 matching cell lines.")
