#@title Cell 9 - Screen 4: same tests with whole-genome doubling (yes/no)
wgd = test_link(expression, "RNA, whole-genome doubling", exposure="wgd")
test_random(wgd)
test_controls(wgd)
print(f"{int(wgd['table'].wgd.sum())} of {len(wgd['table'])} lines have a doubled genome")
display(wgd["result"].loc[wgd["result"].programme.str.startswith("REDOX"),
                          ["programme", "partial_rho", "q_fdr", "p_vs_random", "p_vs_controls"]].round(3))
