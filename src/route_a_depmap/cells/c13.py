#@title Cell 13 - Findings, ranked by strength of evidence
def passes(row):
    return (row["q_fdr"] < FDR_MAX) and (row["p_vs_random"] < RANDOM_P_MAX) and (row["p_vs_controls"] < CONTROL_P_MAX)


results = pd.concat([c["result"] for c in (rna, wgd, prot) if c is not None], ignore_index=True)
redox = results[results.programme.str.startswith("REDOX")].copy()
redox["passes_all_checks"] = redox.apply(passes, axis=1)
findings = []
for name, part in redox.groupby("programme", sort=False):
    get = lambda label: part[part.data == label].iloc[0] if (part.data == label).any() else None
    r, p, w = get("RNA"), get("Protein"), get("RNA, whole-genome doubling")
    passed = [row for row in (r, p) if row is not None and row["passes_all_checks"]]
    signs = {np.sign(row["partial_rho"]) for row in passed}
    if r is not None and p is not None and r["passes_all_checks"] and p["passes_all_checks"] and len(signs) == 1:
        level = "Strong: RNA and protein"
    elif passed:
        level = "Moderate: " + ("RNA" if r is not None and r["passes_all_checks"] else "protein") + " only"
    elif any(row is not None and row["q_fdr"] < FDR_MAX for row in (r, p)):
        level = "Weak: correlates, but not specific"
    else:
        level = "None detected"
    lead = r if r is not None else p
    findings.append({"programme": name[7:], "evidence": level,
                     "direction": "higher with ploidy" if lead["partial_rho"] > 0 else "lower with ploidy",
                     "rho_RNA": None if r is None else round(r["partial_rho"], 3),
                     "rho_protein": None if p is None else round(p["partial_rho"], 3),
                     "whole_genome_doubling": None if w is None else ("passes" if w["passes_all_checks"] else "does not pass"),
                     "closest_control_RNA": None if r is None else r["closest_control"]})
order = {"Strong": 0, "Moderate": 1, "Weak": 2, "None": 3}
findings = pd.DataFrame(findings)
findings = findings.sort_values("evidence", key=lambda s: s.str.split(":").str[0].map(order))
display(findings.reset_index(drop=True))
for _, row in findings.iterrows():
    print(f"{row['programme']}: {row['evidence']}" + ("" if row["evidence"].startswith("None") else f", {row['direction']}"))
