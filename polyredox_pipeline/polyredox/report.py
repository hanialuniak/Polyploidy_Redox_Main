"""Write a plain-language Markdown summary of one run."""
from __future__ import annotations

import os

import pandas as pd


def _table(frame, digits=3):
    if frame.empty:
        return "_none_\n"
    shown = frame.copy()
    for column in shown.columns:
        if pd.api.types.is_float_dtype(shown[column]):
            shown[column] = shown[column].map(lambda v: "" if pd.isna(v) else f"{v:.{digits}g}")
    head = "| " + " | ".join(map(str, shown.columns)) + " |"
    rule = "|" + "---|" * len(shown.columns)
    rows = ["| " + " | ".join(map(str, row)) + " |" for row in shown.to_numpy()]
    return "\n".join([head, rule, *rows]) + "\n"


def write(outdir, config, frame, excess, verdict, interaction, summary, importance, cross, skipped, purity, libraries):
    lines = ["# Polyploidy and redox: run summary\n",
             "All findings are associations in post-mortem nuclei. They do not show direction or causation.\n",
             f"Analysed {len(frame):,} nuclei from {frame['donor'].nunique()} donors.\n",
             "## Is the polyploid signal more than doublets?\n", _table(excess)]
    analysed = excess[excess["cell_class"].isin(frame["cell_class"].unique())]
    if (analysed["excess_ci_low"] <= 0).any():
        lines.append("**Warning: in at least one analysed class, the polyploid signal cannot be separated from "
                     "expected same-class doublets. Treat its results as unresolved.**\n")
    if libraries["n_recovered_estimated"].any():
        lines.append("**Warning: n_recovered was missing for some libraries, so nuclei in the cells table were used. "
                     "That underestimates doublets if the table was filtered.**\n")
    lines.append(f"Estimated share of the signal not explained by doublets, across analysed nuclei: {purity:.2f}\n")
    lines.append("## Redox programmes, within donor\n")
    columns = ["programme", "beta", "ci_low", "ci_high", "p_signflip", "p_adj", "p_matched",
               "p_vs_worst_stress_control", "worst_stress_control", "beta_unadjusted", "beta_depth_iqr", "effective_donors", "verdict",
               "shares_20pct_genes_with"]
    lines.append(_table(verdict[columns]))
    lines.append("Programmes that share many genes are not independent results; read them together.\n")
    lines.append("`beta` is the change in programme score, in within-class standard deviations, for polyploid "
                 "versus diploid.\n")
    if not interaction.empty:
        lines.append("## Polyploidy x pathology interaction\n")
        lines.append(_table(interaction[["programme", "beta", "ci_low", "ci_high", "p_signflip", "p_adj"]]))
    lines.append("## Machine learning, split by donor\n")
    lines.append("This tests whether redox scores predict ploidy beyond technical features. It inherits the same "
                 "depth confounding, so read it only alongside the verdicts above.\n")
    lines.append(_table(pd.DataFrame([summary])))
    lines.append("### Programme importance\n")
    lines.append(_table(importance))
    if not cross.empty:
        lines.append("### Held-out cohort\n")
        lines.append(_table(cross))
    if skipped:
        lines.append(f"Gene sets skipped for too few detected genes: {skipped}\n")
    with open(os.path.join(outdir, "report.md"), "w") as handle:
        handle.write("\n".join(lines))
