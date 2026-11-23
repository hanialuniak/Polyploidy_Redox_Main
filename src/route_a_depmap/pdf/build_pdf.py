import html, os, glob, asyncio
HERE = os.path.dirname(os.path.abspath(__file__))
CELLS = {int(os.path.basename(p)[1:3]): open(p).read().rstrip() for p in glob.glob(os.path.join(HERE, "..", "cells", "c*.py"))}
TITLE = "Route A: ploidy, redox and Geminin in cell lines"

def code(n, title):
    return f'<p class="cellname">Cell {n} - {html.escape(title)}</p><pre>{html.escape(CELLS[n])}</pre>'

def table(head, rows, cls=""):
    h = "".join(f"<th>{c}</th>" for c in head)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<table class="{cls}"><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table>'

def tested(text): return f'<div class="ok">{text}</div>'
def note(text): return f'<div class="warn">{text}</div>'

def glance(question, method, passes, fooled):
    return table(["", ""], [["<b>Question</b>", question], ["<b>How it works</b>", method],
                            ["<b>Passes when</b>", passes], ["<b>What can fool it</b>", fooled]], "glance")

body = f"""
<h1>{TITLE}: an exploratory step-by-step guide</h1>
<p class="sub">Screening which redox programmes and genes differ with ploidy in DepMap cancer cell lines, in RNA and
protein, and where Geminin fits. No hypotheses are set in advance. For Google Colab; written for a beginner.</p>

<h2>1. What this notebook explores</h2>
<p>Each cancer cell line in DepMap has a <b>ploidy</b> estimated from its DNA sequencing (the average number of copies of
its genome), plus gene expression, CRISPR knockout effects and, for about 375 lines, protein levels. The notebook compares
lines from the same tissue and asks, without a prediction, what changes with ploidy. Seven screens build on each other:</p>
{table(["Question", "How it is explored", "Cell"], [
 ["<b>Screen 1. Which redox programmes go with ploidy?</b>", "Correlate each programme's score with ploidy within each tissue, after removing proliferation; p-value by shuffling ploidy within tissue", "6"],
 ["<b>Screen 2. Is a signal stronger than random genes?</b>", "Repeat Screen 1 with 500 random gene sets of the same size and expression level", "7"],
 ["<b>Screen 3. Is it specific to redox?</b>", "Test whether the redox programme still goes with ploidy after subtracting each of four non-redox stress programmes", "8"],
 ["<b>Screen 4. Does it hold for genome doubling?</b>", "Repeat Screens 1-3 with doubled / not doubled instead of ploidy as a number", "9"],
 ["<b>Screen 5. Does it hold for proteins?</b>", "Repeat Screens 1-3 on CCLE mass-spectrometry protein levels", "10"],
 ["<b>Screen 6. Which individual genes drive it?</b>", "Every redox and control gene tested on its own, RNA and protein", "11"],
 ["<b>Screen 7. Where does Geminin fit?</b>", "GMNN and CDT1 expression, their balance and their knockout effects, against ploidy and against redox", "12"]])}

<h3>Starting point: Geminin</h3>
<p>The one thing known from the fly experiments is that <b>Geminin modulates endoreplication</b>. Geminin (human gene
<b>GMNN</b>) blocks <b>CDT1</b>, a protein that licenses DNA for replication, so that the genome is copied only once per
cycle. Losing that brake lets cells copy their DNA again without dividing. Screen 7 asks whether, in human cell lines:</p>
<ul>
<li>GMNN or CDT1 levels, or the balance between them, differ with ploidy;</li>
<li>lines with higher ploidy depend more or less on GMNN (CRISPR knockout effect);</li>
<li>GMNN goes with redox programmes, and whether that link remains once ploidy is taken into account.</li>
</ul>
{note("<b>A fly link may not show here.</b> Cancer cells become polyploid by several routes (failed cell division, skipped mitosis, re-replication), not mainly through endocycles like fly cells. GMNN is also cell-cycle regulated, so its level follows how fast lines divide. A null result in cell lines therefore does not argue against the fly finding.")}

<h3>What stays the same throughout</h3>
<p>These are analysis choices, not predictions. They apply to every screen.</p>
{table(["Item", "Setting"], [
 ["<b>Comparisons</b>", "Always within tissue of origin, adjusted for proliferation (22 cell-cycle genes). Follow-up runs can drop the proliferation adjustment or add TP53 status (section 13)"],
 ["<b>Lines left out</b>", "Lines with KEAP1, CUL3 or NFE2L2 mutations, which switch redox genes on regardless of ploidy; tissues with fewer than 10 lines"],
 ["<b>Redox programmes (7)</b>", "NRF2 targets, glutathione, thioredoxin/peroxiredoxin, antioxidant enzymes, ROS sources, ferroptosis, mitochondrial respiration"],
 ["<b>Stress controls (4)</b>", "Unfolded protein response, heat shock, DNA damage, immediate-early genes"],
 ["<b>Labelling thresholds</b>", "False discovery rate below 0.05; beats random gene sets (p below 0.05); beats every stress control (p below 0.05)"],
 ["<b>Genes of interest</b>", "GMNN and CDT1 (add others in Cell 1, for example other licensing genes)"]])}

<h3>The data</h3>
<p>All seven files sit in your Google Drive folder <b>Colab Notebooks/route_a_data</b>. The notebook reads them from
there; it does not download anything.</p>
{table(["File", "What it is", "Used in", "Source"], [
 ["<b>OmicsSignatures.csv</b> (DepMap 24Q4)", "Ploidy and whole-genome doubling per line, estimated with PureCN", "All screens", "DepMap portal, release 24Q4"],
 ["<b>OmicsExpressionProtein-CodingGenesTPMLogp1.csv</b>", "Gene expression, log TPM, about 19,000 genes", "Screens 1-4, 6, 7", "DepMap portal, release 24Q4"],
 ["<b>CRISPRGeneEffect.csv</b>", "Effect of knocking out each gene on growth; more negative = line depends more on it", "Screen 7", "DepMap portal, release 24Q4"],
 ["<b>Model.csv</b>; mutation matrices", "Tissue of origin, names; mutations", "Adjustment; exclusions", "DepMap portal, release 24Q4 (3 files)"],
 ["<b>protein_quant_current_normalized.csv.gz</b> (Gygi lab)", "Mass-spectrometry protein levels, 375 lines", "Screens 5, 6", "Gygi lab CCLE page, TSV format"]])}
{note("<b>These are cancer cell lines, not neurons.</b> Findings are leads for brain tissue (Routes B and C) and for the fly work, not conclusions. Correlations here do not show which causes which.")}

<h2>2. Before you start: note your settings (5 minutes)</h2>
<p>Exploratory does not need hypotheses, but it does need the analysis fixed before you see results, so that settings are
not changed until something looks interesting.</p>
<ol>
<li><b>Record the date and the Cell 1 settings</b> in your lab notebook or an email to your supervisor before running Cell 6.</li>
<li><b>Report every screen</b>, including those that find nothing. With 7 programmes, about 170 genes and several Geminin measures, some results will pass by chance.</li>
<li><b>If you change a setting after seeing results</b>, keep the first run and report both.</li>
<li><b>Treat any finding as a lead</b> to test in a new dataset (brain Multiome, flies), not as a result in itself.</li>
</ol>
{note("<b>Why the checks exist.</b> In a simulated dataset where ploidy had <i>no</i> effect on redox, but tissue and proliferation affected both, the raw correlation between NRF2 targets and ploidy was 0.56; within tissue and adjusted for proliferation it fell to 0.07. Keeping the KEAP1-mutant lines in produced a false link of 0.20 (false discovery rate 0.002). In the same simulation, GMNN correlated with three redox programmes (up to 0.11, false discovery rate 0.007) only because both follow proliferation; comparing GMNN with genes that follow proliferation equally closely removed every one.")}

<h2>3. Setup</h2>
<ol>
<li>Go to <b>colab.research.google.com</b>, choose <b>File &rarr; Upload notebook</b> and select <i>RouteA_ploidy_redox.ipynb</i>. Save a copy to Google Drive.</li>
<li>Check that DATA_DIR in Cell 1 reads <i>/content/drive/MyDrive/Colab Notebooks/route_a_data</i>. If your folder is elsewhere, change it there: <i>MyDrive</i> is the top of your Drive.</li>
<li>Run <b>Cell 1</b> (the play button on its left). Colab asks to connect to Google Drive: choose your account and click <b>Allow</b>.</li>
<li>Add any extra genes of interest (for example <i>GMNN, CDT1, CDC6, MCM2</i>) to GENES_OF_INTEREST.</li>
</ol>
<p>You can run one cell at a time, as this guide does, or everything with <b>Runtime &rarr; Run all</b>.</p>
{code(1, "Setup, settings and connect Google Drive")}
{tested("<b>Tested</b> on mock data in a local folder. The Drive connection itself runs only in Colab, so I could not test it here.")}
<p><b>What you should see:</b> a line ending <i>data folder: /content/drive/MyDrive/Colab Notebooks/route_a_data</i>. If it says <i>Folder not found</i>, the folder name or location differs: fix DATA_DIR.</p>

<h2>4. Step 1: check and match the data (Cells 2-3, 2-3 minutes)</h2>
<p>Cell 2 checks that all seven files are in your Drive folder and prints each one with its size. Cell 3 lines the files
up by cell line, removes NRF2-pathway mutants and keeps tissues with at least 10 lines. Reading the expression file from
Drive takes a minute or two.</p>
{code(2, "Check the data files in your Drive folder")}
{tested("<b>Tested</b> on mock files, including a missing file and one renamed by the browser to <i>Model (1).csv</i>; Cell 2 named both and stopped.")}
<p><b>What you should see:</b> seven lines starting <i>found</i>, then <i>All files ready.</i></p>
{table(["What you see", "What to do"], [
 ["<b>MISSING ... (similar name in folder: ... - rename it)</b>", "In Google Drive, right-click that file, choose Rename and use the exact name shown after MISSING"],
 ["<b>MISSING CRISPRGeneEffect.csv</b> or the protein file", "Add the file, or untick USE_CRISPR or USE_PROTEINS in Cell 1 to run without it"],
 ["<b>Protein file found as .csv rather than .csv.gz</b>", "Fine: your computer unzipped it; both work"]])}
{code(3, "Load and match the data")}
{tested("<b>Tested</b> on mock files in both the 24Q4 layout and the newer layout (extra columns and duplicate rows per cell line).")}
<p><b>What you should see:</b> lines such as <i>Leaving out 58 lines with KEAP1, CUL3 or NFE2L2 mutations</i>, <i>Cell lines analysed: 1122 from 30 tissues</i>, <i>... protein data: 389</i> and <i>... CRISPR knockout data: 1122</i> (simulated numbers; yours will differ).</p>

<h2>5. Step 2: build the scores (Cells 4-5)</h2>
<p>Cell 4 holds the gene lists. Cell 5 holds the tools every screen uses; it prints only <i>Tools ready</i>. Two ideas matter:</p>
<ul>
<li><b>Programme score.</b> Each gene becomes a z-score (how far a line sits from the average line, in standard deviations). A programme's score is the average z-score of its genes. Genes measured in fewer than 80% of lines are left out.</li>
<li><b>Adjustment.</b> Before any correlation, tissue of origin and proliferation are removed from both ploidy and the score. What is left is how a line differs from other lines of the same tissue that divide equally fast. The correlation of those leftovers is the <b>partial correlation</b> (partial rho), from -1 to 1.</li>
</ul>
{code(4, "Gene programmes")}
{code(5, "Scoring and adjustment tools (no results yet)")}
{tested("<b>Tested:</b> ran correctly on mock data built in the same file formats.")}
{table(["partial rho", "Rough size"], [["below 0.1", "very small; can be significant only because there are many lines"], ["0.1 to 0.3", "small to moderate"], ["0.3 to 0.5", "moderate to strong"], ["above 0.5", "strong"]])}

<h2>6. Screen 1: which programmes go with ploidy? (Cell 6)</h2>
{glance("Within the same tissue, and at equal proliferation, do lines with more genome copies have higher or lower scores, for any programme?",
 "Partial rank correlation between ploidy and each programme score. Ploidy is shuffled within each tissue 2,000 times; the p-value is the share of shuffles giving a correlation as large. The false discovery rate (q) corrects for screening 7 redox programmes.",
 "q below 0.05, in either direction",
 "Anything that varies with ploidy and redox together within a tissue: leftover proliferation, cell size, culture stress. Screens 2 and 3 exist for this.")}
{code(6, "Screen 1: does each programme go with ploidy, within tissue? (RNA)")}
{tested("<b>Tested</b> on two simulated datasets: one with a planted link for NRF2 targets and glutathione (partial rho 0.50 and 0.55; q = 0.001), and one with no link (largest partial rho 0.07, nothing below q = 0.05).")}
<p><b>What you should see:</b> one row per programme (11 rows) with <i>lines</i>, <i>genes_used</i>, <i>partial_rho</i>, <i>p_perm</i> and, for redox programmes, <i>q_fdr</i>. The smallest possible p_perm is 1/2001, printed as 0.0.</p>

<h2>7. Screen 2: stronger than random genes? (Cell 7)</h2>
{glance("Would any set of genes with similar expression levels correlate with ploidy just as much?",
 "For each redox programme, 500 random gene sets are drawn, the same size, each gene matched to a programme gene on average expression. Each goes through Screen 1. A real signal sits outside the range the random sets reach.",
 "p_vs_random below 0.05: fewer than 1 in 20 random sets do as well",
 "A broad shift affecting most genes, such as a change in RNA composition with ploidy. If redox genes move but so does everything else, this screen fails.")}
{code(7, "Screen 2: is it stronger than random genes?")}
{tested("<b>Tested.</b> In the planted simulation, random sets reached -0.04 to 0.08 against 0.50 for NRF2 targets.")}

<h2>8. Screen 3: specific to redox? (Cell 8)</h2>
{glance("Is a signal about redox, or do all stress responses change with ploidy?",
 "For each redox programme and each of the four stress controls, the control score is subtracted from the redox score and the difference goes through Screen 1. The largest of the four p-values is kept.",
 "p_vs_controls below 0.05 (all four differences pass, in the same direction as Screen 1)",
 "Little. This is the strictest screen, and it can miss a real but modest redox-specific link if a control also moves somewhat.")}
{code(8, "Screen 3: is it specific to redox, or general stress?")}
{tested("<b>Tested.</b> In a no-link simulation, Screen 3 rejected every programme (p_vs_controls 0.42 or more).")}
<p><b>What you should see:</b> <i>p_vs_controls</i> and <i>closest_control</i>, the control that moved most like the redox programme. A recurring closest control is itself a lead: unfolded protein response suggests protein stress; immediate-early genes suggest culture or handling stress.</p>
{note("<b>Programmes that share genes are not independent.</b> NRF2 targets share genes with glutathione, thioredoxin and ferroptosis. In the planted simulation, only NRF2 targets and glutathione were changed, yet thioredoxin and ferroptosis also passed in RNA (partial rho 0.14 and 0.17) through their shared genes. Screen 6 shows which genes carry a signal.")}

<h2>9. Screen 4: whole-genome doubling (Cell 9)</h2>
{glance("Does a signal hold if lines are split simply into doubled and not doubled?",
 "Screens 1-3 again, with whole-genome doubling (yes/no, from PureCN) in place of ploidy as a number.",
 "Same three thresholds",
 "Lines near the cut-off. It checks that a signal does not depend on how ploidy is measured.")}
{code(9, "Screen 4: same tests with whole-genome doubling (yes/no)")}
{tested("<b>Tested.</b> In the planted simulation, 434 of 1,122 lines were doubled; NRF2 targets gave partial rho 0.38 and passed all three checks.")}

<h2>10. Screen 5: proteins (Cell 10)</h2>
{glance("Do redox proteins, not only their RNA, change with ploidy?",
 "Screens 1-3 on CCLE mass-spectrometry protein levels. Proliferation is still measured from RNA. Proteins missing in more than 20% of lines are left out.",
 "Same three thresholds; check <i>genes_used</i> is at least 5 and <i>lines</i> at least about 150",
 "Fewer lines (about 375 before filtering) means less power: a real but small link can fail here. A failure is weaker evidence than a pass.")}
{code(10, "Screen 5: same tests on proteins")}
{tested("<b>Tested</b> on mock protein files in the Gygi layout, including a line measured twice (averaged) and missing values. Planted simulation: NRF2 targets 0.52, glutathione 0.57, both passing all checks.")}

<h2>11. Screen 6: gene by gene (Cell 11)</h2>
{glance("Which individual genes change with ploidy, and do RNA and protein agree?",
 "Each redox and control gene, and each gene of interest, goes through Screen 1 on its own. The false discovery rate is computed across all genes tested (about 170).",
 "q below 0.05; strongest when significant in both RNA and protein in the same direction",
 "Single genes are noisier than programmes and are not checked against random genes or controls. Use this to see what drives a programme result, or to find leads outside the programmes.")}
{code(11, "Screen 6: gene by gene")}
{tested("<b>Tested.</b> Planted simulation: 43 of 171 genes passed in RNA, all from the changed programmes plus GMNN (planted). No-link simulation: none passed.")}
<p><b>What you should see:</b> the 15 strongest genes in RNA with the programmes they belong to, then the genes that pass in both RNA and protein.</p>

<h2>12. Screen 7: Geminin and replication licensing (Cell 12)</h2>
{glance("Do GMNN and CDT1 differ with ploidy, and does Geminin go with redox?",
 "For each measure (GMNN RNA, CDT1 RNA, GMNN minus CDT1, and the GMNN and CDT1 knockout effects): (a) partial correlation with ploidy and with doubling; (b) with each redox programme; (c) with each redox programme after also removing ploidy. RNA measures are compared with about 200 genes that have similar expression and follow proliferation equally closely.",
 "q below 0.05 and p_vs_matched_genes below 0.05",
 "Proliferation. GMNN and CDT1 rise and fall with the cell cycle, so they correlate with anything else that follows proliferation. The matched-gene comparison is the guard. Knockout effects have no matched comparison, so read them with more caution.")}
{code(12, "Screen 7: Geminin and replication licensing")}
{tested("<b>Tested.</b> Planted simulation (GMNN lower in higher-ploidy lines): GMNN vs ploidy -0.41, matched genes -0.04 to 0.07; GMNN vs NRF2 targets -0.13, falling to +0.09 and no longer beating matched genes once ploidy was removed. No-link simulation: GMNN vs redox up to 0.11 (q = 0.007), but within the matched-gene range for every programme.")}
<h3>How to read Screen 7</h3>
{table(["Pattern", "What it suggests"], [
 ["<b>GMNN lower (or GMNN minus CDT1 lower) with ploidy</b>", "Weaker licensing brake in high-ploidy lines, consistent with re-replication. Worth comparing with the fly result"],
 ["<b>GMNN knockout more harmful (more negative) in high-ploidy lines</b>", "High-ploidy lines depend more on Geminin, perhaps to prevent further re-replication"],
 ["<b>GMNN goes with redox, but not after removing ploidy</b>", "The Geminin-redox link runs through ploidy. This was the pattern in the planted simulation"],
 ["<b>GMNN goes with redox even after removing ploidy</b>", "A link between Geminin and redox that ploidy does not explain; a separate lead"],
 ["<b>Significant, but within the matched-gene range</b>", "Explained by proliferation, not Geminin specifically"],
 ["<b>Nothing passes</b>", "No Geminin signal in cell lines; does not contradict the fly result (see the box in section 1)"]])}

<h2>13. Follow-up runs: TP53 and proliferation (Cell 1 settings)</h2>
<p>Run the main analysis first, with the default settings. Then rerun with these settings changed, one at a time, and
keep every run. Cell 15 adds the settings to each file name (for example <i>route_a_findings_TP53normal_prolif_alltissues.csv</i>), so runs do not overwrite each other.</p>
{table(["Setting in Cell 1", "What it does", "Why"], [
 ["<b>TP53_MODE = TP53-normal lines only</b>", "Keeps only lines with no damaging or hotspot TP53 mutation", "Cells that double their genome usually need p53 switched off, and p53 also switches on some redox genes (SESN1, TIGAR, GPX1). Mixing TP53-mutant and normal lines can hide a redox difference"],
 ["<b>TP53_MODE = Adjust for TP53 status</b>", "Keeps all lines and removes the TP53-mutant/normal difference, like tissue", "Same reason, with more lines; weaker if the link exists only in TP53-normal lines"],
 ["<b>ADJUST_PROLIFERATION unticked</b>", "Compares within tissue only", "Shows how much the proliferation adjustment changes the results. If ploidy affects redox through how fast cells divide, the adjustment removes that part"]])}
{tested("<b>Tested</b> on a simulation built to look like the first real run: redox rose with ploidy only in TP53-normal lines, and TP53-mutant lines had higher ploidy and lower redox. With all lines, NRF2 targets gave -0.07 and the DNA damage (p53) control -0.28, as in the real data. TP53-normal lines only gave +0.53; adjusting for TP53 status gave +0.29. All 15 cells ran in both modes.")}
{note("<b>The DNA damage control is mostly a p53 readout.</b> Most of its genes (CDKN1A, MDM2, BAX, TP53I3, SESN1, RRM2B, DDB2) are switched on by p53. A negative correlation with ploidy in all lines is expected, since genome-doubled cells usually lack p53 activity. If it moves towards zero in TP53-normal lines, the TP53 split is working as intended.")}
{note("<b>Several runs means more chances of a false positive.</b> Report every run, including the default one, and treat a result that appears only in one follow-up run as a lead to confirm in another dataset.")}

<h2>14. Findings, ranked by evidence (Cell 13)</h2>
{code(13, "Findings, ranked by strength of evidence")}
{tested("<b>Tested.</b> Planted simulation: NRF2 targets and glutathione labelled <i>Strong</i>; thioredoxin and ferroptosis <i>Moderate: RNA only</i> through shared genes. No-link simulation: every programme <i>None detected</i>.")}
{table(["Evidence label", "Meaning"], [
 ["<b>Strong: RNA and protein</b>", "Passes all three checks in RNA and in protein, same direction"],
 ["<b>Moderate: RNA only / protein only</b>", "Passes all three checks in one data type"],
 ["<b>Weak: correlates, but not specific</b>", "Passes Screen 1 only; random genes or stress controls move as much"],
 ["<b>None detected</b>", "No programme-level change with ploidy"]])}
<p>The table also shows the direction, partial rho in RNA and protein, whether whole-genome doubling agrees, and the closest
stress control. Check Screen 6 for programmes that share genes before counting them as separate findings.</p>

<h2>15. Charts and saving (Cells 14-15)</h2>
<p>Results tables and charts are saved into the same Drive folder, so they are kept after the session ends.</p>
<p>The first chart shows one dot per programme (blue RNA, orange protein) with a grey bar for the random-set range. The second
plots each gene's partial correlation against its false discovery rate; blue genes pass, and the 12 strongest are named.</p>
{code(14, "Charts")}
{code(15, "Save and download the results")}
{tested("<b>Tested:</b> both charts drawn and four tables saved on mock data. The download step itself runs only in Colab.")}

<h2>16. Results summary for your supervisors</h2>
{table(["Item", "Your result"], [[f"<b>{i}</b>", ""] for i in [
 "Date and settings recorded", "Cell lines analysed; tissues; NRF2-pathway mutants removed",
 "Lines with protein and CRISPR data", "Programmes labelled Strong or Moderate, with direction (Cell 13)",
 "Closest stress control for programmes that failed Screen 3", "Whole-genome doubling: agrees or not",
 "Genes passing in both RNA and protein (Screen 6)", "GMNN vs ploidy: partial rho, q, beats matched genes?",
 "GMNN minus CDT1 vs ploidy", "GMNN knockout effect vs ploidy",
 "GMNN vs redox: before and after removing ploidy", "TP53-normal lines only: lines kept; partial rho for each redox programme and the DNA damage control",
 "Adjust for TP53 status: same", "Without proliferation adjustment: same",
 "CNS/Brain-only run, if done",
 "Leads to take to brain data or flies"]], "summary")}

<h2>17. Troubleshooting</h2>
{table(["What happens", "Likely cause", "What to do"], [
 ["<b>Cell 1: Folder not found</b>", "Folder name or location differs", "Find the folder in Google Drive and copy its path into DATA_DIR; it starts /content/drive/MyDrive/"],
 ["<b>Cell 1: Drive asks to connect again</b>", "The session restarted", "Normal; click Allow. Files stay in Drive"],
 ["<b>Cell 2: MISSING ...</b>", "File not in the folder, or renamed", "Add or rename it in Google Drive, then rerun Cell 2"],
 ["<b>Cell 3: Expected a column like ...</b>", "DepMap changed a file layout", "Send me the message and the file name"],
 ["<b>Cell 3: protein layout not recognised</b>", "Different protein file", "Untick USE_PROTEINS in Cell 1; send me the first lines of the file"],
 ["<b>Session crashed / out of memory</b>", "Expression and CRISPR files are large", "Runtime &rarr; Restart session and run from Cell 1 (files stay in Drive), or untick USE_CRISPR"],
 ["<b>Cell 12: None of the genes of interest were found</b>", "Typo or alias in GENES_OF_INTEREST", "Use official symbols (GMNN, not Geminin)"],
 ["<b>All p-values exactly 0.0005 or 0.002</b>", "The smallest value the shuffles allow", "Normal; report as below 0.001 or below 0.005"]])}

<h2>18. What I tested and what I did not</h2>
{note("<b>Tested on simulated data:</b> all 15 cells, on mock files built in the DepMap and CCLE layouts (including the newer DepMap expression layout and a duplicated protein sample), in two scenarios: a planted link between ploidy, redox and GMNN (found by every screen), and no link with tissue, proliferation and KEAP1 confounding (nothing labelled; Geminin correlations caught by the matched-gene comparison). <b>Not tested:</b> the Google Drive connection and the real file layouts; DepMap's 25Q2 notes say ploidy is available only for lines released up to 24Q4, which is why the 24Q4 files are used. The TP53 split uses DepMap's damaging and hotspot mutation calls only; lines with TP53 lost by other routes count as normal. <b>Not verified:</b> the gene lists are a curated starter set; have your supervisor check them or replace them with MSigDB sets before results are shared.")}
"""
css = """
@page { size: A4; margin: 18mm 17mm 20mm 17mm; }
body { font-family: 'DejaVu Sans', Arial, sans-serif; font-size: 9.6pt; color: #1d1d1d; line-height: 1.42; }
h1 { color: #1f3b5c; font-size: 22pt; line-height: 1.15; margin: 0 0 6px 0; }
.sub { color: #666; font-size: 11pt; margin-top: 0; }
h2 { color: #1f3b5c; font-size: 15pt; margin: 22px 0 8px 0; page-break-after: avoid; }
h3 { color: #2a7a73; font-size: 11.5pt; margin: 16px 0 6px 0; page-break-after: avoid; }
table { border-collapse: collapse; width: 100%; margin: 6px 0 10px 0; font-size: 8.6pt; page-break-inside: auto; }
tr { page-break-inside: avoid; }
th { background: #1f3b5c; color: white; text-align: left; padding: 6px 7px; font-weight: bold; }
td { padding: 6px 7px; border-bottom: 1px solid #d5dbe3; vertical-align: top; }
tbody tr:nth-child(even) td { background: #eef1f6; }
table.glance thead { display: none; }
table.glance td:first-child { width: 22%; background: #eef1f6; }
table.glance { border-top: 2px solid #1f3b5c; }
table.summary td:first-child { width: 62%; }
pre { font-family: 'DejaVu Sans Mono', monospace; font-size: 6.6pt; line-height: 1.32; background: #f7f8fa;
      border-left: 3px solid #c9d1dc; padding: 7px 9px; white-space: pre-wrap; margin: 4px 0 8px 0; }
.cellname { font-weight: bold; margin: 12px 0 2px 0; page-break-after: avoid; }
.ok { background: #eef7ee; border: 1px solid #9fcb9f; padding: 8px 10px; margin: 6px 0 12px 0; font-size: 8.8pt; page-break-inside: avoid; }
.warn { background: #fdf6e3; border: 1px solid #e3c26a; padding: 8px 10px; margin: 8px 0 12px 0; font-size: 8.8pt; page-break-inside: avoid; }
.small { font-size: 8.4pt; color: #555; }
ol, ul { padding-left: 20px; } li { margin-bottom: 4px; }
"""
doc = f"<!doctype html><html><head><meta charset='utf-8'><title>{TITLE}</title><style>{css}</style></head><body>{body}</body></html>"
open(os.path.join(HERE, "guide.html"), "w").write(doc)

from playwright.async_api import async_playwright
async def main():
    async with async_playwright() as p:
        browser = await p.chromium.launch()
        page = await browser.new_page()
        await page.goto("file://" + os.path.join(HERE, "guide.html"))
        footer = (f"<div style='font-family:DejaVu Sans,Arial;font-size:7.5pt;color:#777;width:100%;padding:0 17mm;"
                  f"display:flex;justify-content:space-between'><span>{TITLE} - exploratory guide - draft to share with "
                  f"supervisors</span><span class='pageNumber'></span></div>")
        await page.pdf(path=os.path.join(HERE, "RouteA_ploidy_redox_guide.pdf"), format="A4", print_background=True,
                       display_header_footer=True, header_template="<div></div>", footer_template=footer,
                       margin={"top": "18mm", "bottom": "20mm", "left": "17mm", "right": "17mm"})
        await browser.close()
asyncio.run(main())
print("pdf written")
