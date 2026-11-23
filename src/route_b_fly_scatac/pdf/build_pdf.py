import html, os, glob, asyncio
HERE = os.path.dirname(os.path.abspath(__file__))
CELLS = {int(os.path.basename(p)[1:3]): open(p).read().rstrip() for p in glob.glob(os.path.join(HERE, "..", "cells", "c*.py"))}
TITLE = "Route B (fly): polyploidy and redox in the fly brain"
BULK = {int(os.path.basename(p)[1:3]): open(p).read().rstrip() for p in glob.glob("/home/claude/route_b_bulk/cells/c*.py")}


def bcode(n, title):
    return f'<p class="cellname">Part 1, Cell {n} - {html.escape(title)}</p><pre>{html.escape(BULK[n])}</pre>'

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
<h1>{TITLE}: a step-by-step guide</h1>
<p class="sub">Two fly analyses: one you can run now with public RNA-seq from optic lobes and central brain at four ages,
and a single-cell ATAC analysis ready for when the fragment files are obtained. Exploratory; for Google Colab; written
for a beginner.</p>

<h2>Data access: what is available now</h2>
{table(["Data", "What it gives", "Status", "Used in"], [
 ["<b>GSE153165</b> (Nandakumar et al. 2020, eLife): bulk RNA-seq of optic lobe, central brain and ventral nerve cord at 0, 2, 7 and 21 days (72 samples)", "Gene expression where polyploidy rises steeply (optic lobe) and where it stays low (central brain)", "<b>Public, on GEO.</b> The notebook downloads it", "Part 1"],
 ["<b>Fly gene annotation</b> (Ensembl release 112, BDGP6.46)", "Gene names and positions", "<b>Already in your Drive</b> (route_b_fly folder)", "Parts 1 and 2"],
 ["<b>Single-cell ATAC fragment files</b> (Janssens et al. 2022, Aerts lab)", "Per-cell ploidy from overlapping fragments", "<b>Available again</b> (since 19 Nov 2026) at flybrain.aertslab.org/downloads/cellRanger/. The helper notebook <i>RouteB_fly_download_Aerts.ipynb</i> lists the runs and downloads them to Drive. GEO (GSE163697) does not have these files", "Part 2"],
 ["<b>Aging fly brain atlas</b> (Davie et al. 2018): single-cell RNA, 0-50 days", "Cell-type expression with age; no ploidy", "Public; not used yet", "Possible extension"]])}
{note("<b>Even with the fragment files, Part 2 has a limit.</b> The Aerts lab single-cell ATAC data come from larvae, pupae and young adults, and their extra adult runs left out the optic lobes. Polyploid cells appear only from day 7 and are mostly in the optic lobes, so those data may contain too few polyploid cells. Part 1 uses older flies and includes the optic lobes.")}

<h1 style="font-size:17pt;margin-top:22px">Part 1. Optic lobe versus central brain with age (run now)</h1>

<h2>P1.1 What it tests</h2>
<p>Nandakumar et al. measured polyploidy by flow cytometry and sequenced RNA from the same brain regions and ages. In the
<b>optic lobe</b> polyploidy rises from about 1% to over 20% by day 21; in the <b>central brain</b> it stays much lower.
If polyploid cells change their redox genes, the optic lobe should change with age <b>more than</b> the central brain.
The notebook measures exactly that difference (an age-by-region interaction) for each gene programme.</p>
{glance("Does a programme change with age differently in the optic lobe than in the central brain?",
 "For each programme score: score = region + age (log2 of days + 1, because the first time point is day 0) + region x age (+ sex). The region x age term is the difference between the optic lobe's and central brain's age slopes. Screens 2 and 3 compare it with random gene sets and stress controls.",
 "False discovery rate below 0.05, beats random gene sets (p below 0.05) and every stress control (p below 0.05)",
 "Anything else that changes with age in the optic lobe only. The paper reports DNA damage and cell cycle genes rising there with age, so a result is a lead, not proof that polyploid cells cause it.")}
{note("<b>Built-in positive control.</b> The paper reported the p53 head DNA damage response genes (Ku80, Irbp, Cht8, CG3344, CG4734, Fanci, lok, rad50, Xrp1) rising with age in the optic lobe specifically. The <i>EXPLORE: p53 head DNA damage response</i> set should therefore pass. If it does not, the samples may be mislabelled: check the table from Cell 3.")}

<h2>P1.2 Setup (Cell 1)</h2>
<p>Upload <i>RouteB_fly_bulk_GSE153165.ipynb</i> to Colab and run Cell 1 (click <b>Allow</b> for Google Drive). It creates a
folder <i>Colab Notebooks/route_b_fly_bulk</i> and uses the annotation you already downloaded.</p>
{bcode(1, "Setup, settings and connect Google Drive")}

<h2>P1.3 Download from GEO (Cell 2)</h2>
<p>Cell 2 downloads the sample information (series matrix) and the expression table from GEO into your Drive folder. If
GEO only has one file per sample, it downloads and unpacks those instead.</p>
{bcode(2, "Download the GEO files (a few minutes the first time)")}
{note("<b>What GEO has for this study:</b> one Excel table, <i>GSE153165_polyploidy_depth_normalized_counts.xlsx</i>, with gene names in an <i>id</i> column and one column per sample. A file of 0.0 MB in the download list is the small sample-information file and is fine.")}

<h2>P1.4 Build the table and label the samples (Cell 3)</h2>
<p>The Excel columns are named with the sequencing facility's codes (<i>Sample_72006</i>, <i>Sample_72007</i>, ...), while GEO
names the samples like <i>0dFemaleOL1</i> (day 0, female, optic lobe, replicate 1). Nothing in the files links the two, so
Cell 3 pairs them <b>in order</b>: the first Excel column with the first GEO sample, and so on. Region, age and sex are
read from the GEO names. The labels are saved as <b>sample_labels.csv</b> in your Drive folder.</p>
<p>Because the pairing is an assumption, Cell 3 then checks it against the data. Samples from the same region resemble
each other much more than samples from different regions, so for each sample it finds the most similar other sample and
asks whether that one carries the same region label. If the order is right, this is close to 100%; if the columns were
in a different order, it falls to about a third (chance). Below 80%, the cell stops.</p>
{bcode(3, "Build the expression table and read sample information")}
<p><b>What you should see:</b> a table with one row per sample and its region, age and sex, counts per region and age
(6 per region and age), then the check line, for example "nearest neighbour ... from the same region for 100%". The
"same group" figure can be low: replicates of neighbouring ages look alike, so this one is for information only.
If any sample is <i>unknown</i>, open sample_labels.csv in Google Drive (it opens in Google Sheets), fill in the missing values,
save it as CSV with the same name, and run Cell 3 again.</p>

<h2>P1.5 Normalise and choose genes (Cells 4-5)</h2>
<p>The table is already depth-normalised, so Cell 4 only log2-transforms it (raw counts, if a table had them, would first be scaled to counts per million). Genes expressed in
fewer than half the samples are left out. Cell 5 lists how many genes of each programme are expressed.</p>
{bcode(4, "Normalise and choose expressed genes")}
{bcode(5, "Gene programmes (fly gene symbols)")}

<h2>P1.6 Screens (Cell 6)</h2>
{bcode(6, "Screens: does the optic lobe change differently with age from the central brain?")}
{tested("<b>Tested</b> on mock GEO data (48 samples: 3 regions x 4 ages x 2 sexes x 2 replicates), in both file layouts. With a cnc (NRF2) change that tracked optic-lobe polyploidy: cnc passed all checks (t = 6.0), as did the p53 positive control (t = 6.1); nothing else passed. With no redox change: no redox programme passed (largest t = 1.1); the p53 positive control still passed.")}
<p><b>How to read it:</b> <i>interaction</i> is the extra change per doubling of age in the optic lobe compared with the
central brain, in standard-deviation units; <i>t</i> is its strength. Verdicts: <i>passes: rises (or falls) more with age in
optic lobe</i>, <i>not specific</i> (first screen only), or <i>no difference</i>.</p>

<h2>P1.7 Gene by gene, charts and saving (Cells 7-8)</h2>
{bcode(7, "Gene by gene, including Geminin, dup, Cdc6 and Orc1")}
{bcode(8, "Charts and saving")}
<p>The chart shows each programme's average score by age, optic lobe in orange and central brain in blue. A programme that
tracks polyploidy should separate after day 7.</p>
{table(["Pattern", "What it suggests"], [
 ["<b>A redox programme passes, rising more in the optic lobe</b>", "Redox genes go up where polyploid cells build up. Lead for FlyPer imaging"],
 ["<b>p53 positive control passes, redox does not</b>", "The data work; no redox shift tracks optic-lobe polyploidy at the tissue level"],
 ["<b>Licensing genes (Gem, dup, Cdc6, Orc1) rise in the optic lobe</b>", "Consistent with the re-entry the paper reported; compare with your Geminin results"],
 ["<b>p53 positive control does not pass</b>", "Check the Cell 3 region check and sample_labels.csv before reading anything else"]])}

<h1 style="font-size:17pt;margin-top:22px">Part 2. Single-cell ATAC (when fragment files are obtained)</h1>
<p>This part uses <i>RouteB_fly_ploidy_redox.ipynb</i>. It is ready, but waits for fragment files from the Aerts lab or
from your own aged optic-lobe experiments. Its sections are numbered B.1 onwards.</p>
<h2>B.1 What this notebook explores</h2>
<p>Single-cell ATAC data record where each cell's DNA is open, as short DNA <b>fragments</b>. One copy of the genome cannot
give two fragments that overlap, so cells whose fragments overlap a lot must hold more copies of their genome. This
notebook uses that to estimate the <b>ploidy of every cell</b>, then compares polyploid and diploid cells of the same type.</p>
{table(["Question", "How it is explored", "Cell"], [
 ["<b>Step A. Which cells are polyploid?</b>", "A probability for each cell of being 2C, 4C, 8C or 16C, from overlapping fragments", "6"],
 ["<b>Step B. Does it work?</b>", "Each cell type's polyploid rate must clearly exceed the rate doublets would give; surface glia are a bonus check", "7"],
 ["<b>Screen 1. Do redox genes open more or less in polyploid cells?</b>", "Correlate each programme's accessibility with P(polyploid) within each cluster; p-value by shuffling", "10"],
 ["<b>Screen 2. Stronger than random genes?</b>", "Same test with 500 random gene sets matched on accessibility level", "10"],
 ["<b>Screen 3. Specific to redox?</b>", "Must beat four non-redox stress programmes", "10"],
 ["<b>Explore. Licensing and DNA damage</b>", "The same screens for a licensing and cell-cycle re-entry set (Gem, dup, Cdc6, Orc1, Mcms, CycE, E2f1, Myc, CycD) and the p53 head DNA damage response", "10"],
 ["<b>Screen 4. Which genes?</b>", "Every programme gene, plus Gem, dup, Cdc6 and Orc1, tested on its own", "11"],
 ["<b>Per cell type</b>", "The same correlation in each cluster that contains polyploid cells", "12"]])}

<h3>What is known (Nandakumar et al. 2020, eLife)</h3>
<ul>
<li><b>Polyploid cells build up after emergence.</b> Newly emerged adult brains are 98-99% diploid. Polyploid cells appear from day 7 and rise until about day 21; most are 4C.</li>
<li><b>Mostly in the optic lobes.</b> By 3 weeks over 20% of optic-lobe cells can be polyploid, far more than in the central brain. Many neuron types, cortex glia and astrocyte-like glia become polyploid.</li>
<li><b>Geminin and Cdc6 control it.</b> Knocking down Geminin in neurons raised optic-lobe polyploidy; knocking down Cdc6 lowered it. Neurons re-enter DNA replication.</li>
<li><b>Oxidative stress drives it.</b> Paraquat and UV raised polyploidy by day 7, and polyploid cells resisted death, independently of p53.</li>
<li><b>Surface glia are highly polyploid but very few</b> (under 5% of polyploid cells), so they may not form their own cluster.</li>
</ul>
<p>So the data to use are <b>optic lobes from flies at least 1-3 weeks old</b>. In young or central-brain samples, polyploid
cells may be as rare as doublets.</p>
{note("<b>Accessibility is not expression.</b> These data measure how open the DNA around a gene is, not how much RNA it makes. Open chromatin usually goes with active genes, but the link is loose for single genes. RNA from the same cells is not available in this dataset. Treat a result as a lead for RNA or protein measurements.")}

<h3>What stays the same throughout</h3>
{table(["Item", "Setting"], [
 ["<b>Comparisons</b>", "Always within cluster (cells of the same type) and sequencing run, adjusted for fragments per genome copy"],
 ["<b>Depth</b>", "Every cell is thinned to the same number of cut sites in genes before scoring, because polyploid cells have more fragments"],
 ["<b>Redox programmes (5)</b>", "cnc (NRF2) pathway, thioredoxin/peroxiredoxin, antioxidant enzymes, iron, mitochondrial respiration"],
 ["<b>Explore sets (2)</b>", "Licensing and cell-cycle re-entry; p53 head DNA damage response (genes reported up in aged optic lobes)"],
 ["<b>Stress controls (4)</b>", "Unfolded protein response, heat shock, DNA damage, immediate-early genes"],
 ["<b>Genes of interest</b>", "Gem (Geminin), dup (double parked, the fly CDT1), Cdc6 and Orc1"],
 ["<b>Labelling thresholds</b>", "False discovery rate below 0.05; beats random gene sets and every stress control at p below 0.05"]])}

<h3>The data</h3>
{table(["File", "What it is", "Where from"], [
 ["<b>fragments.tsv.gz</b> (about 280 MB per run)", "Every DNA fragment, with its cell barcode", "Aerts lab fly brain atlas, Cell Ranger outputs (Janssens et al. 2022, Nature)"],
 ["<b>singlecell.csv</b>", "Which barcodes Cell Ranger called as cells", "Same folder as the fragments"],
 ["<b>Fly gene annotation (GTF)</b>", "Where each fly gene sits on the genome", "Ensembl, Drosophila melanogaster (BDGP6)"]])}

<h2>B.2 Before you start: choose and get the data</h2>
<p>The download site works again. A small helper notebook, <b>RouteB_fly_download_Aerts.ipynb</b>, does the downloading
for you: Cell 2 lists every run on the site with its file sizes, and Cell 3 downloads the runs you choose straight into
<i>Colab Notebooks/route_b_fly</i>, one folder per run, with the file names the main notebook expects. If a download is
cut off, run Cell 3 again and it resumes. Fragment files can be several GB each, so check your Drive's free space.</p>
<ol>
<li>Upload <i>RouteB_fly_download_Aerts.ipynb</i> to Colab and run Cells 1 and 2. Look at the run names. Choose adult runs that include the <b>optic lobes</b> (or whole brains), and the <b>oldest flies</b> available. Avoid larval and pupal runs: dividing cells there pass through 4C and look polyploid. If only a central-brain run exists (<i>FDM__9c338b__10x_scATAC_mix_DGRP_adult_central_brain</i>), it still works, but expect few polyploid cells.</li>
<li>Type the chosen run names into <b>RUNS</b> in the helper's Cell 1, run Cell 1 again, then Cell 3. The last line it prints is the RUNS value to use in the main notebook.</li>
<li>If <i>Drosophila_melanogaster.gtf.gz</i> is not already in <i>route_b_fly</i> (you added it for the bulk notebook): open <b>ftp.ensembl.org/pub/current_gtf/drosophila_melanogaster/</b> and download the file ending <i>.gtf.gz</i> (not <i>.chr.gtf.gz</i> or <i>.abinitio</i>). Put it in <i>route_b_fly</i> and rename it <b>Drosophila_melanogaster.gtf.gz</b>.</li>
<li>In Cell 1, set RUNS to your run folder names, separated by commas.</li>
</ol>
{note("<b>Not checked from here.</b> I cannot reach the Aerts lab site from my side, so the helper was tested on a mock site with the same kind of folder listing, and the main notebook on simulated files. The chromosome names, barcode format and gene symbols are handled flexibly and checked by the notebook as it runs.")}

<h2>B.3 Setup (Cell 1)</h2>
<p>Upload <i>RouteB_fly_ploidy_redox.ipynb</i> to Colab, then run Cell 1. Colab asks to connect to Google Drive: click
<b>Allow</b>. RUNS must match the run folder names, separated by commas.</p>
{code(1, "Setup, settings and connect Google Drive")}
{tested("<b>Tested</b> on simulated data in the same file formats; the Drive connection runs only in Colab.")}

<h2>B.4 Step 1: check the files, read cells and genes (Cells 2-3)</h2>
{code(2, "Check the data files in your Drive folder")}
{code(3, "Cells and genes")}
<p><b>What you should see:</b> each file marked <i>found</i>, then the number of cells Cell Ranger called and the number
of genes read from the annotation (roughly 17,000 to 24,000, depending on the Ensembl release).</p>

<h2>B.5 Step 2: read the fragments (Cell 4, 10-30 minutes per run)</h2>
<p>This is the slow step. For every cell it counts how often fragments overlap (for ploidy) and how many fragment ends
fall in or just upstream of each gene (gene accessibility). The file is read one chromosome at a time to save memory.</p>
{code(4, "Read fragments: overlap counts and gene accessibility (10-30 minutes per run)")}
{tested("<b>Tested</b> on a simulated run of 1,800 cells and 6.7 million fragments, with chromosome names written as chr2L to check they are matched to 2L in the annotation.")}

<h2>B.6 Step 3: clusters and their names (Cell 5)</h2>
<p>Cells are grouped into clusters by which genes are open, then each cluster is named from marker genes: <i>elav, nSyb,
brp</i> (neurons), <i>repo</i> (glia), <i>moody, Mdr65</i> (surface glia), <i>alrm, Eaat1</i> (astrocytes), <i>Eaat2</i>
(ensheathing glia), <i>wrapper</i> (cortex glia). Every cell is then thinned to the same number of cut sites in genes.</p>
{code(5, "Group cells into clusters and name them from marker genes")}
{note("<b>Why thinning matters.</b> Without it, in a simulation where only the cnc pathway was changed, every programme and every stress control looked more open in polyploid cells (up to 0.15 for controls), because polyploid cells simply have more fragments. After thinning, the controls fell to about zero.")}
<p><b>What you should see:</b> a table of clusters with marker scores and a name. Neurons should be the large majority, with
several glia clusters. Surface glia may not appear, because they are few.</p>

<h2>B.7 Step 4: ploidy of every cell (Cell 6)</h2>
{glance("How many genome copies does each cell hold?",
 "For each fragment, count how many fragments from the same cell cover its start. The pattern of those counts depends on the number of copies. A model fitted separately for each cluster gives each cell a probability of being 2C, 4C, 8C or 16C.",
 "Calls are clear (above 0.8) for most cells and surface glia come out polyploid (Cell 7)",
 "Doublets: two cells captured together look like one cell with twice the DNA. Very low-depth cells: calls become uncertain.")}
{code(6, "Ploidy of every cell, from overlapping fragments")}
{tested("<b>Tested.</b> In the simulation, surface glia simulated as 8C and 16C were called 8C or 16C; all 86 neurons simulated as 4C (including doublets) were called 4C, and no diploid neuron was called polyploid.")}

<h2>B.8 Step 5: main check, more polyploidy than doublets? (Cell 7)</h2>
<p>Two cells captured together look exactly like one cell with twice the DNA. About 0.8% of cells per 1,000 recovered
in a run are doublets. Two corrections matter. First, doublets have about twice the fragments, so the fragment filter
in Cell 5 keeps them more often than single cells: when most cells are below MIN_FRAGMENTS, the kept cells can be
several times richer in doublets than the run as a whole. Cell 7 works this out from each run's own fragment counts.
Second, a doublet of two <i>different</i> neurons still lands in a neuron cluster, so doublets are not limited to
same-type pairs. Each cell type's polyploid rate, with a 95% interval, is compared with the expected doublet share of
the kept cells. Cell types whose interval sits above it have more polyploidy than doublets can explain.</p>
<p><b>Clusters with mixed markers</b> (for example neuron markers low and glia markers high, in a cluster named
"neuron") are often made of neuron-glia doublets. A high polyploid rate there is a warning, not a finding.</p>
{code(7, "Main check: more polyploidy than doublets can explain?")}
{tested("<b>Tested.</b> Simulation with polyploid neurons, astrocytes, cortex glia and surface glia: 9 of 12 cell types showed a clear excess over the expected doublet rate; surface glia came out at 1.00 against 0.08 for other cells.")}
{table(["What you see", "What it means"], [
 ["<b>Several cell types with clear excess (optic-lobe neurons, cortex or astrocyte-like glia)</b>", "Ploidy calls are picking up real polyploidy; carry on"],
 ["<b>WARNING: no cell type has more polyploidy than doublets</b>", "Stop. Results from Cell 10 on cannot be trusted. Use runs with optic lobes and older flies"],
 ["<b>Surface glia present and high</b>", "Bonus confirmation"],
 ["<b>A single neuron cluster with very high polyploidy</b>", "Real, or a cluster formed by high-depth cells. Check it has a clear excess and compare with the paper's cell types"]])}

<h2>B.9 Gene programmes and tools (Cells 8-9)</h2>
<p>Gene names are FlyBase symbols. Besides the redox programmes and stress controls, two <b>explore</b> sets come from the
fly papers: licensing and cell-cycle re-entry (Gem, dup, Cdc6, Orc1, Mcm2-7 with dpa as Mcm4, CycE, E2f1, Myc, CycD) and the
p53 head DNA damage response (Ku80, Irbp, Cht8, CG3344, CG4734, Fanci, lok, rad50, Xrp1, p53), reported to rise in aged
optic lobes. Cell 8 prints how many genes of each programme were found in the annotation;
programmes with fewer than 4 are dropped. It also says if Gem or dup are missing, which would mean a different symbol in
your annotation.</p>
{code(8, "Gene programmes (fly gene symbols)")}
{code(9, "Scores and tools")}
<p>A programme's score is the average, over its genes, of each gene's accessibility z-scored within cluster. Tests remove
cluster, run and fragments per genome copy before correlating.</p>

<h2>B.10 Screens 1-3: programmes vs ploidy (Cell 10)</h2>
{glance("Are redox genes more or less open in polyploid cells than in diploid cells of the same type?",
 "Screen 1: partial correlation of each programme score with P(polyploid), p-value from shuffling P(polyploid) within cluster and run. Screen 2: 500 random gene sets matched on accessibility. Screen 3: redox minus each stress control.",
 "q below 0.05, p_vs_random below 0.05 and p_vs_controls below 0.05",
 "Cells, not animals, are the unit here, so p-values are optimistic. Results from one run need a second run, or your own experiments, to count.")}
{code(10, "Screens 1-3: programmes vs ploidy, random genes, stress controls")}
{tested("<b>Tested.</b> Simulation with the cnc pathway and the licensing set opened in polyploid cells: both passed all three screens (partial rho 0.20 and 0.28); nothing else did. Four simulations with no effect: no programme passed all three screens.")}
{note("<b>A small built-in lean.</b> In the no-effect simulations every programme, controls included, was slightly more open in polyploid cells (about +0.03, up to +0.10 by chance). Each genome copy can only give so many fragments at a strongly open gene, so extra copies shift reads towards the most open genes. Random gene sets lean the same way, which is why a result must beat them (Screen 2) and the controls (Screen 3).")}
<p>Cell 10 ends with one line per programme: <i>passes</i> (with direction), <i>not specific</i> (Screen 1 only) or <i>no link</i>.</p>

<h2>B.11 Screen 4: gene by gene, Geminin and dup (Cell 11)</h2>
{glance("Which single genes open or close with ploidy, and what about Gem and dup?",
 "Each programme gene and gene of interest tested on its own; false discovery rate across all genes tested.",
 "q below 0.05",
 "Single genes are noisy in ATAC data; not checked against random genes.")}
{code(11, "Screen 4: gene by gene, including Geminin and dup")}
{tested("<b>Tested.</b> Planted simulation (licensing set opened): Cdc6 0.19 and Gem 0.14 passed (q = 0.017); Orc1 0.08 and dup did not. No-effect simulation: none of the four passed.")}
<p><b>How to read Gem and dup:</b> more open in polyploid cells would fit active licensing or re-replication; less open would
fit a switched-off licensing system after the cells stopped endocycling. Either is a lead to test with RNA or staining.</p>

<h2>B.12 Per cell type, charts and saving (Cell 12)</h2>
{code(12, "Per cell type, charts and saving")}
<p>The table shows the correlation separately for each cluster with enough polyploid cells. A real effect should point the
same way in several clusters. The chart shows mean P(polyploid) per cell type, surface glia in orange. Four tables are
saved to your Drive folder.</p>

<h2>B.13 Next: test it directly in flies</h2>
<p>ATAC data can only suggest a redox difference. The FlyPer lines (van Leeuwen, Aldaz Casanova et al. 2026) express the
H<sub>2</sub>O<sub>2</sub> sensor HyPer7 in the nucleus, mitochondria, cytosol or plasma membrane under UAS control. Because the sensor
reads out a ratio, more sensor in a larger polyploid cell should not change the reading. A direct test: nuclear or
mitochondrial FlyPer with nSyb-GAL4 or repo-GAL4; image optic lobes at 14-21 days; score each nucleus for size or DAPI
intensity and HyPer7 ratio; compare with Geminin RNAi (more polyploidy) and paraquat (more oxidative stress).</p>

<h2>B.14 Results summary for your supervisors</h2>
{table(["Item", "Your result"], [[f"<b>{i}</b>", ""] for i in [
 "Runs used (tissue, fly age); cells called; cells analysed", "Cell types with clear excess polyploidy over expected doublets",
 "Surface glia, if present: polyploid rate", "Licensing and p53 explore sets: verdicts", "Screens 1-3: any redox programme passing, with direction",
 "Screen 4: top genes; Gem, dup, Cdc6, Orc1 partial rho and q", "Per cell type: does any effect agree across clusters?",
 "Leads to test by RNA, staining or genetics in flies"]], "summary")}

<h2>B.15 Troubleshooting</h2>
{table(["What happens", "Likely cause", "What to do"], [
 ["<b>Cell 2: MISSING ...</b>", "Folder or file name differs", "Folder names must match RUNS; the GTF must be renamed or GTF_FILE changed"],
 ["<b>Cell 4 crashes / out of memory</b>", "Large run", "Runtime &rarr; Change runtime type &rarr; High-RAM if available; or use one run at a time"],
 ["<b>Cell 4 reports 0 autosomal fragments</b>", "Barcodes in singlecell.csv do not match the fragments file", "Send me the first lines of both files"],
 ["<b>Cell 5: no surface glia</b>", "They are very few", "Normal; Cell 7's doublet check is the main test. N_CLUSTERS = 30 may separate them"],
 ["<b>Cell 8: few genes found</b>", "Different gene symbols in the annotation", "Send me the printed list; symbols change between FlyBase releases"],
 ["<b>Cell 7: no cell type above the doublet rate</b>", "Young flies or no optic lobes", "Use runs with optic lobes and older flies; send me the Cell 7 table"]])}

<h2>B.16 What I tested and what I did not</h2>
{note("<b>Tested on simulated data:</b> all 12 cells, on fragment, singlecell.csv and GTF files built in the Cell Ranger and Ensembl formats: runs with polyploid surface glia, 4C neurons and glia, doublets and mapping artefacts, once with the cnc pathway and licensing set opened in polyploid cells and four times with no effect. <b>Not tested:</b> the real Aerts lab files, real cluster structure, and the Drive connection. <b>Not verified:</b> the fly gene lists, written from FlyBase symbols as I know them; Cell 8 reports what was found. Cells are the unit of analysis here, so p-values are optimistic; confirm any lead in a second run or in your own experiments.")}

<h2>B.17 References</h2>
<ul class="small">
<li>Nandakumar S, Grushko O, Buttitta LA (2020). Polyploidy in the adult Drosophila brain. eLife 9:e54385.</li>
<li>van Leeuwen LAG, Aldaz Casanova S, et al. (2026). Mapping subcellular H<sub>2</sub>O<sub>2</sub> dynamics reveals tissue specific redox patterns in Drosophila. bioRxiv 2026.08.21.745691 (preprint, not peer reviewed).</li>
<li>Shu Z, Row S, Deng W-M (2018). Endoreplication: the good, the bad, and the ugly. Trends in Cell Biology.</li>
<li>Janssens J, et al. (2022). Decoding gene regulation in the fly brain. Nature 601:630-636. Data: flybrain.aertslab.org.</li>
</ul>
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
                  f"display:flex;justify-content:space-between'><span>{TITLE} - step-by-step guide - draft to share with "
                  f"supervisors</span><span class='pageNumber'></span></div>")
        await page.pdf(path=os.path.join(HERE, "RouteB_fly_guide.pdf"), format="A4", print_background=True,
                       display_header_footer=True, header_template="<div></div>", footer_template=footer,
                       margin={"top": "18mm", "bottom": "20mm", "left": "17mm", "right": "17mm"})
        await browser.close()
asyncio.run(main())
print("pdf written")
