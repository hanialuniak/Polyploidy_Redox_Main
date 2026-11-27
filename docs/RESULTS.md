# Results

Last updated: 24 November 2026.

## Summary

No redox *programme* tracks polyploidy in any dataset tested, and four working positive controls
say the tests would have found one. The association that recurs instead is DNA damage and p53.

**The redox leads did not survive the germline control, and the honest headline is now
negative.** The fly ovary shows a large ER protein-folding rise at the mitotic-to-endocycle
switch, with `Jafrac1` and `Gclc` falling and `Jafrac2` flat. Comparing the germline, which makes
the same switch in the same egg chambers without becoming secretory, shows **every one of those
changes is specific to follicle cells** and therefore tracks secretory differentiation, not
endoreplication. See "The germline control" below for the numbers.

What changes in both lineages is the replication licensing programme, which is expected and serves
as a positive control (*p* = 0.002 as a set), and — once the germline side is calibrated against
shuffled labels — a block of chromatin, satellite-DNA and nuclear-lamina genes, with `mre11` the
one DNA-repair gene among them. `Rrp1`, the fly APEX1, has the largest effect size in the project
and is shared on the face of it, but at *q* = 0.13 it does not clear FDR as a single gene; the same
applies to `geminin` individually (*q* = 0.14) even though its programme passes.

So after eight datasets the position is: **no redox programme and no individual redox gene has
been shown to track polyploidy once differentiation is controlled for.**

What that does and does not license. It rules out a *transcriptional* coupling of the kind these
datasets can see, in the systems tested. It does not distinguish between two remaining
possibilities: that polyploidy and redox state are not causally coupled at all, or that they are
coupled post-translationally, through enzyme activity, oxidation state or localisation, which
every dataset here is blind to. Choosing between those needs a measurement of activity rather than
abundance; it cannot be settled from expression data.

Neither lead is a transcriptional antioxidant response, which is why eight gene-set screens
returned nothing. Taken with the speed of the HyPer7 signal, the nulls argue that the regulation
is post-translational rather than transcriptional — a constraint on the mechanism, not a failure.

## Route A — human cancer cell lines (DepMap)

About 1,200 lines, ploidy and whole-genome doubling from PureCN, expression and CRISPR from
DepMap 24Q4, protein from CCLE proteomics. Adjusted for tissue and proliferation; NRF2-pathway
mutants excluded.

- No redox programme associates with ploidy in RNA (largest effect −0.07), with whole-genome
  doubling, or in protein (266 lines).
- The p53 and DNA damage programme is the clearest signal, at −0.11.
- Gene level: p53 targets plus PRDX3, PRDX1, ABCC1, VDAC2, MSRA, PEBP1 around −0.1, ME1 at
  +0.1. None confirmed at the protein level.
- Geminin: GMNN and CDT1 transcript levels do not track ploidy. GMNN knockout versus NRF2
  activity is −0.107, more negative than 99.2% of genes, but not confirmed in KEAP1-mutant
  lines. CDT1 knockout versus ploidy is +0.14; the essential-gene artefact check was not run.

## Route B — fly brain

### Bulk RNA-seq across ageing (GSE153165)

Optic lobe, central brain and ventral nerve cord at 0, 2, 7 and 21 days, 72 samples. Tests the
age × region interaction, since optic-lobe polyploidy rises from day 7.

- **Positive control works:** p53 DNA damage genes rise more in the optic lobe (+0.22,
  q < 0.001; Xrp1, rad50, Ku80, CG3344, Irbp), reproducing the published result.
- No redox programme beats random gene sets. NRF2 (−0.11) and iron (−0.18) reach q = 0.045 but
  fail the random-set test.
- Replication licensing rises in the optic lobe: Orc1 +0.38 (q = 0.035), Mcm5 +0.28 (q = 0.012),
  dup −0.27 (q = 0.063). Geminin and Cdc6 were below the expression filter.
- Single redox genes falling in the optic lobe with age: SelR (t = −6.2), Tsf1 (t = −6.2),
  GstS1 (−0.44, q = 0.033).

### Single-cell ATAC (Janssens et al. 2022)

#### Adult brain — three FDM runs

30,049 cells called, 9,455 passing the fragment filter.

**Not usable.** Doublets alone can account for all the apparent polyploidy: expected doublet
share among the kept cells is 13–19% per run, against an observed polyploid rate of 6.3%. The
flies are young, the libraries shallow (median 880–1,540 autosomal fragments per cell), and
doublets were not removed. No cell type exceeds the doublet expectation. Apparent redox hits
(ferritins, GstD1, cnc) are consistent with neuron–glia doublets.

#### Larval brain — two DG2 runs, as a positive control

19,363 cells called, 4,699 passing the fragment filter and the depth-thinning step. Larval
brains carry a known polyploid population, the subperineurial (surface) glia, which
endoreplicate to cover the growing brain. The ploidy model is fitted **per run**, not per
cluster: fitting per cluster re-baselines a uniformly polyploid cluster to diploid and hides
exactly the population the control depends on.

**The ploidy caller works; the association screens do not.**

- **Positive control, PARTLY PASSED.** Surface glia 0.174 versus 0.077 in other non-cycling
  cells, a difference of **+0.097 (95% interval +0.056 to +0.141)**. Against all other
  non-dividing cells the difference is +0.088 (+0.051 to +0.131). Cells are scored individually
  for cycling markers and the most cycling 20% dropped, because coarse cluster labels leave some
  S/G2 cells in the comparison group; doing so widens the gap.
- **Second control passing independently:** the dividing clusters are the most polyploid
  population of all (c15 at 0.362). Cells in S/G2 do have more than 2C of DNA, so the method is
  reading DNA content, not noise.
- **Absolute rates are not usable.** Overall polyploid rate 0.110 against a doublet expectation
  of 0.159 among kept cells. Only comparisons *between* cell types in the same run, which share
  the doublet background, can be read.
- **Screens 1–4 are uninterpretable in this data.** Every programme shifts positively, controls
  included — DNA damage +0.053 and immediate-early +0.047 match or beat licensing (+0.054) and
  NRF2 (+0.039) — and `p_vs_controls` is 1.00 for all seven test programmes. This is the doublet
  signature: a neuron–glia doublet looks polyploid *and* carries two cell types' worth of open
  chromatin, so every gene set appears more accessible in the cells scored as polyploid. Gene
  level: nothing survives FDR (best q = 0.092, `kay`, `GstD1`, `Mcm7`, `Hsc70-4`, `Xbp1`,
  `E2f1`, `tefu`); geminin +0.031, q = 0.237.

Doublet-removed or deeper libraries would be needed for the screens. The limitation is in the
data, not the method.

## Route C — mammalian liver

### Mouse hepatocytes sorted by ploidy and size (Lanz et al. 2024)

Proteomics of G1 cells sorted into small 2n, large 2n, small 4n, large 4n, which separates
ploidy from cell size and from genome concentration. 3,719 proteins with ≥2 peptides.

- **Positive control works:** histones track DNA almost exactly (+1.07 per doubling of DNA per
  volume).
- Redox and detox enzymes follow **genome concentration** more strongly than average proteins:
  GSTs −0.64, NADPH supply −0.39, NRF2 targets −0.41, glutathione synthesis −0.38,
  thioredoxin/peroxiredoxin −0.27 per doubling (all q < 0.05).
- At matched DNA per volume, **GSTs are 15–20% lower in 4n cells** (−0.23, q = 0.003), holding
  in the independent large-cell comparison (−0.29, q = 0.004), surviving stricter filters,
  compartment-matched controls, leave-one-out, and removal of size scaling.
- SOD and catalase do not change in any comparison.
- An unbiased screen of 678 GO processes found nothing significant after correction; redox
  terms rank in the top 3% but are not significant.
- Limitation: one TMT experiment, no biological replication.

### Rat hepatocytes sorted 2c/4c/8c (GSE132409)

Five animals, each sorted into 2c, 4c and 8c, compared within animal after quantile
normalisation.

- **Positive control works:** cell-cycle genes rise strongly in 8c (+1.37 per doubling).
- GSTs −0.04 and NADPH supply −0.09 per doubling, both significant against random genes, but
  **both lose significance after adjusting for zonation**. 2c hepatocytes are mostly
  pericentral, where several of these genes are naturally higher (zonation versus ploidy,
  r = −0.79).
- **No agreement with the mouse protein result** gene by gene (r = −0.06 across 53 shared
  enzymes). Gsta4, Gstt1 and Maob barely move in rat mRNA.

**Conclusion across liver:** the GST result does not replicate between mouse protein and rat
mRNA. If real, it is likely post-transcriptional.

## Route D — fly endocycling cells

### Salivary glands versus brain and imaginal discs (GSE121956)

- **Positive control works:** mitotic genes −3.4, E2F1 S-phase genes −0.8 in salivary glands.
- Higher in salivary glands, beyond random genes and all controls: NADPH supply +2.4
  (q = 0.011), glutathione synthesis +2.1 (q = 0.036), cnc/NRF2 targets +1.9 (q = 0.019),
  ER folding load +1.6 (q = 0.019), mitochondrial respiration +1.5 (q = 0.019), GSTs +1.2
  (q = 0.019). GstD1 +3.5, catalase +3.2.
- The NADPH oxidases Duox (−3.0) and Nox (−1.4) are *lower*; Ero1L is higher (+1.6).
- **But this compares different tissues**, so secretory identity is confounded with ploidy.

### S2 cells forced into endocycles (GSE121955)

Cyclin A or Myb knockdown versus GFP control, three replicates each.

- **Positive control works:** mitotic genes fall in both knockdowns.
- **No antioxidant, NRF2 or ER-load programme rises.** Mitochondrial respiration (−0.3 to
  −0.45) and ribosome genes (−0.19) fall in both.
- GSTs, NADPH supply and SOD/catalase fall with Myb knockdown only.

### Compartment-resolved screen across both datasets

The pooled sets average enzymes that work in different places, which hides where H₂O₂ is
handled. Splitting every redox gene by compartment, and adding the nucleus, changes the reading
of Route D. Percentiles are among all expressed genes, so a library-wide drift cannot
manufacture a result.

**Confirms the workload conclusion, compartment by compartment.** Everything that rises in the
gland is flat or falling when endocycles are induced without a change of identity:

| | Salivary gland | CycA KD | Myb KD |
| --- | --- | --- | --- |
| ER folding load | +1.64 (87th pct, p = 0.013) | +0.04 (54th) | −0.11 (34th) |
| `Ero1L` | +1.62 (87th) | −0.24 (28th) | −0.33 (12th) |
| `Jafrac2` (ER peroxiredoxin) | +0.86 (73rd) | −0.46 (14th) | −0.07 (38th) |
| Catalase | +3.16 (95th) | −0.04 (46th) | −0.31 (13th) |
| Glutathione supply | +2.22 (92nd) | −0.01 (49th) | −0.06 (40th) |

Note that `Jafrac2` *does* track folding load in the gland. The failure to track it in the ovary
is therefore specific to the mitotic-to-endocycle switch, not a general property of fly
secretory cells.

**New result, and the first redox finding that survives the identity control: cytosolic
peroxiredoxins fall on entering an endocycle.**

| Cytosolic H₂O₂ removal | log2 | percentile | p vs random |
| --- | --- | --- | --- |
| Myb knockdown | **−0.45** | **6.8th** | **0.009** |
| CycA knockdown | −0.38 | 19th | 0.102 |

Consistent gene by gene in both knockdowns: `Jafrac1` −0.27 / −0.27, `Prx2540-1` −0.51 / −0.48,
`Prx2540-2` −0.48 / −0.65. The direction matches Route E, where `Jafrac1` falls at the switch
(−1.06 and −0.43). **Three systems, same direction**, and S2 cells are the clean one because
identity is held constant.

`Duox` moves the other way: **+0.52 (84th pct) with CycA and +1.11 (98th pct) with Myb**,
against −2.99 in the gland. So in induced endocycles a source rises while cytosolic sinks fall —
the imbalance pattern, for the first time in a design where ploidy rather than tissue identity
is the variable. Caveats: `Duox` is a single gene, `Nox` was not detected, two of the four
cytosolic removal genes are the less-characterised `Prx2540` pair, and Myb knockdown has the
weaker positive control (mitotic genes −0.60 against CycA's −1.44).

**The thioredoxin system does not compensate.** In the gland it is the only removal class that
does not rise (65th percentile, against 86th–95th for catalase, peroxiredoxins, glutathione and
mitochondrial `Prx3`). In S2 iECs it is flat to down (−0.01, −0.21). Since HyPer7 is reduced by
the thioredoxin system, falling peroxiredoxin capacity and unchanged Trx capacity push a HyPer7
ratio the same way, and brightness alone cannot separate them — bolus-and-recovery kinetics can.

**No support for a nuclear H₂O₂ source.** `Su(var)3-3` (fly LSD1, FAD-dependent demethylase,
one H₂O₂ per demethylation) falls in the gland (29th pct) and in both knockdowns (16th, 33rd).
`Rrp1` (fly APEX1) is −2.49 in the gland but −0.85 and −0.01 in the knockdowns, so the gland
value is proliferation, not redox. `dhd` and `TrxT` were not detected, as expected for germline
and testis thioredoxins, so the nuclear thioredoxin arm remains untested.

**Conclusion:** the salivary gland pattern reflects secretory workload, not ploidy. What does
track endocycling, in the design that controls for identity, is a **fall in cytosolic
peroxiredoxin capacity** alongside a rise in `Duox`.

## Route E — fly ovary follicle cells

The cleanest comparison available: follicle cells divide to stage 6, then endocycle. Same cell
type, same tissue, same egg chamber. Two independent atlases, with groups matched on
dissociation-stress score before comparison.

| Programme | Jevitt GSE146040 (1 library) | Rust GSE136162 (3 libraries) | Replicated |
| --- | --- | --- | --- |
| ER folding load | +1.35, q = 0.005 | +0.72, q = 0.005, same direction in all three | Yes |
| Licensing and re-replication | −0.65, q = 0.014 | −0.47, q = 0.014 | Yes (expected) |
| Glutathione S-transferases | +0.40, fails control test | +0.18, not significant | No |
| cnc/NRF2 targets, SOD and catalase, NADPH supply, thioredoxin/peroxiredoxin, mitochondrial respiration | not significant | not significant | Consistent null |

Gene level, consistent in both atlases:

- **Up:** Sec61beta +1.59, Sec61alpha +1.28, Hsc70-3 (BiP) +1.16, Pdi +0.81, ERp60 +0.76,
  Xbp1 +0.69, crc +0.63, Calr +0.56, **Ero1L +0.50**.
- **Down:** **Jafrac1** (cytosolic peroxiredoxin) and **Gclc** (rate-limiting for glutathione
  synthesis). Direction replicates, magnitude does not: Jafrac1 is −1.06 in Rust (2nd percentile
  of all genes; −1.26, −1.05, −0.86 across the three libraries) but only −0.43 in Jevitt (18th
  percentile). Gclc −0.88 (4th percentile) and −0.59 (11th percentile).
- **Flat:** GstD1 (−0.04 and +0.16). The GstD-GFP reporter would not register this.
- Not replicated: Sod3 (+1.91 in one atlas, +0.30 and inconsistent in the other).

The first reading of this was "more H₂O₂ produced, less capacity to remove it". **The absolute
abundance check does not support the first half.** In linear CPM, `Ero1L` sits at 2.6–3.15 log2
CPM against `Jafrac1` at 8.6–9.2 — Ero1L is 60–90× less abundant, so a +0.50 rise in it is a
small absolute change against a large quenching capacity, and several other quenchers rise.
What the data support is a **fall in specific reducing capacity** (`Jafrac1`, `Gclc`) alongside
a large rise in ER folding load, not a demonstrated net rise in H₂O₂.

### The thioredoxin arm, enzyme by enzyme

Pulled from the saved gene tables of both atlases (`RouteE_thioredoxin_arm.ipynb`) to ask whether
reducing capacity falls *as a group* or whether `Jafrac1` moves alone.

| Enzyme | Compartment | Rust GSE136162 | Jevitt GSE146040 |
| --- | --- | --- | --- |
| `Jafrac1` | cytosol | −1.06 (2nd pct) | −0.43 (18th pct) |
| `Jafrac2` | **ER** | **+0.08 (56th pct)** | **+0.04 (54th pct)** |
| `Trxr-1` | thioredoxin reductase | −0.25 (24th pct) | −0.18 (35th pct) |
| `Trx-2` | thioredoxin | −0.27 (23rd pct) | −0.42 (19th pct) |
| `Prx3` | mitochondria | −0.37 (16th pct) | −0.14 (38th pct) |
| `Prx2540-1`, `Prx2540-2` | peroxiredoxin | not detected | **+0.95, +0.83 (90th+ pct)** |
| `Sod1` | cytosol | +0.38 (80th pct) | +0.54 (82nd pct) |
| `Cat` | peroxisome | −0.28 (22nd pct) | +0.58 (83rd pct) |

- **The recycling machinery falls in both, consistently but weakly.** `Trxr-1` and `Trx-2` are
  down in both atlases and in all three Rust libraries, but at 19th–35th percentile they sit
  within the general downward drift rather than standing out.
- **As a programme the arm does not pass:** −0.27 against a random-gene median of +0.11
  (q = 0.26) in Rust, −0.16 against +0.06 (q = 0.36) in Jevitt. In Jevitt two peroxiredoxins
  rise strongly, so the group mean is +0.05.
- **`Jafrac2` is flat in both.** The ER peroxiredoxin does not move while ER folding load rises
  to the 98th percentile (`Hsc70-3` +1.16 and +2.21) and `Ero1L` rises to the 82nd–87th. The ER
  increases its folding load and its oxidase without increasing its peroxidase. This is the
  clearest compartment-specific imbalance in the dataset.
- **`Sod1` rises in both** (80th–82nd percentile), which adds a second H₂O₂ source, since SOD
  converts superoxide to H₂O₂.
- **Power is limited.** Dropout leaves 5 of 12 arm enzymes detected in Rust and 8 of 12 in
  Jevitt; `Prx5`, `Trxr-2`, `TrxT` and `Srx` are missing from both, and `Gr1` from either.

**Consequence for the HyPer7 confound.** HyPer7 is reduced by the thioredoxin system, so a fall
in that system raises the probe's basal ratio without any change in H₂O₂. The fall here is real
in direction but modest in size, which **weakens that confound without eliminating it** — it is
unlikely to account for a large HyPer7 rise on its own. The bolus-and-recovery control is still
needed to separate the two.

### Compartment-resolved screen, both atlases

| | Rust GSE136162 | Jevitt GSE146040 |
| --- | --- | --- |
| ER folding load | +0.76, **94th pct**, p = 0.001 | +1.48, **~97th pct** |
| `Ero1L` (makes one H₂O₂ per disulfide) | +0.50, 87th | +0.53, 82nd |
| **`Jafrac2`** (ER peroxiredoxin) | **+0.08, 56th** | **+0.04, 54th** |
| `Trx-2` | −0.27, 23rd | −0.42, 19th |
| `Trxr-1` | −0.25, 24th | −0.18, 35th |
| **`Rrp1`** (fly APEX1) | **−2.04, 0.4th pct** | **−1.73, 0.6th pct** |
| `Prx3` (mitochondrial) | −0.37, 16th | −0.14, 38th |
| `Sod2` | −0.22, 26th | −0.44, 17th |
| `Cat` | −0.28, 22nd | +0.58, 83rd |

- **The ER imbalance replicates.** `Jafrac2` sits at the 56th and 54th percentile — the middle of
  the genome — while the folding machinery it serves is in the top 3–6% (`Sec61beta` +1.59 and
  +2.36, `Hsc70-3` +1.16 and +2.21). Since the salivary gland *does* scale `Jafrac2` up with its
  load, this failure to scale is specific to the switch, not to secretory cells.
- **`Rrp1` has the largest consistent effect size anywhere in the project** (but see the
  calibration note in the germline section: as a single gene it is *q* = 0.13 on the germline side,
  so the effect size is large and the evidence that it is shared is not conclusive): bottom 0.5% of all
  genes in both atlases, and consistent across all three Rust libraries (−2.34, −1.82, −1.97).
  It also falls hard in the salivary gland (−2.49). Caveat: this is a replication-coupled repair
  enzyme and every one of these comparisons involves cells leaving the mitotic cycle, so
  proliferation is the obvious confound. Against that, S2 cells also leave the mitotic cycle and
  `Rrp1` barely moves there (−0.85 with CycA, −0.01 with Myb). Worth checking whether the
  Ref-1 redox-activating function of human APEX1 (which depends on Cys65) is conserved in
  *Drosophila* Rrp1 before reading anything redox into it — in many invertebrates it is not.
- **The thioredoxin arm falls consistently but weakly** in both, as the enzyme-by-enzyme table
  above records.
- `Su(var)3-3` was not detected in either atlas, so the nuclear-source question is untestable in
  ovary. It was testable in Route D, where it failed.

### The germline control (GSE136162)

Follicle cells endocycle and differentiate into an eggshell-secreting epithelium over the same
hours, so Route E alone cannot separate those. The germline can: cystocytes divide in the
germarium and become endocycling nurse cells — the same switch, the same egg chambers, the same
libraries, with a completely different differentiation programme and no eggshell.

GSE136162 resolves a single germline cluster of 499 cells, so the split was made by ranking those
cells on M-phase score and comparing the top and bottom thirds (167 each). That is weaker than an
independent cluster split, but both controls passed and the comparison demonstrably has power:

- **M-phase genes fall −1.08** between the groups, past the −0.8 threshold.
- **`bam` is enriched 56.3% against 24.6%**, a 2.3-fold difference, so the groups really do differ
  in germarial character.
- **Licensing falls in both lineages** (−0.49 follicle, −0.59 germline, 12th percentile,
  p = 0.002) — the expected positive control.
- The germline contrast detects large changes where they exist: GSTs +0.95 (95th percentile,
  p < 0.001) and ribosomal genes +0.81 (93rd percentile, p < 0.001). A null here is a real null.

| Gene | Follicle | Germline | Reading |
| --- | --- | --- | --- |
| `Sec61beta` | +1.59 (99th pct) | +0.11 (60th) | follicle only |
| `Hsc70-3` (BiP) | +1.16 (98th) | +0.03 (54th) | follicle only |
| `Ero1L` | +0.50 (87th) | −0.21 (36th) | follicle only |
| `Jafrac1` | −1.06 (2nd) | −0.06 (46th) | follicle only |
| `Gclc` | −0.88 (4th) | +0.01 (52nd) | follicle only |
| `Trx-2` / `Trxr-1` | −0.27 / −0.25 | −0.08 / +0.13 | follicle only |
| `Cat` | −0.28 (22nd) | −0.13 (41st) | follicle only |
| **`Rrp1`** | **−2.04 (0.4th)** | **−0.67 (9th)** | **both lineages, q = 0.13** |
| **`geminin`** | **−1.11 (1.7th)** | **−0.53 (15th)** | **both lineages, q = 0.14** |

**Downgraded on 27 November 2026, after calibration.** The germline side of this table was later
given a permutation null (Cell 13: germline labels shuffled within library, 200 times, per-gene
*p* against a null pooled within expression deciles). Two things changed.

First, the shared list as a whole is only partly real: 726 genes were called shared with the
germline, against a median of 226 from shuffled labels, 95th percentile 339. The count is above
chance (*p* = 0.005) but roughly a third of it is noise, and 242 genes survive at *q* < 0.10.

Second, and this corrects the claim made above and in the summary PDF and slides, **neither `Rrp1`
nor `geminin` clears FDR 0.10 as an individual gene**: `Rrp1` is *q* = 0.130 and `geminin` is
*q* = 0.140. Both are suggestive and neither is significant, on 197 germline cells after stress
matching. Note the distinction for `geminin`: the *licensing programme* as a set still falls in
both lineages at *p* = 0.002, which is a better-powered test than any single gene, so the
programme-level statement stands and the single-gene statement does not. For `Rrp1` there is no
programme-level result to fall back on — the p53/DNA-damage set is *p* = 0.37.

What does clear calibration is a block of chromatin, satellite-DNA and nuclear-lamina genes falling
in both lineages: `HmgD`, `D1`, `Lam`, `ball`, `Det` (*q* = 0.076), `Df31`, `HipHop`, `HP4`,
`mod(mdg4)`, `CtBP`, `Aos1`, `Bub3` (*q* = 0.08–0.09), with replication and repair machinery
alongside (`RPA2`, `RfC4`, `dUTPase`, `Cdk1`, and `mre11` at *q* = 0.092, the one repair gene that
does clear it). A few metal- and NADPH-handling genes rise in both (`MtnA` +0.94/+2.17, `Men`
+1.02/+1.89, `GstE1` +0.80/+0.90). **These themes were read out of a ranked list after the fact,
so Cell 14 tests them as named sets against expression-matched random sets before any of them is
treated as a finding.** Heterochromatin under-replication is well described in *Drosophila*
polyploid tissue, which makes the chromatin block plausible — and plausibility is exactly what
made the earlier redox leads survive longer than they deserved to.

**Consequence.** The ER folding rise, the `Jafrac1` fall and the `Gclc` fall are all absent from
the germline. They are therefore properties of follicle cells becoming a secretory epithelium, not
of entering an endocycle. This is the same conclusion the salivary gland and S2 pair reached, now
reached inside the one tissue that was supposed to be free of that confound.

It also withdraws the cross-system `Jafrac1` claim above: it falls in follicle cells and in S2
cells, but not in the germline, which is the comparison that holds identity constant in vivo.

**Caveats.** The germline captured in ovary scRNA-seq is biased toward small cells, because mature
nurse cells reach 512C and are too large for droplet capture, so the "nurse cell" group may be
earlier than its name suggests. The split is by score rather than by cluster. 167 cells per group.

### Is `Jafrac2` flat, or unmeasurable?

The ER imbalance claim needed `Jafrac2` to have been able to rise. It could: expression at the
63rd percentile of all genes, detected in 52.4% of mitotic and 39.7% of endocycling cells —
neither at ceiling nor too sparse to judge. So the flat value is real. But the per-library values
disagree in sign (+0.64, −0.13, −0.27), so it is better described as "no consistent change" than
as a measured zero.

`Ero1L` is the weak link on the other side: expression in the **bottom 5%** of all genes, detected
in only 11–13% of cells, which barely clears the sparsity threshold. The rise in ER folding load
is the robust part of that story — and it does not need `Ero1L` transcript to rise, since more
substrate through existing Ero1 generates more H₂O₂ regardless.

**Cross-reference to Route D:** the `Jafrac1` fall here is the same direction as the significant
fall in cytosolic peroxiredoxins when S2 cells are pushed into endocycles, which is the one
comparison that holds cell identity constant. `Jafrac2` behaves differently in the two systems:
flat here, rising with folding load in the salivary gland.

**Caveat:** in this lineage ploidy and developmental stage rise together, so this is
endocycling versus mitotic, not ploidy alone.

## Why the public data show so little

- **Wrong readout.** These datasets measure transcripts, proteins or chromatin. H₂O₂ level is
  set by enzyme activity, oxidation state and metabolic flux, none of which appear here.
- **Dilution.** Polyploid cells are at most about 20% of the fly optic lobe.
- **Different polyploidy.** Cancer lines and cultured hepatocytes double by failed division,
  not endocycling.
- **Timing.** A redox change that triggers endoreplication may be gone once cells are polyploid.
- **Confounders.** Cell size, zonation, doublets and dissociation stress each mimic or mask a
  ploidy effect.

## Open questions

1. Does H₂O₂ rise before or after entry into endoreplication? Nothing here can date it; every
   dataset is cross-sectional.
2. **How fast does the HyPer7 signal actually rise?** This has not been measured. A rate faster
   than a round of endoreplication would argue against a replication-dependent cause; a slower one
   would be consistent with it. Until it is measured, neither inference is available, and earlier
   drafts of this document and of the PDF summaries wrongly treated a fast rise as established.
3. Is Ero1L the source of the H₂O₂ in polyploid cells?
4. Does raised H₂O₂ act through licensing factors, and is Geminin protein redox-sensitive?
5. What explains `Rrp1`? It is the most extreme consistent change in the project and the only
   non-licensing gene that survives the germline control.

## Supporting literature

- Zhou et al. 2026, *J Cell Biol* 225(5):e202507096 — polyploid wing disc cells carry more ROS
  than diploid neighbours; catalase, Sod1 or NAC lowers it; proposed route is protein synthesis
  → ER stress → ROS.
- Tien et al. 2008, *J Cell Biol* 182:1113 — Ero1L is required for Notch signalling, including
  follicle cell differentiation; Ero1L clones give giant compound egg chambers.
- Deshpande et al. 2024, *J Cell Biol* 223(12):e202401084 — raising ROS restores Hindsight, a
  Notch-driven differentiation marker, in follicle cells.
- Li, Young & Sun 2018, *PNAS* — Nox-generated superoxide, converted by SOD3 to H₂O₂, is
  required for *Drosophila* ovulation.
- Nandakumar et al. 2020, *eLife* — fly brains are 98–99% diploid at eclosion; polyploidy
  appears from day 7 and the optic lobe exceeds 20% polyploid cells at three weeks.
