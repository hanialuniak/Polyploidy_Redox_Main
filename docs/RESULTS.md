# Results

Last updated: 24 November 2026.

## Summary

No redox *programme* tracks polyploidy in any dataset tested, and four working positive controls
say the tests would have found one. The association that recurs instead is DNA damage and p53.

Two mechanistic leads survive. The first is from the fly ovary, where the ER protein-folding
machinery rises sharply at the mitotic-to-endocycle switch — including `Ero1L`, which generates
one H₂O₂ per disulfide — while `Jafrac1` (cytosolic peroxiredoxin) and `Gclc` (glutathione
synthesis) fall, and `Jafrac2` (the ER peroxiredoxin) stays flat.

The second is newer. **`Jafrac1`, the cytosolic peroxiredoxin, falls in all four
mitotic-to-endocycle comparisons available**: −0.27 and −0.27 in S2 cells pushed into endocycles
by CycA or Myb knockdown, −1.06 (2nd percentile) in the Rust ovary atlas and −0.43 (18th) in
Jevitt. S2 cells are the clean case, because identity is held constant. In the same S2
comparison the cytosolic removal class as a whole falls (−0.45, 6.8th percentile, p = 0.009 with
Myb) and `Duox` rises to the 98th percentile.

**That class-level claim does not generalise, and `Jafrac1` alone should be quoted.** The class
falls in three of four comparisons but rises in Jevitt (+0.68), because `Prx2540-1` (+0.95) and
`Prx2540-2` (+0.83) move opposite to the way they move in S2 (−0.51/−0.48 and −0.48/−0.65);
neither was detected in Rust.

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
- **`Rrp1` is the most extreme consistent change anywhere in the project:** bottom 0.5% of all
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

1. Does H₂O₂ rise before or after entry into endoreplication?
2. Is Ero1L the source of the H₂O₂ in polyploid cells?
3. Does raised H₂O₂ act through licensing factors, and is Geminin protein redox-sensitive?

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
