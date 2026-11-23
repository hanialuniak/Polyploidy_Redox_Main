# Results

Last updated: 23 November 2026.

## Summary

No redox programme tracks polyploidy in any dataset tested. The association that recurs is
with DNA damage and p53. The one mechanistic lead comes from the fly ovary, where the ER
protein-folding machinery rises at the mitotic-to-endocycle switch — including `Ero1L`, which
generates H₂O₂ — while two enzymes that remove H₂O₂ or supply glutathione fall. That pattern
would raise H₂O₂ without ever appearing as an antioxidant response in a gene set.

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

Three adult runs, 30,049 cells called, 9,455 passing the fragment filter.

**Not usable.** Doublets alone can account for all the apparent polyploidy: expected doublet
share among the kept cells is 13–19% per run, against an observed polyploid rate of 6.3%. The
flies are young, the libraries shallow (median 880–1,540 autosomal fragments per cell), and
doublets were not removed. No cell type exceeds the doublet expectation. Apparent redox hits
(ferritins, GstD1, cnc) are consistent with neuron–glia doublets.

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

**Conclusion:** the salivary gland pattern reflects secretory workload, not ploidy.

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
- **Down:** **Jafrac1 −1.06** (peroxiredoxin; −1.26, −1.05, −0.86 across the three libraries)
  and **Gclc −0.88** (rate-limiting for glutathione synthesis).
- **Flat:** GstD1 (−0.04 and +0.16). The GstD-GFP reporter would not register this.
- Not replicated: Sod3 (+1.91 in one atlas, +0.30 and inconsistent in the other).

More H₂O₂ produced, less capacity to remove it — a combination that raises H₂O₂ without
looking like an antioxidant response.

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
