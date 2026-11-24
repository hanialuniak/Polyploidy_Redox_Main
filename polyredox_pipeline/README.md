# polyredox

Tests whether neurons with inferred polyploid genomes differ in redox gene programmes from diploid
neurons of the same type in the same donor. It uses 10x Multiome data: ploidy comes from the
ATAC fragments and redox state from the RNA of the same nucleus.

This pipeline finds associations only. It cannot show whether oxidative stress causes polyploidy,
whether polyploidy causes oxidative stress, or neither.

## Steps

1. **Fragment overlaps** (`fragments.py`). For each fragment, count the same-nucleus fragments
   covering its left end. One genome copy cannot produce overlapping fragments, so depth cannot
   exceed the copy number except at mapping artefacts. The idea comes from scPloidy
   (Takeuchi & Kato 2024). Fragment files are streamed one chromosome at a time.
2. **Ploidy calls** (`ploidy.py`). Each nucleus gets P(2C), P(4C) and P(8C) from a beta-binomial
   model, fitted by EM separately for each `cell_class`. Unlike scPloidy's unconstrained mixture,
   the polyploid components keep their binomial form, so rare polyploid neurons are not absorbed
   into a component that only separates diploid nuclei by depth.
3. **Doublet check** (`doublets.py`). Two same-class diploid nuclei in one droplet look exactly
   like one tetraploid. The expected same-class doublet rate is computed for each library from
   `n_recovered`, then compared with the observed polyploid rate.
4. **Redox scores** (`scoring.py`). RNA is first thinned to the same UMI count in every nucleus,
   because polyploid nuclei carry more RNA and rank scores depend on library size. Programmes are
   then scored by rank (UCell-style), with cell-cycle genes removed.
5. **Within-donor tests** (`stats.py`). Score is regressed on P(polyploid), with donor × region ×
   cell-type fixed effects and per-copy ATAC depth. P values come from donor-level sign flips.
6. **Specificity**. Each redox programme must beat expression-matched random gene sets and
   non-redox stress programmes (unfolded protein response, heat shock, DNA damage, immediate-early
   genes). It must also hold without depth adjustment and within the interquartile range of
   per-copy depth.
7. **Machine learning** (`ml.py`). Gradient boosting, split by donor, tests whether redox scores
   predict ploidy beyond technical features, compared with matched random programmes. It also
   reports leave-one-cohort-out performance and programme importance.

## Running

```bash
pip install -r requirements.txt            # add anndata for h5ad input
python -m polyredox run config.yaml
python -m pytest tests                     # unit tests plus a file-based end-to-end run
bash validation/run_scenarios.sh /tmp/sim  # simulated scenarios with known answers
```

See `config.example.yaml` for every option. For fly data, set `autosome_regex` and supply
Drosophila gene sets.

## Outputs (in `outdir`)

| File | Contents |
|---|---|
| `report.md` | Plain-language summary of everything below |
| `redox_verdicts.tsv` | One row per redox programme: effect, CIs, all checks, verdict |
| `association.tsv` | All programmes, all analysed nuclei and each class separately |
| `doublet_excess.tsv` | Observed vs doublet-expected polyploid rate, with donor-bootstrap CI |
| `ploidy_calls.tsv.gz` | Per-nucleus posteriors, calls, per-copy overlap rate, expected doublet rate |
| `ploidy_by_depth.tsv` | Polyploid rate by depth decile; a strong rise only in the bottom deciles suggests miscalibration |
| `ploidy_params.json` | Fitted parameters for each class |
| `interaction_pathology.tsv` | Whether the association changes with pathology |
| `ml_summary.json`, `ml_importance.tsv`, `ml_cross_cohort.tsv` | Machine-learning results |

### Verdicts

| Verdict | Meaning |
|---|---|
| `supported_association` | Passed all checks |
| `depth_sensitive` | Sign or significance changes with the depth adjustment |
| `not_specific_to_redox` | Non-redox stress programmes move about as much |
| `not_above_matched_gene_sets` | Random genes of similar expression move as much |
| `no_association` | FDR ≥ 0.05 |

Even a supported association means little if `doublet_excess.tsv` shows no excess for that class.

## Using SEA-AD

- **Fragments are not in the open AWS release.** The open ATAC AnnData there holds peak counts,
  which cannot be used to infer ploidy. Fragment files are in the controlled-access raw data on the
  AD Knowledge Portal (Synapse syn26223298), which needs a data use agreement.
- Use **Multiome** nuclei only, so ATAC and RNA come from the same nucleus. Barcodes in your cells
  table must match each library's fragment file. A forum thread reports that some barcodes do not
  match, so check the match rate for each library.
- Set `rna.layer` to the raw UMI layer. The pipeline stops if counts are not integers.
- I could not confirm SEA-AD's metadata column names, so map them with `columns` after checking
  `adata.obs.columns`.
- **Check how the atlas filtered doublets.** Overlap-based doublet callers such as AMULET flag
  nuclei with more than two overlapping fragments, which removes tetraploids by design. Filters on
  high counts also remove them. If such filters were applied, start from unfiltered barcodes.
- Validate on ROSMAP Multiome data if you can get access.

## Main limitations

- Ploidy calls were validated only on simulated data. Benchmark them on real nuclei of known DNA
  content (FACS-sorted 2C/4C) before trusting the rates.
- In simulation, calls below about 2,000 fragments overstated polyploidy, hence the default
  `min_fragments: 3000`. Calibration relies on mapping artefacts scaling with fragment count.
  Applying `exclude_bed` (simple repeats) makes this assumption safer.
- Every depth measure is partly caused by ploidy, so no single adjustment is unbiased. The depth
  sensitivity check exists for this reason. Add ratio QC metrics through `extra_covariates`.
- `kappa` (site heterogeneity) is fixed, not estimated.
- The packaged gene sets are a curated starter list. Replace them with MSigDB or GO sets for any
  real analysis. Sets sharing many genes are flagged in the report and are not independent results.
- Post-mortem data cannot show direction. Test that in flies or cultured neurons (for example,
  induced oxidative stress versus forced endocycling, measured with a redox sensor).
