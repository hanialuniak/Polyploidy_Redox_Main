# Polyploidy and redox state: an analysis of public datasets

Code for testing whether polyploid cells differ in redox state, across human cancer cell
lines, mouse and rat liver, and four *Drosophila* systems.

Every analysis runs in Google Colab from a single notebook. Each notebook downloads its own
data, checks its own assumptions, and stops with a plain-language message when something is
wrong.

**Short answer from the analyses:** no redox *programme* tracks polyploidy in any dataset
tested. The signal that recurs is DNA damage and p53. In the fly ovary, where mitotic and
endocycling cells of the same type can be compared directly, the ER protein-folding machinery
rises sharply at the switch — including `Ero1L`, which produces H₂O₂ — while `Jafrac1`
(peroxiredoxin) and `Gclc` (glutathione synthesis) fall. See [docs/RESULTS.md](docs/RESULTS.md).

## What is here

| Route | System | Question | Notebook |
| --- | --- | --- | --- |
| A | Human cancer cell lines (DepMap, CCLE proteomics) | Do redox programmes track ploidy or whole-genome doubling? | `RouteA_ploidy_redox.ipynb` |
| B | Fly brain, bulk RNA-seq (GSE153165) | Does the optic lobe shift as polyploidy rises with age? | `RouteB_fly_bulk_GSE153165.ipynb` |
| B | Fly brain, single-cell ATAC (Aerts lab atlas) | Do polyploid cells differ from diploid neighbours? | `RouteB_fly_ploidy_redox.ipynb` (+ `RouteB_fly_download_Aerts.ipynb`) |
| C | Mouse hepatocytes sorted by ploidy and size (Lanz 2024) | Does ploidy change redox enzymes at the protein level? | `RouteC_liver_redox_enzymes.ipynb` |
| C | Rat hepatocytes sorted 2c/4c/8c (GSE132409) | Does it reproduce at the mRNA level in vivo? | `RouteC_rat_GSE132409.ipynb` |
| D | Fly salivary glands and S2 cells (GSE121955/6) | Is the salivary gland pattern ploidy or secretory workload? | `RouteD_salivary_endocycle_redox.ipynb` |
| E | Fly ovary follicle cells (GSE146040, GSE136162) | What changes at the mitotic-to-endocycle switch, within one cell type? | `RouteE_ovary_endocycle_redox.ipynb` |

`polyredox_pipeline/` is a standalone Python package for inferring ploidy from single-cell
ATAC fragment overlaps and testing gene-set associations. It is independent of the notebooks.

## Running a notebook

1. Open [Google Colab](https://colab.research.google.com) and upload the notebook.
2. Run Cell 1. It mounts Google Drive and creates a folder under `Colab Notebooks/`.
3. Choose **Runtime → Run all**.

Most notebooks download their own data. Route A and Route C (mouse) need files placed in the
Drive folder by hand, because the hosts block automated downloads; each notebook's Cell 2 says
which files and where to get them.

No installation is needed: the notebooks use only packages Colab already has
(see `requirements.txt` for running them outside Colab).

## How the statistics work

The same structure is used throughout, because a gene set moving is not by itself evidence.

1. **Association** between the programme score and ploidy, within strata (donor, tissue,
   cell type, cluster or animal), by permutation or regression.
2. **Against random gene sets** matched on expression level or measurement depth, usually 500
   to 5,000 of them. This is the main defence: it catches anything that moves because the whole
   transcriptome shifted.
3. **Against stress controls** — heat shock, unfolded protein response, DNA damage,
   immediate-early genes, housekeeping. A redox programme has to beat these to count as specific.
4. **Benjamini–Hochberg FDR** across the programmes tested.
5. **A positive control that must move**, so a null result means something: mitotic genes in
   endocycling cells, cell-cycle genes in 8c hepatocytes, p53 targets in the ageing optic lobe,
   histones tracking DNA content.

Confounders handled explicitly, each of which produced a false positive at some point during
development:

- **Sequencing depth**: polyploid cells have more DNA, so raw depth is a collider. Handled by
  per-copy depth, binomial thinning, or rank-based scores.
- **Doublets**: two cells captured together look like one polyploid cell. Expected rates are
  computed per run, allowing for the fragment filter keeping doublets preferentially.
- **Cell size**: separated from ploidy using the size-sorted hepatocyte design.
- **Liver zonation**: 2c hepatocytes are mostly pericentral, where several redox genes are
  naturally higher.
- **Dissociation stress**: tissue dissociation switches on heat shock and immediate-early genes
  and hits large fragile cells hardest. Groups are matched on a stress score before comparison.
- **Composition**: yolk and chorion genes dominate late follicle cell libraries, so
  normalisation is median-of-ratios rather than simple scaling.

## Testing

Every notebook is built from numbered cell files in `src/` and tested against simulated data
before use:

```bash
cd src/route_e_ovary
python3 build.py                       # regenerate the notebook from cells/
python3 simulate.py /tmp/mock effect   # mock data with a planted effect
python3 run.py /tmp/out http://...     # run every cell against it
```

Each route is tested in two scenarios: one with a planted effect, which must be recovered at
roughly the right size, and one without, which must come back empty. Several routes are also
tested against deliberately broken input — shuffled sample order, mislabelled columns, missing
files, the wrong file format — to check that the notebook stops rather than reporting nonsense.

## Repository layout

```
notebooks/   ready-to-run Colab notebooks (generated; this is what you upload)
src/         numbered cell sources, notebook builders, simulators and test runners
docs/        step-by-step PDF guides and the results summary
results/     output tables worth keeping
polyredox_pipeline/   standalone ploidy-inference package with its own tests
```

Notebooks are generated from `src/`, so edits belong in the cell files, followed by
`python3 build.py`.

## Data sources

| Dataset | Accession or source |
| --- | --- |
| DepMap Public 24Q4 | [depmap.org](https://depmap.org/portal/download/) |
| CCLE proteomics (Nusinow et al. 2020) | [gygi.hms.harvard.edu](https://gygi.hms.harvard.edu/publications/ccle.html) |
| Fly brain ageing, bulk RNA-seq (Nandakumar et al. 2020, *eLife*) | GSE153165 |
| Fly brain single-cell ATAC (Janssens et al. 2022, *Nature*) | GSE163697; fragments from flybrain.aertslab.org |
| Mouse hepatocytes by ploidy and size (Lanz et al. 2024, *NSMB*) | PXD052786, Supplementary Table S1 |
| Rat hepatocytes 2c/4c/8c (Katsuda et al. 2020, *CMGH*) | GSE132409 |
| Fly endocycling cells (Rotelli et al. 2019, *PLoS Genetics*) | GSE121955, GSE121956 |
| Fly ovary single cells (Jevitt et al. 2020, *PLoS Biology*) | GSE146040 |
| Fly ovary single cells (Rust et al. 2020, *Nat Commun*) | GSE136162 |

## Limitations

- Transcript and protein abundance are not redox state. A change in enzyme activity,
  oxidation state or localisation would be invisible in all of these datasets.
- No dataset here measures ploidy and redox in the same cells. Ploidy is either inferred,
  sorted in bulk, or stands in for developmental stage.
- In the follicle cell lineage, ploidy and developmental stage rise together and cannot be
  separated.
- Several results rest on a single experiment with no biological replication, and are labelled
  as leads rather than findings.

## Licence

MIT. See [LICENSE](LICENSE).
