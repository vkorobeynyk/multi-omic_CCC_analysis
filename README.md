# Multi-omic cell–cell communication in the ageing mouse hippocampus

**🔗 Reports and evidence graph online: https://vkorobeynyk.github.io/multi-omic_CCC_analysis/**

Which ligand–receptor (LR) interactions between hippocampal cell types change with age, and how consistent is that change across technologies?

This repository combines four independent data types:

| Modality | Data | What it contributes |
|---|---|---|
| **Visium** (spot-based spatial transcriptomics) | In-house, young / middle / old mouse hippocampus (16 sections) | Cell-type colocalisation (RCTD) and LR co-expression across age |
| **scRNA-seq** | In-house hippocampus reference (Jessberger lab) | Cell-type-resolved LR scores per age (LIANA+) |
| **STARmap PLUS** (single-cell spatial) | Zeng et al., *Nat Neurosci* 2023 (SCP1375), control mice, 8 and 13 months | Independent spatial test of the candidate pairs (liana+ bivariate Moran's I, stLearn) |
| **Proteomics** | Huang et al. 2024 (mouse, 3 vs 20 months); Mol et al. 2022 (human hippocampal subfields, 1–10 to 81–90 y) | Protein-level age change of the LR genes |

Candidate pairs are first selected from Visium + scRNA-seq, then checked in STARmap PLUS and proteomics. The results for each candidate pair are collected in an interactive **evidence graph** (see [Evidence graph](#evidence-graph)).

---

## Contact

Vladyslav Korobeynyk, HIFO / DMLS, University of Zurich (Jessberger lab / Mark D. Robinson lab)

## Generative AI statement
Generative AI was used throughout this entire benchmark to make code nicer to read and more efficient. The entire logic of the benchmark was created by myself and I assume responsability of the content within this repo.


## Pipeline

All analyses are Quarto documents (`.qmd`, R) or Jupyter notebooks (`.ipynb`, Python) that render to self-contained HTML reports. Each script writes to `output/<script name>/`.

| # | Script | Language | What it does |
|---|---|---|---|
| 01 | `01_processing_subset_to_HPC.qmd` | R | Loads the Visium slides, computes QC metrics and keeps the manually selected hippocampus spots (Loupe Browser exports). |
| 02 | `02_Integration_clustering.qmd` | R | SCTransform per sample, merge, Harmony integration, clustering; regions ML, CA1_2, CA3, DG. |
| 03a | `03a_RCTD_spotdeconvolution.qmd` | R | RCTD deconvolution per region (DG, ML) with the scRNA-seq reference; marker and colocalisation QC; LR co-expression vs colocalisation across age for 4 cell-type pairs (Astrocyte–Vascular, Vascular–Microglia, OPC–Cajal.Retzius, Neurogenic_lineage–GC). Exports data for 03b. |
| 03b | `03b_runLiana_scRNAseq.ipynb` | Python | LIANA+ (CellPhoneDB + geometric mean) on the scRNA-seq reference, per age, for the candidate pairs from 03a. |
| 03c | `03c_crossModality_plots.qmd` | R | Joins spatial (03a) and scRNA-seq (03b) results; keeps pairs whose spatial slope, spatial detection enrichment and scRNA-seq score agree in direction with age → `significant_lr_pairs.csv`. |
| 04 | `04_KEGG_PathwayConvergence.qmd` | R | KEGG pathway membership of the LR genes; pathways hit by several cell-type pairs (with a pathview diagram). |
| 05 | `05_STRING_NetworkAnalysis.qmd` | R | STRING network of the LR genes, k-means clusters and GO enrichment. |
| 06 | `06_Synthesis_FunctionBySpatialPair.qmd` | R | Combines 04 and 05 into one cell-type-pair × function view (redundant GO terms merged with rrvgo). |
| 07 | `07_STARmapPLUS_ControlSpatialObject.qmd` | R | Builds a Seurat object of the 4 STARmap PLUS control sections and a hippocampus-only subset. |
| 08 | `08_STARmapPLUS_CCC_GenePanelCheck.qmd` | R | Keeps the 03c pairs whose genes are on the STARmap PLUS panel; normalises and exports data for Python. |
| 09 | `09_STARmapPLUS_BivariateMoransI.ipynb` | Python | liana+ bivariate Moran's I (global and local) in one section; cell-type enrichment of the local scores. |
| 10 | `10_STARmapPLUS_BivariateMoransI_AgeComparison.ipynb` | Python | Moran's I in all 4 sections, 8 vs 13 months; null model of detection-matched random gene pairs. |
| 09b | `09b_STARmapPLUS_stlearn.ipynb` | Python | stLearn CCI per cell in one section: LR significance against expression-matched random pairs, comparison with the connectomeDB2020 database, cell-type interactions. |
| 10b | `10b_STARmapPLUS_stlearn_Age_comparison.ipynb` | Python | stLearn in all 4 sections, 8 vs 13 months; comparison with detection-matched database pairs. |
| 11 | `11_Proteomics_LR_AgeTrends.qmd` | R | Protein-level age change of the LR genes (mouse and human), the subset on the STARmap PLUS panel, and the log2FC distribution of all proteins. |

`helper_functions.R` holds the shared functions for 03a and 03c.

### Flow

```
01 → 02 → 03a ─┬─────────────→ 03c → 04 → 05 → 06
               └→ 03b (scRNA) ─┘     │
                                     ├→ 08 → 09 / 10 (Moran's I)
                         07 ─────────┘    └→ 09b / 10b (stLearn)
                                     └→ 11 (proteomics)
```

---

## Data

Raw data are **not** in the repository (see `.gitignore`).

| Data | Location | Notes |
|---|---|---|
| Visium slides | outside the repo (paths in `01_processing_subset_to_HPC.qmd`) | 4 groups (Y, M, O_1, O_2) × 4 slides; hippocampus spots selected manually as `HPC_coords_*.csv`. |
| scRNA-seq reference | outside the repo (path in `03a_RCTD_spotdeconvolution.qmd`) | Annotated Seurat object (`Final_Annotation.rds`). |
| STARmap PLUS | `STARmap_PLUS_data/` | Download from the Broad Single Cell Portal (SCP1375); see `STARmap_PLUS_data/README.md`. |
| Proteomics | `Proteomics_data/Huang_2024/`, `Proteomics_data/Mol_2022/` | Supplementary tables of the two papers; each folder has a README. |

---

## Requirements

**R** (tested with R 4.5.2). Core packages: Seurat, harmony, spacexr, Matrix, hdf5r, dplyr, tidyr, purrr, tibble, ggplot2, patchwork, ggrepel, forcats, ppcor, liana (R), babelgene, readxl. For 04–06: KEGGREST, pathview, org.Mm.eg.db, STRINGdb, rrvgo, GOSemSim, ggraph, tidygraph, igraph.

**Python**, two Jupyter kernels:

- `lianaPlus_localenv` (03b, 09, 10): liana, scanpy, anndata, numpy, pandas, scipy, matplotlib.
- `stlearn` (09b, 10b): stLearn 1.4.1 (Python ≥ 3.12). `setup_stlearn_env.sh` creates this conda environment and registers the kernel.

**Quarto** to render everything (≥ 1.5), and `nbclient` in the Python environment that Quarto uses for notebooks.

---

## Running

Render the whole pipeline in order, stopping at the first failure:

```bash
bash render_all.sh          # everything, 01 → 11
bash render_all.sh 08       # start at 08 (skip 01–07)
```

Logs go to `render_logs/`. The script sets `QUARTO_R` to the R installation used in RStudio (`/opt/R/4.5.2/bin/R`); change it there or export `QUARTO_R` if your R lives elsewhere.

Single files:

```bash
quarto render 03a_RCTD_spotdeconvolution.qmd
quarto render 09b_STARmapPLUS_stlearn.ipynb --execute -M embed-resources:true
```

Notes:

- **03a:** with `eval_RCTD = TRUE`, RCTD runs on every sample × region and needs a lot of memory; lower `max_cores` if the render is killed. Set `eval_RCTD = FALSE` to reuse saved results.
- **05** downloads STRING files (~100 MB) once into `output/05_STRING_NetworkAnalysis/string_cache/`.
- **Seeds** are fixed in all stochastic steps; stLearn's cell-type permutations run on one thread so their p-values are reproducible.

---

## Evidence graph

`evidence_graph/` collects the conclusions for each candidate LR pair:

- **`Evidence_graph.docx`:** the write-up per pair: technology, report, samples, method, description, overall verdict, comments, literature.
- **`evidence.yaml`:** the same content in structured form (the "database").
- **`build_graph.R`:** builds `evidence_graph.html` from the YAML.

```bash
cd evidence_graph
Rscript build_graph.R      # needs yaml, jsonlite
xdg-open evidence_graph.html
```

There is one network per LR pair, chosen from a dropdown. Ligand and receptor are linked through an *interaction* node. Each technology points to the interaction or to a single gene, coloured by the direction of the age effect. Technology nodes are coloured by modality (transcriptomics / proteomics). Clicking a node or an arrow shows the details, the links to the HTML reports and papers, the verdict and the literature. To add a pair or a new result, edit `evidence.yaml` and re-run the script.

---

## Current candidate pairs (summary)

- **Igfbpl1 → Dcc** (Neurogenic_lineage → GC): the spatial methods disagree on the age direction, and Igfbpl1 is barely detected (STARmap PLUS) or not detected (proteomics). The most robust result is a moderate decrease of DCC protein in old mouse hippocampus.
- **Adam10 → Notch2** (Neurogenic_lineage → GC): the spatial evidence is weak and inconsistent, but ADAM10 protein increases with age in mouse (significant) and in all human subfields (descriptive). NOTCH2 protein is not detected.

The details, numbers and caveats are in `evidence_graph/`.

---

## Caveats

- Visium spots (~55 µm) contain several cells, so LR co-expression is not cell-resolved, and spots are treated as independent observations.
- STARmap PLUS "old" means 13 months, which is closer to middle age in the other datasets; the ages differ between all four modalities.
- liana's Moran's I p-value tests against shuffled expression, which spatially patterned genes pass easily. The detection-matched random-pair controls in 10 and 10b are the stricter test.
- The human proteomics has one pooled value per age group, so it is descriptive only.
- The Cajal.Retzius RCTD weights correlate poorly with Cajal–Retzius markers, so the OPC–Cajal.Retzius results should be read with care.

---

## References

- Zeng H, et al. Integrative in situ mapping of single-cell transcriptional states and tissue histopathology in a mouse model of Alzheimer's disease. *Nat Neurosci* 2023. https://doi.org/10.1038/s41593-022-01251-x
- Huang H, et al. Global quantitative proteomic analysis of aged mouse hippocampus. *Proteomics* 2024. https://doi.org/10.1002/pmic.202300276
- Mol P, et al. Age-associated molecular changes in human hippocampus subfields as determined by quantitative proteomics. *OMICS* 2022. https://doi.org/10.1089/omi.2022.0053
- Cable DM, et al. Robust decomposition of cell type mixtures in spatial transcriptomics (RCTD). *Nat Biotechnol* 2022. https://doi.org/10.1038/s41587-021-00830-w
- Dimitrov D, et al. LIANA+ provides an all-in-one framework for cell–cell communication inference. *Nat Cell Biol* 2024. https://doi.org/10.1038/s41556-024-01469-w
- Pham D, et al. Robust mapping of spatiotemporal trajectories and cell–cell interactions in healthy and diseased tissues (stLearn). *Nat Commun* 2023. https://doi.org/10.1038/s41467-023-43120-6

---

## License

The code in this repository is licensed under the GNU General Public License v3.0; see [LICENSE](LICENSE). The third-party datasets (STARmap PLUS, Huang et al. 2024, Mol et al. 2022) are covered by their original terms.
