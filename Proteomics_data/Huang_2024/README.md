# Hippocampal proteomics, young vs old mice (Huang et al. 2024)

Bulk proteomics of the mouse hippocampus, used as protein-level, orthogonal evidence for the LR candidates in this repository. The data has no spatial or cell-type resolution.

## Source

- **Paper:** Huang H, van Waardenberg AJ, Graham ME, Anggono V, Widagdo J. *Global quantitative proteomic analysis of aged mouse hippocampus.* **Proteomics** 24(7):2300276 (2024). DOI: [10.1002/pmic.202300276](https://doi.org/10.1002/pmic.202300276)
- **File here:** `pmic13786-sup-0003-tables1-s3.xlsx` (the paper's Supporting Information, Tables S1–S3)

## Design

- **Mice:** C57BL/6, males only. Young = 3 months, old = 20 months, n = 5 per group.
- **Tissue:** whole hippocampus.
- **Method:** TMT-labelled LC-MS/MS. Differential abundance uses limma-style statistics.
- **Result:** 5,768 proteins quantified. 324 are differentially abundant (FDR < 0.05 and |log2FC| > 0.2): 236 up and 88 down in old.

## Excel sheets

| Sheet | Contents |
|---|---|
| Table S1 | All 5,768 proteins: `Protein ID`, `Genename`, `logFC` (old vs young, positive = higher in old), `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`, plus log2 intensities per sample (`old_1`–`old_5`, `young_1`–`young_5`). |
| Table S2 | GO enrichment of the down- and up-regulated proteins. |
| Table S3 | The 94 genes changed in the same direction at both RNA and protein level. RNA comes from the authors' earlier RNA-seq of a similar cohort (Huang et al. 2023, *Aging Cell* 22:e13755). |

## Caveats

- **Different ages from our data:** the comparison is 3 vs 20 months, with no middle-age group. Our Visium/scRNA data has young, middle and old.
- **Males only:** no female mice.
- **Bulk tissue:** a protein change can't be assigned to a cell type or a location.
