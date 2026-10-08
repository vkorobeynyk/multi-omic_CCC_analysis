# Human hippocampal subfield proteomics across age (Mol et al. 2022)

Human proteomics of individual hippocampal subfields from post-mortem tissue. We use it as protein-level, cross-species evidence for the LR candidates in this repository; it has subregion resolution but no cell-type resolution.

## Source

- **Paper:** Mol P, Chatterjee O, Gopalakrishnan L, et al., Prasad TSK. *Age-associated molecular changes in human hippocampus subfields as determined by quantitative proteomics.* **OMICS** 26(7):382–391 (2022). DOI: [10.1089/omi.2022.0053](https://doi.org/10.1089/omi.2022.0053)
- **File here:** `sj-xlsx-2-omi-10.1089_omi.2022.0053.xlsx`, which is the paper's Supplementary Table S2.

## Design

- **Donors:** 12 histopathologically normal donors from the NIMHANS Human Brain Tissue Repository (Bangalore), 3 per age group:

  | Age group | Ages |
  |---|---|
  | AG1 | 1–10 years |
  | AG2 | 21–30 years |
  | AG3 | 31–40 years |
  | AG4 | 81–90 years |

  There are no donors between 41 and 79 years.
- **Tissue:** formalin-fixed, paraffin-embedded (FFPE) sections of posterior hippocampus. Five subfields were manually microdissected: DG, CA4, CA3, CA2 and CA1.
- **Method:** TMT 10-plex, Orbitrap Fusion, with one TMT run per subfield and one channel per age group.
- **Statistics (in the paper):** AG4 is the reference group. Hits are called with t-tests across technical replicates (Perseus), using fold change ≥ 1.5 and p ≤ 0.05.
- **Proteins identified:** 5,698 in DG, 5,508 in CA4, 5,428 in CA3, 4,271 in CA2 and 4,773 in CA1.

## Excel sheets

The first sheet, `S2 Table Information`, lists the captions. Then there is one sheet per subfield (`DG`, `CA1`, `CA2`, `CA3`, `CA4`), each with these columns:

- RefSeq accession (`NP_Accession`)
- `Gene Symbol` and `Gene ID`
- `Protein name`
- `# PSMs` and `# Unique Peptides`
- four `Abundances (Grouped)` columns, one per age group

Things to watch when reading it:

- **Header row:** there are 2 title rows above the header, so read with `skip = 2`.
- **Column order differs between sheets:** `Gene Symbol` and `Gene ID` swap places in CA1, so select columns by name, not position.
- **Age-group column names differ between sheets:**
  - CA1 and CA3 use `A1`–`A4`.
  - CA4 uses `AG1`–`AG4`.
  - CA2 uses channel plus group, e.g. `130C_AG1`.
  - DG uses only TMT channel names (`126`, `127N`, `127C`, `128N`), with no age labels. The channel-to-age mapping is in Supplementary Table S1, which isn't included here.
- **Values:** raw grouped reporter-ion abundances, not log-transformed. There are no fold-change or p-value columns.
- **Isoforms:** one gene can appear on several rows (RefSeq isoforms, e.g. PLEC).

## Caveats

- **No biological replicates:** there is one abundance value per age group per subfield, so the 3 donors per group can't be compared with each other.
- **Very small sample:** n = 3 per group, with an age gap from 41 to 79 years.
- **Donor and tissue differences:** the samples are human FFPE autopsy tissue, so expect donor variability.
