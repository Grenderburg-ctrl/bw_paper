# QC: figures/source_data/FigVG_ext_source_data.csv

Generated 2026-10-12 07:45 UTC. Source: the data objects that `R/figVG_forest.R` passes to ggplot (`make_rows()` + `prep_rows()` for panels A, B, C of FigVG_ext); statistics from `results/deseq2_R_results_no_exclusion.csv`. Nothing recomputed, no shrinkage, no TPM.

| check | expected | obtained | OK |
|---|---|---|---|
| Rows in table | 110 | 110 | ✅ |
| plotted + not_tested_padjNA + not_detected | 110 | 110 | ✅ |
| Gene set = scripts/voltage_gated_genes.csv, no duplicates | TRUE | TRUE | ✅ |
| Significant (padj < 0.05) | 21 | 21 | ✅ |
| ... down / up | 19 / 2 | 19 / 2 | ✅ |
| direction consistent with padj and sign | TRUE | TRUE | ✅ |
| padj_stars consistent with padj | TRUE | TRUE | ✅ |
| Reference Scn5a | -0.96; 0.20; 3.3e-05 | -0.96; 0.20; 3.3e-05 | ✅ |
| Reference Scn3b | -2.06; 0.32; 4.1e-09 | -2.06; 0.32; 4.1e-09 | ✅ |
| Reference Hcn4 | -1.90; 0.33; 1.9e-07 | -1.90; 0.33; 1.9e-07 | ✅ |
| Reference Cacna2d4 | 12.66; 1.03; 1.2e-31 | 12.66; 1.03; 1.2e-31 | ✅ |
| Reference Kcnh7 | 6.13; 1.30; 4.0e-05 | 6.13; 1.30; 4.0e-05 | ✅ |
| Max |diff| vs results CSV (plotted): log2FC / lfcSE / padj | 0 / 0 / 0 | 0 / 0 / 0 | ✅ |
| Excel (ru) copy: max |diff| vs results CSV | 0 | 0 | ✅ |
| Excel copy identical to main CSV in text columns | TRUE | TRUE | ✅ |
| whisker ends = log2FC -/+ lfcSE | TRUE | TRUE | ✅ |
| Non-plotted rows have empty numeric/derived columns | TRUE | TRUE | ✅ |
| off_axis genes | Cacna2d4, Kcnh7, Kcnt1 | Cacna2d4, Kcnh7, Kcnt1 | ✅ |
| Gene symbol + protein label found in FigVG_ext.pdf text | all 70 | 70 | ✅ |
| Group (strip) labels found in FigVG_ext.pdf text | all 27 | 27 | ✅ |
| Row order top-to-bottom in PDF = table order (A then C; B) | TRUE | TRUE | ✅ |
| Panels A/C in left column, B in right column (x of labels) | TRUE | TRUE | ✅ |

Status counts: plotted 70, not_tested_padjNA 6, not_detected 34. Genes per panel: A 40, B 54, C 16 (plotted: A 26, B 33, C 11).
low_expr = TRUE for 21 plotted genes (definition in the figure script: mean DESeq2-normalised counts < 5 in at least one group): Scn2a, Scn2b, Cacnb4, Cacna2d3, Cacna2d4, Cacng7, Kcna3, Kcnc1, Kcnc4, Kcnd1, Kcng1, Kcnv1, Kcnh1, Kcnh7, Kcnab3, Kcne2, Kcne3, Kcnip4, Kcnma1, Kcnt1, Catsper4.

Visual comparison with `figures/FigVG_ext.tiff` (opened as an image): all channel labels and group labels of panels A, B and C are written identically in the table (gene in italics + protein in the figure; group labels are wrapped in two lines on the figure, e.g. "Cav1 / (L-type)", "Kv5/6/8/9" under "KvS"; the table keeps them on one line). No spelling differences found. The two grey "auxiliary subunits" sub-header rows on the figure are layout elements, not genes, and are not in the table.

Not plotted (end of table): panel A — not_tested_padjNA: Cacna1e, Cacng5; not_detected: Scn8a, Scn9a, Scn10a, Scn11a, Cacna1s, Cacna1f, Cacna1b, Cacna1i, Cacng1, Cacng2, Cacng3, Cacng8. Panel B and C: see rows with `status` ≠ `plotted`.

Formats: `log2FC`, `lfcSE`, whisker ends and `padj` are written with 15 significant digits (as stored in `results/deseq2_R_results_no_exclusion.csv`); `padj` in scientific notation; only `log2FC_pm_lfcSE` is rounded (2 decimals, hyphen-minus). Main file: comma-separated, UTF-8, NA = empty. Excel copy: `;` separator, decimal comma, UTF-8 with BOM.
