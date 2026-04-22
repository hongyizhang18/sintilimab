# Code for LungMate-009

Code for reproducing figures from:

Tao G, Zhang H, Zhu X, et al. *Efficiency and safety of sintilimab-based induction therapy in potentially resectable stage III non-small cell lung cancer (LungMate-009): an open-label, phase 2 trial*. **EClinicalMedicine**. 2025;90:103669. doi:10.1016/j.eclinm.2025.103669

The repository contains R scripts for main and supplementary figures, plus one notebook for Scissor-based single-cell integration.

## Contents

| File | Focus | Output |
| --- | --- | --- |
| `Figure1A.R` | Differential expression in tumor samples | Volcano plots |
| `Figure1B.R` | Immune features, ESTIMATE, marker genes | Heatmap |
| `Figure1D.R` | Scissor-labeled single-cell visualization | UMAP, circlize plot |
| `Figure1E.R` | Plasma/immune signatures and pathway scoring | Signature comparisons |
| `Figure1F.R` | Plasma signature in baseline samples | Boxplot |
| `Figure2A.R` | Broad pathway scoring and entropy-like metrics | Pathway-level analysis |
| `Figure2B.R` | Metabolic pathway activity across groups | Heatmap |
| `Figure2C.R` | Metabolic GSEA | Enrichment plot |
| `Figure2D.R` | AKR1C signature analysis | Violin plot |
| `Figure2E.R` | OAK/POPLAR validation of AKR1C signature | Survival analysis |
| `SFigure1A.R` | Cohort/sample annotation summary | Annotated heatmap |
| `SFigure1C.R` | PCA of cohort subsets | PCA plots |
| `SFigure1D.R` | Differential expression in lymph node samples | Volcano plots |
| `SFigure1G.R` | OAK/POPLAR validation of plasma-related signatures | Survival and forest plots |
| `sicssor.ipynb` | Scissor workflow and export of `subset_data.h5` | Notebook record |

## Dependencies

Main R packages include:

- `tidyverse`
- `ggplot2`
- `ggrepel`
- `edgeR`
- `DESeq2`
- `GSVA`
- `clusterProfiler`
- `enrichplot`
- `pheatmap`
- `ComplexHeatmap`
- `circlize`
- `Seurat`
- `dior`
- `plot1cell`
- `survival`
- `survminer`
- `openxlsx`
- `msigdbr`

The notebook additionally uses `Scissor`, `SeuratObject`, `Matrix`, `SingleCellExperiment`, and related packages.

## Notes

- This is a script collection, not a packaged pipeline.
- Input data and helper files are not bundled.
- Several scripts depend on local paths or precomputed intermediate files.
- `sicssor.ipynb` contains exploratory steps and recorded errors; it should be treated as a workflow record, not a finalized pipeline.

## Citation

Article DOI: `10.1016/j.eclinm.2025.103669`
