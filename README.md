# PDAC multi-cohort transcriptome analysis

R scripts for a multi-cohort transcriptome analysis of pancreatic ductal adenocarcinoma (PDAC):
tumor-versus-normal meta-analysis of six microarray cohorts, tissue-composition annotation, survival,
TCGA/GTEx RNA-seq confirmation, DepMap dependencies, published-subtype overlap, and single-cell localisation
in two public atlases.

Author: Hassan Dastsooz, PhD (IIGM-Epigenetics Unit, University of Turin).

License: MIT.

Code archive (v1.0): https://doi.org/10.5281/zenodo.23459619

Code archive (all versions): https://doi.org/10.5281/zenodo.23459618

Result tables (data): https://doi.org/10.5281/zenodo.23468919

## Data sources (not included; download from the original sources)
- Microarray cohorts (GEO): GSE28735, GSE62452, GSE15471, GSE16515, GSE62165 (core); GSE32676 (sensitivity only)
- RNA-seq: TCGA-PAAD and GTEx pancreas via recount3 (saved as `results/rnaseq/rse_tcga_PAAD.rds`)
- DepMap CRISPR gene effect, release 26Q3
- Single cell (GEO): GSE155698 (atlas 1), GSE212966 (atlas 2)
- Published subtype lists: pdacR repository (Moffitt, Bailey, Puleo)

## Requirements
R with GEOquery, limma, Biobase, survival, edgeR, recount3, SummarizedExperiment, Seurat (v5), metafor-style
random effects coded in base R. Record `sessionInfo()` with your run.
Run everything from the project folder (`~/pdac_geo`); each script reads and writes under `results/`.

## Run order
1. Per-cohort limma: `geo_limma_pdac.R <GSE_ID> results` (paired, duplicateCorrelation or unpaired design set in its CONFIG block).
   Helpers: `check_pairs.R`, `patch_remap.R` (GSE62452 pairing), `sens_GSE28735.R` (sensitivity), `compare_datasets.R`, `compare_three.R`.
2. Meta-analysis: `meta_three.R`, `meta_tiers.R` (core and extended tiers, robust DEGs), `enrich_meta.R` (GO enrichment).
3. Survival: `survival_meta.R`, `robust_survival.R`, `survival_robust_arrays.R`, `survival_tcga.R`, `fdr_survival.R`, `fdr_survival2.R`, `epi_prolif_cox.R`; ARNTL2 curves: `arntl2_km3.R`.
   Refit with true standard errors (replaces the p-value-derived SE of `fdr_survival2.R` for the reported pooled results): `survival_refit_se.R`.
4. RNA-seq arm: `rnaseq_inspect.R`, `rnaseq_de.R` (TCGA tumors vs GTEx pancreas, limma-voom).
5. Prioritisation and dependencies: `prioritize.R`, `depmap_check.R`, `depmap_epi.R`, `epigenetic_relaxed.R`, `epi_check.R`, `funnel_tcga.R`, `funnel_merge.R`.
   DepMap sensitivity analysis (all other lines vs epithelial-type lines, strict rule): `depmap_sensitivity.R`.
6. Subtype overlap and modules: `get_signatures.R`, `modules_subtype.R`, `module_overlap_stats.R`, `module_overlap_stats_v2.R` (hypergeometric overlap against the 280-gene clustering background and the 720-gene background; used for the reported statistics), `novel18_evidence.R`.
   Compartment validation: `compartment_vs_singlecell.R` (marker-based compartment vs single-cell assignment for each gene in both atlases).
7. Literature counts: `pubmed_counts.R`, `alias_counts.R`, `epmc_check.R`, `epmc_fix.R`, `epmc_fix2.R`.
8. Single cell: atlas 1 `sc_analysis.R`, `sc_module_check.R`, `sc_basal_check.R`, `sc_mucin_check.R`, `sc_mucin_within.R`, `sc_mucin_depth.R`;
   atlas 2 `sc2_build.R`, `sc2_tests.R`, `sc2_adj5.R`, `sc2_target_table.R`, `sc2_figure.R`;
   sample-level comparisons for both atlases: `sc_sample_level.R`.
9. Figures: `make_figures.R`, `sc2_figure.R`.

## utils/ (metadata inspection, no results written)
`inspect_candidates.R` and `check_candidates2.R` (vetting candidate microarray cohorts and GSE16515 pairing),
`sc_candidates.R`, `sc_samples2.R`, `sc_inspect.R`, `sc_inspect2.R` (choosing and inspecting the single-cell atlases).
They print GEO metadata to the screen and are provided for transparency of cohort selection.

## Notes
- DEG rule: FDR < 0.05, |pooled log2FC| >= 1, same direction in every cohort.
- Survival hazard ratios are per SD of expression; pooling is inverse-variance, BH-adjusted.
- Compartment marker sets and the assignment rule (Spearman rho >= 0.4 with marker scores in TCGA-PAAD) are given in Supplementary Table S4 of the manuscript.
- Results tables used in the manuscript are deposited with the archived release (see the result-tables DOI above).