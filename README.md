# PDAC multi-cohort transcriptome analysis

R scripts for a multi-cohort transcriptome analysis of pancreatic ductal adenocarcinoma (PDAC):
tumor-versus-normal meta-analysis of six microarray cohorts, tissue-composition annotation, survival,
TCGA/GTEx RNA-seq confirmation, DepMap dependencies, published-subtype overlap, and single-cell localisation
in two public atlases.

Author: Hassan Dastsooz, PhD (IIGM-Epigenetics Unit, University of Turin).

License: MIT.

Code archive (v1.0, first analysis): https://doi.org/10.5281/zenodo.23459619

Code archive (v2.0, final analysis with corrected clustering): https://doi.org/10.5281/zenodo.24136446

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

## Final analysis with corrected clustering (v2.0)
The module analysis was rerun with hierarchical clustering by Ward.D2 on sqrt(2(1 - rho)) (k = 6). Scripts in `ward_corrected/` reproduce
every module-dependent result in the manuscript (modules of 69/65/59/49/22/16 genes; 22 Module 1 candidates). Steps that do not depend on
the clustering (meta-analysis, tissue composition, pooled survival, RNA-seq, DepMap) are unchanged from v1.0. Run from `~/pdac_geo`, in this order:
1. `modules_subtype_ward_corrected.R` - module assignments, module-subtype correlations, module and prognostic-gene Cox models.
2. `module_overlap_stats_ward_corrected.R` - hypergeometric overlap with Moffitt, Bailey and Puleo lists; BH over the unique gene sets (280-gene primary background, 720-gene sensitivity background).
3. `module1_candidate_evidence_ward_corrected.R` - the 22 Module 1 candidates not in published lists: bulk, pooled survival and TCGA evidence.
   The Europe PMC publication counts in the manuscript were obtained in a separate lookup (see `epmc_fix.R` and `epmc_fix2.R`); the count column of this script's output is not used.
4. `sc_sample_level_ward_corrected.R` - Module 1 score and key genes at sample level in both atlases (ductal cells, exact Wilcoxon tests).
5. `sc_sample_level_candidates_ward_corrected.R` - the 22 candidates at sample level, BH within each atlas.
6. `sc_mucin_within_ward_corrected.R` and `sc_mucin_depth_ward_corrected.R` - mucin correlations within tumours, and after adjusting for sequencing depth.
7. `make_fig6_ward_corrected.R` - Figure 6. (`sc_module_check_ward_corrected.R` is a preliminary check and is not used for the reported results.)

Supplementary table numbering in the final manuscript follows citation order: S1 leave-one-out and REML/Hartung-Knapp sensitivity, S2 compartment markers,
S3 module membership, S4 candidate evidence (S4a 22 candidates, S4b literature counts, S4c protein evidence), S5 signature, S6 compartment vs single cell,
S7 survival, S8 module-signature overlap, S9 DepMap, S10 single-cell sample-level results (S10a-S10e).
Result tables for this version are deposited with the v2.0 data archive.

## utils/ (metadata inspection, no results written)
`inspect_candidates.R` and `check_candidates2.R` (vetting candidate microarray cohorts and GSE16515 pairing),
`sc_candidates.R`, `sc_samples2.R`, `sc_inspect.R`, `sc_inspect2.R` (choosing and inspecting the single-cell atlases).
They print GEO metadata to the screen and are provided for transparency of cohort selection.

## Notes
- DEG rule: FDR < 0.05, |pooled log2FC| >= 1, same direction in every cohort.
- Survival hazard ratios are per SD of expression; pooling is inverse-variance, BH-adjusted.
- Compartment marker sets and the assignment rule (Spearman rho >= 0.4 with marker scores in TCGA-PAAD) are given in Supplementary Table S4 of the manuscript.
- Results tables used in the manuscript are deposited with the archived release (see the result-tables DOI above).
