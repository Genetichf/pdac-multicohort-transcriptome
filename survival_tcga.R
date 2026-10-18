suppressPackageStartupMessages({ library(SummarizedExperiment); library(recount3); library(edgeR); library(survival) })
rt <- readRDS("results/rnaseq/rse_tcga_PAAD.rds")
assay(rt, "counts") <- compute_read_counts(rt)
sel <- rt$tcga.cgc_sample_sample_type %in% "Primary Tumor" & rt$tcga.cgc_case_histological_diagnosis %in% "Pancreas-Adenocarcinoma Ductal Type"
rt <- rt[, sel]
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
dead <- rt$tcga.gdc_cases.diagnoses.vital_status == "dead"
time <- ifelse(dead, num(rt$tcga.gdc_cases.diagnoses.days_to_death), num(rt$tcga.gdc_cases.diagnoses.days_to_last_follow_up)) / 30.44
stage <- tolower(as.character(rt$tcga.gdc_cases.diagnoses.tumor_stage))
adv <- ifelse(stage %in% c("stage iii", "stage iv"), 1, ifelse(is.na(stage) | stage == "not reported", NA, 0))
age <- num(rt$tcga.gdc_cases.diagnoses.age_at_diagnosis) / 365.25
ok <- is.finite(time) & time > 0 & !is.na(dead)
cat("TCGA ductal primary tumors with survival data:", sum(ok), "| deaths:", sum(dead[ok]), "\n")
rt <- rt[, ok]; time <- time[ok]; event <- as.numeric(dead[ok]); adv <- adv[ok]; age <- age[ok]
S <- read.csv("results/meta/meta_three_datasets_DEGs.csv", stringsAsFactors = FALSE)
M <- read.csv("results/survival/cox_GSE28735_GSE62452.csv", stringsAsFactors = FALSE)
rd <- rowData(rt)
lcpm <- cpm(calcNormFactors(DGEList(assay(rt, "counts"))), log = TRUE, prior.count = 1)
pc <- rd$gene_type == "protein_coding" & rd$gene_name %in% S$gene
lcpm <- lcpm[pc, ]; sym <- rd$gene_name[pc]
mu <- rowMeans(lcpm); best <- tapply(seq_along(sym), sym, function(i) i[which.max(mu[i])])
X <- lcpm[unlist(best), ]; rownames(X) <- names(best)
res <- t(sapply(rownames(X), function(g) {
  z <- as.numeric(scale(X[g, ]))
  s0 <- summary(coxph(Surv(time, event) ~ z))$coefficients
  d <- data.frame(time, event, z, adv, age)
  f1 <- try(coxph(Surv(time, event) ~ z + adv + age, data = d), silent = TRUE)
  s1 <- if (inherits(f1, "try-error")) c(NA, NA) else summary(f1)$coefficients["z", c("coef", "Pr(>|z|)")]
  c(s0[1, "coef"], s0[1, "Pr(>|z|)"], s1)
}))
colnames(res) <- c("logHR", "p", "logHR_adj", "p_adj")
R <- data.frame(gene = rownames(res), res, row.names = NULL)
R$HR <- exp(R$logHR); R$HR_adj <- exp(R$logHR_adj)
write.csv(R, "results/survival/cox_TCGA_PAAD.csv", row.names = FALSE)
cat("Genes tested:", nrow(R), "| nominal p<0.05:", sum(R$p < 0.05, na.rm = TRUE), "(chance ~", round(0.05 * nrow(R)), ")\n")
J <- merge(M, R[, c("gene", "logHR", "p", "HR", "HR_adj", "p_adj")], by = "gene", suffixes = c("", "_T"))
cat("Correlation of TCGA logHR with mean array logHR across", nrow(J), "genes (Spearman):", round(cor((J$logHR_A + J$logHR_B) / 2, J$logHR, method = "spearman", use = "complete.obs"), 3), "\n")
C <- J[J$both, c("gene", "HR_A", "HR_B", "HR", "p", "HR_adj", "p_adj")]
C$same_dir <- sign(log(C$HR_A)) == sign(log(C$HR))
cat("\nThe 15 candidates in TCGA-PAAD (HR per SD):\n"); print(C[order(C$p), ], row.names = FALSE, digits = 3)
cat("Same direction as arrays:", sum(C$same_dir), "of", nrow(C), "| p<0.05:", sum(C$p < 0.05), "| adjusted p<0.05:", sum(C$p_adj < 0.05, na.rm = TRUE), "\n")
