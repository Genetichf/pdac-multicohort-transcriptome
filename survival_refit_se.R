suppressPackageStartupMessages({ library(survival); library(SummarizedExperiment); library(recount3); library(edgeR) })
old <- read.csv("results/survival/survival_meta_fdr.csv", stringsAsFactors = FALSE)
genes <- old$gene
cat("Genes in current pooled table:", length(genes), "\n")

## ---- arrays (same filters as survival_meta.R) ----
get_surv <- function(id) {
  pd <- read.csv(sprintf("results/%s/metadata_preview.csv", id), check.names = FALSE, stringsAsFactors = FALSE)
  ch <- grep("characteristics", names(pd), value = TRUE)
  pick <- function(rx) {
    for (cn in ch) if (any(grepl(rx, pd[[cn]], ignore.case = TRUE)))
      return(suppressWarnings(as.numeric(sub("^[^:]*: *", "", pd[[cn]]))))
    rep(NA_real_, nrow(pd))
  }
  data.frame(gsm = pd$geo_accession, time = pick("^survival[ _]month"), event = pick("^(survival status|cancer_death)"))
}
fit_arr <- function(id, tag) {
  ex <- read.csv(sprintf("results/%s/expression_gene_level.csv", id), row.names = 1, check.names = FALSE)
  an <- read.csv(sprintf("results/%s/sample_annotation_used.csv", id), stringsAsFactors = FALSE)
  sv <- get_surv(id)
  sv <- sv[sv$gsm %in% an$gsm[an$tissue == "Tumor"] & is.finite(sv$time) & is.finite(sv$event) & sv$time > 0, ]
  cat(id, ": tumors =", nrow(sv), "| events =", sum(sv$event == 1), "\n")
  g <- intersect(genes, rownames(ex))
  r <- t(sapply(g, function(gn) {
    z <- as.numeric(scale(as.numeric(ex[gn, sv$gsm])))
    f <- try(coxph(Surv(sv$time, sv$event) ~ z), silent = TRUE)
    if (inherits(f, "try-error")) return(c(NA, NA, NA))
    s <- summary(f)$coefficients
    c(s[1, "coef"], s[1, "se(coef)"], s[1, "Pr(>|z|)"])
  }))
  out <- data.frame(gene = g, r, row.names = NULL)
  names(out) <- c("gene", paste0(c("b_", "se_", "p_"), tag))
  out
}
A <- fit_arr("GSE28735", "A")
B <- fit_arr("GSE62452", "B")

## ---- TCGA-PAAD (same filters as survival_tcga.R) ----
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
cat("TCGA tumors =", sum(ok), "| deaths =", sum(dead[ok]), "\n")
rt <- rt[, ok]; time <- time[ok]; event <- as.numeric(dead[ok]); adv <- adv[ok]; age <- age[ok]
rd <- rowData(rt)
lcpm <- cpm(calcNormFactors(DGEList(assay(rt, "counts"))), log = TRUE, prior.count = 1)
pc <- rd$gene_type == "protein_coding" & rd$gene_name %in% genes
lcpm <- lcpm[pc, ]; sym <- rd$gene_name[pc]
mu <- rowMeans(lcpm); best <- tapply(seq_along(sym), sym, function(i) i[which.max(mu[i])])
X <- lcpm[unlist(best), ]; rownames(X) <- names(best)
r <- t(sapply(rownames(X), function(g) {
  z <- as.numeric(scale(X[g, ]))
  s0 <- summary(coxph(Surv(time, event) ~ z))$coefficients
  d <- data.frame(time, event, z, adv, age)
  f1 <- try(coxph(Surv(time, event) ~ z + adv + age, data = d), silent = TRUE)
  s1 <- if (inherits(f1, "try-error")) c(NA, NA, NA) else summary(f1)$coefficients["z", c("coef", "se(coef)", "Pr(>|z|)")]
  c(s0[1, "coef"], s0[1, "se(coef)"], s0[1, "Pr(>|z|)"], s1)
}))
Tt <- data.frame(gene = rownames(r), r, row.names = NULL)
names(Tt) <- c("gene", "b_T", "se_T", "p_T", "b_Tadj", "se_Tadj", "p_Tadj")

## ---- pool directly (inverse variance, fixed effect) ----
N <- Reduce(function(x, y) merge(x, y, by = "gene"), list(A, B, Tt))
N <- N[complete.cases(N[, c("b_A", "se_A", "b_B", "se_B", "b_T", "se_T")]), ]
b <- cbind(N$b_A, N$b_B, N$b_T); w <- 1 / cbind(N$se_A, N$se_B, N$se_T)^2
N$b_pool <- rowSums(w * b) / rowSums(w); N$se_pool <- 1 / sqrt(rowSums(w))
N$HR_pool <- exp(N$b_pool)
N$p_pool <- 2 * pnorm(-abs(N$b_pool / N$se_pool))
N$fdr_pool <- p.adjust(N$p_pool, "BH")
Q <- rowSums(w * (b - N$b_pool)^2); N$I2 <- round(pmax(0, (Q - 2) / Q) * 100)
N$same_dir <- sign(N$b_A) == sign(N$b_B) & sign(N$b_B) == sign(N$b_T)
write.csv(N, "results/survival/survival_meta_refit_se.csv", row.names = FALSE)

## ---- compare with the current (p-value-derived) table ----
M <- merge(old[, c("gene", "HR_A", "HR_B", "HR_T", "HR_pooled", "fdr_pooled", "same_dir")], N, by = "gene", suffixes = c("_old", ""))
cat("\n--- consistency of refit with saved HRs (max |difference| in log HR) ---\n")
cat("GSE28735:", signif(max(abs(log(M$HR_A) - M$b_A), na.rm = TRUE), 3),
    "| GSE62452:", signif(max(abs(log(M$HR_B) - M$b_B), na.rm = TRUE), 3),
    "| TCGA:", signif(max(abs(log(M$HR_T) - M$b_T), na.rm = TRUE), 3), "\n")
cat("Genes pooled (new):", nrow(N), "\n")
cat("p-values below the 5e-5 clamp: A =", sum(N$p_A < 5e-5), "B =", sum(N$p_B < 5e-5), "T =", sum(N$p_T < 5e-5), "\n")
cat("FDR<0.05 old:", sum(old$fdr_pooled < 0.05, na.rm = TRUE), "| new:", sum(N$fdr_pool < 0.05),
    "| new and same direction in all 3:", sum(N$fdr_pool < 0.05 & N$same_dir), "\n")
cat("In both:", sum(M$fdr_pooled < 0.05 & M$fdr_pool < 0.05, na.rm = TRUE),
    "| lost:", paste(M$gene[M$fdr_pooled < 0.05 & M$fdr_pool >= 0.05], collapse = ", "),
    "| gained:", paste(M$gene[M$fdr_pooled >= 0.05 & M$fdr_pool < 0.05], collapse = ", "), "\n")
cat("Spearman of pooled p (old vs new):", round(cor(M$fdr_pooled, M$fdr_pool, method = "spearman", use = "complete.obs"), 4), "\n")
