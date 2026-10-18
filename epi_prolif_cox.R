suppressPackageStartupMessages({ library(SummarizedExperiment); library(recount3); library(edgeR); library(survival) })
genes <- c("ATAD2","NCAPG2","SDR16C5","PARP9","HELLS","TOP2A","MKI67")
pm <- c("MKI67","CCNB1","CDK1","MCM2","PCNA","AURKA")
rt <- readRDS("results/rnaseq/rse_tcga_PAAD.rds")
assay(rt, "counts") <- compute_read_counts(rt)
rt <- rt[, rt$tcga.cgc_sample_sample_type %in% "Primary Tumor" & rt$tcga.cgc_case_histological_diagnosis %in% "Pancreas-Adenocarcinoma Ductal Type"]
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
dead <- rt$tcga.gdc_cases.diagnoses.vital_status == "dead"
time <- ifelse(dead, num(rt$tcga.gdc_cases.diagnoses.days_to_death), num(rt$tcga.gdc_cases.diagnoses.days_to_last_follow_up)) / 30.44
stage <- tolower(as.character(rt$tcga.gdc_cases.diagnoses.tumor_stage))
adv <- ifelse(stage %in% c("stage iii", "stage iv"), 1, ifelse(is.na(stage) | stage == "not reported", NA, 0))
age <- num(rt$tcga.gdc_cases.diagnoses.age_at_diagnosis) / 365.25
ok <- is.finite(time) & time > 0 & !is.na(dead)
rt <- rt[, ok]; time <- time[ok]; event <- as.numeric(dead[ok]); adv <- adv[ok]; age <- age[ok]
cat("Tumors:", ncol(rt), "| deaths:", sum(event), "\n")
rd <- rowData(rt)
lcpm <- cpm(calcNormFactors(DGEList(assay(rt, "counts"))), log = TRUE, prior.count = 1)
pc <- rd$gene_type %in% "protein_coding" & rd$gene_name %in% unique(c(genes, pm))
lcpm <- lcpm[pc, ]; sym <- rd$gene_name[pc]
mu <- rowMeans(lcpm); best <- tapply(seq_along(sym), sym, function(i) i[which.max(mu[i])])
X <- lcpm[unlist(best), ]; rownames(X) <- names(best)
pz <- function(excl) { g <- setdiff(intersect(pm, rownames(X)), excl); colMeans(t(scale(t(X[g, , drop = FALSE])))) }
fit <- function(d, f) { m <- try(coxph(f, data = d), silent = TRUE); if (inherits(m, "try-error")) return(c(NA, NA)); s <- summary(m)$coefficients; c(exp(s["z", "coef"]), s["z", "Pr(>|z|)"]) }
out <- do.call(rbind, lapply(intersect(genes, rownames(X)), function(g) {
  d <- data.frame(time, event, z = as.numeric(scale(X[g, ])), pro = as.numeric(scale(pz(g))), adv, age)
  a <- fit(d, Surv(time, event) ~ z); b <- fit(d, Surv(time, event) ~ z + pro); c3 <- fit(d, Surv(time, event) ~ z + pro + adv + age)
  data.frame(gene = g, r_prolif = round(cor(d$z, d$pro, method = "spearman"), 2), HR = a[1], p = a[2], HR_prolif = b[1], p_prolif = b[2], HR_full = c3[1], p_full = c3[2])
}))
write.csv(out, "results/survival/epi_prolif_adjusted_cox.csv", row.names = FALSE)
print(out, row.names = FALSE, digits = 3)
