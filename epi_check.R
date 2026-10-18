suppressPackageStartupMessages({ library(SummarizedExperiment); library(recount3); library(edgeR); library(survival) })
genes <- c("ATAD2","HELLS","KDM2A","L3MBTL3","SMYD3","JADE1","KAT2B","SDR16C5","ZNFX1","PARP9","NCAPG2","SMC6")
mk <- list(stroma = c("COL1A1","COL1A2","COL3A1","DCN","LUM","FAP"), immune = c("PTPRC","CD53","LAPTM5","CORO1A","CD3E","CD68"),
           epithelial = c("KRT19","KRT7","EPCAM","KRT8","KRT18","CLDN4"), acinar = c("PRSS1","CPA1","PNLIP","CTRB1"),
           prolif = c("MKI67","CCNB1","CDK1","MCM2","PCNA","AURKA"), classical = c("GATA6","TFF1","AGR2","LGALS4"), basal = c("KRT5","KRT6A","TP63","S100A2"))
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
cat("Tumors with survival data:", ncol(rt), "| deaths:", sum(event), "\n")
rd <- rowData(rt)
lcpm <- cpm(calcNormFactors(DGEList(assay(rt, "counts"))), log = TRUE, prior.count = 1)
pc <- rd$gene_type %in% "protein_coding" & rd$gene_name %in% unique(c(genes, unlist(mk)))
lcpm <- lcpm[pc, ]; sym <- rd$gene_name[pc]
mu <- rowMeans(lcpm); best <- tapply(seq_along(sym), sym, function(i) i[which.max(mu[i])])
X <- lcpm[unlist(best), ]; rownames(X) <- names(best)
sc <- sapply(mk, function(g) { g <- intersect(g, rownames(X)); colMeans(t(scale(t(X[g, , drop = FALSE])))) })
G <- intersect(genes, rownames(X)); cat("Genes found in TCGA data:", length(G), "of", length(genes), "\n")
r <- t(sapply(G, function(g) as.numeric(cor(X[g, ], sc, method = "spearman"))))
colnames(r) <- paste0("r_", colnames(sc))
main <- c("stroma", "immune", "epithelial", "acinar")
rm4 <- r[, paste0("r_", main)]
comp <- ifelse(apply(rm4, 1, max) >= 0.4, main[apply(rm4, 1, which.max)], "unassigned")
cx <- t(sapply(G, function(g) {
  z <- as.numeric(scale(X[g, ])); d <- data.frame(time, event, z, adv, age)
  s0 <- summary(coxph(Surv(time, event) ~ z, data = d))$coefficients
  f1 <- try(coxph(Surv(time, event) ~ z + adv + age, data = d), silent = TRUE)
  s1 <- if (inherits(f1, "try-error")) c(NA, NA) else summary(f1)$coefficients["z", c("coef", "Pr(>|z|)")]
  c(HR = exp(s0[1, "coef"]), p = s0[1, "Pr(>|z|)"], HR_adj = exp(s1[1]), p_adj = s1[2])
}))
E <- read.csv("results/meta/epigenetic_regulators_relaxed.csv", stringsAsFactors = FALSE)[, c("gene", "lfc", "I2", "pdac")]
D <- merge(E, data.frame(gene = G, comp = comp, round(r, 2), round(cx, 3), row.names = NULL), by = "gene")
D$prolif_like <- D$r_prolif >= 0.5
D <- D[order(-D$r_prolif), ]
write.csv(D, "results/meta/epigenetic_shortlist_check.csv", row.names = FALSE)
print(D, row.names = FALSE, digits = 3)
