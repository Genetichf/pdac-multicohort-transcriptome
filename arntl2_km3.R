suppressPackageStartupMessages({library(GEOquery); library(survival); library(SummarizedExperiment); library(recount3); library(edgeR)})
gene <- "ARNTL2"
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
get_surv <- function(id) {
  g <- getGEO(id, GSEMatrix = TRUE, getGPL = FALSE)[[1]]; pd <- pData(g)
  pick <- function(rx) pd[[grep(rx, colnames(pd), value = TRUE, ignore.case = TRUE)[1]]]
  data.frame(gsm = pd$geo_accession, time = num(pick("^survival[ _]month")),
             event = as.numeric(grepl("^(1|dead|deceased|death|yes)", trimws(tolower(as.character(pick("^(survival status|cancer_death)")))))))
}
dat <- list()
for (id in c("GSE28735", "GSE62452")) {
  sv <- get_surv(id)
  ex <- read.csv(file.path("results", id, "expression_gene_level.csv"), check.names = FALSE, row.names = 1)
  an <- read.csv(file.path("results", id, "sample_annotation_used.csv"))
  sv <- sv[sv$gsm %in% an$gsm[an$tissue == "Tumor"], ]
  sv$expr <- as.numeric(ex[gene, sv$gsm]); sv <- sv[complete.cases(sv), ]
  dat[[id]] <- sv
}
rt <- readRDS("results/rnaseq/rse_tcga_PAAD.rds")
assay(rt, "counts") <- compute_read_counts(rt)
rt <- rt[, rt$tcga.cgc_sample_sample_type %in% "Primary Tumor" & rt$tcga.cgc_case_histological_diagnosis %in% "Pancreas-Adenocarcinoma Ductal Type"]
dead <- rt$tcga.gdc_cases.diagnoses.vital_status == "dead"
time <- ifelse(dead, num(rt$tcga.gdc_cases.diagnoses.days_to_death), num(rt$tcga.gdc_cases.diagnoses.days_to_last_follow_up)) / 30.44
ok <- is.finite(time) & time > 0 & !is.na(dead)
rt <- rt[, ok]; time <- time[ok]; event <- as.numeric(dead[ok])
lcpm <- cpm(calcNormFactors(DGEList(assay(rt, "counts"))), log = TRUE, prior.count = 1)
rd <- rowData(rt); idx <- which(rd$gene_type == "protein_coding" & rd$gene_name == gene)
idx <- idx[which.max(rowMeans(lcpm[idx, , drop = FALSE]))]
dat[["TCGA-PAAD"]] <- data.frame(time = time, event = event, expr = lcpm[idx, ])
png("FigS3_KM_ARNTL2.png", width = 3300, height = 1100, res = 300)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))
for (id in names(dat)) {
  sv <- dat[[id]]; cat(id, ": n =", nrow(sv), "events =", sum(sv$event), "\n")
  sv$grp <- factor(ifelse(sv$expr > median(sv$expr), "High", "Low"), levels = c("Low", "High"))
  p <- 1 - pchisq(survdiff(Surv(time, event) ~ grp, sv)$chisq, 1)
  cx <- summary(coxph(Surv(time, event) ~ scale(expr), sv))$coefficients
  cat(id, ": logrank p =", signif(p, 3), "| HR per SD =", round(exp(cx[1, 1]), 2), "p =", signif(cx[1, 5], 3), "\n")
  plot(survfit(Surv(time, event) ~ grp, sv), col = c("#2c7bb6", "#d7191c"), lwd = 2, xlab = "Months", ylab = "Survival probability",
       main = paste0(gene, " - ", id, " (n=", nrow(sv), ", events=", sum(sv$event), ")"), cex.main = 0.9)
  legend("topright", c("Low (<= median)", "High (> median)"), col = c("#2c7bb6", "#d7191c"), lwd = 2, bty = "n", cex = 0.8)
  text(0, 0.05, sprintf("log-rank p = %.3g\nHR/SD = %.2f (p = %.3g)", p, exp(cx[1, 1]), cx[1, 5]), adj = 0, cex = 0.8)
}
dev.off()
