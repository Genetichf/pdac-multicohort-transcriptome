library(survival)
M <- read.csv("results/survival/cox_GSE28735_GSE62452.csv", stringsAsFactors = FALSE)
G <- M$gene[M$both]
load_cohort <- function(id) {
  pd <- read.csv(sprintf("results/%s/metadata_preview.csv", id), check.names = FALSE, stringsAsFactors = FALSE)
  ch <- grep("characteristics", names(pd), value = TRUE)
  pick <- function(rx) {
    for (cn in ch) if (any(grepl(rx, pd[[cn]], ignore.case = TRUE))) return(trimws(sub("^[^:]*: *", "", pd[[cn]])))
    rep(NA_character_, nrow(pd))
  }
  sv <- data.frame(gsm = pd$geo_accession, time = suppressWarnings(as.numeric(pick("^survival[ _]month"))),
                   event = suppressWarnings(as.numeric(pick("^(survival status|cancer_death)"))),
                   grade = pick("^grading"), stage = pick("^stage"), stringsAsFactors = FALSE)
  an <- read.csv(sprintf("results/%s/sample_annotation_used.csv", id), stringsAsFactors = FALSE)
  sv <- sv[sv$gsm %in% an$gsm[an$tissue == "Tumor"] & is.finite(sv$time) & is.finite(sv$event) & sv$time > 0, ]
  ex <- read.csv(sprintf("results/%s/expression_gene_level.csv", id), row.names = 1, check.names = FALSE)
  list(sv = sv, ex = as.matrix(ex[G, sv$gsm]))
}
D <- list(GSE28735 = load_cohort("GSE28735"), GSE62452 = load_cohort("GSE62452"))
d <- D$GSE62452$sv
d$adv <- d$stage %in% c("III", "IV", ">IIB")
d$g3 <- d$grade == "G3"
ok <- d$stage %in% c("IA", "IB", "IIA", "IIB", "III", "IV", ">IIB") & d$grade %in% c("G1", "G2", "G3")
cat("GSE62452 tumors with complete stage and grade:", sum(ok), "of", nrow(d), "\n")
out <- do.call(rbind, lapply(G, function(gn) {
  dd <- d; dd$z <- as.numeric(scale(D$GSE62452$ex[gn, ]))
  f0 <- coxph(Surv(time, event) ~ z, data = dd)
  f1 <- coxph(Surv(time, event) ~ z + adv + g3, data = dd[ok, ])
  s1 <- summary(f1)$coefficients
  data.frame(gene = gn, HR_unadj = exp(coef(f0)[["z"]]), HR_adj = exp(s1["z", "coef"]), p_adj = s1["z", "Pr(>|z|)"])
}))
write.csv(out, "results/survival/GSE62452_adjusted_cox.csv", row.names = FALSE)
cat("\nGSE62452: unadjusted vs stage+grade-adjusted HR per SD\n"); print(out, row.names = FALSE, digits = 3)
cols <- c("#3B7DD8", "#D8503B")
pdf("results/survival/KM_top_genes.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
for (gn in head(G, 6)) for (id in names(D)) {
  s <- D[[id]]$sv; x <- D[[id]]$ex[gn, ]
  grp <- factor(ifelse(x > median(x), "High", "Low"), levels = c("Low", "High"))
  lr <- survdiff(Surv(s$time, s$event) ~ grp)
  plot(survfit(Surv(s$time, s$event) ~ grp), col = cols, lwd = 2, xlab = "Months", ylab = "Survival probability",
       main = sprintf("%s, %s (log-rank p = %.3f)", gn, id, 1 - pchisq(lr$chisq, 1)))
  legend("topright", c("Low", "High"), col = cols, lwd = 2, bty = "n")
}
dev.off()
