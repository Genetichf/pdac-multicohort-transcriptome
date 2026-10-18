suppressPackageStartupMessages(library(survival))
D <- read.csv("results/meta/robust_funnel_table.csv", stringsAsFactors = FALSE)
get_surv <- function(id) {
  pd <- read.csv(sprintf("results/%s/metadata_preview.csv", id), check.names = FALSE, stringsAsFactors = FALSE)
  ch <- grep("characteristics", names(pd), value = TRUE)
  pick <- function(rx) { for (cn in ch) if (any(grepl(rx, pd[[cn]], ignore.case = TRUE))) return(suppressWarnings(as.numeric(sub("^[^:]*: *", "", pd[[cn]])))); rep(NA_real_, nrow(pd)) }
  data.frame(gsm = pd$geo_accession, time = pick("^survival[ _]month"), event = pick("^(survival status|cancer_death)"))
}
res <- list()
for (id in c("GSE28735", "GSE62452")) {
  ex <- read.csv(sprintf("results/%s/expression_gene_level.csv", id), row.names = 1, check.names = FALSE)
  an <- read.csv(sprintf("results/%s/sample_annotation_used.csv", id), stringsAsFactors = FALSE)
  sv <- get_surv(id); sv <- sv[sv$gsm %in% an$gsm[an$tissue == "Tumor"] & is.finite(sv$time) & is.finite(sv$event) & sv$time > 0, ]
  cat(id, ": tumors with survival =", nrow(sv), "| events =", sum(sv$event == 1), "\n")
  g <- intersect(D$gene, rownames(ex))
  r <- t(sapply(g, function(gn) { z <- as.numeric(scale(as.numeric(ex[gn, sv$gsm])))
    f <- try(coxph(Surv(sv$time, sv$event) ~ z), silent = TRUE); if (inherits(f, "try-error")) return(c(NA, NA))
    s <- summary(f)$coefficients; c(exp(s[1, "coef"]), s[1, "Pr(>|z|)"]) }))
  res[[id]] <- data.frame(gene = g, HR = r[, 1], p = r[, 2])
}
M <- merge(res[[1]], res[[2]], by = "gene", suffixes = c("_A", "_B"))
M <- merge(M, D[, c("gene", "new", "comp", "HR", "p", "p_pro", "lfc_core", "I2_core")], by = "gene")
names(M)[names(M) == "HR"] <- "HR_T"; names(M)[names(M) == "p"] <- "p_T"
M$n_sig <- (M$p_A < 0.05) + (M$p_B < 0.05) + (M$p_T < 0.05)
M$same <- sign(log(M$HR_A)) == sign(log(M$HR_B)) & sign(log(M$HR_B)) == sign(log(M$HR_T))
write.csv(M[order(-M$n_sig, M$p_T), ], "results/survival/cox_arrays_robust.csv", row.names = FALSE)
cat("\nGenes tested in both arrays + TCGA:", nrow(M), "| same direction in all 3:", sum(M$same, na.rm = TRUE), "| p<0.05 in >=2 cohorts (same direction):", sum(M$n_sig >= 2 & M$same, na.rm = TRUE), "| in all 3:", sum(M$n_sig == 3 & M$same, na.rm = TRUE), "\n")
g18 <- c("ARNTL2","MET","LAMA3","RARRES3","HOXB2","ITGA2","TPX2","ITGA3","ITGB6","GPR87","LAMC2","NT5E","GNA15","DMBT1","LYZ","F5","COL17A1","LPCAT2")
cols <- c("gene", "new", "HR_A", "p_A", "HR_B", "p_B", "HR_T", "p_T", "n_sig", "same")
cat("\nThe 18 TCGA-prognostic genes in all three cohorts:\n"); print(M[match(g18[g18 %in% M$gene], M$gene), cols], row.names = FALSE, digits = 3)
cat("\nAny gene with p<0.05 in >=2 cohorts, same direction, and p_pro<0.05 in TCGA:\n")
print(M[M$n_sig >= 2 & M$same & !is.na(M$p_pro) & M$p_pro < 0.05, c(cols, "comp", "p_pro")][order(-M$n_sig[M$n_sig >= 2 & M$same & !is.na(M$p_pro) & M$p_pro < 0.05]), ], row.names = FALSE, digits = 3)
