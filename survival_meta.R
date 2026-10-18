library(survival)
S <- read.csv("results/meta/meta_three_datasets_DEGs.csv", stringsAsFactors = FALSE)
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
res <- list()
for (id in c("GSE28735", "GSE62452")) {
  ex <- read.csv(sprintf("results/%s/expression_gene_level.csv", id), row.names = 1, check.names = FALSE)
  an <- read.csv(sprintf("results/%s/sample_annotation_used.csv", id), stringsAsFactors = FALSE)
  sv <- get_surv(id)
  sv <- sv[sv$gsm %in% an$gsm[an$tissue == "Tumor"] & is.finite(sv$time) & is.finite(sv$event) & sv$time > 0, ]
  cat(id, ": tumors with survival data =", nrow(sv), "| events coded 1 =", sum(sv$event == 1), "| coded 0 =", sum(sv$event == 0), "\n")
  g <- intersect(S$gene, rownames(ex))
  r <- t(sapply(g, function(gn) {
    z <- as.numeric(scale(as.numeric(ex[gn, sv$gsm])))
    f <- try(coxph(Surv(sv$time, sv$event) ~ z), silent = TRUE)
    if (inherits(f, "try-error")) return(c(NA, NA))
    s <- summary(f)$coefficients
    c(s[1, "coef"], s[1, "Pr(>|z|)"])
  }))
  res[[id]] <- data.frame(gene = g, logHR = r[, 1], p = r[, 2])
  res[[id]]$fdr <- p.adjust(res[[id]]$p, "BH")
  cat(id, ": genes tested =", length(g), "| nominal p<0.05 =", sum(res[[id]]$p < 0.05, na.rm = TRUE),
      "(chance expectation ~", round(0.05 * length(g)), ") | FDR<0.1 =", sum(res[[id]]$fdr < 0.1, na.rm = TRUE), "\n")
}
M <- merge(res[[1]], res[[2]], by = "gene", suffixes = c("_A", "_B"))
M$HR_A <- exp(M$logHR_A); M$HR_B <- exp(M$logHR_B)
M$both <- M$p_A < 0.05 & M$p_B < 0.05 & sign(M$logHR_A) == sign(M$logHR_B)
dir.create("results/survival", showWarnings = FALSE)
write.csv(M[order(pmax(M$p_A, M$p_B)), ], "results/survival/cox_GSE28735_GSE62452.csv", row.names = FALSE)
cat("\nGenes nominally significant in BOTH cohorts, same direction:", sum(M$both, na.rm = TRUE), "\n")
print(head(M[order(pmax(M$p_A, M$p_B)), c("gene", "HR_A", "p_A", "HR_B", "p_B")], 15), row.names = FALSE)
