f <- read.csv("results/meta/robust_funnel_table.csv")
cat("Funnel table columns:\n"); print(colnames(f))
if ("p" %in% colnames(f)) {
  ok <- !is.na(f$p); f$fdr_tcga <- NA
  f$fdr_tcga[ok] <- p.adjust(f$p[ok], "BH")
  cat("\nGenes with TCGA Cox p:", sum(ok), "| p<0.05:", sum(f$p < 0.05, na.rm=TRUE),
      "| BH FDR<0.05:", sum(f$fdr_tcga < 0.05, na.rm=TRUE), "\n")
  o <- order(f$p); cols <- intersect(c("gene","HR","p","fdr_tcga","HR_adj","p_adj"), colnames(f))
  print(head(f[o, cols], 15), row.names=FALSE)
}
a <- read.csv("results/survival/cox_arrays_robust.csv")
cat("\nArray Cox file columns:\n"); print(colnames(a)); print(head(a, 3))
