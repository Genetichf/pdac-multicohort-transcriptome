core_ids <- c("GSE28735", "GSE62452", "GSE15471", "GSE16515", "GSE62165")
ext_ids  <- c(core_ids, "GSE32676")
rd <- function(id) {
  x <- read.csv(sprintf("results/%s/%s_limma_Tumor_vs_Normal.csv", id, id), stringsAsFactors = FALSE)
  x$se <- abs(x$logFC / x$t); x <- x[, c("gene", "logFC", "se")]
  names(x)[2:3] <- paste0(c("lfc_", "se_"), id); x
}
meta <- function(ids) {
  M <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), lapply(ids, rd))
  lf <- paste0("lfc_", ids); sf <- paste0("se_", ids)
  res <- t(apply(M[, -1], 1, function(r) {
    y <- r[lf]; s <- r[sf]; ok <- is.finite(y) & is.finite(s) & s > 0; k <- sum(ok)
    if (k < 2) return(c(k = k, lfc = NA, se = NA, z = NA, p = NA, I2 = NA, agree = NA, frac_same = NA))
    y <- y[ok]; s <- s[ok]; w <- 1 / s^2
    mu_f <- sum(w * y) / sum(w); Q <- sum(w * (y - mu_f)^2)
    tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
    ws <- 1 / (s^2 + tau2); mu <- sum(ws * y) / sum(ws); se <- sqrt(1 / sum(ws)); z <- mu / se
    c(k = k, lfc = mu, se = se, z = z, p = 2 * pnorm(-abs(z)), I2 = if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0,
      agree = as.numeric(all(sign(y) == sign(mu))), frac_same = mean(sign(y) == sign(mu)))
  }))
  R <- data.frame(gene = M$gene, M[, -1], res, row.names = NULL); R$fdr <- p.adjust(R$p, "BH"); R
}
core <- meta(core_ids); ext <- meta(ext_ids)
write.csv(core[order(core$fdr), ], "results/meta/meta_core5_all_genes.csv", row.names = FALSE)
write.csv(ext[order(ext$fdr), ], "results/meta/meta_ext6_all_genes.csv", row.names = FALSE)
deg <- function(R) R[!is.na(R$fdr) & R$fdr < 0.05 & abs(R$lfc) >= 1 & R$agree == 1, ]
dc <- deg(core); de <- deg(ext)
old <- read.csv("results/meta/meta_three_datasets_DEGs.csv", stringsAsFactors = FALSE)
rob <- intersect(dc$gene, de$gene)
cat("Core (5 cohorts): genes tested", sum(!is.na(core$p)), "| in all 5:", sum(core$k == 5, na.rm = TRUE), "| DEGs:", nrow(dc), "(up", sum(dc$lfc > 0), "/ down", sum(dc$lfc < 0), ") | median I2 of DEGs:", round(median(dc$I2), 1), "\n")
cat("Extended (6 cohorts): DEGs:", nrow(de), "(up", sum(de$lfc > 0), "/ down", sum(de$lfc < 0), ")\n")
cat("Robust (DEG in both tiers):", length(rob), "| overlap with earlier 402-gene list:", length(intersect(rob, old$gene)), "| earlier genes lost from core:", length(setdiff(old$gene, dc$gene)), "\n")
cat("Core genes with FDR<0.05 & |lfc|>=1:", sum(core$fdr < 0.05 & abs(core$lfc) >= 1, na.rm = TRUE), "| of which same sign in every cohort:", nrow(dc), "\n")
R <- core[core$gene %in% rob, c("gene", "lfc", "I2", "k", "fdr")]; names(R)[2:5] <- paste0(names(R)[2:5], "_core")
R <- merge(R, ext[, c("gene", "lfc", "I2", "fdr")], by = "gene", suffixes = c("", "_ext"))
write.csv(R[order(-abs(R$lfc_core)), ], "results/meta/meta_robust_DEGs.csv", row.names = FALSE)
write.csv(dc[order(-dc$lfc), ], "results/meta/meta_core5_DEGs.csv", row.names = FALSE)
cat("\nSpearman correlation of logFC between cohorts (genes in both):\n")
print(round(cor(ext[, paste0("lfc_", ext_ids)], method = "spearman", use = "pairwise.complete.obs"), 2))
cat("\nTop 10 up / down in the robust set:\n")
print(head(R[order(-R$lfc_core), c("gene", "lfc_core", "I2_core", "k_core")], 10), row.names = FALSE)
print(head(R[order(R$lfc_core), c("gene", "lfc_core", "I2_core", "k_core")], 10), row.names = FALSE)
