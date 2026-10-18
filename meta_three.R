ids <- c("GSE28735", "GSE62452", "GSE15471")
rd <- function(id) {
  x <- read.csv(sprintf("results/%s/%s_limma_Tumor_vs_Normal.csv", id, id), stringsAsFactors = FALSE)
  x$se <- abs(x$logFC / x$t)
  x[, c("gene", "logFC", "se")]
}
L <- lapply(ids, rd)
for (i in seq_along(L)) names(L[[i]])[2:3] <- paste0(c("lfc_", "se_"), ids[i])
M <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), L)
lf <- paste0("lfc_", ids); sf <- paste0("se_", ids)
res <- t(apply(M[, -1], 1, function(r) {
  y <- r[lf]; s <- r[sf]
  ok <- is.finite(y) & is.finite(s) & s > 0
  k <- sum(ok)
  if (k < 2) return(c(k = k, lfc = NA, se = NA, z = NA, p = NA, I2 = NA, agree = NA))
  y <- y[ok]; s <- s[ok]; w <- 1 / s^2
  mu_f <- sum(w * y) / sum(w); Q <- sum(w * (y - mu_f)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (s^2 + tau2); mu <- sum(ws * y) / sum(ws); se <- sqrt(1 / sum(ws))
  z <- mu / se
  I2 <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
  c(k = k, lfc = mu, se = se, z = z, p = 2 * pnorm(-abs(z)), I2 = I2, agree = as.numeric(all(sign(y) == sign(mu))))
}))
R <- data.frame(gene = M$gene, M[, -1], res, row.names = NULL)
R$fdr <- p.adjust(R$p, "BH")
dir.create("results/meta", showWarnings = FALSE)
write.csv(R[order(R$fdr), ], "results/meta/meta_three_datasets_all_genes.csv", row.names = FALSE)
S <- R[!is.na(R$fdr) & R$fdr < 0.05 & abs(R$lfc) >= 1 & R$agree == 1, ]
S <- S[order(-S$lfc), ]
write.csv(S, "results/meta/meta_three_datasets_DEGs.csv", row.names = FALSE)
cat("Genes meta-analysed:", sum(!is.na(R$p)), "| present in all 3 datasets:", sum(R$k == 3, na.rm = TRUE), "\n")
cat("Meta DEGs (FDR<0.05, |pooled log2FC|>=1, same direction in every dataset):", nrow(S),
    "| up:", sum(S$lfc > 0), "| down:", sum(S$lfc < 0), "\n")
cat("Of these, with high heterogeneity (I2 > 75%):", sum(S$I2 > 75), "\n")
core <- read.csv("results/replicated/shared_DEGs_3datasets.csv")$gene
cat("Overlap with the 235-gene core list:", length(intersect(S$gene, core)), "\n")
cat("\nTop 10 UP:\n"); print(head(S[, c("gene", "lfc", "fdr", "I2", "k")], 10), row.names = FALSE)
cat("\nTop 10 DOWN:\n"); print(head(S[order(S$lfc), c("gene", "lfc", "fdr", "I2", "k")], 10), row.names = FALSE)
