rd <- function(id) read.csv(sprintf("results/%s/%s_limma_Tumor_vs_Normal.csv", id, id), stringsAsFactors = FALSE)
S <- read.csv("results/replicated/shared_DEGs_GSE28735_GSE62452.csv")
C <- rd("GSE15471")
k <- merge(S[, c("gene","logFC_A","logFC_B")], C[, c("gene","logFC","adj.P.Val")], by = "gene")
k$mean_logFC_AB <- (k$logFC_A + k$logFC_B) / 2
same <- sign(k$logFC) == sign(k$mean_logFC_AB)
cat("Shared genes (GSE28735 & GSE62452):", nrow(S), "| present in GSE15471:", nrow(k), "\n")
cat("Same direction in GSE15471:", sum(same), "/", nrow(k), "\n")
cat("Same direction and FDR<0.05:", sum(same & k$adj.P.Val < 0.05), "/", nrow(k), "\n")
cat("Same direction, FDR<0.05 and |log2FC|>=1:", sum(same & k$adj.P.Val < 0.05 & abs(k$logFC) >= 1), "/", nrow(k), "\n")
for (id in c("GSE28735", "GSE62452")) {
  j <- merge(rd(id)[, c("gene","logFC")], C[, c("gene","logFC")], by = "gene")
  cat(id, "vs GSE15471 logFC Pearson/Spearman:", round(cor(j$logFC.x, j$logFC.y), 3), "/",
      round(cor(j$logFC.x, j$logFC.y, method = "spearman"), 3), "(n =", nrow(j), ")\n")
}
ok <- k[same & k$adj.P.Val < 0.05, ]
write.csv(ok[order(-ok$mean_logFC_AB), ], "results/replicated/shared_DEGs_3datasets.csv", row.names = FALSE)
bad <- k[!(same & k$adj.P.Val < 0.05), ]
cat("\nShared genes NOT supported in GSE15471 (top 15 by effect size):\n")
print(head(bad[order(-abs(bad$mean_logFC_AB)), c("gene","logFC_A","logFC_B","logFC","adj.P.Val")], 15), row.names = FALSE)
