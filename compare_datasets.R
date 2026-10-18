a_id <- "GSE28735"; b_id <- "GSE62452"
rd <- function(id) read.csv(sprintf("results/%s/%s_limma_Tumor_vs_Normal.csv", id, id), stringsAsFactors = FALSE)
A <- rd(a_id); B <- rd(b_id)
j <- merge(A[, c("gene","logFC","adj.P.Val")], B[, c("gene","logFC","adj.P.Val")], by = "gene", suffixes = c("_A","_B"))
cat("Genes tested in both:", nrow(j), "\n")
cat("logFC correlation (Pearson):", round(cor(j$logFC_A, j$logFC_B), 3), "| Spearman:", round(cor(j$logFC_A, j$logFC_B, method = "spearman"), 3), "\n")
sigA <- j$adj.P.Val_A < 0.05 & abs(j$logFC_A) >= 1
sigB <- j$adj.P.Val_B < 0.05 & abs(j$logFC_B) >= 1
both <- sigA & sigB
same <- both & sign(j$logFC_A) == sign(j$logFC_B)
cat("DEGs in", a_id, ":", sum(sigA), "| in", b_id, ":", sum(sigB), "\n")
cat("Shared DEGs:", sum(both), "| same direction:", sum(same), "\n")
cat("Share of", b_id, "DEGs also found in", a_id, ":", round(100 * sum(both) / sum(sigB), 1), "%\n")
p <- phyper(sum(both) - 1, sum(sigA), nrow(j) - sum(sigA), sum(sigB), lower.tail = FALSE)
cat("Overlap enrichment (hypergeometric) p =", signif(p, 3), "\n")
out <- j[same, ]; out$mean_logFC <- (out$logFC_A + out$logFC_B) / 2
out <- out[order(-out$mean_logFC), ]
dir.create("results/replicated", showWarnings = FALSE)
write.csv(out, "results/replicated/shared_DEGs_GSE28735_GSE62452.csv", row.names = FALSE)
cat("\nTop 10 shared UP:\n"); print(head(out[, c("gene","logFC_A","logFC_B")], 10), row.names = FALSE)
cat("\nTop 10 shared DOWN:\n"); print(head(out[order(out$mean_logFC), c("gene","logFC_A","logFC_B")], 10), row.names = FALSE)
