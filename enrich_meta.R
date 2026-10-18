if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
for (p in c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot) })
R <- read.csv("results/meta/meta_three_datasets_all_genes.csv", stringsAsFactors = FALSE)
S <- read.csv("results/meta/meta_three_datasets_DEGs.csv", stringsAsFactors = FALSE)
universe <- R$gene[!is.na(R$p)]
dir.create("results/enrichment", showWarnings = FALSE)
run <- function(genes, tag, show = TRUE) {
  e <- enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500)
  if (is.null(e) || nrow(as.data.frame(e)) == 0) { cat("\n", tag, ": no significant terms\n"); return(invisible(NULL)) }
  e <- simplify(e, cutoff = 0.7, by = "p.adjust", select_fun = min)
  d <- as.data.frame(e)
  write.csv(d, sprintf("results/enrichment/GO_BP_%s.csv", tag), row.names = FALSE)
  pdf(sprintf("results/enrichment/GO_BP_%s.pdf", tag), width = 9, height = 7)
  print(dotplot(e, showCategory = 15, title = paste("GO BP -", tag)))
  dev.off()
  cat("\n", tag, ":", length(genes), "genes,", nrow(d), "terms after simplification\n")
  if (show) print(head(d[, c("Description", "Count", "GeneRatio", "p.adjust")], 10), row.names = FALSE)
  invisible(e)
}
up <- S$gene[S$lfc > 0]; dn <- S$gene[S$lfc < 0]
run(up, "up"); run(dn, "down")
run(S$gene[S$lfc > 0 & S$I2 <= 50], "up_consistent", show = FALSE)
run(S$gene[S$lfc < 0 & S$I2 <= 50], "down_consistent", show = FALSE)
