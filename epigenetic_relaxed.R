suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi); library(rentrez) })
A <- read.csv("results/meta/meta_three_datasets_all_genes.csv", stringsAsFactors = FALSE)
P <- read.csv("results/meta/gene_prioritization.csv", stringsAsFactors = FALSE)[, c("gene", "comp", "surv_n")]
rel <- A[!is.na(A$fdr) & A$fdr < 0.05 & abs(A$lfc) >= 0.5 & A$agree == 1 & A$k >= 2, ]
go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = c("GO:0006325", "GO:0016569", "GO:0006306", "GO:0080111", "GO:0042393", "GO:0042054", "GO:0004407", "GO:0004402", "GO:0032452"), keytype = "GOALL", columns = "SYMBOL"))
epi <- unique(na.omit(go$SYMBOL))
tested <- A$gene[!is.na(A$p)]
m <- sum(tested %in% epi); x <- sum(rel$gene %in% epi)
cat("Relaxed DEGs:", nrow(rel), "| epigenetic-annotated among them:", x, "(", round(100 * x / nrow(rel), 1), "% vs", round(100 * m / length(tested), 1), "% of all tested genes ) hypergeometric p =", signif(phyper(x - 1, m, length(tested) - m, nrow(rel), lower.tail = FALSE), 3), "\n")
E <- rel[rel$gene %in% epi, c("gene", "lfc", "I2", "k", "fdr")]
cnt <- function(q) { for (i in 1:3) { r <- try(entrez_search(db = "pubmed", term = q, retmax = 0)$count, silent = TRUE); if (!inherits(r, "try-error")) return(r); Sys.sleep(1) }; NA }
pd <- "(pancreatic cancer[tiab] OR pancreatic ductal adenocarcinoma[tiab] OR PDAC[tiab])"
E$pdac <- NA
for (i in seq_len(nrow(E))) { E$pdac[i] <- cnt(sprintf("%s[tiab] AND %s", E$gene[i], pd)); Sys.sleep(0.4) }
E <- merge(E, P, by = "gene", all.x = TRUE)
E <- E[order(-abs(E$lfc)), ]
write.csv(E, "results/meta/epigenetic_regulators_relaxed.csv", row.names = FALSE)
cat("\nDirection:\n"); print(table(ifelse(E$lfc > 0, "up", "down")))
cat("\nTop 40 epigenetic-annotated genes by effect size:\n"); print(head(E, 40), row.names = FALSE, digits = 3)
cat("\nWith 5 or fewer PDAC papers:\n"); print(E[!is.na(E$pdac) & E$pdac <= 5, ], row.names = FALSE, digits = 3)
