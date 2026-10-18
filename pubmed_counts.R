if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez", repos = "https://cloud.r-project.org")
library(rentrez)
S <- read.csv("results/meta/meta_three_datasets_DEGs.csv", stringsAsFactors = FALSE)
cnt <- function(q) {
  for (i in 1:3) { r <- try(entrez_search(db = "pubmed", term = q, retmax = 0)$count, silent = TRUE); if (!inherits(r, "try-error")) return(r); Sys.sleep(1) }
  NA
}
pd <- "(pancreatic cancer[tiab] OR pancreatic ductal adenocarcinoma[tiab] OR PDAC[tiab])"
S$pdac <- NA; S$pdac_prog <- NA
for (i in seq_len(nrow(S))) {
  S$pdac[i] <- cnt(sprintf("%s[tiab] AND %s", S$gene[i], pd)); Sys.sleep(0.4)
  S$pdac_prog[i] <- cnt(sprintf("%s[tiab] AND %s AND (prognosis[tiab] OR survival[tiab])", S$gene[i], pd)); Sys.sleep(0.4)
  if (i %% 50 == 0) cat("done", i, "of", nrow(S), "\n")
}
S <- S[order(S$pdac), ]
write.csv(S, "results/meta/meta_DEGs_pubmed_counts.csv", row.names = FALSE)
cat("\nGenes with 0 PDAC papers:", sum(S$pdac == 0, na.rm = TRUE), "| 1-5 papers:", sum(S$pdac >= 1 & S$pdac <= 5, na.rm = TRUE), "\n")
cat("\nLeast-studied genes (lowest PDAC counts):\n")
print(head(S[, c("gene", "lfc", "I2", "k", "pdac", "pdac_prog")], 30), row.names = FALSE)
