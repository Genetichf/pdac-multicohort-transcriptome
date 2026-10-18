suppressPackageStartupMessages({ library(rentrez); library(org.Hs.eg.db); library(AnnotationDbi) })
S <- read.csv("results/meta/meta_DEGs_pubmed_counts.csv", stringsAsFactors = FALSE)
low <- S$gene[!is.na(S$pdac) & S$pdac <= 5]
eg <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = low, keytype = "ALIAS", columns = c("ENTREZID", "SYMBOL")))
al <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(na.omit(eg$ENTREZID)), keytype = "ENTREZID", columns = c("SYMBOL", "ALIAS")))
cnt <- function(q) {
  for (i in 1:3) { r <- try(entrez_search(db = "pubmed", term = q, retmax = 0)$count, silent = TRUE); if (!inherits(r, "try-error")) return(r); Sys.sleep(1) }
  NA
}
pd <- "(pancreatic cancer[tiab] OR pancreatic ductal adenocarcinoma[tiab] OR PDAC[tiab])"
S$pdac_alias <- NA
S$aliases <- NA
for (g in low) {
  ids <- unique(na.omit(eg$ENTREZID[eg$ALIAS == g]))
  a <- unique(c(g, al$SYMBOL[al$ENTREZID %in% ids], al$ALIAS[al$ENTREZID %in% ids]))
  a <- a[!is.na(a) & nchar(a) >= 4 & grepl("[A-Za-z]", a) & !grepl("^[0-9]+$", a)]
  a <- head(unique(c(g, a)), 8)
  i <- which(S$gene == g)
  S$aliases[i] <- paste(a, collapse = ";")
  S$pdac_alias[i] <- cnt(sprintf("(%s) AND %s", paste(sprintf("\"%s\"[tiab]", a), collapse = " OR "), pd))
  Sys.sleep(0.4)
}
S$pdac_best <- pmax(S$pdac, S$pdac_alias, na.rm = TRUE)
write.csv(S, "results/meta/meta_DEGs_pubmed_counts_alias.csv", row.names = FALSE)
cat("Genes with 0 PDAC papers after alias expansion:", sum(S$pdac_best == 0, na.rm = TRUE), "| 1-5:", sum(S$pdac_best >= 1 & S$pdac_best <= 5, na.rm = TRUE), "\n")
T <- S[!is.na(S$pdac_best) & S$pdac_best <= 2 & S$I2 < 50 & S$k >= 2, ]
cat("\nConsistent (I2<50) genes with <=2 PDAC papers, by effect size:\n")
print(head(T[order(-abs(T$lfc)), c("gene", "lfc", "I2", "k", "pdac", "pdac_best", "aliases")], 30), row.names = FALSE)
