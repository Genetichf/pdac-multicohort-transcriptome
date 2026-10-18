suppressPackageStartupMessages(library(GEOquery))
ids <- c("GSE71989", "GSE16515", "GSE62165", "GSE41368", "GSE32676")
for (id in ids) {
  cat("\n====", id, "====\n")
  g <- try(getGEO(id, GSEMatrix = TRUE, getGPL = FALSE), silent = TRUE)
  if (inherits(g, "try-error")) { cat("Download failed\n"); next }
  for (i in seq_along(g)) {
    e <- g[[i]]; p <- pData(e)
    cat("Platform:", annotation(e), "| samples:", ncol(e), "| expression range:", paste(round(range(exprs(e), na.rm = TRUE), 1), collapse = " to "), "\n")
    write.csv(p, sprintf("results/candidates/%s_%s_metadata.csv", id, annotation(e)), row.names = FALSE)
    cols <- grep("^(title|source_name_ch1|characteristics_ch1)", names(p), value = TRUE)
    for (cn in cols) {
      tb <- sort(table(p[[cn]]), decreasing = TRUE)
      cat(" [", cn, "] ", length(tb), " distinct values\n", sep = "")
      print(head(tb, 6))
    }
  }
}
