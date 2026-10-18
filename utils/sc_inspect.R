suppressMessages(library(GEOquery))
gse <- getGEO("GSE155698", GSEMatrix = FALSE)
gsms <- GSMList(gse)
cat("Samples:", length(gsms), "\n\n")
for (g in names(gsms)) {
  m <- Meta(gsms[[g]])
  cat(g, "|", m$title, "|",
      paste(m$characteristics_ch1, collapse = "; "), "\n")
  cat("   files:", paste(basename(m$supplementary_file), collapse = ", "), "\n")
}
cat("\nSeries-level files:\n")
print(Meta(gse)$supplementary_file)
