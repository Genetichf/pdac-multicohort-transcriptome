suppressMessages(library(GEOquery))
gse <- getGEO("GSE155698", GSEMatrix = FALSE)
gsms <- GSMList(gse)
for (g in names(gsms)) {
  m <- Meta(gsms[[g]])
  cat(g, m$title, paste(m$characteristics_ch1, collapse=";"), "\n")
}
