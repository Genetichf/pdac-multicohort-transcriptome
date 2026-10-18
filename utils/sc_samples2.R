suppressMessages(library(GEOquery))
for (id in c("GSE212966","GSE154778")) {
  cat("\n=====", id, "=====\n")
  gs <- GSMList(getGEO(id, GSEMatrix=FALSE))
  for (n in names(gs)) {
    m <- Meta(gs[[n]])
    f <- if (is.null(m$supplementary_file)) "none" else paste(basename(m$supplementary_file), collapse=",")
    cat(n, "|", m$title, "|", paste(m$characteristics_ch1, collapse="; "), "| files:", f, "\n")
  }
}
