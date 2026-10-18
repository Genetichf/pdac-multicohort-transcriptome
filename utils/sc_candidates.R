suppressMessages(library(GEOquery))
for (id in c("GSE202051","GSE154778","GSE165399","GSE212966")) {
  cat("\n=====", id, "=====\n")
  g <- tryCatch(getGEO(id, GSEMatrix=FALSE), error=function(e) NULL)
  if (is.null(g)) { cat("could not fetch\n"); next }
  m <- Meta(g)
  cat("Title:", m$title, "\n")
  cat("Summary:", substr(paste(m$summary, collapse=" "), 1, 250), "\n")
  gs <- GSMList(g)
  cat("Samples:", length(gs), "\n")
  txt <- sapply(gs, function(x) tolower(paste(Meta(x)$title,
                  paste(Meta(x)$characteristics_ch1, collapse=" "))))
  for (w in c("normal","adjacent","tumor|tumour|pdac","metast","blood|pbmc","nucle"))
    cat(sprintf("  %-18s %d\n", w, sum(grepl(w, txt))))
  cat("First titles:", paste(head(sapply(gs, function(x) Meta(x)$title), 6), collapse=" | "), "\n")
  cat("Series files:", paste(basename(m$supplementary_file), collapse=", "), "\n")
}
