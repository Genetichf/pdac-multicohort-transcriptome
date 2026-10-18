suppressMessages(library(jsonlite))
g <- list(
  VSIG2=c("VSIG2"),
  ST6GALNAC1=c("ST6GALNAC1","\"ST6GalNAc-I\"","SIAT7A","STYN"),
  VSIG1=c("VSIG1","GPA34"),
  CTSE=c("CTSE","\"cathepsin E\""),
  GPX2=c("GPX2","\"glutathione peroxidase 2\""),
  ANXA10=c("ANXA10","\"annexin A10\""),
  TSPAN8=c("TSPAN8","\"tetraspanin 8\"","CO-029"),
  ZG16B=c("ZG16B","PAUF","\"zymogen granule protein 16B\"")
)
pan <- "(pancreatic OR pancreas OR PDAC)"
base <- "https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=lite&pageSize=8&query="
get <- function(q) tryCatch(fromJSON(paste0(base, URLencode(q, reserved=TRUE))),
                            error=function(e) NULL)
for (n in names(g)) {
  a <- paste0("(", paste(g[[n]], collapse=" OR "), ")")
  ta <- get(paste0("TITLE_ABS:", a, " AND TITLE_ABS:", pan)); Sys.sleep(1)
  fu <- get(paste0(a, " AND ", pan)); Sys.sleep(1)
  cat("\n==", n, "| title/abstract hits:", if (is.null(ta)) "ERR" else ta$hitCount,
      "| all-text hits:", if (is.null(fu)) "ERR" else fu$hitCount, "\n")
  r <- if (!is.null(ta)) ta$resultList$result else NULL
  if (!is.null(r) && nrow(r) > 0)
    for (i in seq_len(min(5, nrow(r))))
      cat("  ", r$pubYear[i], r$pmid[i], substr(r$title[i], 1, 110), "\n")
}
