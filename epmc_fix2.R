suppressMessages(library(jsonlite))
run <- function(label, query, n=5){ Sys.sleep(0.3)
  u <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=", URLencode(query, reserved=TRUE),
              "&format=json&resultType=lite&pageSize=", n)
  r <- tryCatch(fromJSON(u), error=function(e) NULL)
  cat("\n##", label, "- hits:", if (is.null(r)) NA else r$hitCount, "\n")
  if (!is.null(r) && r$hitCount > 0) print(substr(r$resultList$result$title, 1, 130)) }
pdac <- '(pancreatic OR pancreas)'
# full text + aliases for the low-count candidates
run("FAR2 title/abs + alias",  paste0('(TITLE_ABS:"FAR2" OR TITLE_ABS:"fatty acyl-CoA reductase 2" OR TITLE_ABS:"MLSTD2") AND TITLE_ABS:', pdac))
run("FAR2 full text",          paste0('BODY:"FAR2" AND BODY:"pancreatic cancer"'))
run("SLC9B2 + alias",          paste0('(TITLE_ABS:"SLC9B2" OR TITLE_ABS:"NHA2" OR TITLE_ABS:"NHEDC2") AND TITLE_ABS:', pdac))
run("SLC9B2 full text",        paste0('(BODY:"SLC9B2" OR BODY:"NHA2") AND BODY:"pancreatic cancer"'))
run("ELOVL6 titles",           paste0('TITLE_ABS:"ELOVL6" AND TITLE_ABS:', pdac), 12)
# disambiguated versions of the ambiguous symbols
run("NHS gene",   paste0('(TITLE_ABS:"NHS actin remodeling regulator" OR TITLE_ABS:"Nance-Horan") AND TITLE_ABS:', pdac))
run("PI3 elafin", paste0('(TITLE_ABS:"elafin" OR TITLE_ABS:"peptidase inhibitor 3" OR TITLE_ABS:"skin-derived antileukoproteinase") AND TITLE_ABS:', pdac))
run("F5 factor V",paste0('(TITLE_ABS:"coagulation factor V" OR TITLE_ABS:"factor V" ) AND TITLE_ABS:"pancreatic cancer" AND TITLE_ABS:"F5"'))
