suppressMessages(library(jsonlite))
x <- read.csv("results/meta/module3_novel18_evidence.csv")
q <- function(g) sprintf('TITLE_ABS:"%s" AND (TITLE_ABS:pancreatic OR TITLE_ABS:pancreas)', g)
get <- function(g, n) { Sys.sleep(0.3)
  u <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=", URLencode(q(g), reserved=TRUE),
              "&format=json&resultType=lite&pageSize=", n)
  tryCatch(fromJSON(u), error=function(e) NULL) }
x$epmc_pdac_papers <- sapply(x$gene, function(g){ r <- get(g,1); if (is.null(r)) NA else r$hitCount })
write.csv(x, "results/meta/module3_novel18_evidence.csv", row.names=FALSE)
print(x[, c("gene","epmc_pdac_papers")], row.names=FALSE)
for (g in x$gene[!is.na(x$epmc_pdac_papers) & x$epmc_pdac_papers > 0 & x$epmc_pdac_papers <= 4]) {
  r <- get(g, 5); cat("\n##", g, "\n"); print(substr(r$resultList$result$title, 1, 130)) }
