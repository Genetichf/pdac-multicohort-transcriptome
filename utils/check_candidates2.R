suppressPackageStartupMessages(library(GEOquery))
options(timeout = 900)
cat("==== GSE16515 pairing ====\n")
g <- getGEO("GSE16515", GSEMatrix = TRUE, getGPL = FALSE)[[1]]; p <- pData(g)
p$grp <- ifelse(grepl("Normal", p$title), "Normal", "Tumor")
p$num <- sub("^Pancreatic Sample[[:space:]]*([0-9]*)-.*$", "\\1", p$title)
cat("Samples with no number in title:", sum(p$num == ""), "\n")
tb <- table(p$num[p$num != ""], p$grp[p$num != ""])
cat("Numbered IDs:", nrow(tb), "| IDs with both N and T:", sum(tb[, "Normal"] > 0 & tb[, "Tumor"] > 0), "| normals without a numbered tumor:", sum(tb[, "Normal"] > 0 & tb[, "Tumor"] == 0), "\n")
cat("Duplicate-ID check (max samples per ID and group):", max(tb), "\n")
cat("\n==== GSE71989 RNA vendor vs group ====\n")
q <- pData(getGEO("GSE71989", GSEMatrix = TRUE, getGPL = FALSE)[[1]])
print(table(group = q$characteristics_ch1.2, vendor = q$characteristics_ch1.3))
cat("\n==== GSE41368 vs GSE28735 duplicate check ====\n")
a <- getGEO("GSE41368", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
b <- getGEO("GSE28735", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
ids <- intersect(rownames(exprs(a)), rownames(exprs(b)))
cat("Shared probes:", length(ids), "\n")
cm <- cor(exprs(a)[ids, ], exprs(b)[ids, ])
o <- t(apply(cm, 1, function(r) { s <- sort(r, decreasing = TRUE); c(best = s[1], second = s[2], median = median(r)) }))
res <- data.frame(GSE41368 = pData(a)$title, best_match_in_GSE28735 = pData(b)$title[apply(cm, 1, which.max)], round(o, 4))
print(res, row.names = FALSE)
cat("\n==== GSE62165 retry ====\n")
r <- try(getGEO("GSE62165", GSEMatrix = TRUE, getGPL = FALSE), silent = TRUE)
if (inherits(r, "try-error")) cat("Error:", conditionMessage(attr(r, "condition")), "\n") else {
  e <- r[[1]]; pp <- pData(e); cat("Platform:", annotation(e), "| samples:", ncol(e), "| range:", paste(round(range(exprs(e), na.rm = TRUE), 1), collapse = " to "), "\n")
  for (cn in grep("^(title|source_name_ch1|characteristics_ch1)", names(pp), value = TRUE)) { cat("[", cn, "]\n", sep = ""); print(head(sort(table(pp[[cn]]), decreasing = TRUE), 6)) }
}
