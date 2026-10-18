if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("recount3", quietly = TRUE)) BiocManager::install("recount3", update = FALSE, ask = FALSE)
suppressPackageStartupMessages({ library(recount3); library(SummarizedExperiment) })
hp <- available_projects()
sel <- subset(hp, project_type == "data_sources" & ((file_source == "tcga" & project == "PAAD") | (file_source == "gtex" & project == "PANCREAS")))
print(sel[, c("project", "file_source", "project_type", "n_samples")])
dir.create("results/rnaseq", recursive = TRUE, showWarnings = FALSE)
for (i in seq_len(nrow(sel))) {
  rse <- create_rse(sel[i, ])
  saveRDS(rse, sprintf("results/rnaseq/rse_%s_%s.rds", sel$file_source[i], sel$project[i]))
  cat("\n==", sel$file_source[i], sel$project[i], ":", ncol(rse), "samples,", nrow(rse), "genes\n")
  cd <- as.data.frame(colData(rse))
  keep <- grep("type|histolog|diagnos|purity|tissue|smts|smtsd", names(cd), value = TRUE, ignore.case = TRUE)
  cat("Candidate annotation columns:\n"); print(head(keep, 40))
  for (cn in keep) {
    u <- unique(cd[[cn]])
    if (length(u) <= 12) { cat("\n--", cn, "\n"); print(table(cd[[cn]], useNA = "ifany")) }
  }
}
