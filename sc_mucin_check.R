suppressMessages(library(Seurat))
obj <- readRDS("sc/seurat_obj.rds")
ma <- read.csv("results/meta/module_assignments.csv")
mucins <- c("MUC5AC","MUC6","MUC1","MUC13","MUC4","MUC2")
g3 <- intersect(ma$gene[as.character(ma$module) == "3"], rownames(obj))
cat("Module 3 genes:", length(g3), "| mucin genes inside module:",
    paste(intersect(g3, mucins), collapse=",") , "\n")
g3 <- setdiff(g3, mucins)
obj <- AddModuleScore(obj, features=list(g3), name="gast", seed=1)
d <- subset(obj, subset = celltype == "Ductal" & group == "PDAC")
mucins <- intersect(mucins, rownames(d))
ex <- GetAssayData(d, layer="data")[mucins, , drop=FALSE]
s <- d$gast1
rho <- apply(ex, 1, function(x) suppressWarnings(cor(x, s, method="spearman")))
pct <- 100*rowMeans(ex > 0)
q <- cut(s, quantile(s, c(0, .25, .75, 1)), labels=c("low25","mid","high25"), include.lowest=TRUE)
hi <- 100*rowMeans(ex[, q=="high25", drop=FALSE] > 0)
lo <- 100*rowMeans(ex[, q=="low25", drop=FALSE] > 0)
res <- data.frame(gene=mucins, pct_cells=round(pct,1), spearman_rho=round(rho,3),
                  pct_in_high_score=round(hi,1), pct_in_low_score=round(lo,1))
cat("\nPDAC ductal cells:", ncol(d), "\n"); print(res, row.names=FALSE)
ps <- do.call(rbind, lapply(split(seq_len(ncol(d)), d$sample), function(i)
  data.frame(cells=length(i), module=round(mean(s[i]),3),
             t(round(100*rowMeans(ex[, i, drop=FALSE] > 0), 1)))))
cat("\nPer tumor:\n"); print(ps[order(-ps$module), ])
write.csv(res, "results/sc/mucin_vs_module3.csv", row.names=FALSE)
