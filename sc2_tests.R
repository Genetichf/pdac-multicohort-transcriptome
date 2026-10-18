suppressMessages(library(Seurat))
obj <- readRDS("sc2/GSE212966/seurat_obj.rds")
ma <- read.csv("results/meta/module_assignments.csv")
mucins <- c("MUC5AC","MUC13","MUC4","MUC1","MUC6","MUC2")
g3 <- setdiff(intersect(ma$gene[as.character(ma$module) == "3"], rownames(obj)), mucins)
cat("Module 3 genes found:", length(g3), "of 56\n")
obj <- AddModuleScore(obj, features=list(g3), name="gast", seed=1)
key <- intersect(c("TFF1","CTSE","VSIG2","CAPN8","ARNTL2","VSIG1","ST6GALNAC1",
                   "EPS8L3","TMEM45B","TFF2","ACACB", mucins), rownames(obj))
X <- GetAssayData(obj, layer="data")[key, , drop=FALSE]
grp <- paste(obj$celltype, obj$group, sep="|")
pct <- sapply(split(seq_len(ncol(obj)), grp), function(i) 100*rowMeans(X[, i, drop=FALSE] > 0))
oth <- pct[, !grepl("^Ductal", colnames(pct)), drop=FALSE]
cat("\n% of cells expressing: ductal PDAC vs ductal adjacent vs best other cell type\n")
print(data.frame(PDAC_ductal=round(pct[, "Ductal|PDAC"], 1),
                 Adj_ductal=round(pct[, "Ductal|Adjacent"], 1),
                 max_other=round(apply(oth, 1, max), 1),
                 in_group=colnames(oth)[apply(oth, 1, which.max)]))
d <- subset(obj, subset = celltype == "Ductal")
sub6 <- intersect(c("TFF1","CTSE","VSIG2","ARNTL2","MUC5AC","MUC13"), rownames(d))
Xd <- as.matrix(GetAssayData(d, layer="data")[sub6, , drop=FALSE])
ps <- do.call(rbind, lapply(split(seq_len(ncol(d)), d$sample), function(i)
  data.frame(group=d$group[i[1]], ductal=length(i), module=round(mean(d$gast1[i]), 3),
             t(round(100*rowMeans(Xd[, i, drop=FALSE] > 0), 1)))))
cat("\nPer sample (ductal cells):\n"); print(ps[order(ps$group, -ps$module), ])
t <- subset(d, subset = group == "PDAC"); mu <- intersect(mucins, rownames(t))
Xt <- as.matrix(GetAssayData(t, layer="data")[mu, , drop=FALSE])
dep <- cbind(log(t$nCount_RNA), log(t$nFeature_RNA)); s <- t$gast1
part <- function(x, i) { if (length(i) < 30) return(NA)
  D <- apply(dep[i, ], 2, rank)
  suppressWarnings(cor(resid(lm(rank(x[i]) ~ D)), resid(lm(rank(s[i]) ~ D)))) }
res <- sapply(mu, function(m) sapply(split(seq_along(s), t$sample),
                                      function(i) part(Xt[m, ], i)))
cat("\nWithin-tumor partial rho (module score vs mucin, depth removed):\n"); print(round(res, 2))
print(rbind(median=round(apply(res, 2, median, na.rm=TRUE), 2), n_pos=colSums(res > 0.1, na.rm=TRUE)))
