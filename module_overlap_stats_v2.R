mod  <- read.csv("results/meta/module_assignments.csv", stringsAsFactors = FALSE)
mod3 <- unique(mod$gene[mod$module == 3])
U_cl <- unique(mod$gene)                                   # genes that entered the module clustering
U_72 <- unique(read.csv("results/meta/robust_tcga_features.csv", stringsAsFactors = FALSE)$gene)  # tested genes
cat("Module 3 genes:", length(mod3), "| clustering universe:", length(U_cl), "| tested-gene universe:", length(U_72), "\n")
cat("Module 3 genes inside the tested-gene universe:", length(intersect(mod3, U_72)), "\n")

mof <- read.csv("ref/moffitt_lists.csv", stringsAsFactors = FALSE)
b <- read.delim("ref/pdacR/data/Bailey_readable_list.tsv", header = FALSE, col.names = c("gene", "class"))
p <- read.delim("ref/pdacR/data/Puleo_readable_list.tsv",  header = FALSE, col.names = c("gene", "class"))
sig <- rbind(mof, data.frame(signature = paste0("Bailey_", b$class), gene = b$gene),
                  data.frame(signature = paste0("Puleo_",  p$class), gene = p$gene))

test <- function(U) do.call(rbind, lapply(split(sig$gene, sig$signature), function(g) {
  g <- intersect(unique(g), U); m <- intersect(mod3, U); k <- length(intersect(g, m))
  data.frame(n_list = length(g), overlap = k,
             p = if (length(g) == 0) NA else phyper(k - 1, length(g), length(U) - length(g), length(m), lower.tail = FALSE))
}))
rc <- test(U_cl); r7 <- test(U_72)
old <- do.call(rbind, lapply(split(sig$gene, sig$signature), function(g) {
  g <- unique(g); k <- length(intersect(g, mod3))
  data.frame(list_size = length(g), overlap_old = k, p_old15000 = phyper(k - 1, length(g), 15000 - length(g), length(mod3), lower.tail = FALSE))
}))
out <- data.frame(signature = rownames(old), list_size = old$list_size, overlap = old$overlap_old,
                  n_in_clustering_set = rc$n_list, p_clustering_bg = rc$p, fdr_clustering_bg = p.adjust(rc$p, "BH"),
                  n_in_720 = r7$n_list, p_720_bg = r7$p, fdr_720_bg = p.adjust(r7$p, "BH"), p_old15000 = old$p_old15000)
out <- out[order(out$p_clustering_bg), ]
num <- sapply(out, is.numeric); out[num] <- lapply(out[num], signif, 3)
print(head(out[out$overlap > 0, ], 12), row.names = FALSE)
write.csv(out, "results/meta/module3_overlap_stats_refined.csv", row.names = FALSE)
