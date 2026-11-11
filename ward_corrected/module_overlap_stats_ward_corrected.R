mod  <- read.csv("results/meta/module_assignments_ward_corrected.csv", stringsAsFactors = FALSE)
mod1 <- unique(mod$gene[mod$module == 1])
U_cl <- unique(mod$gene)                                   # genes that entered the module clustering
U_72 <- unique(read.csv("results/meta/robust_tcga_features.csv", stringsAsFactors = FALSE)$gene)  # tested genes
cat("Corrected Module 1 genes:", length(mod1), "| clustering universe:", length(U_cl), "| tested-gene universe:", length(U_72), "\n")
cat("Corrected Module 1 genes inside the tested-gene universe:", length(intersect(mod1, U_72)), "\n")

mof <- read.csv("ref/moffitt_lists.csv", stringsAsFactors = FALSE)
b <- read.delim("ref/pdacR/data/Bailey_readable_list.tsv", header = FALSE, col.names = c("gene", "class"))
p <- read.delim("ref/pdacR/data/Puleo_readable_list.tsv",  header = FALSE, col.names = c("gene", "class"))
sig <- rbind(mof, data.frame(signature = paste0("Bailey_", b$class), gene = b$gene),
                  data.frame(signature = paste0("Puleo_",  p$class), gene = p$gene))
sigs <- split(sig$gene, sig$signature)

test <- function(U) do.call(rbind, lapply(sigs, function(g) {
  g <- intersect(unique(g), U); m <- intersect(mod1, U); k <- length(intersect(g, m))
  data.frame(n_list = length(g), overlap = k,
             p = if (length(g) == 0) NA else phyper(k - 1, length(g), length(U) - length(g), length(m), lower.tail = FALSE))
}))

# BH across unique effective gene sets: lists with identical gene content in the universe are counted once
adj_unique <- function(r, U) {
  key <- sapply(sigs, function(g) paste(sort(intersect(unique(g), U)), collapse = "|"))
  ok <- !is.na(r$p); k <- key[ok]; first <- !duplicated(k)
  q <- p.adjust(r$p[ok][first], "BH")
  fdr <- rep(NA_real_, nrow(r)); fdr[ok] <- q[match(k, k[first])]
  cat("Unique effective gene sets tested:", sum(first), "\n")
  fdr
}
rc <- test(U_cl); r7 <- test(U_72)
rc$fdr <- adj_unique(rc, U_cl); r7$fdr <- adj_unique(r7, U_72)
old <- do.call(rbind, lapply(sigs, function(g) {
  g <- unique(g); k <- length(intersect(g, mod1))
  data.frame(list_size = length(g), overlap_old = k, p_old15000 = phyper(k - 1, length(g), 15000 - length(g), length(mod1), lower.tail = FALSE))
}))
out <- data.frame(signature = rownames(old), list_size = old$list_size, overlap = old$overlap_old,
                  n_in_clustering_set = rc$n_list, p_clustering_bg = rc$p, fdr_clustering_bg = rc$fdr,
                  n_in_720 = r7$n_list, p_720_bg = r7$p, fdr_720_bg = r7$fdr, p_old15000 = old$p_old15000)
out <- out[order(out$p_clustering_bg), ]
num <- sapply(out, is.numeric); out[num] <- lapply(out[num], signif, 3)
print(head(out[out$overlap > 0, ], 12), row.names = FALSE)
write.csv(out, "results/meta/module1_ward_corrected_overlap_stats_refined.csv", row.names = FALSE)
