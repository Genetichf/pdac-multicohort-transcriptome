gl <- readRDS("ref/pdacR/data/gene_lists.rds")
cat("Lists found:\n"); print(names(gl))

# Moffitt basal/classical (25-gene lists)
mof <- do.call(rbind, lapply(grep("Moffitt", names(gl), value = TRUE), function(n)
  data.frame(signature = n, gene = as.character(gl[[n]]))))
write.csv(mof, "ref/moffitt_lists.csv", row.names = FALSE)
print(table(mof$signature))

b <- read.delim("ref/pdacR/data/Bailey_readable_list.tsv", header=FALSE, col.names=c("gene","class")); print(table(b$class))
p <- read.delim("ref/pdacR/data/Puleo_readable_list.tsv", header=FALSE, col.names=c("gene","class"));  print(table(p$class))

# Overlap with our 56-gene module
m <- read.csv("results/meta/module_assignments.csv")
mod3 <- m$gene[m$module == 3]
allsig <- rbind(
  mof,
  data.frame(signature = paste0("Bailey_", b$class), gene = b$gene),
  data.frame(signature = paste0("Puleo_",  p$class), gene = p$gene))
ov <- aggregate(gene ~ signature, allsig, function(g) {
  h <- intersect(g, mod3); paste0(length(h), ": ", paste(h, collapse = ","))})
print(ov)
write.csv(ov, "results/meta/module3_signature_overlap.csv", row.names = FALSE)
