mod3 <- read.csv("results/meta/module_assignments.csv"); mod3 <- mod3$gene[mod3$module==3]
mof <- read.csv("ref/moffitt_lists.csv")
b <- read.delim("ref/pdacR/data/Bailey_readable_list.tsv", header=FALSE, col.names=c("gene","class"))
p <- read.delim("ref/pdacR/data/Puleo_readable_list.tsv", header=FALSE, col.names=c("gene","class"))
sig <- rbind(mof, data.frame(signature=paste0("Bailey_",b$class), gene=b$gene),
                  data.frame(signature=paste0("Puleo_",p$class),  gene=p$gene))
N <- 15000  # nominal expressed-gene universe (approximation)
res <- do.call(rbind, lapply(split(sig$gene, sig$signature), function(g){
  g <- unique(g); k <- length(intersect(g, mod3))
  data.frame(list_size=length(g), overlap=k,
    p_hyper=phyper(k-1, length(g), N-length(g), length(mod3), lower.tail=FALSE))}))
res$signature <- rownames(res); res <- res[order(res$p_hyper),]
print(head(res[res$overlap>0, c("signature","list_size","overlap","p_hyper")], 12), row.names=FALSE)
# genes in the module not in any classical/progenitor/basal list
pub <- unique(sig$gene[grepl("Classical|Basal|Progenitor|Squamous|ADEX|Immunogenic|Pure|pure|Stroma|Desno", sig$signature)])
novel <- setdiff(mod3, pub)
cat("\nModule genes in NO published list (", length(novel), "):\n", paste(novel, collapse=", "), "\n", sep="")
write.csv(res, "results/meta/module3_overlap_stats.csv", row.names=FALSE)
write.csv(data.frame(gene=novel), "results/meta/module3_not_in_published_lists.csv", row.names=FALSE)
