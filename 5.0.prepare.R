#-------------------20210927
setwd("F:/gxe_data/YQ180-8/5.machine learning/0.prepare")

#-------------------input data
immune_gene <- read.table("immune_gene.txt", sep = "\t", header = T)
immune_gene <- unique(immune_gene$Symbol)
wgcna_gene <- read.table("wgcna_gene.txt", sep = "\t", header = T)
DEGs <- read.table("diff_gene.txt", sep = "\t")

gene <- Reduce(intersect,list(immune_gene, wgcna_gene$gene, DEGs$V1))  #三者取交集

write.table(gene, "gene.txt", quote = F, row.names = F, col.names = "gene")
write.table(immune_gene, "immune_gene_new.txt", quote = F, row.names = F, col.names = "gene")

