#---------------------20210924
setwd("F:/gxe_data/YQ180-8/4.DEG/4.plot")

#---------------------loading packages
library(reshape2)

#---------------------input data
kegg <- read.table("kegg_new.txt", sep = "\t", header = T)
interaction <- read.table("interactions.txt", sep = "\t", header = T)
id <- read.table("id.txt", sep = "\t", header = T)

###melt
kegg_melt <- melt(kegg, id = "Description")
kegg_melt <- na.omit(kegg_melt)

#ID转换
kegg_melt$symbol <- id[match(kegg_melt$value, id$ENTREZID),1]
kegg_melt <- kegg_melt[c(1,4)]
#write.table(kegg_melt, "pathway.txt", sep = "\t", row.names = F, quote = F)

#---------------------网络输入文件整理
interaction_new <- interaction[interaction$node1 %in% kegg_melt$symbol & interaction$node2 %in% kegg_melt$symbol,]
colnames(kegg_melt) <- colnames(interaction_new)
result <- rbind(interaction_new, kegg_melt)
write.table(result, "mRNA-Pathway.txt", sep = "\t", quote = F, row.names = F)  #网络互作文件

#网络节点属性文件
node_mRNA <- unique(kegg_melt$node2)
node_pathway <- unique(kegg_melt$node1)

node_mRNA <- as.data.frame(node_mRNA)
node_pathway <- as.data.frame(node_pathway)

node_mRNA$text <- rep("mRNA", length(rownames(node_mRNA))) #注释其属于哪种RNA
node_pathway$text <- rep("pathway", length(rownames(node_pathway)))

write.table(node_mRNA, "node_mRNA.txt", sep = "\t", row.names = F, quote = F)
write.table(node_pathway, "node_pathway.txt", sep = "\t", row.names = F, quote = F)
