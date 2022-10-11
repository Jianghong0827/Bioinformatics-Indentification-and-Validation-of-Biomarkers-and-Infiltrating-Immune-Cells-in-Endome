#----------------------20210924
setwd("F:/gxe_data/YQ180-8/4.DEG/2.pathway")

#----------------------loading packages
library(clusterProfiler)
library(org.Hs.eg.db)

#----------------------input data
genes <- read.table('diff_gene.txt')
eg <- bitr(genes$V1, fromType="SYMBOL", toType=c("ENTREZID"), 
           OrgDb="org.Hs.eg.db")   #将symbol转换为ENTREZID和ENSEMBL

#write.table(eg, "id.txt", sep = "\t", quote = F, row.names = F)

#-------------------GO富集分析
enrich.go <- enrichGO(gene = eg$ENTREZID, OrgDb = "org.Hs.eg.db", 
                      keyType = "ENTREZID", ont = "ALL", 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2)    #富集分析

write.table(as.data.frame(enrich.go), 'go.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#--------------------KEGG富集分析
kegg <- enrichKEGG(gene = eg$ENTREZID, keyType = "kegg", organism = "human",
                   pAdjustMethod = "none", pvalueCutoff = 0.05, qvalueCutoff = 0.4)

write.table(kegg, 'kegg.txt', sep = '\t', quote = FALSE, row.names = FALSE)
