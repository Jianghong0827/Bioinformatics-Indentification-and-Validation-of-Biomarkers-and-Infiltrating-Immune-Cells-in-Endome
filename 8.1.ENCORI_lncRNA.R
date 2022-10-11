###starbase批量下载lncRNA—miRNA互作关系
###curl 'http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=lncRNA&miRNA=all&clipExpNum=0&degraExpNum=0&pancancerNum=0&programNum=1&program=None&target=all&cellType=all'
#-------------------20210910
setwd("F:/gxe_data/YQ180-8/8.ceRNA/2.miRNA_lncRNA")

#-------------------input data
lncRNA_miRNA_all <- read.table("lncRNA_miRNA(all).txt", sep = "\t", header = T)
lncRNA_miRNA_all <- lncRNA_miRNA_all[lncRNA_miRNA_all$clipExpNum > 1 &lncRNA_miRNA_all$degraExpNum > 1,]
miRNA <- read.table("miRNA.txt", header = T)

#-------------------获取lncRNA-miRNA关系对
lncRNA_miRNA <- lncRNA_miRNA_all[lncRNA_miRNA_all$miRNAname %in% miRNA$miRNA,]
lncRNA_miRNA <- cbind(lncRNA_miRNA$miRNAname, lncRNA_miRNA$geneName)
lncRNA_miRNA <- lncRNA_miRNA[!duplicated(lncRNA_miRNA),]
lncRNA_miRNA <- data.frame(lncRNA_miRNA)
colnames(lncRNA_miRNA) <- c("miRNA", "lncRNA")

#-------------------output data
write.table(lncRNA_miRNA, "lncRNA_miRNA.txt", sep = "\t", quote = F, row.names = F)
write.table(unique(lncRNA_miRNA$lncRNA), "lncRNA.txt", quote = F, row.names = F, col.names = "lncRNA")
