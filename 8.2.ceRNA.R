#--------------------20211020
setwd("F:/gxe_data/YQ180-8/8.ceRNA/3.ceRNA")

#--------------------input data
mRNA_miRNA <- read.table("miRNA_mRNA.txt", sep = "\t", header = T)
miRNA_lncRNA <- read.table("lncRNA_miRNA.txt", sep = "\t", header = T)

#--------------------ceRNA网络数据准备
miRNA <- data.frame(unique(mRNA_miRNA$miRNA))
miRNA$text <- rep("miRNA", length(rownames(miRNA)))
colnames(miRNA)[1] <- "miRNA"

mRNA <- data.frame(unique(mRNA_miRNA$mRNA))
mRNA$text <- rep("mRNA", length(rownames(mRNA)))
colnames(mRNA)[1] <- "mRNA"

lncRNA <- data.frame(unique(miRNA_lncRNA$lncRNA))
lncRNA$text <- rep("lncRNA", length(rownames(lncRNA)))
colnames(lncRNA)[1] <- "lncRNA"

colnames(mRNA_miRNA)[2] <- "others"
colnames(miRNA_lncRNA)[2] <- "others"
ceRNA <- rbind(mRNA_miRNA, miRNA_lncRNA)

#---------------------output data
write.table(mRNA, "mRNA.txt", quote = F, sep = "\t", row.names = F)
write.table(miRNA, "miRNA.txt", quote = F, sep = "\t", row.names = F)
write.table(lncRNA, "lncRNA.txt", quote = F, sep = "\t", row.names = F)
write.table(ceRNA, "ceRNA.txt", quote = F, sep = "\t", row.names = F)
