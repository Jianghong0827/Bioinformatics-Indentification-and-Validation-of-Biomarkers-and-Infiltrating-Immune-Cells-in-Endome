#------------------20210922
setwd("F:/gxe_data/YQ180-8/1.batch effect")

#------------------loading packages
library(sva)

#------------------input data
data1 <- read.table("mRNA_expr_GSE11691.txt", sep = "\t", header = T, row.names = 1)
data1 <- log2(data1 + 1)
data2 <- read.table("mRNA_expr_GSE25628.txt", sep = "\t", header = T, row.names = 1)
data3 <- read.table("mRNA_expr_GSE86534.txt", sep = "\t", header = T, row.names = 1)
group1 <- read.table("group_GSE11691.txt", sep = "\t", header = T)
group2 <- read.table("group_GSE25628.txt", sep = "\t", header = T)
group3 <- read.table("group_GSE86534.txt", sep = "\t", header = T)

#------------------对三个文件的基因取交集
gene <- intersect(intersect(rownames(data1), rownames(data2)), rownames(data3))
data <- cbind(data1[gene,], data2[gene,], data3[gene,])  #表达矩阵

#------------------批次矫正
batchType <- c(rep(1,ncol(data1)), rep(2,ncol(data2)),rep(3,ncol(data3)))
outTab <- ComBat(data, batchType, par.prior = T)
outTab <- rbind(gene = colnames(outTab), outTab)

#------------------output data
group <- rbind(group1, group2, group3)
group$batch <- c(rep("GSE11691",ncol(data1)), rep("GSE25628",ncol(data2)),rep("GSE86534",ncol(data3)))
data <- rbind(gene = colnames(data), data)
write.table(data, "non-normalize.txt", quote = F, row.names = T, sep = "\t", col.names = F)
write.table(outTab, "normalize.txt", quote = F, row.names = T, sep = "\t", col.names = F)
write.table(group, "group.txt", quote = F, row.names = F, sep = "\t")
