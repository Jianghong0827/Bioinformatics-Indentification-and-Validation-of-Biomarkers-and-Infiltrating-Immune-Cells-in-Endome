#-----------------20210927
setwd("F:/gxe_data/YQ180-8/5.machine learning/1.ROC")

#-----------------loading packages
library(pROC)

#-----------------input data
roc_gene <- read.table("roc_gene.txt", header = T)
expr <- read.csv("mRNA_diff_expr_matrix.csv", header = T, row.names = 1)
expr_gene <- expr[roc_gene$gene,]
expr_gene <- as.matrix(expr_gene)
group <- read.table("group_last.txt", sep = "\t", header = T)

#-----------------ROCÇúÏß»æÖÆ
group <- group$group
group <- factor(group, levels = c("EU", "EC"))

par(mfrow = c(2,2))
for(i in 1:length(roc_gene$gene)){
  plot.roc(group, expr_gene[i,], main = roc_gene$gene[i], print.auc = T, percent = T, cex.lab = 1.5, print.auc.cex = 1.5)
}
