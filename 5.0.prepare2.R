#------------------20211012
setwd("F:/gxe_data/YQ180-8/5.machine learning/1.random_forest")

#------------------input data
expr <- read.table("exprmatrix.txt", sep = "\t", row.names = 1, header = T)
group <- read.table("group_last.txt", sep = "\t", row.names = 1, header = T)
gene <- read.table("gene.txt", header = T)

#-----------------��ȡ�ؼ�����ı������
expr <- expr[gene$gene,]
expr <- as.data.frame(t(expr))
expr$group <- group[match(rownames(expr), group$accession), 2]

#-----------------output data
write.table(expr, "expr_kyegene.txt", sep = "\t", row.names = F, quote = F)