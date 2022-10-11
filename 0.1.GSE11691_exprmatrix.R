#---------------------20210918
setwd("F:/gxe_data/YQ180-8/0.download data/GSE11691")

#---------------------loading packages
library(GEOquery)

#---------------------download data
gse <- getGEO("GSE11691", destdir = ".", AnnotGPL = F, getGPL = F)  #表达矩阵
gset <- exprs(gse[[1]])
ann <- read.table("ann.txt", sep = "\t", header = T)  #平台文件：探针-gene symbol

#---------------------ID转化
gset <- gset[rownames(gset) %in% ann$ID,]  
rownames(gset) <- ann[match(rownames(gset), ann$ID),2]  #将表达矩阵的行名转换为mRNA的名字
gset <- aggregate(gset, by = list(rownames(gset)), FUN = mean)  #将行名相同的取平均值
colnames(gset)[1] <- "mRNA"

#-------------------分组文件
pdata <- pData(gse[[1]])  #分组数据  #本实验中Endometrium--EU， Endometriosis--EC
pdata$group <- rep(c("EU", "EC"), each = 9)
group <- cbind(title = pdata$title, accession = pdata$geo_accession, group = pdata$group)

#-------------------output data
write.table(gset, "mRNA_expr_GSE11691.txt", quote = F, row.names = F, sep = "\t")
write.table(group, "group.txt", quote = F, row.names = F, sep = "\t")
