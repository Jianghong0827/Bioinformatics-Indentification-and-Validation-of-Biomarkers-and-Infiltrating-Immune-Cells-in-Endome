#---------------------20210918
setwd("F:/gxe_data/YQ180-8/0.download data/GSE11691")

#---------------------loading packages
library(GEOquery)

#---------------------download data
gse <- getGEO("GSE11691", destdir = ".", AnnotGPL = F, getGPL = F)  #�������
gset <- exprs(gse[[1]])
ann <- read.table("ann.txt", sep = "\t", header = T)  #ƽ̨�ļ���̽��-gene symbol

#---------------------IDת��
gset <- gset[rownames(gset) %in% ann$ID,]  
rownames(gset) <- ann[match(rownames(gset), ann$ID),2]  #��������������ת��ΪmRNA������
gset <- aggregate(gset, by = list(rownames(gset)), FUN = mean)  #��������ͬ��ȡƽ��ֵ
colnames(gset)[1] <- "mRNA"

#-------------------�����ļ�
pdata <- pData(gse[[1]])  #��������  #��ʵ����Endometrium--EU�� Endometriosis--EC
pdata$group <- rep(c("EU", "EC"), each = 9)
group <- cbind(title = pdata$title, accession = pdata$geo_accession, group = pdata$group)

#-------------------output data
write.table(gset, "mRNA_expr_GSE11691.txt", quote = F, row.names = F, sep = "\t")
write.table(group, "group.txt", quote = F, row.names = F, sep = "\t")