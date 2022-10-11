#-------------------20210927
setwd("F:/gxe_data/YQ180-8/3.WGCNA/table")

#-------------------input data
data <- read.table("GS_MM.xls", sep = "\t", header = T, row.names = 1)

#-------------------��ѡģ�������ϸ��
module <- "brown"
column <- c("GS.T.cells.follicular.helper", "GS.NK.cells.activated", "GS.Macrophages.M2", "MMbrown")
data_new <- data[data$moduleColor == module, column]

#--------------------��ѡ����Ҫ���ģ�顢��������ϸ��
cell <- list()
for(i in 1:(length(column)-1)){
  cell[[i]] <- data_new[abs(data_new[i]) > 0.4 & abs(data_new$MMbrown) > 0.8, c(i,length(colnames(data_new)))]
  names(cell)[i] <- sapply(strsplit(column[i], "GS."), "[", 2)
}

#--------------------��ȡ������ϸ����صĻ���(����)
gene <- Reduce(intersect,list(rownames(cell[[1]]),rownames(cell[[2]]),rownames(cell[[3]])))
#gene <- intersect(rownames(cell[[1]]), intersect(rownames(cell[[2]]), rownames(cell[[3]])))

#--------------------output data
write.table(gene, "gene.txt", quote = F, row.names = F, col.names = "gene")
