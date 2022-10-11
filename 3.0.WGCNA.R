#-----------------20210922
setwd("F:/gxe_data/YQ180-8/3.WGCNA")

#-----------------loading packages
library(WGCNA)
enableWGCNAThreads()  #打开多线程

#-----------------input data
data <- read.table("normalize.txt", sep = "\t", row.names = 1, header = T)
group <- read.table("group_new.txt", sep = "\t", header = T)
data <- data[group$accession]  #去除免疫浸润预测结果差的样本之后的表达矩阵，用于后续分析

data <- t(data)

#-----------------样本聚类树
gsg <- goodSamplesGenes(data, verbose = 3)  #对有着过多（minNSamples、minNGenes）缺失的样本和基因进行迭代过滤
gsg$allOK  #查看是否存在过多的缺失值，没有则TRUE
#绘制样本聚类树
sampleTree <- hclust(dist(data), method = "average") #hclust层次聚类，average，计算两个cluster各自数据点的两两距离的平均值
pdf("sampleClustering.pdf", width = 12, height = 9)
sizeGrWindow(12, 9);par(cex = 1);par(mar = c(2, 4, 2, 0)) #打开图形窗口，设置图形参数
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5)
#abline(h = 105, col = "red")
dev.off()

#剔除离群样本
#clust <- cutreeStatic(sampleTree, cutHeight = 105, minSize = 10)
#table(clust)
#keepSamples <- (clust == 1)
#data <- data[keepSamples, ]
#group <- group[group$accession %in% rownames(data),]
write.table(t(data), "exprmatrix.txt", quote = F, row.names = T, sep = "\t")  #最终的表达矩阵
write.table(group, "group_last.txt", quote = F, row.names = F, sep = "\t")  #最终的分组信息

#-----------------软阈值筛选
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))  #为了使网络更符合无尺度分布
sft <- pickSoftThreshold(data, powerVector = powers, verbose = 5, RsquaredCut = 0.85) 
#相关系数的平方越高，网络越接近无尺度网络分布，R^2 > 0.85 的power值用于网络构建
pdf("soft-thresholding.power.pdf", width = 9, height = 5)
par(mfrow = c(1,2))
plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit,signed R^2",type="n",
  main = paste("Scale independence")  
)

text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  labels=powers,
  cex = 0.85,
  col = "red"
)

abline(h = 0.85, col = "red")

plot(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)

text(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  labels = powers,
  cex = 0.85,
  col = "red"
)

power = sft$powerEstimate
power
dev.off()

#---------------构建共表达网络，划分模块
net <- blockwiseModules(data,power = sft$powerEstimate, maxBlockSize = 20000, TOMType = "unsigned", 
                       minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.3, numericLabels = TRUE,
                       pamRespectsDendro = FALSE, verbose = 3)  
#核心步骤，将成千上万个基因划分为几十个模块
#maxBlockSize要大于基因数目确保一次性计算，mergeCutHeight合并模块的阈值(模块间相关系数大于1-0.3合并为一个模块)，minModuleSize模块内最少基因数量
mergedColors <- labels2colors(net$colors)  #每个基因对应哪个颜色
table(mergedColors)  #查看模块划分数量及每个模块内基因数量

#绘制模块聚类树
pdf("Gene_ClusterDendrogram.pdf")  
plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off() 
#上方为基因的聚类树，下方为根据聚类结果分类的modules,
#不同的modules标示不同的颜色，其中灰色代表那些没有归入任何modules的基因

#----------------绘制模块和性状相关性图
nSamples <- nrow(data)
MEList <- moduleEigengenes(data, colors = mergedColors) #计算模块特征值（PC1）
MEs <- MEList$eigengenes
MET <- orderMEs(MEs)  #列为模块，行为样本，数值代表特征值

###免疫细胞（性状）
dataTrait <- read.table("CIBERSORT-Results.txt", sep = "\t", header = T, row.names = 1)
signif_cell <- read.table("signif_cell.txt", header = T)
dataTrait <- dataTrait[group$accession,signif_cell$cell]  #去除了预测结果不好、离群的样本和差异不显著的免疫细胞

###相关性--筛选与表型相关的模块
moduleTraitCor <- cor(MET, dataTrait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
pdf("module_trait_cor.pdf")
textMatrix <- paste(signif(moduleTraitCor, 2),
                   "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(8, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(dataTrait),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix,setStdMargins = FALSE,
               cex.text = 0.4, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

###----------模块基因导出
modNames  <- substring(names(MEs), 3) #模块名称
AllGenes <- colnames(data) #基因名
for (module in modNames)
{modGenes = (mergedColors==module)
modLLIDs = AllGenes[modGenes]; 
fileName = paste("table/LinkID-", module, ".txt", sep="");
write.table(as.data.frame(modLLIDs), file = fileName, row.names = FALSE, col.names = FALSE, quote = F)
}

###导出Cytoscape文件
moduleColors <- labels2colors(net$colors)
#TOM <- TOMsimilarityFromExpr(data, power = sft$powerEstimate)
# Select modules--与大多数免疫细胞显著相关的模块
#module <- "brown"
# Select module genes---brown模块内的基因
#Genes <- colnames(data)
#inModule <- (moduleColors == module)
#modGenes <- Genes[inModule]
# Select the corresponding Topological Overlap
#modTOM <- TOM[inModule, inModule]
#dimnames(modTOM) <- list(modGenes, modGenes)
# cyt <- exportNetworkToCytoscape(modTOM, 
#                                 edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse = "-"), ".txt", sep = ""),
#                                 nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse = "-"), ".txt", sep = ""),
#                                 weighted = TRUE,
#                                 threshold = 0.02,
#                                 nodeNames = modGenes,
#                                 nodeAttr = moduleColors[inModule])

###GS和MM
geneModuleMembership <- as.data.frame(cor(data, MET, use = "p")) #基因表达矩阵与模块特征值相关性
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")
traitNames <- names(dataTrait)
geneTraitSignificance <- as.data.frame(cor(data, dataTrait, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", traitNames, sep="")
names(GSPvalue) <- paste("p.GS.", traitNames, sep="")

###散点图
for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors == module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf = paste("GS_MM/", trait, "_", module, ".pdf", sep="")
      pdf(file = outPdf, width = 7, height = 7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      abline(v = 0.8, h = 0.4, col = "red")
      dev.off()
    }
  }
}

###GS_MM数据
probes = colnames(data)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "table/GS_MM.xls",sep = "\t",row.names = F)
