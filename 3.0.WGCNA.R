#-----------------20210922
setwd("F:/gxe_data/YQ180-8/3.WGCNA")

#-----------------loading packages
library(WGCNA)
enableWGCNAThreads()  #�򿪶��߳�

#-----------------input data
data <- read.table("normalize.txt", sep = "\t", row.names = 1, header = T)
group <- read.table("group_new.txt", sep = "\t", header = T)
data <- data[group$accession]  #ȥ�����߽���Ԥ�����������֮��ı���������ں�������

data <- t(data)

#-----------------����������
gsg <- goodSamplesGenes(data, verbose = 3)  #�����Ź��ࣨminNSamples��minNGenes��ȱʧ�������ͻ�����е�������
gsg$allOK  #�鿴�Ƿ���ڹ����ȱʧֵ��û����TRUE
#��������������
sampleTree <- hclust(dist(data), method = "average") #hclust��ξ��࣬average����������cluster�������ݵ�����������ƽ��ֵ
pdf("sampleClustering.pdf", width = 12, height = 9)
sizeGrWindow(12, 9);par(cex = 1);par(mar = c(2, 4, 2, 0)) #��ͼ�δ��ڣ�����ͼ�β���
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5)
#abline(h = 105, col = "red")
dev.off()

#�޳���Ⱥ����
#clust <- cutreeStatic(sampleTree, cutHeight = 105, minSize = 10)
#table(clust)
#keepSamples <- (clust == 1)
#data <- data[keepSamples, ]
#group <- group[group$accession %in% rownames(data),]
write.table(t(data), "exprmatrix.txt", quote = F, row.names = T, sep = "\t")  #���յı������
write.table(group, "group_last.txt", quote = F, row.names = F, sep = "\t")  #���յķ�����Ϣ

#-----------------����ֵɸѡ
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))  #Ϊ��ʹ����������޳߶ȷֲ�
sft <- pickSoftThreshold(data, powerVector = powers, verbose = 5, RsquaredCut = 0.85) 
#���ϵ����ƽ��Խ�ߣ�����Խ�ӽ��޳߶�����ֲ���R^2 > 0.85 ��powerֵ�������繹��
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

#---------------�������������磬����ģ��
net <- blockwiseModules(data,power = sft$powerEstimate, maxBlockSize = 20000, TOMType = "unsigned", 
                       minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.3, numericLabels = TRUE,
                       pamRespectsDendro = FALSE, verbose = 3)  
#���Ĳ��裬����ǧ��������򻮷�Ϊ��ʮ��ģ��
#maxBlockSizeҪ���ڻ�����Ŀȷ��һ���Լ��㣬mergeCutHeight�ϲ�ģ�����ֵ(ģ������ϵ������1-0.3�ϲ�Ϊһ��ģ��)��minModuleSizeģ�������ٻ�������
mergedColors <- labels2colors(net$colors)  #ÿ�������Ӧ�ĸ���ɫ
table(mergedColors)  #�鿴ģ�黮��������ÿ��ģ���ڻ�������

#����ģ�������
pdf("Gene_ClusterDendrogram.pdf")  
plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off() 
#�Ϸ�Ϊ����ľ��������·�Ϊ���ݾ����������modules,
#��ͬ��modules��ʾ��ͬ����ɫ�����л�ɫ������Щû�й����κ�modules�Ļ���

#----------------����ģ�����״�����ͼ
nSamples <- nrow(data)
MEList <- moduleEigengenes(data, colors = mergedColors) #����ģ������ֵ��PC1��
MEs <- MEList$eigengenes
MET <- orderMEs(MEs)  #��Ϊģ�飬��Ϊ��������ֵ��������ֵ

###����ϸ������״��
dataTrait <- read.table("CIBERSORT-Results.txt", sep = "\t", header = T, row.names = 1)
signif_cell <- read.table("signif_cell.txt", header = T)
dataTrait <- dataTrait[group$accession,signif_cell$cell]  #ȥ����Ԥ�������á���Ⱥ�������Ͳ��첻����������ϸ��

###�����--ɸѡ�������ص�ģ��
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

###----------ģ����򵼳�
modNames  <- substring(names(MEs), 3) #ģ������
AllGenes <- colnames(data) #������
for (module in modNames)
{modGenes = (mergedColors==module)
modLLIDs = AllGenes[modGenes]; 
fileName = paste("table/LinkID-", module, ".txt", sep="");
write.table(as.data.frame(modLLIDs), file = fileName, row.names = FALSE, col.names = FALSE, quote = F)
}

###����Cytoscape�ļ�
moduleColors <- labels2colors(net$colors)
#TOM <- TOMsimilarityFromExpr(data, power = sft$powerEstimate)
# Select modules--����������ϸ��������ص�ģ��
#module <- "brown"
# Select module genes---brownģ���ڵĻ���
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

###GS��MM
geneModuleMembership <- as.data.frame(cor(data, MET, use = "p")) #������������ģ������ֵ�����
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")
traitNames <- names(dataTrait)
geneTraitSignificance <- as.data.frame(cor(data, dataTrait, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", traitNames, sep="")
names(GSPvalue) <- paste("p.GS.", traitNames, sep="")

###ɢ��ͼ
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

###GS_MM����
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