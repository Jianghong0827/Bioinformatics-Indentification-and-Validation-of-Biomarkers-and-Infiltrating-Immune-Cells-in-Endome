#--------------------20210924
setwd("F:/gxe_data/YQ180-8/4.DEG/1.heatmap")

#--------------------loading packages
library(limma)
library(ComplexHeatmap)

#--------------------input data
expr <- read.table("exprmatrix.txt", sep = "\t", row.names = 1, header = T)
group <- read.table("group_last.txt", sep = "\t", header = T)
group_list <- group$group

#--------------------差异分析
design <- model.matrix(~0+factor(group_list, levels = c("EC", "EU")))
colnames(design) <- c("EC", "EU")
rownames(design) <- colnames(expr) #分组矩阵
contrast.matrix <- makeContrasts(EC-EU, levels = design) #差异比较矩阵

#-------------------差异分析
fit <- lmFit(expr,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput <- topTable(fit2, coef = 1, n = Inf)
nrDEG <- na.omit(tempOutput)  
diff_expr <- nrDEG   #差异表达矩阵
#write.csv(diff_expr, "diff_expr.csv", row.names = T, quote = F)

#-------------------筛选差异表达mRNA
#diff_expr <- read.csv("diff_expr.csv", header = T, row.names = 1)
diff_expr$DEG <- ifelse(diff_expr$adj.P.Val < 0.05 &abs(diff_expr$logFC) >= 1,  
                        ifelse(diff_expr$logFC >= 1, "Up", "Down"), "no diff")  #定义上下调基因

num_up_gene <- length(rownames(diff_expr[diff_expr$DEG == "Up",]))      #上调基因
diff_up_gene <- rownames(diff_expr[diff_expr$DEG == "Up",])  
#write.table(diff_up_gene, "diff_up_gene.txt", quote = F, row.names = F, col.names = F)

num_down_gene <- length(rownames(diff_expr[diff_expr$DEG == "Down",]))  #下调基因
diff_down_gene <- rownames(diff_expr[diff_expr$DEG == "Down",])  
#write.table(diff_down_gene, "diff_down_gene.txt", quote = F, row.names = F, col.names = F)

diff_gene <- rownames(diff_expr[diff_expr$DEG != "no diff",])           #所有差异基因
#write.table(diff_gene, "diff_gene.txt", quote = F, row.names = F, col.names = F) 

#--------------------热图
diff_expr_matrix <- expr[rownames(expr) %in% diff_gene,]  #筛选差异基因的表达矩阵
#write.csv(diff_expr_matrix, "mRNA_diff_expr_matrix.csv", quote = F, row.names = T)
#colnames(diff_expr_matrix) <- group$title
diff_expr_matrix <- as.matrix(diff_expr_matrix)
pdf("heatmap.pdf", width = 10.0, height = 10.0)
p2 <- Heatmap(diff_expr_matrix,  
              show_row_names = F,
              column_km = 2,column_gap = unit(3, "mm"),
              column_names_rot = 45,
              heatmap_legend_param = list(title = "expr"),
              row_dend_width = unit(1.5, "cm"))
print(p2)
dev.off()
