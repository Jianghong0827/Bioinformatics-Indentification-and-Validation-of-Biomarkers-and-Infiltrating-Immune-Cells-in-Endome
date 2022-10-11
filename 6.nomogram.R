#----------------------20210928
setwd("F:/gxe_data/YQ180-8/6.Nomogram")

#----------------------loading packages
library(rms)
library(rmda)

#----------------------input data
gene <- read.table("gene.txt", header = T)
expr <- read.table("exprmatrix.txt", header =  T, row.names = 1)
group <- read.table("group_last.txt", header = T, sep = "\t")

#----------------------数据整理
expr_gene <- expr[gene$gene,]
expr_gene <- data.frame(t(expr_gene))
expr_gene$group <- group[match(rownames(expr_gene), group$accession),3]
expr_gene$group <- factor(expr_gene$group, levels = c("EU", "EC"))
str(expr_gene$group)
#expr_gene$group <- ifelse(expr_gene$group == "EU", 0, 1)
#expr_gene$group <- factor(expr_gene$group, levels = c(0,1), labels = c("EU", "EC"))

#-----------------------列线图绘制
ddist <- datadist(expr_gene); options(datadist='ddist')
fit1 <- lrm(group ~ CXCL12+PDGFRL+AGTR1+PTGER3+S1PR1, data = expr_gene)  #模型拟合
nom <- nomogram(fit1, fun = function(x)plogis(x), funlabel = "EC", lp = F)  #lp = F不显示Linear Predictor;plogis(x)改为1/(1+expr(-x))也可，因为因变量是二分类的logistic回归
pdf("nomogram.pdf", width = 8, height = 8)
plot(nom, col.grid = gray(c(0.8, 0.95)))
dev.off()

#-----------------------校准曲线
fit2 <- lrm(group ~ CXCL12+PDGFRL+AGTR1+PTGER3+S1PR1, data = expr_gene, x = T, y = T) #模型拟合
call <- calibrate(fit2, cmethod = "KM", method = "boot", B = 1000)
pdf("calibration.pdf", width = 8, height = 8)
plot(call)  #Apparent代表参考线（理论预测与实际相吻合），Ideal代表未校正的，Bias代表校正后的，看这两条曲线是否接近于Apparent 
dev.off()

#------------------------DCA（决策）曲线绘制
expr_gene$group <- ifelse(expr_gene$group == "EU", 0, 1)  #和上面两个分析不同，这里要用数字表示
# model1 <- decision_curve(group ~ VCAM1, expr_gene[c(1,9)],
#                         family = binomial(link = "logit"),
#                         thresholds = seq(0, 1, by = 0.01),
#                         confidence.intervals = 0.95)
model2 <- decision_curve(group ~ CXCL12, data = expr_gene[c(2,9)],
                         family = binomial(link = "logit"),
                         thresholds = seq(0, 1, by = 0.01),
                         confidence.intervals = 0.95)
# model3 <- decision_curve(group ~ NGF, data = expr_gene[c(3,9)],
#                          family = binomial(link = "logit"),
#                          thresholds = seq(0, 1, by = 0.01),
#                          confidence.intervals = 0.95)
model4 <- decision_curve(group ~ PDGFRL, data = expr_gene[c(4,9)],
                         family = binomial(link = "logit"),
                         thresholds = seq(0, 1, by = 0.01),
                         confidence.intervals = 0.95)
model5 <- decision_curve(group ~ AGTR1, data = expr_gene[c(5,9)],
                         family = binomial(link = "logit"),
                         thresholds = seq(0, 1, by = 0.01),
                         confidence.intervals = 0.95)
model6 <- decision_curve(group ~ PTGER3, data = expr_gene[c(6,9)],
                         family = binomial(link = "logit"),
                         thresholds = seq(0, 1, by = 0.01),
                         confidence.intervals = 0.95)
# model7 <- decision_curve(group ~ PTGFR, data = expr_gene[c(7,9)],
#                          family = binomial(link = "logit"),
#                          thresholds = seq(0, 1, by = 0.01),
#                          confidence.intervals = 0.95)
model8 <- decision_curve(group ~ S1PR1, data = expr_gene[c(8,9)],
                         family = binomial(link = "logit"),
                         thresholds = seq(0, 1, by = 0.01),
                         confidence.intervals = 0.95)
model9 <- decision_curve(group ~ CXCL12+PDGFRL+AGTR1+PTGER3+S1PR1, expr_gene,
                         family = binomial(link = "logit"),
                         thresholds = seq(0, 1, by = 0.01),
                         confidence.intervals = 0.95)  #所有基因，即Nonmogram
model_all <- list(model2, model4, model5, model6, model8, model9)
pdf("decision_curve.pdf", width = 8, height = 8)
plot_decision_curve(model_all,
                    confidence.intervals = F,
                    curve.names = c("CXCL12", "PDGFRL", "AGTR1", "PTGER3", "S1PR1", "Nomogram"),
                    legend.position = "bottomright",
                    standardize = F)
dev.off()

#------------------------临床影响曲线（Clinical Impact Curve）
pdf("clinical_impact_curve.pdf", width = 8, height = 8)
plot_clinical_impact(model9, population.size = 1000, cost.benefit.axis = T,
                     n.cost.benefits = 8, col = c("red", "blue"),
                     confidence.intervals = T)  #confidence.intervals = T展示置信区间
dev.off()
