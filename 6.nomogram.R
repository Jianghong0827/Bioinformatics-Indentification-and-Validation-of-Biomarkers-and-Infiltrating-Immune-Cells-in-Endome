#----------------------20210928
setwd("F:/gxe_data/YQ180-8/6.Nomogram")

#----------------------loading packages
library(rms)
library(rmda)

#----------------------input data
gene <- read.table("gene.txt", header = T)
expr <- read.table("exprmatrix.txt", header =  T, row.names = 1)
group <- read.table("group_last.txt", header = T, sep = "\t")

#----------------------��������
expr_gene <- expr[gene$gene,]
expr_gene <- data.frame(t(expr_gene))
expr_gene$group <- group[match(rownames(expr_gene), group$accession),3]
expr_gene$group <- factor(expr_gene$group, levels = c("EU", "EC"))
str(expr_gene$group)
#expr_gene$group <- ifelse(expr_gene$group == "EU", 0, 1)
#expr_gene$group <- factor(expr_gene$group, levels = c(0,1), labels = c("EU", "EC"))

#-----------------------����ͼ����
ddist <- datadist(expr_gene); options(datadist='ddist')
fit1 <- lrm(group ~ CXCL12+PDGFRL+AGTR1+PTGER3+S1PR1, data = expr_gene)  #ģ�����
nom <- nomogram(fit1, fun = function(x)plogis(x), funlabel = "EC", lp = F)  #lp = F����ʾLinear Predictor;plogis(x)��Ϊ1/(1+expr(-x))Ҳ�ɣ���Ϊ������Ƕ������logistic�ع�
pdf("nomogram.pdf", width = 8, height = 8)
plot(nom, col.grid = gray(c(0.8, 0.95)))
dev.off()

#-----------------------У׼����
fit2 <- lrm(group ~ CXCL12+PDGFRL+AGTR1+PTGER3+S1PR1, data = expr_gene, x = T, y = T) #ģ�����
call <- calibrate(fit2, cmethod = "KM", method = "boot", B = 1000)
pdf("calibration.pdf", width = 8, height = 8)
plot(call)  #Apparent�����ο��ߣ�����Ԥ����ʵ�����Ǻϣ���Ideal����δУ���ģ�Bias����У����ģ��������������Ƿ�ӽ���Apparent 
dev.off()

#------------------------DCA�����ߣ����߻���
expr_gene$group <- ifelse(expr_gene$group == "EU", 0, 1)  #����������������ͬ������Ҫ�����ֱ�ʾ
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
                         confidence.intervals = 0.95)  #���л��򣬼�Nonmogram
model_all <- list(model2, model4, model5, model6, model8, model9)
pdf("decision_curve.pdf", width = 8, height = 8)
plot_decision_curve(model_all,
                    confidence.intervals = F,
                    curve.names = c("CXCL12", "PDGFRL", "AGTR1", "PTGER3", "S1PR1", "Nomogram"),
                    legend.position = "bottomright",
                    standardize = F)
dev.off()

#------------------------�ٴ�Ӱ�����ߣ�Clinical Impact Curve��
pdf("clinical_impact_curve.pdf", width = 8, height = 8)
plot_clinical_impact(model9, population.size = 1000, cost.benefit.axis = T,
                     n.cost.benefits = 8, col = c("red", "blue"),
                     confidence.intervals = T)  #confidence.intervals = Tչʾ��������
dev.off()