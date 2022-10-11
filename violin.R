#20211125
setwd('F:/gxe_data/YQ180-8/陆康补充')
#-------------------loading packages
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggpubr)

#-------------------input data
data <- read.table("exprmatrix.txt", sep = "\t", row.names = 1, header = T)
gene <- c("AGTR1","CXCL12","PDGFRL","PTGER3","S1PR1")
data <- data[gene,]
group <- read.table("group_last.txt", sep = "\t", header = T)

#--------------------小????图
data$gene <- rownames(data)
data_violin <- melt(data)
data_violin$group <- group[match(data_violin$variable, group$accession),3]
#--统?萍???????????
data_signif <- list() 
wilcox <- list()
result <- data.frame()
gene <- unique(data_violin$gene)
for(i in 1:length(gene)){
  data_signif[[i]] <- data_violin[data_violin$gene == gene[i],]
  wilcox[[i]] <- wilcox.test(value~group, data = data_signif[[i]])
  a <- cbind(gene[i], wilcox[[i]]$p.value)
  result <- rbind(result, a)
}   #result????每??????细????EC??EU之??????????p值

num <- data.frame()  #????小????图?斜?签位??
for(i in 1:length(data_signif)){
  b <- max(data_signif[[i]]$value)
  num <- rbind(num, b)
}

result$V2 <- format(as.numeric(result$V2), scientific = F) #??学????????为???郑???然??识??????

label <- ifelse(result$V2 < 0.001, "p<0.001", paste0("p=",sprintf("%.3f", as.numeric(result$V2))))
label_new <- cbind(result$V1, label, num)
colnames(label_new) <- c("gene", "label", "num")
data_violin_new <- cbind(data_violin, label_new[2:3])  #????签???????系???图??????

#--------plot
pdf("violin.pdf", height = 8, width = 14)
p3 <- ggplot(data_violin_new, aes(x = gene, y = value, fill = group)) +
  geom_violin(scale = "width", width = 0.6, position = position_dodge(0.6)) +
  stat_summary(fun = median, geom = "point", size = 1.5, color = "white", position = position_dodge(0.6)) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "expression", title = "") +
  geom_text(data = data_violin_new[1:7,], aes(label = label, x = gene, y = num + 0.18)) +
  geom_errorbar(data = data_violin_new[1:7,], aes(y = num + 0.04, ymin = num + 0.04, ymax = num + 0.04, width = 0.3))

print(p3)
dev.off()

#--------plot
pdf("boxplot.pdf", height = 8, width = 10)
p7 <- ggplot(data_violin_new, aes(x = gene, y = value, fill = group)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6)) +
  stat_summary(fun = median, geom = "point", size = 1.5, color = "white", position = position_dodge(0.6)) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "expression", title = "") +
  geom_text(data = data_violin_new[1:7,], aes(label = label, x = gene, y = num + 0.18)) +
  geom_errorbar(data = data_violin_new[1:7,], aes(y = num + 0.04, ymin = num + 0.04, ymax = num + 0.04, width = 0.3))

print(p7)
dev.off()
