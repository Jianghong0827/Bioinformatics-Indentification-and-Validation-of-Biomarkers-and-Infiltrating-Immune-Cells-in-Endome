#-------------------20210922
setwd("F:/gxe_data/YQ180-8/2.immune")

#-------------------loading packages
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggpubr)

#-------------------input data
data <- read.table("CIBERSORT-Results.txt", sep = "\t", row.names = 1, header = T)
group <- read.table("group.txt", sep = "\t", header = T)

#-------------------删除p值大于0.05的预测
data <- data[data$P.value < 0.05, 1:22]
group <- group[group$accession %in% rownames(data),]

#-------------------输出预测结果好的样本
#write.table(group, "group_new.txt", sep = "\t", quote = F, row.names = F)

#-------------------将EU和EC分类，并取均值
data <- as.data.frame(t(data))
group_EC <- subset(group, group == "EC");group_EU <- subset(group, group == "EU")  #将EU和EC分开
data_EC <- data[,group_EC$accession];data_EU <- data[,group_EU$accession]
data_EC$mean <- rowMeans(data_EC);data_EU$mean <- rowMeans(data_EU)  #免疫细胞在所有样本中占比的均值

#-------------------环状条形图
data_cir_EC <- data_EC[ncol(data_EC)];data_cir_EU <- data_EU[ncol(data_EU)] #取作图数值列
data_cir_EC$cell <- rownames(data_cir_EC);data_cir_EU$cell <- rownames(data_cir_EU)
data_cir_EC <- data_cir_EC[order(data_cir_EC$mean),];data_cir_EU <- data_cir_EU[order(data_cir_EU$mean),]  #按比例均值的大小排序

#--------EC
data_cir_EC$id <- seq(1,nrow(data_cir_EC),1)
label_EC <- data_cir_EC  #添加标签
number_of_bar_EC <- nrow(label_EC)
angle_EC <-  90 - 360 * (label_EC$id-0.5) /number_of_bar_EC  #计算标签旋转角度
label_EC$hjust<-ifelse(angle_EC < -90, 1, 0)
label_EC$angle<-ifelse(angle_EC < -90, angle_EC + 180, angle_EC)
col_fun1 <- brewer.pal(9, "Set1");col_fun2 <- brewer.pal(8, "Set2");col_fun3 <- brewer.pal(5, "Set3")
col_fun <- c(col_fun1, col_fun2, col_fun3)  #颜色设置
pdf("EC_circle.pdf", width = 8,height = 8)
p1 <- ggplot(data_cir_EC, aes(x = as.factor(id), y = mean)) +
  geom_bar(stat="identity", fill = col_fun) +
  theme_minimal() +
  ylim(-0.05,max(data_cir_EC$mean) + 0.05) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm")) +
  coord_polar(start = 0) +
  geom_text(data = label_EC, aes(x = id, y = 0.01, label = cell, hjust = hjust), 
            color = "black", fontface ="bold", alpha = 0.6, size = 4.0, 
            angle = label_EC$angle, inherit.aes = FALSE)
print(p1)
dev.off()

#--------EU
data_cir_EU$id <- seq(1,nrow(data_cir_EU),1)
label_EU <- data_cir_EU  #添加标签
number_of_bar_EU <- nrow(label_EU)
angle_EU <-  90 - 360 * (label_EU$id-0.5) /number_of_bar_EU  #计算标签旋转角度
label_EU$hjust<-ifelse(angle_EU < -90, 1, 0)
label_EU$angle<-ifelse(angle_EU < -90, angle_EU + 180, angle_EU)
pdf("EU_circle.pdf", width = 8,height = 8)
p1 <- ggplot(data_cir_EU, aes(x = as.factor(id), y = mean)) +
  geom_bar(stat="identity", fill = col_fun) +
  theme_minimal() +
  ylim(-0.065,max(data_cir_EU$mean) + 0.05) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm")) +
  coord_polar(start = 0) +
  geom_text(data = label_EU, aes(x = id, y = 0.01, label = cell, hjust = hjust), 
            color = "black", fontface ="bold", alpha = 0.6, size = 4.0, 
            angle = label_EU$angle, inherit.aes = FALSE)
print(p1)
dev.off()

#--------------------小提琴图
data$cell <- rownames(data)
data_violin <- melt(data)
data_violin$group <- group[match(data_violin$variable, group$accession),3] #作图数据
#--统计检验组间差异
data_signif <- list()   #指定免疫细胞在所有样本中的比例
wilcox <- list()
result <- data.frame()
cell <- unique(data_violin$cell)
for(i in 1:length(cell)){
  data_signif[[i]] <- data_violin[data_violin$cell == cell[i],]
  wilcox[[i]] <- wilcox.test(value~group, data = data_signif[[i]])
  a <- cbind(cell[i], wilcox[[i]]$p.value)
  result <- rbind(result, a)
}   #result中是每个免疫细胞在EC和EU之间的显著性p值

num <- data.frame()  #定义小提琴图中标签位置
for(i in 1:length(data_signif)){
  b <- max(data_signif[[i]]$value)
  num <- rbind(num, b)
}

result$V2 <- format(as.numeric(result$V2), scientific = F) #科学计数法改为数字，不然会识别错误

label <- ifelse(result$V2 < 0.05,
                ifelse(result$V2 < 0.001,
                       "p<0.001", "p<0.05"), paste0("p=",sprintf("%.3f", as.numeric(result$V2))))
label_new <- cbind(result$V1, label, num)
colnames(label_new) <- c("cell", "label", "num")
label_new$label[5] <- ""  #在所有样本中均为0的免疫细胞，p值不显示
data_violin_new <- cbind(data_violin, label_new[2:3])  #将标签数据整合到作图数据中

#--------plot
pdf("violin.pdf", height = 8, width = 14)
p3 <- ggplot(data_violin_new, aes(x = cell, y = value, fill = group)) +
  geom_violin(scale = "width", width = 0.6, position = position_dodge(0.6)) +
  stat_summary(fun = median, geom = "point", size = 1.5, color = "white", position = position_dodge(0.6)) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "ratio", title = "Immune cell infiltration in ectopic and eutopic endometrium") +
  geom_text(data = data_violin_new[1:22,], aes(label = label, x = cell, y = num + 0.02)) +
  geom_errorbar(data = data_violin_new[1:22,], aes(y = num + 0.01, ymin = num + 0.01, ymax = num + 0.01, width = 0.4))

print(p3)
dev.off()

#-------------将有显著差异的免疫细胞输出
result <- result[result$V2 != "          NaN",]
out <- result[result$V2 < 0.05,1]
write.table(out, "signif_cell.txt", quote = F, row.names = F, col.names = "cell")
