#--------------20210922
setwd("F:/gxe_data/YQ180-8/1.batch effect")

#--------------loading packages
library(ggplot2)

#--------------input data
non_normalize <- read.table("non-normalize.txt", sep = "\t", header = T, row.names = 1)
normalize <- read.table("normalize.txt", sep = "\t", header = T, row.names = 1)
group <- read.table("group.txt", sep = "\t", header = T)

#--------------PCA--批次矫正前
non_normalize <- as.data.frame(t(non_normalize))
pca_non <- prcomp(non_normalize, scale. = T)  #计算PC1.PC2...
pcaPredict_non <- predict(pca_non)  #提取PC1.PC2...
pdf("PCA_before batch-effect removal.pdf", width = 5.0, height = 4.8)
p1 <- ggplot(as.data.frame(pcaPredict_non), aes(PC1, PC2)) +
  geom_point(aes(color = group$batch), size = 2.5) +
  scale_color_manual(values = c("#ee4266", "#be0aff", "#0d47a1")) +
  theme(panel.grid.major = element_line(color = "#e9ecef", size = 0.1),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.title = element_blank(),
        legend.position = c(.15,.9),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Before batch-effect removal")
print(p1)
dev.off()

#--------------PCA--批次矫正后
normalize <- as.data.frame(t(normalize))
pca <- prcomp(normalize, scale. = T)
pcaPredict <- predict(pca)
pdf("PCA_after batch-effect removal.pdf", width = 5.0, height = 4.8)
p1 <- ggplot(as.data.frame(pcaPredict), aes(PC1, PC2)) +
  geom_point(aes(color = group$batch), size = 2.5) +
  scale_color_manual(values = c("#ee4266", "#be0aff", "#0d47a1")) +
  theme(panel.grid.major = element_line(color = "#e9ecef", size = 0.1),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.title = element_blank(),
        legend.position = c(.85,.9),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "After batch-effect removal")
print(p1)
dev.off()