#20211117
setwd("F:/gxe_data/YQ180-8/10.��������/1.RBP-lncRNA")

#����RBP-lncRNA������ϵ��
AGTR1 <- read.table("AGTR1.txt", sep = "\t", header = T)
CXCL12 <- read.table("CXCL12.txt", sep = "\t", header = T)
PDGFRL <- read.table("PDGFRL.txt", sep = "\t", header = T)
PTGER3 <- read.table("PTGER3.txt", sep = "\t", header = T)
S1PR1 <- read.table("S1PR1.txt", sep = "\t", header = T)

#��ȡ��������
AGTR1 <- cbind(data.frame(lncRNA = sapply(strsplit(AGTR1$RNA_ID, "_"), "[", 2)), RBP = rep("AGTR1", nrow(AGTR1)))[1:20,]
CXCL12 <- cbind(data.frame(lncRNA = sapply(strsplit(CXCL12$RNA_ID, "_"), "[", 2)), RBP = rep("CXCL12", nrow(CXCL12)))[1:20,]
PDGFRL <- cbind(data.frame(lncRNA = sapply(strsplit(PDGFRL$RNA_ID, "_"), "[", 2)), RBP = rep("PDGFRL", nrow(PDGFRL)))[1:20,]
PTGER3 <- cbind(data.frame(lncRNA = sapply(strsplit(PTGER3$RNA_ID, "_"), "[", 2)), RBP = rep("PTGER3", nrow(PTGER3)))[1:20,]
S1PR1 <- cbind(data.frame(lncRNA = sapply(strsplit(S1PR1$RNA_ID, "_"), "[", 2)), RBP = rep("S1PR1", nrow(S1PR1)))[1:20,]
result <- rbind(AGTR1, CXCL12, PDGFRL, PTGER3, S1PR1)  #�����ļ�

#�����ļ�
RBP <- unique(result$RBP)
lncRNA <- unique(result$lncRNA)
ann <- data.frame(rbind(cbind(RBP, rep("RBP", length(RBP))), cbind(lncRNA, rep("lncRNA", length(lncRNA)))))
colnames(ann)[2] <- "text"

#������
write.table(result, "network.txt", quote = F, row.names = F, sep = "\t")
write.table(ann, "ann.txt", quote = F, row.names = F, sep = "\t")