#----------------20211020
setwd("F:/gxe_data/YQ180-8/8.ceRNA/1.mRNA_miRNA")
library(dplyr)

#----------------input data
ENCORI_AGTR1 <- read.table("ENCORI_hg19_CLIP-seq_miRNA-target_all_AGTR1.txt", sep = "\t", header = T)
ENCORI_CXCL12 <- read.table("ENCORI_hg19_CLIP-seq_miRNA-target_all_CXCL12.txt", sep = "\t", header = T)
ENCORI_PDGFRL <- read.table("ENCORI_hg19_CLIP-seq_miRNA-target_all_PDGFRL.txt", sep = "\t", header = T)
ENCORI_PTGER3 <- read.table("ENCORI_hg19_CLIP-seq_miRNA-target_all_PTGER3.txt", sep = "\t", header = T)
ENCORI_S1PR1 <- read.table("ENCORI_hg19_CLIP-seq_miRNA-target_all_S1PR1.txt", sep = "\t", header = T)
miRWalk <- read.csv("miRWalk_miRNA_Targets.csv", header = T)

#-----------------分别获取两个数据库的miRNA-mRNA pairs
ENCORI <- rbind(ENCORI_AGTR1, ENCORI_CXCL12, ENCORI_PDGFRL, ENCORI_PTGER3, ENCORI_S1PR1)
ENCORI <- ENCORI[c(2,4)]
ENCORI <- ENCORI[!duplicated(ENCORI),]
colnames(ENCORI) <- c("miRNA", "mRNA")   #ENCORI数据库的miRNA-mRNA pairs

miRWalk <- miRWalk[c(1,3)]    
miRWalk <- miRWalk[!duplicated(miRWalk),]
colnames(miRWalk) <- c("miRNA", "mRNA")    #miRWalk数据库的miRNA-mRNA pairs

#------------------miRWalk和ENCORI数据库取交集
miRNA_mRNA <- inner_join(ENCORI, miRWalk)

#------------------output data
write.table(unique(miRNA_mRNA$miRNA), "miRNA.txt", quote = F, row.names = F, col.names = "miRNA")
write.table(ENCORI, "ENCORI.txt", quote = F, row.names = F, sep = "\t")
write.table(miRWalk, "miRWalk.txt", quote = F, row.names = F, sep = "\t")
write.table(miRNA_mRNA, "miRNA_mRNA.txt", quote = F, row.names = F, sep = "\t")
