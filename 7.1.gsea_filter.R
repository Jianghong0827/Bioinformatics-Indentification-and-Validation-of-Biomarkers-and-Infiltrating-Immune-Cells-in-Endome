setwd("F:/gxe_data/YQ180-8/7.GSEA/table")

AGTR1_neg <- read.table("gsea_report_for_AGTR1_neg_1634281314878.xls", sep = "\t", header = T)
AGTR1_pos <- read.table("gsea_report_for_AGTR1_pos_1634281314878.xls", sep = "\t", header = T)

CXCL12_neg <- read.table("gsea_report_for_CXCL12_neg_1634282016895.xls", sep = "\t", header = T)
CXCL12_pos <- read.table("gsea_report_for_CXCL12_pos_1634282016895.xls", sep = "\t", header = T)

PDGFRL_neg <- read.table("gsea_report_for_PDGFRL_neg_1634282096919.xls", sep = "\t", header = T)
PDGFRL_pos <- read.table("gsea_report_for_PDGFRL_pos_1634282096919.xls", sep = "\t", header = T)

PTGER3_neg <- read.table("gsea_report_for_PTGER3_neg_1634282578807.xls", sep = "\t", header = T)
PTGER3_pos <- read.table("gsea_report_for_PTGER3_pos_1634282578807.xls", sep = "\t", header = T)

S1PR1_neg <- read.table("gsea_report_for_S1PR1_neg_1634282627975.xls", sep = "\t", header = T)
S1PR1_pos <- read.table("gsea_report_for_S1PR1_pos_1634282627975.xls", sep = "\t", header = T)

gene <- c("AGTR1", "CXCL12", "PDGFRL", "PTGER3", "S1PR1")
data <- list()
data[[1]] <- rbind(AGTR1_neg, AGTR1_pos)
data[[2]] <- rbind(CXCL12_neg, CXCL12_pos)
data[[3]] <- rbind(PDGFRL_neg, PDGFRL_pos)
data[[4]] <- rbind(PTGER3_neg, PTGER3_pos)
data[[5]] <- rbind(S1PR1_neg, S1PR1_pos)

for(i in 1:length(data)){
  data[[i]] <- data[[i]][data[[i]]$NOM.p.val < 0.05,]
  write.table(data[[i]], paste0(gene[i], ".txt"), sep = "\t", quote = F, row.names = F)
}
