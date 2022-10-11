#------------------20211020
setwd("F:/gxe_data/YQ180-8/9.drug")

#------------------input data
drug <- read.csv("dgidb_export_2021-10-20.tsv", sep = "\t", header = T)
ceRNA <- read.table("ceRNA.txt", sep = "\t", header = T)

#------------------网络属性文件
drug <- drug[c(1,5)]
drug <- drug[!duplicated(drug),]
write.table(drug, "drug_mRNA.txt", sep = "\t", row.names = F, quote = F)
CXCL12 <- drug[drug$search_term == "CXCL12",]
AGTR1 <- drug[drug$search_term == "AGTR1",]
PTGER3 <- drug[drug$search_term == "PTGER3",]
S1PR1 <- drug[drug$search_term == "S1PR1",]
write.table(CXCL12, "CXCL12-drug.txt", sep = "\t", row.names = F, quote = F)
write.table(AGTR1, "AGTR1-drug.txt", sep = "\t", row.names = F, quote = F)
write.table(PTGER3, "PTGER3-drug.txt", sep = "\t", row.names = F, quote = F)
write.table(S1PR1, "S1PR1-drug.txt", sep = "\t", row.names = F, quote = F)


ceRNA <- ceRNA[ceRNA$others != "PDGFRL",]

colnames(drug) <- colnames(ceRNA)
network <- rbind(ceRNA, drug)   #网络文件

#网络属性文件
colnames(drug) <- c("mRNA", "drug")
drug <- data.frame(unique(drug$drug))
drug$text <- rep("drug", length(rownames(drug)))

#------------------output data
write.table(ceRNA, "ceRNA_new.txt", sep = "\t", row.names = F, quote = F)
write.table(network, "network.txt", sep = "\t", row.names = F, quote = F)
write.table(drug, "drug.txt", sep = "\t", row.names = F, quote = F)


