setwd("F:/gxe_data/YQ180-8/7.GSEA/GO/CC/table")

gene <- c("AGTR1", "CXCL12", "PDGFRL", "PTGER3", "S1PR1")
data <- list()
for(i in 1:length(gene)){
  data[[i]] <- read.table(paste0(gene[i], ".txt"), sep = "\t", header = T)
}

for(i in 1:length(data)){
  data[[i]] <- data[[i]][1]
}

result <- data.frame()
for(i in 1:length(data)){
  a <- data[[i]]
  result <- rbind(result, a)
}
result <- as.data.frame(unique(result$NAME))
colnames(result) <- "CC"

for(j in 1:length(rownames(result))){
  result$AGTR1[j] <- ifelse(result$CC[j] %in% data[[1]][,1], 1, 0)
  result$CXCL12[j] <- ifelse(result$CC[j] %in% data[[2]][,1], 1, 0)
  result$PDGFRL[j] <- ifelse(result$CC[j] %in% data[[3]][,1], 1, 0)
  result$PTGER3[j] <- ifelse(result$CC[j] %in% data[[4]][,1], 1, 0)
  result$S1PR1[j] <- ifelse(result$CC[j] %in% data[[5]][,1], 1, 0)
}

result$sum <- rowSums(result[2:6])
result <- result[order(result$sum, decreasing = T),]
result <- result[result$sum > 2, 1:6]

write.table(result, "result.txt", sep = "\t", row.names = F, quote = F)
