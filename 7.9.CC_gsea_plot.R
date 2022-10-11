setwd("F:/gxe_data/YQ180-8/7.GSEA/GO/CC/table")

result <- read.table("result.txt", sep = "\t", row.names = 1, header = T)
result <- result[rowSums(result) == 5,]


library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(1, 0), c("darkblue", "gray"))
pdf("heatmap.pdf", width = 10,height = 7)
Heatmap(as.matrix(result[1:10,]), 
        cluster_rows = F,
        cluster_columns = F,
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 48),
        column_names_rot = 0,
        row_names_side = "left",
        show_heatmap_legend = F,
        row_names_gp = gpar(cex = 0.5))
dev.off()
