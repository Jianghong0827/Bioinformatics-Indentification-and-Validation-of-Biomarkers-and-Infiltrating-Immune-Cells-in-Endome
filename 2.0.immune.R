#--------------------20210922
setwd("F:/gxe_data/YQ180-8/2.immune")

#--------------------loading packages
source("CIBERSORT.R")
library("e1071")

#--------------------CIBERSORT
result <- CIBERSORT("LM22.txt", "normalize.txt", perm = 1000, QN = T)
###本地获得results文件，即每个样本中22种免疫细胞的比例
