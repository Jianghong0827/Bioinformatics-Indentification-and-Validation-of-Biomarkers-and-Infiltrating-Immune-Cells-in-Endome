#--------------------20210922
setwd("F:/gxe_data/YQ180-8/2.immune")

#--------------------loading packages
source("CIBERSORT.R")
library("e1071")

#--------------------CIBERSORT
result <- CIBERSORT("LM22.txt", "normalize.txt", perm = 1000, QN = T)
###���ػ��results�ļ�����ÿ��������22������ϸ���ı���