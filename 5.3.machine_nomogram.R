#----------------20211014
setwd("F:/gxe_data/YQ180-8/5.machine learning/1.random_forest")

#----------------loading packages
library(tidyverse)
library(caret)
library(DALEX)
library(tibble)
library(tidyr)

#----------------input data
disease <- read.csv("expr_keygene.csv",header = T)
disease$group <- factor(disease$group)

set.seed(678)
#训练集和验证集
trainIndex <- createDataPartition(disease$group, p = 0.6, list = FALSE, times = 1)
diseaseTrain <- disease[ trainIndex,]
diseaseTest  <- disease[-trainIndex,]

#建立模型
classif_rf <- train(group~., data = diseaseTrain, method="rf", ntree = 100, tuneLength = 1)
classif_glm <- train(group~., data = diseaseTrain, method="glm", family="binomial")
classif_svm <- train(group~., data = diseaseTrain, method="svmRadial", prob.model = TRUE, tuneLength = 1)

#定义求得似然估计的函数
p_fun <- function(object, newdata){predict(object, newdata=newdata, type="prob")[,2]} 

#把测试集的响应变量转化为数值
yTest <- as.numeric(as.character(diseaseTest$group))

#对模型进行解释
explainer_classif_rf <- DALEX::explain(classif_rf, label = "rf",
                                       data = diseaseTest, y = yTest,
                                       predict_function = p_fun)


explainer_classif_glm <- DALEX::explain(classif_glm, label = "glm",  
                                        data = diseaseTest, y = yTest, 
                                        predict_function = p_fun)


explainer_classif_svm <- DALEX::explain(classif_svm,  label = "svm", 
                                        data = diseaseTest, y = yTest,
                                        predict_function = p_fun)

#模型表现
mp_classif_rf <- model_performance(explainer_classif_rf)
mp_classif_glm <- model_performance(explainer_classif_glm)
mp_classif_svm <- model_performance(explainer_classif_svm)

#累计残差分布图
pdf(file = "Cumulative residual distribution.pdf", width = 9, height = 8)
plot(mp_classif_rf, mp_classif_glm, mp_classif_svm)
dev.off()

#累计残差箱线图
pdf(file = "boxplot.pdf",width = 9,height = 7)
plot(mp_classif_rf, mp_classif_glm, mp_classif_svm, geom = "boxplot")
dev.off()

#变量重要性分析
vi_classif_rf <- variable_importance(explainer_classif_rf, loss_function = loss_root_mean_square)
vi_classif_glm <- variable_importance(explainer_classif_glm, loss_function = loss_root_mean_square)
vi_classif_svm <- variable_importance(explainer_classif_svm, loss_function = loss_root_mean_square)

pdf(file = "Importance of variables.pdf",width = 9,height = 8.5)
plot(vi_classif_rf, vi_classif_glm, vi_classif_svm)
dev.off()

pdf(file = "svm.pdf",width = 9,height = 8.5)
plot(vi_classif_svm)
dev.off()
