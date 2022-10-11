#----------------20211014
setwd("F:/gxe_data/YQ180-8/5.machine learning/1.random_forest")

#----------------loading packages
library(tidyverse)
library(caret)
library(DALEX)
library(tibble)
library(tidyr)

##�������ݲ�����ѵ�����Ͳ��Լ�

#airquality <- read.table("expr_keygene.txt",header = T)
#airquality$group <- ifelse(airquality$group == "EU", "0", "1")
#write.csv(airquality, "expr_keygene.csv", quote = F, row.names = F)

airquality <- read.csv("expr_keygene.csv",header = T)
airquality <- as_tibble(airquality)  #���ݸ�ʽת��
airquality

set.seed(114)
index <- createDataPartition(
  airquality$group,
  p = 0.7,
  list = FALSE
)
airqualitytrain <- airquality[index, ]
airqualitytest <- airquality[-index, ]

##ʹ��caret������ģ��
ozone_rf <- train(group~.,data = airqualitytrain,
                method = "rf",
                ntree = 20)
ozone_glm <- train(group~.,data = airqualitytrain,
                 method = "glm")
ozone_svm <- train(group~.,data = airqualitytrain,
                 method = "svmLinear")
##��ģ�ͽ��н���

explainer_rf <- explain(ozone_rf,label = "rf",
                      data = airqualitytest,
                      y = airqualitytest$group)

explainer_glm <- explain(ozone_glm,label = "glm",
                       data = airqualitytest,
                       y = airqualitytest$group)
explainer_svm <- explain(ozone_svm,label = "svmLinear",
                       data = airqualitytest,
                       y = airqualitytest$group)

##ģ�ͱ���

per_rf <- model_performance(explainer_rf)
per_glm <- model_performance(explainer_glm)
per_svm <- model_performance(explainer_svm)

#�ۻ��в�ֲ�ͼ

pdf(file = "Cumulative residual distribution3.pdf",width = 9,height = 8)

plot(per_rf, per_glm, per_svm)

# a <- per_rf
# a$residuals$diff <- a$residuals$diff * 20
# 
# b <- per_svm
# b$residuals$diff <- b$residuals$diff * 20
# 
# c <- per_glm
# c$residuals$diff <- c$residuals$diff * 20
# 
# plot(a,b,c)

dev.off()


#����ͼ�ֲ�ͼ
pdf(file = "boxplot3.pdf",width = 9,height = 7)
plot(per_rf, per_glm, per_svm, geom = "boxplot")
dev.off()

#������Ҫ�Է���
importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)


pdf(file = "Importance of variables3.pdf",width = 9,height = 8.5)
plot(importance_rf,importance_glm,importance_svm)
dev.off()