rm(list = ls())
setwd("~/help_for_others/DYQ/output/1_preclean/data//")
library(stringr)
library(gridExtra)
library(future)
library(sva)
# library(e1071)
library(pROC)
library(ROCit)
library(caret)
library(doParallel)
library(cancerclass)




data <- read.csv("./C1_df.csv",check.names = FALSE)
set.seed(13942)
trainIndex <- createDataPartition(data$R, p = .7,  ## 80% training set; 20% validation set
                                  list = FALSE, 
                                  times = 1)

training <- data[trainIndex,]
validation  <- data[-trainIndex,]

table(training$R)
table(validation$R)

write.csv(training,file = "~/help_for_others/DYQ/output/3_baseline_model/csv/C1D1_response_df_training.csv",row.names = FALSE)
write.csv(validation,file = "~/help_for_others/DYQ/output/3_baseline_model/csv//C1D1_response_df_validation.csv",row.names = FALSE)
