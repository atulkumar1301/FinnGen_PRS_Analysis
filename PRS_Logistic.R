#! /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/bin/Rscript
library(ROCR)
library(pROC)
library(bdpv)
library(data.table)
library(dplyr)
library(caret)
library(Rcpp)
args <- commandArgs(trailingOnly = TRUE)
TABLE<-as.data.frame(matrix(ncol=13, nrow=7)) 
names(TABLE)<-c("Model", "Effect", "SE", "R2", "L95", "U95", "AIC", "BIC", "Sensitivity","Specificity", "AUC", "Number_of_Variants", "P")
Clinical_data <- fread (file = paste0 ("/Volumes/ATUL_6TB/Work/Projects/Cell_Specific_PRSs/Full_Data/", args[1])) #Give the path for clinical data file
Full_File <- fread ("/Volumes/ATUL_6TB/Work/Projects/Cell_Specific_PRSs/Full_Data/PCA_FILE.eigenvec") %>% merge(Clinical_data, by = 'IID')
PRS_1 <- fread ("p_value_0.05.sscore") %>% merge(Full_File, by = 'IID')
PRS_2 <- fread ("p_value_0.005.sscore") %>% merge(Full_File, by = 'IID')
PRS_3 <- fread ("p_value_0.0005.sscore") %>% merge(Full_File, by = 'IID')
PRS_4 <- fread ("p_value_5e-05.sscore") %>% merge(Full_File, by = 'IID')
PRS_5 <- fread ("p_value_5e-06.sscore") %>% merge(Full_File, by = 'IID')
PRS_6 <- fread ("p_value_5e-07.sscore") %>% merge(Full_File, by = 'IID')
PRS_7 <- fread ("p_value_5e-08.sscore") %>% merge(Full_File, by = 'IID')

i = 1

data2 <- PRS_1
#data2 <- subset (PRS_1, PRS_1$E4! = 0)
data2$PHENO <- data2$Diagnosis
TP = nrow (subset (data2, data2$PHENO == "1"))
TN = nrow (subset (data2, data2$PHENO == "0"))
TABLE [i, 1] <- "PRS 1"

m1 <- mean (data2$SCORE1_AVG)
sd1 <- sd (data2$SCORE1_AVG)
data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1

modeldata20 <- glm (data2$PHENO ~ 1, family = binomial, data = data2)
model_PRS1 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10, family = binomial, data = data2)


l0 <- deviance (modeldata20)
df0 <- df.residual (modeldata20) #degree of freedom
l1 <- deviance (model_PRS1)
df1 <- df.residual (model_PRS1)

TABLE[i, 2] <- summary (model_PRS1)$coefficients[2, "Estimate"]
TABLE[i, 3] <- summary(model_PRS1)$coefficients[2, "Std. Error"]

TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
TABLE[i, 5] <- confint (model_PRS1) [2, 1] 
TABLE[i, 6] <- confint (model_PRS1) [2, 2]
TABLE[i, 7] <- AIC (model_PRS1)
TABLE[i, 8] <- BIC (model_PRS1)

obs = model_PRS1$y
pred = as.numeric (model_PRS1$fitted.values > 0.5)
#mat <- confusionMatrix(data = as.factor(pred), reference = as.factor (obs))
TABLE[i, 9] <- sensitivity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 10] <- specificity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 11] <- auc (obs, pred)
TABLE[i, 12] <- nrow (fread ("p_value_0.05.txt"))
TABLE[i, 13] <- summary (model_PRS1)$coefficients [2, "Pr(>|z|)"]

#=========================================================================================================

i = 2

data2 <- PRS_2
#data2<- subset (PRS_2, PRS_2$E4! = 0)
data2$PHENO <- data2$Diagnosis
TP = nrow (subset (data2, data2$PHENO == "1"))
TN = nrow (subset (data2, data2$PHENO == "0"))
TABLE [i, 1] <- "PRS 2"

m1 <- mean (data2$SCORE1_AVG)
sd1 <- sd (data2$SCORE1_AVG)
data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1

modeldata20 <- glm (data2$PHENO ~ 1, family = binomial, data = data2)
model_PRS2 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10, family = binomial, data = data2)


l0 <- deviance (modeldata20)
df0 <- df.residual (modeldata20) #degree of freedom
l1 <- deviance (model_PRS2)
df1 <- df.residual (model_PRS2)

TABLE[i, 2] <- summary (model_PRS2)$coefficients[2, "Estimate"]
TABLE[i, 3] <- summary(model_PRS2)$coefficients[2, "Std. Error"]

TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
TABLE[i, 5] <- confint (model_PRS2) [2, 1] 
TABLE[i, 6] <- confint (model_PRS2) [2, 2]
TABLE[i, 7] <- AIC (model_PRS2)
TABLE[i, 8] <- BIC (model_PRS2)

obs = model_PRS2$y
pred = as.numeric (model_PRS2$fitted.values > 0.5)
#mat <- confusionMatrix(data = as.factor(pred), reference = as.factor (obs))
TABLE[i, 9] <- sensitivity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 10] <- specificity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 11] <- auc (obs, pred)
TABLE[i, 12] <- nrow (fread ("p_value_0.005.txt"))
TABLE[i, 13] <- summary (model_PRS2)$coefficients [2, "Pr(>|z|)"]

#================================================================================================

i = 3

data2 <- PRS_3
#data2<- subset (PRS_2, PRS_2$E4! = 0)
data2$PHENO <- data2$Diagnosis
TP = nrow (subset (data2, data2$PHENO == "1"))
TN = nrow (subset (data2, data2$PHENO == "0"))
TABLE [i, 1] <- "PRS 3"

m1 <- mean (data2$SCORE1_AVG)
sd1 <- sd (data2$SCORE1_AVG)
data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1

modeldata20 <- glm (data2$PHENO ~ 1, family = binomial, data = data2)
model_PRS3 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10, family = binomial, data = data2)


l0 <- deviance (modeldata20)
df0 <- df.residual (modeldata20) #degree of freedom
l1 <- deviance (model_PRS3)
df1 <- df.residual (model_PRS3)

TABLE[i, 2] <- summary (model_PRS3)$coefficients[2, "Estimate"]
TABLE[i, 3] <- summary(model_PRS3)$coefficients[2, "Std. Error"]

TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
TABLE[i, 5] <- confint (model_PRS3) [2, 1] 
TABLE[i, 6] <- confint (model_PRS3) [2, 2]
TABLE[i, 7] <- AIC (model_PRS3)
TABLE[i, 8] <- BIC (model_PRS3)

obs = model_PRS3$y
pred = as.numeric (model_PRS3$fitted.values > 0.5)
#mat <- confusionMatrix(data = as.factor(pred), reference = as.factor (obs))
TABLE[i, 9] <- sensitivity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 10] <- specificity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 11] <- auc (obs, pred)
TABLE[i, 12] <- nrow (fread ("p_value_0.0005.txt"))
TABLE[i, 13] <- summary (model_PRS3)$coefficients [2, "Pr(>|z|)"]

#=============================================================================================

i = 4

data2 <- PRS_4
#data2<- subset (PRS_2, PRS_2$E4! = 0)
data2$PHENO <- data2$Diagnosis
TP = nrow (subset (data2, data2$PHENO == "1"))
TN = nrow (subset (data2, data2$PHENO == "0"))
TABLE [i, 1] <- "PRS 4"

m1 <- mean (data2$SCORE1_AVG)
sd1 <- sd (data2$SCORE1_AVG)
data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1

modeldata20 <- glm (data2$PHENO ~ 1, family = binomial, data = data2)
model_PRS4 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10, family = binomial, data = data2)


l0 <- deviance (modeldata20)
df0 <- df.residual (modeldata20) #degree of freedom
l1 <- deviance (model_PRS4)
df1 <- df.residual (model_PRS4)

TABLE[i, 2] <- summary (model_PRS4)$coefficients[2, "Estimate"]
TABLE[i, 3] <- summary(model_PRS4)$coefficients[2, "Std. Error"]

TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
TABLE[i, 5] <- confint (model_PRS4) [2, 1] 
TABLE[i, 6] <- confint (model_PRS4) [2, 2]
TABLE[i, 7] <- AIC (model_PRS4)
TABLE[i, 8] <- BIC (model_PRS4)

obs = model_PRS4$y
pred = as.numeric (model_PRS4$fitted.values > 0.5)
#mat <- confusionMatrix(data = as.factor(pred), reference = as.factor (obs))
TABLE[i, 9] <- sensitivity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 10] <- specificity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 11] <- auc (obs, pred)
TABLE[i, 12] <- nrow (fread ("p_value_5e-05.txt"))
TABLE[i, 13] <- summary (model_PRS4)$coefficients [2, "Pr(>|z|)"]

#====================================================================================================

i = 5

data2 <- PRS_5
#data2<- subset (PRS_2, PRS_2$E4! = 0)
data2$PHENO <- data2$Diagnosis
TP = nrow (subset (data2, data2$PHENO == "1"))
TN = nrow (subset (data2, data2$PHENO == "0"))
TABLE [i, 1] <- "PRS 5"

m1 <- mean (data2$SCORE1_AVG)
sd1 <- sd (data2$SCORE1_AVG)
data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1

modeldata20 <- glm (data2$PHENO ~ 1, family = binomial, data = data2)
model_PRS5 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10, family = binomial, data = data2)


l0 <- deviance (modeldata20)
df0 <- df.residual (modeldata20) #degree of freedom
l1 <- deviance (model_PRS5)
df1 <- df.residual (model_PRS5)

TABLE[i, 2] <- summary (model_PRS5)$coefficients[2, "Estimate"]
TABLE[i, 3] <- summary(model_PRS5)$coefficients[2, "Std. Error"]

TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
TABLE[i, 5] <- confint (model_PRS5) [2, 1] 
TABLE[i, 6] <- confint (model_PRS5) [2, 2]
TABLE[i, 7] <- AIC (model_PRS5)
TABLE[i, 8] <- BIC (model_PRS5)

obs = model_PRS5$y
pred = as.numeric (model_PRS5$fitted.values > 0.5)
#mat <- confusionMatrix(data = as.factor(pred), reference = as.factor (obs))
TABLE[i, 9] <- sensitivity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 10] <- specificity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 11] <- auc (obs, pred)
TABLE[i, 12] <- nrow (fread ("p_value_5e-06.txt"))
TABLE[i, 13] <- summary (model_PRS5)$coefficients [2, "Pr(>|z|)"]

#===================================================================================================

i = 6

data2 <- PRS_6
#data2<- subset (PRS_2, PRS_2$E4! = 0)
data2$PHENO <- data2$Diagnosis
TP = nrow (subset (data2, data2$PHENO == "1"))
TN = nrow (subset (data2, data2$PHENO == "0"))
TABLE [i, 1] <- "PRS 6"

m1 <- mean (data2$SCORE1_AVG)
sd1 <- sd (data2$SCORE1_AVG)
data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1

modeldata20 <- glm (data2$PHENO ~ 1, family = binomial, data = data2)
model_PRS6 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10, family = binomial, data = data2)


l0 <- deviance (modeldata20)
df0 <- df.residual (modeldata20) #degree of freedom
l1 <- deviance (model_PRS6)
df1 <- df.residual (model_PRS6)

TABLE[i, 2] <- summary (model_PRS6)$coefficients[2, "Estimate"]
TABLE[i, 3] <- summary(model_PRS6)$coefficients[2, "Std. Error"]

TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
TABLE[i, 5] <- confint (model_PRS6) [2, 1] 
TABLE[i, 6] <- confint (model_PRS6) [2, 2]
TABLE[i, 7] <- AIC (model_PRS6)
TABLE[i, 8] <- BIC (model_PRS6)

obs = model_PRS6$y
pred = as.numeric (model_PRS6$fitted.values > 0.5)
#mat <- confusionMatrix(data = as.factor(pred), reference = as.factor (obs))
TABLE[i, 9] <- sensitivity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 10] <- specificity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 11] <- auc (obs, pred)
TABLE[i, 12] <- nrow (fread ("p_value_5e-07.txt"))
TABLE[i, 13] <- summary (model_PRS6)$coefficients [2, "Pr(>|z|)"]

#=========================================================================================================

i = 7

data2 <- PRS_7
#data2<- subset (PRS_2, PRS_2$E4! = 0)
data2$PHENO <- data2$Diagnosis
TP = nrow (subset (data2, data2$PHENO == "1"))
TN = nrow (subset (data2, data2$PHENO == "0"))
TABLE [i, 1] <- "PRS 7"

m1 <- mean (data2$SCORE1_AVG)
sd1 <- sd (data2$SCORE1_AVG)
data2$NORMSCORE <- (data2$SCORE1_AVG - m1) / sd1

modeldata20 <- glm (data2$PHENO ~ 1, family = binomial, data = data2)
model_PRS6 <- glm (data2$PHENO ~ data2$NORMSCORE + data2$Age + data2$Gender + data2$PC1 + data2$PC2 + data2$PC3 + data2$PC4 + data2$PC5 + data2$PC6 + data2$PC7 + data2$PC8 + data2$PC9 + data2$PC10, family = binomial, data = data2)


l0 <- deviance (modeldata20)
df0 <- df.residual (modeldata20) #degree of freedom
l1 <- deviance (model_PRS6)
df1 <- df.residual (model_PRS6)

TABLE[i, 2] <- summary (model_PRS6)$coefficients[2, "Estimate"]
TABLE[i, 3] <- summary(model_PRS6)$coefficients[2, "Std. Error"]

TABLE[i, 4] <- (1 - exp ((l1 - l0)/nrow(data2)))/(1 - exp (-l0/nrow (data2)))        #Nagelkerke
TABLE[i, 5] <- confint (model_PRS6) [2, 1] 
TABLE[i, 6] <- confint (model_PRS6) [2, 2]
TABLE[i, 7] <- AIC (model_PRS6)
TABLE[i, 8] <- BIC (model_PRS6)

obs = model_PRS6$y
pred = as.numeric (model_PRS6$fitted.values > 0.5)
#mat <- confusionMatrix(data = as.factor(pred), reference = as.factor (obs))
TABLE[i, 9] <- sensitivity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 10] <- specificity (data = as.factor (pred), reference = as.factor (obs))
TABLE[i, 11] <- auc (obs, pred)
TABLE[i, 12] <- nrow (fread ("p_value_5e-08.txt"))
TABLE[i, 13] <- summary (model_PRS6)$coefficients [2, "Pr(>|z|)"]

#=======================================================================================================

TABLE$Bonferroni <- p.adjust(TABLE$P, method = "bonferroni", n = length(TABLE$P))
TABLE$FDR <- p.adjust(TABLE$P, method = "fdr", n = length(TABLE$P))
write.table (TABLE, file = paste0 ("PRS_Result_FinnGen.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
