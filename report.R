#' ---
#' title: "analysis process"
#' author: '1'
#' date: "`r Sys.Date()`"
#' output:
#'   word_document:
#'     toc: yes
#'     # reference_docx: wordstyle.docx
#' ---
#' 
## ----setup, include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(echo = T,warning = F,message = F,fig.width = 6,fig.height = 4,dpi = 250)

#' 
#' - discription of variables:
#' 
#' - age: age of patient anaemia: diagnosis of anaemie (a decrease in total
#' red blood cells or hemoglobin in blood) creatinine_phosphokinase: blood
#' concentration of enzyme CPK (mcg/L) diabetes: diagnosis of diabetes
#' - ejection_fraction: proportion of blood leaving the heart on each
#' contraction (%) high_blood_pressure: diagnosis of hypertension
#' - platelets: blood concentration of platelets (kiloplatelets/mL)
#' - serum_creatinine: blood concentration of serum creatinine (mg/dL).
#' Normal levels approx 0.6 to 1.35. Elevation may imply poor kidney
#' function. 
#' - serum_sodium: blood concentration of serum sodium (mEq/L).
#' Normal levels approx 135 to 145. Hyponatremia risk at low levels. 
#' - sex:male=1, female=0 smoking: patient is a smoker 
#' - time: follow-up period
#' (days) 
#' - fatal_mi: whether the patient suffered a fatal myocardial
#' infarction during the follow-up period
#' 
#' ## Problem description
#' 
#' ### read data
#' 
## -----------------------------------------------------------------------------------
library(tidyverse)
library(flextable)
library(caret)
rawdata=read_csv("heart_failure.csv")
heart=rawdata %>% select(-time)#Time is not a good predictor

#' 
#' ### characterisation of the distribution of variables：
#' 
#' 
## -----------------------------------------------------------------------------------
heart %>% pivot_longer(col=everything()) %>% ggplot(aes(x=value,fill=name,col=1))+
  geom_histogram(show.legend = F)+facet_wrap(vars(name),scales = "free")+labs(x="")+theme_minimal()

#' 
#' ### Variable correlation
#' 
## -----------------------------------------------------------------------------------
library(ggcorrplot)
heart%>% cor() %>%  ggcorrplot(type = "lower",
  insig = "blank",
  lab = TRUE,
  digits = 2,tl.cex = 8,lab_size = 2,  title = "Correlation coefficient plot")

#' 
#' ## Variable screening (analysis of influencing factors)
#' 
#' Because of the presence of highly correlated and non-significant variables, directly fitted logistic regression is not suitable.
#' 
#' ### glm
#' Logistic regression has good interpretability and is suitable for this problem study
#' Direct fit, correlation of variables affects interpretation of results, consider screening variables
#' 
## -----------------------------------------------------------------------------------
X=as.matrix(heart[,-12])
y=heart$fatal_mi
mod_log=glm(fatal_mi~.,family = "binomial",heart)

#' 
#' Logistic regression direct fit results
## -----------------------------------------------------------------------------------
mod_log %>% as_flextable()
prob<-predict(object =mod_log,newdata=heart,type = "response")
pred<-ifelse(prob>=0.5,1,0)

#' Confusion matrix in raw data
## -----------------------------------------------------------------------------------
table(heart$fatal_mi,pred)

#' 
#' 
#' ### PCA
## -----------------------------------------------------------------------------------
set.seed(123)
library(factoextra)
heart.pr <- prcomp(heart[,-12],scale = TRUE)
fviz_eig(heart.pr)
summary(heart.pr)

#' So we're looking at eight principal components with 80 per cent of the variance explained.
#' 
## -----------------------------------------------------------------------------------
pcdat=as.data.frame(cbind(heart.pr$x[,1:8],fatal_mi=heart$fatal_mi))

#' 
#' A plot of the coefficients of the principal components, which allows you to see the relationship between each component and the original variable
## -----------------------------------------------------------------------------------
library(pheatmap)
heart.pr$rotation %>% pheatmap(cluster_row=F,cluster_col=F)

#' 
#' 
#' Principal component regression results
## -----------------------------------------------------------------------------------
mod_log=glm(fatal_mi~.,family = "binomial",pcdat)
mod_log%>% as_flextable()

#' PC1 is a significant component but not significant. The principal components are considered to be generally effective and missing explanatory and significant variables.
#' 
## -----------------------------------------------------------------------------------
prob<-predict(object =mod_log,newdata=pcdat,type = "response")
pred<-ifelse(prob>=0.5,1,0)
table(heart$fatal_mi,pred)

#' 
#' ### lasso
#' 
#' Try lasso screening variables so that the significant variables of the model are interpretable and effective in disease prevention.。
#' 
## -----------------------------------------------------------------------------------
library(glmnet)
fit.lasso= glmnet(X,y,alpha=1)
plot(fit.lasso, xvar="lambda", cex.axis=2, cex.lab=1.1, cex=1.1)

#' 
#' Cross-validation to find the most suitable lambda
## -----------------------------------------------------------------------------------
set.seed(2)
cv.lasso = cv.glmnet(X,y,alpha=1)
plot(cv.lasso)
lambda <- cv.lasso$lambda.min 
log(lambda)#min mean error

#' 
#' Screened model coefficients
## -----------------------------------------------------------------------------------
coef(cv.lasso, s="lambda.min")
non_0_ind <- which(coef(cv.lasso, s = lambda) != 0) - 1
v=colnames(X)[non_0_ind]#find coef !=0 var

#' 
#' 
#' ## Cross-validation and model optimisation
#' 
#' Further, a logistic regression model and a new Gradient Boosting model are fitted on the basis of screening variables .
#' To improve the generalisation and prediction ability of the model, 3/4 samples are used for the training set and 1/4 for the prediction set.
#' And a 5-fold cross-validation with 2 repetitions was performed to fit the model parameters during training.
#' 
#' ### glm
## -----------------------------------------------------------------------------------
library(caret)
y=as.factor(y)
set.seed(111)
tr=createDataPartition(y,p=3/4)
xtrain=X[tr$Resample1,v]
ytrain=y[tr$Resample1]
xtest=X[-tr$Resample1,v]
ytest=y[-tr$Resample1]
#交叉验证等
ctrl=trainControl(method = "repeatedcv", number = 5, repeats=2)
cv.glm=train(x=xtrain,y=ytrain,method='glm',
             trControl=ctrl)
pred.glm=predict(cv.glm,newdata = xtest,type="prob")
pred<-ifelse(pred.glm[,2]>=0.5,1,0)
table(ytest,pred)

#' 
#' Final glm result:
## -----------------------------------------------------------------------------------
cv.glm$finalModel %>% as_flextable()

#' 
#' ### Gradient Boosting Model parameter selection and optimisation
#' 
#' Tuning parameter 'n.minobsinnode' was held constant at a value of 10
#' Accuracy was used to select the optimal model using the largest value.
#' The final values used for the model were n.trees = 150, interaction.depth
#'  = 1, shrinkage = 0.1 and n.minobsinnode = 10.
#'  
## -----------------------------------------------------------------------------------
cv.gbm=train(x=xtrain,y=ytrain,method='gbm',trControl=ctrl,verbose=F)
plot(cv.gbm)
cv.gbm

#' 
#' Forecast results and confusion matrix
## -----------------------------------------------------------------------------------
pred.gbm=predict(cv.gbm,newdata = xtest,type="prob")
pred<-ifelse(pred.gbm[,2]>=0.5,1,0)
library(gt)
table(ytest,pred)

#' 
#' 
#' ### Model ROC Comparison
## -----------------------------------------------------------------------------------
ROC=data_frame(outcome=ytest,glm=pred.glm[,2],gbm=pred.gbm[,2])
library(pROC)
## Calculation of roc, can batch calculate a, b, c three sets of data at a time
res<-roc(outcome~glm+gbm,data=ROC,aur=TRUE,
         ci=TRUE, 
         )
ggroc(res, legacy.axes = TRUE)+
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
    theme_bw()+ggtitle('ROC plot')+ggsci::scale_color_lancet()+
    annotate("text",x=0.75,y=0.125,label=paste("glm-AUC = ", round(res$glm$auc,3)))+
    annotate("text",x=0.75,y=0.25,label=paste("gbm-AUC = ", round(res$gbm$auc,3)))

#' 
#' 
#' AUC is better in comparison to gbm. We are more concerned with identifying the deaths of patients so that we can do better prevention, so we choose a higher threshold when applying the model.
#' 
#' 
#' 
#' <!-- knitr::purl("report.Rmd","report.R", documentation = 2) -->
#' 
#' 
#' 
