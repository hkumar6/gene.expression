library(glmnet)
library(Hmisc)

sample(rownames(mixedSpecies100), 0.5*nrow(mixedSpecies100)) -> selectedGenes
simData <- mixedSpecies100[selectedGenes,]

#use simData for analysis of cells
#use t(simData) for analysis of genes

lasso.mixed.data <- function(ID,simData){
  i = which(colnames(simData)==ID)
  simData <- as.matrix(simData)
  simData.test <- simData[which(simData[,i]==0),,drop=FALSE]
  simData.learn <- simData[-c(which(simData[,i]==0)),,drop=FALSE]
  simData_learn <- simData.learn[,-c(which(colnames(simData.learn)==ID)),drop = FALSE]
  simData_test <- simData.test[,-c(which(colnames(simData.test)==ID)),drop = FALSE]
  vec.learn <- simData.learn[,which(colnames(simData.learn)==ID)]
  vec.test <- simData.test[,which(colnames(simData.test)==ID)]
  n = dim(simData_learn)[1]
  if (dim(simData_learn)[1]>1){
    fit.expr.lin <- glmnet(x = log(1+simData_learn), y = log(1+vec.learn),
                           family = "gaussian", standardize = TRUE, alpha = 1)
  if (sum(simData_learn>0)/(dim(simData_learn)[1]*dim(simData_learn)[2]) > 0.5){
    cv.expr.lin <- cv.glmnet(x = log(1+simData_learn), y = log(1+vec.learn), type.measure="mse")
    prediction <- predict(fit.expr.lin, newx = log(1+simData_test), s = c(cv.expr.lin$lambda.min))
    cov <- sum(coef(fit.expr.lin, cv.expr.lin$lambda.min)>0)
  } else {
    prediction.lin <- predict(fit.expr.lin, newx = log(1+simData_test))
    test <- apply(prediction.lin,2,function(x) mean((vec.test - (exp(x)-1))^2))
    prediction <- prediction.lin[,which(test == min(test))[1]]
    if (is.na(fit.expr.lin$lambda)[1]){
      cov = NaN
    } else {
    cov <- sum(coef(fit.expr.lin, s = which(test == min(test))[1])>0) }
  
  mse.lin <- mean((vec.test - (exp(prediction)-1))^2)
  c <- rcorr((exp(prediction)-1),vec.test, type = "spearman")
  Spear_corr <- c$r[1,2]}
  }else{
      mse.lin = NaN
      Spear_corr = NaN
      cov = NaN
    }
  
  result <- data.frame(mse = mse.lin, Spear_corr = Spear_corr, number.cov = cov, nonzero = n)
}

ID <- rownames(simData)[700:710]

result <- mapply(lasso.mixed.data,ID,MoreArgs = list(t(simData)))





