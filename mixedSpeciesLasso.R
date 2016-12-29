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
  n1 = length(vec.test)
  if (dim(simData_learn)[1]>1){
    x=simData_learn
    fit.expr.lin <- glmnet(x, y = log(1+vec.learn),
                           family = "gaussian", standardize = TRUE, alpha = 1)
  }else{
    mse.lin = NaN
    Spear_corr = NaN
    cov = NaN
  }
  if (n>=9){
    cv.expr.lin <- cv.glmnet(x, y = log(1+vec.learn), nfold=3, type.measure="mse")
    prediction <- predict(fit.expr.lin, newx = simData_test, s = c(cv.expr.lin$lambda.min))
    prediction2 <- exp(prediction)-1
    cov <- sum(coef(fit.expr.lin, cv.expr.lin$lambda.min)!=0)
  } else {
    prediction.lin <- predict(fit.expr.lin, newx = simData_test)
    test <- apply(prediction.lin,2,function(x) mean(((exp(x)-1))^2))
    prediction <- prediction.lin[,which(test == min(test))[1]]}
    if (is.na(fit.expr.lin$lambda)[1]){
      cov = NaN
    } else {
      j  <- which(test == min(test))[1]
    cov <- sum(coef(fit.expr.lin, s = fit.expr.lin$lambda[j])!=0) }
  
  mse.lin <- mean(((exp(prediction)-1))^2)
  
  
  result <- data.frame(mse = mse.lin, number.cov = cov, nonzero = n,length.test = n1)
  return(result)
}

ID <- colnames(simData)[700:710]

result <- mapply(lasso.mixed.data,ID,MoreArgs = list(simData))
result




