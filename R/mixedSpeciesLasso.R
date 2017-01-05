library(glmnet)
library(Hmisc)

# reduce randomly gene dimension by 50%
sample(rownames(mixedSpecies100), 0.5*nrow(mixedSpecies100)) -> selectedGenes
simData <- mixedSpecies100[selectedGenes,]

#use simData for analysis of cells
#use t(simData) for analysis of genes

lasso.mixed.data <- function(ID,simData){
  i = which(colnames(simData)==ID)
  
  #formate for lasso function
  simData <- as.matrix(simData) 
  
  #define test and learn data
  simData.test <- simData[which(simData[,i]==0),,drop=FALSE] # define all zero entries of ID as test
  simData.learn <- simData[-c(which(simData[,i]==0)),,drop=FALSE] # define all non-zero entries as learn
  
  # take out ID from the data set for analysis
  simData_learn <- simData.learn[,-c(which(colnames(simData.learn)==ID)),drop = FALSE] 
  simData_test <- simData.test[,-c(which(colnames(simData.test)==ID)),drop = FALSE]
  
  #define the vector ID for analysis
  vec.learn <- simData.learn[,which(colnames(simData.learn)==ID)]
  vec.test <- simData.test[,which(colnames(simData.test)==ID)]
  n = dim(simData_learn)[1] # number of non-zero entries
  if (dim(simData_learn)[1]>1){ #if there is more than one entry non zero do regression
    x=simData_learn
    fit.expr.lin <- glmnet(x, y = log(1+vec.learn),
                           family = "gaussian", standardize = TRUE, alpha = 1)
  }else{
    mse.lin = NaN
    Spear_corr = NaN
    cov = NaN
  }
  if (n>=9){ # if there are more than 9 non zero entries cross-validation is possible
    cv.expr.lin <- cv.glmnet(x, y = log(1+vec.learn), nfold=3, type.measure="mse")
    prediction <- predict(fit.expr.lin, newx = simData_test, s = c(cv.expr.lin$lambda.min))
    prediction2 <- exp(prediction)-1
    cov <- sum(coef(fit.expr.lin, cv.expr.lin$lambda.min)!=0)
  } else { # if there are less than 9 non zero entries compare results with mse
    prediction.lin <- predict(fit.expr.lin, newx = simData_test)
    test <- apply(prediction.lin,2,function(x) mean(((exp(x)-1))^2))
    prediction <- prediction.lin[,which(test == min(test))[1]]
    if (is.na(fit.expr.lin$lambda)[1]){
      cov = NaN
    } else {
      j  <- which(test == min(test))[1]
    cov <- sum(coef(fit.expr.lin, s = fit.expr.lin$lambda[j])!=0) }}
  
  mse.lin <- mean(((exp(prediction)-1))^2)
  
  
  result <- data.frame(mse = mse.lin, number.cov = cov, nonzero = n)
  return(result)
}


