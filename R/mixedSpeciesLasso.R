#library(glmnet)
#library(Hmisc)
#library(mpath)
#library(kknn)
#library(randomForest)

#' lasso imputation - single column
#' 
#' Impute results for a specific gene or cell using lasso
#' 
#' @param ID string for the cell or gene to impute, must be a column name in the data
#' @param simData data
#' @param genes TRUE/FALSE depending on wether genes are analysed
#' 
#' @return data.frame containing results of the simulation
#'      mean error, percentage of humans
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @export



#use simData for analysis of cells
#use t(simData) for analysis of genes

lasso.mixed.data <- function(ID,simData, genes = TRUE){
  i = which(colnames(simData)==ID)
  
  #formate for lasso function
  simData <- as.matrix((simData)) 
  
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
  if (n > 2){ #if there is more than one entry non zero do regression
    fit.expr.lin <- glmnet(x = log(1+simData_learn), y = log(1+vec.learn),
                          family = "gaussian", standardize = TRUE, alpha = 1)
  
    # if there are more than 9 non zero entries cross-validation is possible
    #cv.expr.lin <- cv.glmnet(x=log(1+simData_learn), y = log(1+vec.learn), nfold=5, type.measure="mse")
    if (genes){
      lambda = 0.0058
    }else{
      lambda = 0.0185  
    }
    prediction <- predict(fit.expr.lin, newx = log(1+simData_test), s = lambda) 
    prediction <- exp(prediction)-1
    coef <- sum(coef(fit.expr.lin, lambda)!=0)

    if (length(grep("HUMAN",colnames(simData)[i]))>0){
      test <- colnames(simData)[which(coef(fit.expr.lin, lambda)!=0)]
      H.coef <- grep("HUMAN",test)
      coef.perc <- length(H.coef)/coef
      M.test <- grep("MOUSE",rownames(prediction))
      prediction.new <- prediction[M.test]
      counts <- sum(prediction.new == 0)
      mse.lin <- mean(prediction.new)
    }else{
      test <- colnames(simData)[which(coef(fit.expr.lin, lambda)!=0)]
      H.coef <- grep("HUMAN",test)
      coef.perc <- length(H.coef)/coef
      H.test <- grep("HUMAN",rownames(prediction))
      prediction.new <- prediction[H.test,drop=FALSE]
      counts <- sum(prediction.new == 0)
      mse.lin <- mean(prediction.new)
    }
  
  }else{
    mse.lin = NaN
    coef = NaN
    coef.perc = NaN
    counts = NaN
  }
  
  result <- data.frame(mean.error = mse.lin, number.coef = coef, nonzero = n, perc.human = coef.perc,
                       zero.predicted = counts)
  return(result)
}

#resultgenes <- mapply(lasso.mixed.data, sample(colnames(simData),50), MoreArgs = list((simData)))
#resultgenes <- mapply(lasso.mixed.data,colnames(simData)[1:10], MoreArgs = list((simData), genes=FALSE))
#0.0058 for genes
#0.0185 for cells

