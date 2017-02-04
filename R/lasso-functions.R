#' lasso imputation - single column
#' 
#' Impute results for a specific gene or cell using lasso
#' 
#' @param ID string for the cell or gene to impute, must be a column name in the data
#' @param simData.learn training data, as matrix
#' @param simData.test testing data, for predictions, as matrix
#' 
#' @return data.frame containing results of the simulation
#'      mean squared error, Spearman correlation, optimal parameters for kknn
#'
#' @importFrom stats coef predict
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom Hmisc rcorr
#' @export
lassoImpute <- function(ID, simData.learn, simData.test, returnPredict = FALSE){
  simData_learn <- simData.learn[,-c(which(colnames(simData.learn) == ID))]
  simData_test <- simData.test[,-c(which(colnames(simData.test) == ID))]
  vec.learn <- simData.learn[,which(colnames(simData.learn) == ID)]
  vec.test <- simData.test[,which(colnames(simData.test) == ID)]
  fit.expr <- glmnet(x = simData_learn, y = vec.learn, family = "gaussian", standardize = TRUE)
  CV.expr <- cv.glmnet(simData_learn,vec.learn, type.measure="mse")
  predicted <- predict(fit.expr, newx = simData_test, s = c(CV.expr$lambda.min))
  mse <- mean((predicted-vec.test)^2)
  c <- rcorr(predicted,vec.test, type = "spearman")
  Spear_corr <- c$r[1,2]
  cov <- sum(coef(fit.expr, CV.expr$lambda.min)>0)
  result <- data.frame(mse = mse, Spear_corr = Spear_corr, number.cov = cov)
  if(returnPredict) {
    return(predicted)
  }
  return(result)
}

#' lasso imputation - multiple columns
#' 
#' Impute results for multiple genes or cells using lasso
#' 
#' @param selectedList vector of strings for the cells or genes to impute, must be a column name in the data
#' @param simData.learn training data, as matrix
#' @param simData.test testing data, for predictions, as matrix
#' 
#' @return data.frame containing results of the simulation
#'      mean squared error, Spearman correlation, optimal parameters for kknn
#'
#' @importFrom stats coef predict
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom Hmisc rcorr
#' @export
lassoImputeAll <- function(selectedList, simData.learn, simData.test) {
  simulationResult <- mapply(lassoImpute, selectedList,
                             MoreArgs = list(
                               simData.learn,
                               simData.test
                             ))
  return(simulationResult)
}