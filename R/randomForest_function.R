#' kknn imputation - single column
#' 
#' Impute results for a specific gene or cell using randomForest
#' 
#' @param id string for the cell or gene to impute, must be a column name in the data
#' @param simData.learn training data as data.frame
#' @param simData.test testing data, for predictions, as data.frame
#' 
#' @return data.frame containing results of the simulation
#'      mean squared error, Spearman correlation
#'
#' @importFrom stats as.formula
#' @importFrom randomForest randomForest
#' @importFrom Hmisc rcorr
#' @export
randomForestImpute <- function(id, simData.learn, simData.test, mixedSpeciesData = FALSE) {
  takeout<-which(colnames(simData.learn)==id)
  # r <- randomForest(as.formula(paste(id, ".", sep = "~")), simData.learn[,-takeout],simData.learn[,id])
  r <- randomForest(simData.learn[,-takeout],simData.learn[,id])
  genePredict=predict(r,simData.test[,-takeout])
  mse <- mean((simData.test[,id]-genePredict)^2)
  x <- data.frame(simData.test[,id])
  x$predicted = genePredict
  
  if (dim(x)[1]>4){
    c <- rcorr(as.matrix(x), type = "spearman")
    outputList <- data.frame(mse = mse, Spear_corr = c$r[1,2])}
  else{
    outputList <- data.frame(mse = mse)
    }

  if(TRUE == mixedSpeciesData) {
    if(length(grep("HUMAN", id)) > 0) {
      zeroPredictions = sum(x$predicted[grep("MOUSE", rownames(simData.test))] == 0)
      outputList$mse = mean(x$predicted[grep("MOUSE", rownames(simData.test))])
    } else {
      zeroPredictions = sum(x$predicted[grep("HUMAN", rownames(simData.test))] == 0)
      outputList$mse = mean(x$predicted[grep("HUMAN", rownames(simData.test))])
    }
    outputList$zero.predicted = zeroPredictions
  }
  
  return(outputList)
}
#' Random Forest imputation - multiple columns
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
randomForestImputeAll <- function(selectedList, simData.learn, simData.test) {
  simulationResult <- mapply(randomForestImpute, selectedList,
                             MoreArgs = list(
                               simData.learn,
                               simData.test
                             ))
  return(simulationResult)
}