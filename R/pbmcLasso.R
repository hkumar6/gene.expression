#' kknn imputation for mixed species data
#' 
#' Impute results for a specific gene or cell using kknn
#' 
#' @param id string for the cell or gene to impute, must be a column name in the data
#' @param simData a matrix containing all the gene expression counts
#' 
#' @return data.frame containing results of the simulation
#'
#' @export


lasso.pbmc <- function(id,simData){
  
  simData <- as.matrix(simData)
  # define all zero entries of ID as test
  simData.test <- simData[which(simData[,id]==0),,drop=FALSE]
  # define all non-zero entries as learn
  simData.learn <- simData[-c(which(simData[,id]==0)),,drop=FALSE]
  
  outputVector <- simData[,id]
  outputVector[-c(which(simData[,id]==0))] <- simData.learn[,id]
  
  simData_learn <- log(simData.learn[,-which(colnames(simData.learn)==id), drop = FALSE]+1)
  simData_test <- log(simData.test[,-which(colnames(simData.learn)==id), drop = FALSE]+1)
  
  vec.learn <- log(1+simData.learn[,id, drop = FALSE])
  ve.test <- log(1+simData.test[,id, drop = FALSE])
  
  if(nrow(simData.learn) > 3) {
    #fit.expr <- glmnet(x = simData_learn, y = vec.learn, family = "gaussian", standardize = TRUE)
    CV.expr <- cv.glmnet(x = simData_learn,y = vec.learn, type.measure="mse", family = "gaussian",standardize = TRUE)
    prediction <- predict(CV.expr$glmnet.fit, newx = simData_test, s = CV.expr$lambda.min )
    prediction <- exp(prediction)-1
    outputVector[which(simData[,id]==0)] <- prediction
  } else {
    outputVector[which(simData[,id]==0)] <- rep(0, nrow(simData.test))
  }
  return(outputVector)
}
