#' kknn imputation for mixed species data
#' 
#' Impute results for a specific gene or cell using kknn
#' 
#' @param id string for the cell or gene to impute, must be a column name in the data
#' @param simData a matrix containing all the gene expression counts
#' 
#' @return data.frame containing results of the simulation
#'      mean squared error, Spearman correlation, optimal parameters for kknn
#'
#' @export


kknn.pbmc <- function(id,simData){
  
  # define all zero entries of ID as test
  simData.test <- simData[which(simData[,id]==0),,drop=FALSE]
  # define all non-zero entries as learn
  simData.learn <- simData[-c(which(simData[,id]==0)),,drop=FALSE]
  
  outputVector <- t(rep(0, nrow(simData)))
  outputVector[-c(which(simData[,id]==0))] <- simData.learn[-c(which(simData[,id]==0)),id,drop=FALSE]
  
  if(nrow(simData.learn) > 3) {
    trainingInfo <- train.kknn(as.formula(paste(id, ".", sep = "~")), simData.learn, kernel = c("triangular", "rectangular", "epanechnikov", "optimal"), distance = 2)
    r <- kknn(as.formula(paste(id, ".", sep = "~")), simData.learn, simData.test, distance = 2, kernel = trainingInfo$best.parameters$kernel, k = trainingInfo$best.parameters$k)
    outputVector[which(simData[,id]==0)] <- r$fittedValues
  } else {
    outputVector[which(simData[,id]==0)] <- rep(0, nrow(simData.test))
  }
  return(outputVector)
}