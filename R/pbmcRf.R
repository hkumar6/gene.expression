#' random forest imputation for mixed species data
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


randomForest.pbmc <- function(id,simData){
  
  # define all zero entries of ID as test
  simData.test <- simData[which(simData[,id]==0),,drop=FALSE]
  # define all non-zero entries as learn
  simData.learn <- simData[-c(which(simData[,id]==0)),,drop=FALSE]
  
  outputVector <- simData[,id]
  
  if(nrow(simData.learn) > 3) {
    tryCatch({
      takeout<-which(colnames(simData.learn)==id)
      r <- randomForest(simData.learn[,-takeout],simData.learn[,id])
      outputVector[which(simData[,id]==0)] <- predict(r,simData.test[,-takeout])
    }, error = function(err){
      return(outputVector)
    })
    
  } else {
    outputVector[which(simData[,id]==0)] <- rep(0, nrow(simData.test))
  }
  return(outputVector)
}
