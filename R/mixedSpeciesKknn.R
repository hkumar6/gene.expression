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


kknn.mixed.data <- function(id,simData){
  
  # define all zero entries of ID as test
  simData.test <- simData[which(simData[,id]==0),,drop=FALSE]
  # define all non-zero entries as learn
  simData.learn <- simData[-c(which(simData[,id]==0)),,drop=FALSE]
  
  if(nrow(simData.learn) > 3) {
    kmax <- sum(apply(simData.learn, 2, sum) != 0)
    kmax <- if(kmax > 20) 20 else kmax
    output <- kknnImpute(id, simData.learn = as.data.frame(simData.learn), simData.test = as.data.frame(simData.test), kmax, mixedSpeciesData = TRUE)
    # output <- estimateKmax(id, simData.learn = simData_learn, simData.test = simData_test, kmax)
    output$nonzero = nrow(simData.learn)
    return(output)
  } else {
    noImputationOutput <- data.frame(mse = NA, Spear_corr = NA,
                             optimalK = NA,
                             optimalKernel = NA,
                             zero.predicted = NA,
                             nonzero = nrow(simData.learn))
    return(noImputationOutput)
  }
}