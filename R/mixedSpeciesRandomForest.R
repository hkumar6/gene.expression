#' randomForest imputation - single column
#' 
#' Impute results for a specific gene or cell using randomForest
#' 
#' @param id string for the cell or gene to impute, must be a column name in the data
#' @param simData count data
#' 
#' @return data.frame containing results of the simulation
#'      mean error
#'
#' @importFrom stats as.formula
#' @importFrom randomForest randomForest
#' @importFrom Hmisc rcorr
#' @export

randomForest.mixed.data <- function(id,simData){

  i = which(colnames(simData)==id)
  simData.test <- simData[which(simData[,i]==0),,drop=FALSE] # define all zero entries of id as test
  simData.learn <- simData[-c(which(simData[,i]==0)),,drop=FALSE] # define all non-zero entries as learn
  if(nrow(simData.learn) > 4) {
    output <- randomForestImpute(id, simData.learn, simData.test, mixedSpeciesData = TRUE)
    # output <- estimateKmax(id, simData.learn = simData_learn, simData.test = simData_test, kmax)
    return(output)
  } else {
    noImputationOutput <- data.frame(mse = NA, Spear_corr = NA,
                                     optimalK = NA,
                                     optimalKernel = NA)
    return(noImputationOutput)
  }
}
