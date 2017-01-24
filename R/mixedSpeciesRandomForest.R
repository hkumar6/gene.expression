

randomForest.mixed.data <- function(id,simData){

  i = which(colnames(simData)==id)
  simData.test <- simData[which(simData[,i]==0),,drop=FALSE] # define all zero entries of ID as test
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