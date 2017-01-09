kknn.mixed.data <- function(id,simData){
  
  # define all zero entries of ID as test
  simData.test <- simData[which(simData[,id]==0),,drop=FALSE]
  # define all non-zero entries as learn
  simData.learn <- simData[-c(which(simData[,id]==0)),,drop=FALSE]
  print(nrow(simData.learn))
  if(nrow(simData.learn) > 3) {
    kmax <- sum(apply(simData.learn, 2, sum) != 0)
    kmax <- if(kmax > 20) 20 else kmax
    # output <- kknnImpute(id, simData.learn = simData.learn, simData.test = simData.test, 10)
    output <- estimateKmax(id, simData.learn = simData.learn, simData.test = simData.test, kmax)
    print(output)
  } else {
    noImputationOutput <- data.frame(mse = NA, Spear_corr = NA,
                             optimalK = NA,
                             optimalKernel = NA)
    print(noImputationOutput)
  }
}


estimateKmax <- function(id, simData.learn, simData.test, kmax) {
  result <- tryCatch(
    {
      kknnImpute(id, simData.learn = simData.learn, simData.test = simData.test, kmax)
    },
    error = function(err) {
      warning(paste("Kmax", kmax, "generated error."))
      if(kmax <= 2) {
        noImputationOutput <- data.frame(mse = NA, Spear_corr = NA,
                                         optimalK = NA,
                                         optimalKernel = NA)
        return(noImputationOutput)
      }
      return(estimateKmax(id, simData.learn = simData.learn, simData.test = simData.test, kmax-3))
    }
  )
  return(result)
}