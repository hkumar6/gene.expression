randomForest.mixed.data <- function(id,simData){

    
  i = which(colnames(simData)==id)
  simData.test <- simData[which(simData[,i]==0),,drop=FALSE] # define all zero entries of ID as test
  simData.learn <- simData[-c(which(simData[,i]==0)),,drop=FALSE] # define all non-zero entries as learn    
  # define all zero entries of ID as test
 # simData.test <- simData[which(simData[,id]==0),,drop=FALSE]
  # define all non-zero entries as learn
#  simData.learn <- simData[-c(which(simData[,id]==0)),,drop=FALSE]
  
  
  
  # define all zero entries of ID as test
 # simData.test <- simData[which(simData[,id]==0),,drop=FALSE]
  # define all non-zero entries as learn
#  simData.learn <- simData[-c(which(simData[,id]==0)),,drop=FALSE]
  
#        if (length(grep("HUMAN",colnames(simData)[i]))>0){
 #         
  #        M.test <- grep("MOUSE",rownames(prediction))
   #       prediction.new <- prediction[M.test]
    #      counts <- sum(prediction.new == 0)
     #     mse.lin <- mean(prediction.new)
      #  }else{
       #   H.test <- grep("HUMAN",rownames(prediction))
        #  prediction.new <- prediction[H.test,drop=FALSE]
         # counts <- sum(prediction.new == 0)
          #mse.lin <- mean(prediction.new)
        #}
  
  
  
  if(nrow(simData.learn) > 4) {
    output <- randomForestImpute(id, simData.learn, simData.test)
    # output <- estimateKmax(id, simData.learn = simData_learn, simData.test = simData_test, kmax)
    return(output)
  } else {
    noImputationOutput <- data.frame(mse = NA, Spear_corr = NA,
                                     optimalK = NA,
                                     optimalKernel = NA)
    return(noImputationOutput)
  }
}