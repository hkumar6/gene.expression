randomForest.mixed.data <- function(id,simData){

    
  i = which(colnames(simData)==id)
  simData.test <- simData[which(simData[,i]==0),,drop=FALSE] # define all zero entries of id as test
  simData.learn <- simData[-c(which(simData[,i]==0)),,drop=FALSE] # define all non-zero entries as learn    
  # define all zero entries of id as test
 # simData.test <- simData[which(simData[,id]==0),,drop=FALSE]
  # define all non-zero entries as learn
#  simData.learn <- simData[-c(which(simData[,id]==0)),,drop=FALSE]
  
  
  
  # define all zero entries of id as test
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
  
  vec.learn <- simData.learn[,which(colnames(simData.learn)==id)]
  vec.test <- simData.test[,which(colnames(simData.test)==id)]
  n = dim(simData_learn)[1] # number of non-zero entries
  if (dim(simData_learn)[1]>2){ #if there is more than one entry non zero do regression
    fit.expr.lin <- glmnet(x = simData_learn, y = log(1+vec.learn),
                           family = "gaussian", standardize = TRUE, alpha = 1)
    
    if (n>=12){ # if there are more than 9 non zero entries cross-validation is possible
      cv.expr.lin <- cv.glmnet(x=simData_learn, y = log(1+vec.learn), nfold=3, type.measure="mse")
      prediction <- predict(fit.expr.lin, newx = simData_test, s = c(cv.expr.lin$lambda.min))
      prediction2 <- exp(prediction)-1
      cov <- sum(coef(fit.expr.lin, cv.expr.lin$lambda.min)!=0)
    } else { # if there are less than 9 non zero entries compare results with mse
      prediction.lin <- predict(fit.expr.lin, newx = simData_test)
      test <- apply(prediction.lin,2,function(x) mean(((exp(x)-1))^2))
      prediction <- prediction.lin[,which(test == min(test))[1]]
      if (is.na(fit.expr.lin$lambda)[1]){
        cov = NaN
      } else {
        j  <- which(test == min(test))[1]
        cov <- sum(coef(fit.expr.lin, s = fit.expr.lin$lambda[j])!=0) }}
    
    mse.lin <- mean(((exp(prediction)-1))^2)
  }else{
    mse.lin = NaN
    Spear_corr = NaN
    cov = NaN
  }
  
  result <- data.frame(mse = mse.lin, number.cov = cov, nonzero = n)
  return(result)
}
  
  
  
  
  
  
  if(nrow(simData.learn) > 3) {
    output <- randomForestImpute(id, simData.learn, simData.test)
    # output <- estimateKmax(id, simData.learn = simData_learn, simData.test = simData_test, kmax)
    return(output)
  } else {
    noImputationOutput <- data.frame(mse = NA, Spear_corr = NA)
    return(noImputationOutput)
  }
}