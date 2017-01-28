
lasso.negbin.mixed.data <- function(ID,simData, genes = TRUE){
  i = which(colnames(simData)==ID)
  
  #formate for lasso function
  simData <- as.data.frame((simData)) 
  
  #define test and learn data
  simData.test <- simData[which(simData[,i]==0),,drop=FALSE] # define all zero entries of ID as test
  simData.learn <- simData[-c(which(simData[,i]==0)),,drop=FALSE] # define all non-zero entries as learn
  
  # take out ID from the data set for analysis
  simData_learn <- simData.learn[,-c(which(colnames(simData.learn)==ID)),drop = FALSE] 
  simData_test <- simData.test[,-c(which(colnames(simData.test)==ID)),drop = FALSE]
  
  #define the vector ID for analysis
  vec.learn <- simData.learn[,which(colnames(simData.learn)==ID)]
  vec.test <- simData.test[,which(colnames(simData.test)==ID)]
  n = dim(simData_learn)[1] # number of non-zero entries
  if (dim(simData_learn)[1]>2 && dim(simData_test)[1]>2){ #if there is more than one entry non zero do regression
    fit.expr.bin <- glmreg(vec.learn ~ ., data = simData_learn, family = "negbin",  theta=1, alpha = 1)
    #cv.bin <- cv.glmreg(vec.learn ~ ., data = simData_learn, family = "negbin",  theta=1, nfold = 5, alpha = 1)
    
    if (genes){
      lambda = 0.0065
    }else{
      lambda = 0.2
    }
    
    lambda.pos <- which(abs(fit.expr.bin$lambda-lambda) == min(abs(fit.expr.bin$lambda-lambda)))
    prediction <- predict(fit.expr.bin, newx = simData_test, which = lambda.pos, type = "response")
    
    coef = sum(coef(fit.expr.bin, which = lambda.pos) !=0)
    if (length(grep("HUMAN",colnames(simData)[i]))>0){
      test <- colnames(simData)[which(coef(fit.expr.bin, which = lambda.pos)!=0)]
      H.coef <- grep("HUMAN",test)
      coef.perc <- length(H.coef)/coef
      M.test <- grep("MOUSE",rownames(simData_test))
      prediction.new <- prediction[M.test]
      counts <- sum(prediction.new == 0)
      mse.bin <- mean(prediction.new)
    }else{
      test <- colnames(simData)[which(coef(fit.expr.bin, which = lambda.pos)!=0)]
      H.coef <- grep("HUMAN",test)
      coef.perc <- length(H.coef)/coef
      H.test <- grep("HUMAN",rownames(simData_test))
      prediction.new <- prediction[H.test,drop=FALSE]
      counts <- sum(prediction.new == 0)
      mse.bin <- mean(prediction.new)
    }
    
  }else{
    mse.bin = NaN
    coef = NaN
    coef.perc = NaN
    counts = NaN
  }
  result <- data.frame(mean.error = mse.bin, number.coef = coef, nonzero = n, perc.human= coef.perc,
                       zero.predicted = counts)
  return(result)
}

#system.time(
 # test <- mapply(lasso.negbin.mixed.data, sample(rownames(simData),50), MoreArgs = list(t(simData))))

#system.time(
#test <- mapply(lasso.negbin.mixed.data, rownames(simData)[22], MoreArgs = list(t(simData))))

#cells lambda 0.2
#genes 0.0065


