
lasso.negbin.mixed.data <- function(ID,simData){
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
    fit.expr.bin <- glmreg(vec.learn ~ ., data = simData_learn, family = "negbin",  theta=1)
    cv.bin <- cv.glmreg(vec.learn ~ ., data = simData_learn, family = "negbin",  theta=1)
    prediction <- predict(fit.expr.bin, newx = simData_test, which =cv.bin$lambda.which, type = "response")
    
    lambda <- cv.bin$lambda.optim
    lambda.w <- cv.bin$lambda.which
    
    coef = sum(coef(fit.expr.bin, which = 30) !=0)
    if (length(grep("HUMAN",colnames(simData)[i]))>0){
      test <- colnames(simData)[which(coef(fit.expr.bin, which = 30)!=0)]
      H.coef <- grep("HUMAN",test)
      coef.perc <- length(H.coef)/coef
      M.test <- grep("MOUSE",rownames(simData_test))
      prediction.new <- prediction[M.test]
      counts <- sum(prediction.new == 0)
      mse.bin <- mean(prediction.new)
    }else{
      test <- colnames(simData)[which(coef(fit.expr.bin, which = 30)!=0)]
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
  result <- data.frame(mse = mse.bin, number.coef = coef, nonzero = n, perc.human= coef.perc,
                       zero.predicted = counts, lambda.which = lambda.w, lambda = lambda)
  return(result)
}

system.time(

#system.time(
 # test <- mapply(lasso.mixed.data, colnames(simData)[c(9)], MoreArgs = list(simData)))

#system.time(
 # test <- mapply(randomForest.mixed.data, colnames(simData)[c(9)], MoreArgs = list(simData)))

