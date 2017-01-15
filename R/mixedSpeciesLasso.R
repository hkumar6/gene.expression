#library(glmnet)
#library(Hmisc)
#mixedSpeciesdata <- mixedSpecies100

#human.simData <- mixedSpeciesdata[grep("HUMAN",rownames(mixedSpeciesdata)),]
#mouse.simData <- mixedSpeciesdata[grep("MOUSE",rownames(mixedSpeciesdata)),]

#local.means.human <- data.frame(rowMeans(human.simData))
#names(local.means.human) <- c("rowMeans")
#rownames(local.means.human)[order(local.means.human$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.human
#local.means.mouse <- data.frame(rowMeans(mouse.simData))
#names(local.means.mouse) <- c("rowMeans")
#rownames(local.means.mouse)[order(local.means.mouse$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.mouse

#simData.human <- mixedSpeciesdata[selectedGenes.human,]
#simData.mouse <- mixedSpeciesdata[selectedGenes.mouse,]

#simData <- rbind(simData.human,simData.mouse)


#use simData for analysis of cells
#use t(simData) for analysis of genes

lasso.mixed.data <- function(ID,simData){
  i = which(colnames(simData)==ID)
  
  #formate for lasso function
  simData <- as.matrix(simData) 
  
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

#resultgenes <- mapply(lasso.mixed.data, rownames(simData), MoreArgs = list(t(simData)))