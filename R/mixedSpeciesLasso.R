library(glmnet)
library(Hmisc)
mixedSpeciesdata <- mixedSpecies100

human.simData <- mixedSpeciesdata[grep("HUMAN",rownames(mixedSpeciesdata)),]

mouse.simData <- mixedSpeciesdata[grep("MOUSE",rownames(mixedSpeciesdata)),]

local.means.human <- data.frame(rowMeans(human.simData))
names(local.means.human) <- c("rowMeans")
rownames(local.means.human)[order(local.means.human$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.human
local.means.mouse <- data.frame(rowMeans(mouse.simData))
names(local.means.mouse) <- c("rowMeans")
rownames(local.means.mouse)[order(local.means.mouse$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.mouse

simData.human <- mixedSpeciesdata[selectedGenes.human,]
simData.mouse <- mixedSpeciesdata[selectedGenes.mouse,]

simData <- rbind(simData.human,simData.mouse)

human.cell <- colnames(simData)[which(colSums(simData.human) > colSums(simData.mouse))]
human.cell.simData <- simData[,which(colSums(simData.human) > colSums(simData.mouse))]
human.cell <- unlist(lapply(human.cell, function(x){ toString(c("HUMAN",x))}))
colnames(human.cell.simData) <- human.cell

mouse.cell <- colnames(simData)[which(colSums(simData.human) < colSums(simData.mouse))]
mouse.cell.simData <- simData[,which(colSums(simData.human) < colSums(simData.mouse))]
mouse.cell <- unlist(lapply(mouse.cell, function(x){ toString(c("MOUSE",x))}))
colnames(mouse.cell.simData) <- mouse.cell

simData <- cbind(human.cell.simData, mouse.cell.simData)

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
    prediction <- exp(prediction)-1
    coef <- sum(coef(fit.expr.lin, cv.expr.lin$lambda.min)!=0)
    
    if (length(grep("HUMAN",colnames(simData)[i]))>0){
      test <- colnames(simData)[which(coef(fit.expr.lin, cv.expr.lin$lambda.min)!=0)]
      H.coef <- grep("HUMAN",test)
      coef.perc <- length(H.coef)/coef
      M.test <- grep("MOUSE",rownames(prediction))
      prediction.new <- prediction[M.test]
      counts <- sum(prediction.new == 0)
      mse.lin <- mean(prediction.new)
    }else{
      test <- colnames(simData)[which(coef(fit.expr.lin, cv.expr.lin$lambda.min)!=0)]
      H.coef <- grep("HUMAN",test)
      coef.perc <- length(H.coef)/coef
      H.test <- grep("HUMAN",rownames(prediction))
      prediction.new <- prediction[H.test,drop=FALSE]
      counts <- sum(prediction.new == 0)
      mse.lin <- mean(prediction.new)
    }
    
  } else { # if there are less than 9 non zero entries compare results with mse
    prediction.lin <- predict(fit.expr.lin, newx = simData_test)
    test <- apply(prediction.lin,2,function(x) mean(((exp(x)-1))^2))
    prediction <- prediction.lin[,which(test == min(test))[1]]
    if (is.na(fit.expr.lin$lambda)[1]){
      coef= NaN
    } else {
      j  <- which(test == min(test))[1]
    coef <- sum(coef(fit.expr.lin, s = fit.expr.lin$lambda[j])!=0) 
    
    if (length(grep("HUMAN",colnames(simData)[i]))>0){
      test <- colnames(simData)[which(coef(fit.expr.lin, s = fit.expr.lin$lambda[j])!=0)]
      H.coef <- grep("HUMAN",test)
      coef.perc <- length(H.coef)/coef
      M.test <- grep("MOUSE",rownames(simData_test))
      prediction.new <- prediction[M.test]
      counts <- sum(prediction.new == 0)
      mse.lin <- mean(prediction.new)
    }else{
      test <- colnames(simData)[which(coef(fit.expr.lin, s = fit.expr.lin$lambda[j])!=0)]
      H.coef <- grep("HUMAN",test)
      coef.perc <- length(H.coef)/coef
      H.test <- grep("HUMAN",rownames(simData_test))
      prediction.new <- prediction[H.test,drop=FALSE]
      counts <- sum(prediction.new == 0)
      mse.lin <- mean(prediction.new)
    }}}
  
  }else{
    mse.lin = NaN
    coef = NaN
    coef.perc = NaN
    counts = NaN
  }
  
  result <- data.frame(mse = mse.lin, number.coef = coef, nonzero = n, perc.human = coef.perc,
                       zero.predicted = counts)
  return(result)
}

#resultgenes <- mapply(lasso.mixed.data, rownames(simData)[1:15], MoreArgs = list(t(simData)))
