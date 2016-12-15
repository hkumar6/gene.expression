LASSO <- function(ID){
  simData_learn <- simData.learn[,-c(which(colnames(simData.learn) == ID))]
  simData_test <- simData.test[,-c(which(colnames(simData.test) == ID))]
  vec.learn <- simData.learn[,which(colnames(simData.learn) == ID)]
  vec.test <- simData.test[,which(colnames(simData.test) == ID)]
  fit.expr <- glmnet(x = simData_learn, y = vec.learn, family = "gaussian", standardize = TRUE)
  CV.expr <- cv.glmnet(simData_learn,vec.learn, type.measure="mse")
  predicted <- predict(fit.expr, newx = simData_test, s = c(CV.expr$lambda.min))
  mse <- mean((predicted-vec.test)^2)
  c <- rcorr(predicted,vec.test, type = "spearman")
  Spear_corr <- c$r[1,2]
  cov <- sum(coef(fit.expr, CV.expr$lambda.min)>0)
  result <- data.frame(mse = mse, Spear_corr = Spear_corr, number.cov = cov)
  return(result)
}


LASSOImputeAllGenes <- function(){
  genes <- mapply(LASSO,selectedGenes)
  return(genes)
}