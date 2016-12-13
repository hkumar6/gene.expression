LASSO <- function(geneID){
  simData.gene.learn <- simData.learn[,-c(which(colnames(simData.learn) == geneID))]
  simData.gene.test <- simData.test[,-c(which(colnames(simData.test) == geneID))]
  gene.learn <- simData.learn[,which(colnames(simData.learn) == geneID)]
  gene.test <- simData.test[,which(colnames(simData.test) == geneID)]
  fit.expr <- glmnet(x = simData.gene.learn, y = gene.learn, family = "gaussian", standardize = TRUE)
  CV.expr <- cv.glmnet(simData.gene.learn,t(gene.learn), type.measure="mse")
  gene.predicted <- predict(fit.expr, newx = simData.gene.test, s = c(CV.expr$lambda.min))
  mse <- sum((gene.predicted-gene.test)^2)/length(selectedCells)
  c <- rcorr(gene.predicted,gene.test, type = "spearman")
  Spear_corr <- c$r[1,2]
  cov <- sum(coef(fit.expr, CV.expr$lambda.min)>0)
  result <- list(mse = mse, Spear_corr = Spear_corr, real = gene.test, predicted = gene.predicted,
                 number.cov = cov)
  return(result)
}


LASSOImputeAllGenes <- function(){
  genes <- mapply(LASSO,selectedGenes)
  return(genes)
}