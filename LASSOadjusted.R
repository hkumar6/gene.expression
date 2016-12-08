# select 500 most abundantly expressed genes
means <- data.frame(rowMeans(expr))
names(means) <- c("rowMeans")
rownames(means)[order(means$rowMeans, decreasing = TRUE)[1:500]] -> selectedGenes
sample(colnames(expr), 0.05*ncol(expr)) -> selectedCells

# prepare data for simulation
expr[selectedGenes, ] -> simData
simData[, -which(colnames(simData) %in% selectedCells)] -> simData.learn
simData[, selectedCells] -> simData.test

# transpose the data for glmnet library
t(simData.test) -> simData.test
t(simData.learn) -> simData.learn

library(glmnet)
library(Hmisc)

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
  Spear_corr <- c$r
  cov <- sum(coef(fit.expr, CV.expr$lambda.min)>0)
  result <- list(mse = mse, Spear_corr = Spear_corr, real = gene.test, predicted = gene.predicted,
                 number.cov = cov)
  return(result)
}

mapply(LASSO,selectedGenes[1:2])

LASSOImputeAllGenes <- function(){
  genes <- mapply(LASSO,selectedGenes)
  return(genes)
}














