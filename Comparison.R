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

# transpose as data.frame for kknn library
data.frame(simData.test) -> kknn.simData.test
data.frame(simData.learn) -> kknn.simData.learn

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


# function to compute results for a specific gene
kknnImputeGene <- function(geneID) {
  trainingInfo <- train.kknn(as.formula(paste(geneID, ".", sep = "~")), kknn.simData.learn, kmax=40, kernel = c("triangular", "rectangular", "epanechnikov", "optimal"), distance = 2)
  r <- kknn(as.formula(paste(geneID, ".", sep = "~")), kknn.simData.learn, kknn.simData.test, distance = 2, kernel = trainingInfo$best.parameters$kernel, k = trainingInfo$best.parameters$k)
  mse <- sum((kknn.simData.test[geneID]-r$fitted.values)^2)/dim(simData.test[geneID])[1]
  
  x <- kknn.simData.test[geneID]
  x$predicted = r$fitted.values
  c <- rcorr(as.matrix(x), type = "spearman")
  outputList <- list(mse = mse, Spear_corr = c$r[1,2],
                     optimalK = trainingInfo$best.parameters$k,
                     optimalKernel = trainingInfo$best.parameter$kernel, real = simData.test[geneID],
                     predicted = r$fitted.values)
  return(outputList)
}

# function to compute results for all genes
kknnImputeAllGenes <- function() {
  genes <- mapply(kknnImputeGene, selectedGenes)
  return(genes)
}

genes.perc <- c()
genes.perc2 <- c()
for (i in 5:7){
  sample(colnames(expr), i/100*ncol(expr)) -> selectedCells
  # prepare data for simulation
  expr[selectedGenes, ] -> simData
  simData[, -which(colnames(simData) %in% selectedCells)] -> simData.learn
  simData[, selectedCells] -> simData.test
  # transpose the data for glmnet library
  t(simData.test) -> simData.test
  t(simData.learn) -> simData.learn
  genes1 <- mapply(LASSO,selectedGenes[1:2])
  genes.perc <- rbind(genes.perc,genes1)
  
  # transpose as data.frame for kknn library
  data.frame(simData.test) -> kknn.simData.test
  data.frame(simData.learn) -> kknn.simData.learn
  
  genes2 <- mapply(kknnImputeGene,selectedGenes[1:2])
  genes.perc2 <- rbind(genes.perc2,genes2)
  
}









