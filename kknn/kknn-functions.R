# function to compute results for a specific gene
kknnImputeGene <- function(geneID, simData.learn, simData.test) {
  trainingInfo <- train.kknn(as.formula(paste(geneID, ".", sep = "~")), simData.learn, kmax=40, kernel = c("triangular", "rectangular", "epanechnikov", "optimal"), distance = 2)
  r <- kknn(as.formula(paste(geneID, ".", sep = "~")), simData.learn, simData.test, distance = 2, kernel = trainingInfo$best.parameters$kernel, k = trainingInfo$best.parameters$k)
  mse <- sum((simData.test[geneID]-r$fitted.values)^2)/(dim(simData.test[geneID])[1])
  
  x <- simData.test[geneID]
  x$predicted = r$fitted.values
  c <- rcorr(as.matrix(x), type = "spearman")
  outputList <- list(mse = mse, Spear_corr = c$r[1,2],
                     optimalK = trainingInfo$best.parameters$k,
                     optimalKernel = trainingInfo$best.parameter$kernel)
  return(outputList)
}

# function to compute results for all genes
kknnImputeAllGenes <- function(simData.learn, simData.test) {
  genes <- mapply(kknnImputeGene, selectedGenes)
  return(genes)
}