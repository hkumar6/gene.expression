# function to compute results for a specific gene
kknnImputeGene <- function(geneID) {
  trainingInfo <- train.kknn(as.formula(paste(geneID, ".", sep = "~")), simData.learn, kmax=40, kernel = c("triangular", "rectangular", "epanechnikov", "optimal"), distance = 2)
  r <- kknn(as.formula(paste(geneID, ".", sep = "~")), simData.learn, simData.test, distance = 2, kernel = trainingInfo$best.parameters$kernel, k = trainingInfo$best.parameters$k)
  outputList <- list(optimalK = trainingInfo$best.parameters$k, optimalKernel = trainingInfo$best.parameter$kernel, predicted = r$fitted.values)
  return(outputList)
}

# function to compute results for all genes
kknnImputeAllGenes <- function() {
  genes <- mapply(kknnImputeGene, selectedGenes)
  return(genes)
}