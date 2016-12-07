# function to compute results for a specific gene
kknnImputeGene <- function(geneID) {
  trainingInfo <- train.kknn(as.formula(paste(geneID, ".", sep = "~")), simData.learn, kmax=40, kernel = c("triangular", "rectangular", "epanechnikov", "optimal"), distance = 2)
  return(c(geneID, trainingInfo$best.parameters$k, trainingInfo$best.parameter$kernel))
}

# function to compute results for all genes
kknnImputeAllGenes <- function() {
  t <- apply(as.array(selectedGenes), 1, kknnImputeGene)
  return(t)
}