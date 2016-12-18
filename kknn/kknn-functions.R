# function to compute results for a specific gene/cell
kknnImpute <- function(id, simData.learn, simData.test) {
  trainingInfo <- train.kknn(as.formula(paste(id, ".", sep = "~")), simData.learn, kmax=40, kernel = c("triangular", "rectangular", "epanechnikov", "optimal"), distance = 2)
  r <- kknn(as.formula(paste(id, ".", sep = "~")), simData.learn, simData.test, distance = 2, kernel = trainingInfo$best.parameters$kernel, k = trainingInfo$best.parameters$k)
  mse <- sum((simData.test[id]-r$fitted.values)^2)/(dim(simData.test[id])[1])
  
  x <- simData.test[id]
  x$predicted = r$fitted.values
  c <- rcorr(as.matrix(x), type = "spearman")
  outputList <- data.frame(mse = mse, Spear_corr = c$r[1,2],
                     optimalK = trainingInfo$best.parameters$k,
                     optimalKernel = trainingInfo$best.parameter$kernel)
  return(outputList)
}

# TODO
# function to compute results for all genes/cells
kknnImputeAll <- function(selectedList, simData.learn, simData.test) {
  simulationResult <- mapply(kknnImputeGene, selectedList)
  return(simulationResult)
}