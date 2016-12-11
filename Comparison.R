genes.perc <- c()
genes.perc2 <- c()
for (i in 5:25){
  sample(colnames(expr), i/100*ncol(expr)) -> selectedCells
  # prepare data for simulation
  expr[selectedGenes, ] -> simData
  simData[, -which(colnames(simData) %in% selectedCells)] -> simData.learn
  simData[, selectedCells] -> simData.test
  # transpose the data for glmnet library
  t(simData.test) -> simData.test
  t(simData.learn) -> simData.learn
  
  genes1 <- LASSOImputeAllGenes
  genes.perc <- rbind(genes.perc,genes1)
  
  # transpose as data.frame for kknn library
  data.frame(simData.test) -> kknn.simData.test
  data.frame(simData.learn) -> kknn.simData.learn
  
  genes2 <- kknnImputeAllGenes
  genes.perc2 <- rbind(genes.perc2,genes2)
  
}









