genes.perc <- c()
genes.perc2 <- c()
for (i in 5:25){
  expr[selectedGenesforSim, ] -> simData
  sample(rownames(simData), i/100*nrow(simData)) -> selectedGenes
  Cells <- rownames(simData)
  
  # prepare data for simulation
  simData[-which(rownames(simData) %in% selectedGenes),] -> simData.learn
  simData[selectedGenes,] -> simData.test
  
  genes1 <- mapply(LASSO,Cells[1:50])
  genes.perc <- rbind(genes.perc,genes1)
  
  # transform to data.frame for kknn library
  data.frame(simData.test) -> kknn.simData.test
  data.frame(simData.learn) -> kknn.simData.learn
  
  genes2 <- mapply(kknnImputeGene,Cells[1:50])
  genes.perc2 <- rbind(genes.perc2,genes2)
  
}