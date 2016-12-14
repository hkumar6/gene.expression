# select 500 most abundantly expressed genes
means <- data.frame(rowMeans(expr))
names(means) <- c("rowMeans")
rownames(means)[order(means$rowMeans, decreasing = TRUE)[1:500]] -> selectedGenesforSim
expr[selectedGenesforSim, ] -> simData
sample(rownames(simData), 0.05*nrow(simData)) -> selectedGenes

# prepare data for simulation
simData[-which(rownames(simData) %in% selectedGenes),] -> simData.learn
simData[selectedGenes,] -> simData.test

# transform to data.frame for kknn library
data.frame(simData.test) -> kknn.simData.test
data.frame(simData.learn) -> kknn.simData.learn