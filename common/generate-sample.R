# select 500 most abundantly expressed genes
means <- data.frame(rowMeans(expr))
names(means) <- c("rowMeans")
rownames(means)[order(means$rowMeans, decreasing = TRUE)[1:500]] -> selectedGenes
sample(colnames(expr), 0.05*ncol(expr)) -> selectedCells

# prepare data for simulation
expr[selectedGenes, ] -> simData
simData[, -which(colnames(simData) %in% selectedCells)] -> simData.learn
simData[, selectedCells] -> simData.test

# transpose the data for kknn library
data.frame(t(simData.test)) -> simData.test
data.frame(t(simData.learn)) -> simData.learn