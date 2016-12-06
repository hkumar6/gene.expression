# configure environment
setwd("~/Code/comp-bio")
load("~/Code/comp-bio/treutlein2016.rdata")
library("kknn")
library("ggplot2")

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

results <- train.kknn(Sod1~., simData.learn, kmax=25, kernel = c("triangular", "rectangular", "epanechnikov", "optimal"), distance = 2)
plot(results)

# simulation
simData.test[,"Sod1"] = 0
simData.kknn <- kknn(Sod1~., simData.learn, simData.test, distance = 2, kernel = results$best.parameters$kernel, k = results$best.parameters$k)
plotData <- data.frame(as.data.frame(expr["Sod1", selectedCells], row.names = selectedCells), as.data.frame(simData.kknn$fitted.values, row.names = selectedCells))
names(plotData) <- c("originalExpression", "fittedExpression")
ggplot(plotData, aes(y=fittedExpression, x=originalExpression)) + geom_point(aes(color=fittedExpression)) + geom_line(aes(y=originalExpression))
sum((plotData$originalExpression - plotData$fittedExpression)^2)/dim(plotData)[1]
