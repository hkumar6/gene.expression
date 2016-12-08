source("common/env-kh.R")
source("common/generate-sample.R")
source("kknn/kknn-functions.R")

allGenes <- kknnImputeAllGenes()

plotData <- data.frame(unlist(allGenes["optimalK",]), unlist(allGenes["optimalKernel",]))
names(plotData) <- c("optimalK", "optimalKernel")
ggplot(plotData, aes(x=optimalK)) + geom_bar() + ggtitle("Variation of number of neighbours for various genes")
ggplot(plotData, aes(x=optimalKernel)) + geom_bar() + coord_flip() + ggtitle("Distribution over optimal kernel for kknn")
