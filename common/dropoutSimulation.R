# Class definition to store results for simulation

DropoutSimulation <- setClass("DropoutSimulation",
                              slots = 
                                c(
                                  data.good = "ANY",
                                  dropout.percentage = "vector",
                                  n = "numeric",
                                  simulation.result = "list"
                                ),
                              prototype = list(
                                simulation.result = list()
                              ))

# method to assign data for the gene expression
setGeneric(name = "test",
           def = function(theObject, expressionData, dropoutPercentage, n) {
             standardGeneric("test")
           })

setMethod(f = "test",
          signature = "DropoutSimulation",
          definition = function(theObject, expressionData, dropoutPercentage, n) {
            theObject@data.good <- expressionData
            theObject@dropout.percentage <- dropoutPercentage
            theObject@n <- n
            for(i in 1:n) {
              theObject@simulation.result[[length(theObject@simulation.result)+1]] <- c(1,2,3)
            }
            return(theObject)
          })



# method to assign data for the gene expression
setGeneric(name = "simulateDropoutGene",
           def = function(theObject, expressionData, dropoutPercentage, n) {
             standardGeneric("simulateDropoutGene")
           })

setMethod(f = "simulateDropoutGene",
          signature = "DropoutSimulation",
          definition = function(theObject, expressionData, dropoutPercentage, n) {
            theObject@data.good <- expressionData
            theObject@dropout.percentage <- dropoutPercentage
            theObject@n <- n
            
            # select 500 most abundantly expressed genes
            local.means <- data.frame(rowMeans(theObject@data.good))
            names(local.means) <- c("rowMeans")
            rownames(local.means)[order(local.means$rowMeans, decreasing = TRUE)[1:500]] -> local.selectedGenes
            expr[local.selectedGenes, ] -> local.simData
            for (p in theObject@dropout.percentage) {
              for (simID in 1:n) {
                sample(colnames(theObject@data.good), p*ncol(theObject@data.good)) -> local.selectedCells
                # t(simData[, -which(colnames(local.simData) %in% local.selectedCells)]) -> local.simData.learn
                # t(simData[, local.selectedCells]) -> local.simData.test
                # print(dim(local.simData.learn))
                # print(dim(local.simData.test))
                td <- mapply(kknnImputeGene, local.selectedGenes[1:3],
                             MoreArgs = list(as.data.frame(t(local.simData[, -which(colnames(local.simData) %in% local.selectedCells)])),
                                             as.data.frame(t(local.simData[, local.selectedCells]))))
                theObject@simulation.result[[length(theObject@simulation.result)+1]] <- td
              }
            }
            return(theObject)
          })