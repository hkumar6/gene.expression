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
                
                # simulation for lasso
                td <- mapply(lassoImpute, local.selectedGenes[1:3],
                             MoreArgs = list(as.matrix(t(local.simData[, -which(colnames(local.simData) %in% local.selectedCells)])), 
                                             as.matrix(t(local.simData[, local.selectedCells]))))
                theObject@simulation.result[[length(theObject@simulation.result)+1]] <- td
                attr(theObject@simulation.result[[length(theObject@simulation.result)]], "drop-percentage") <- p
                attr(theObject@simulation.result[[length(theObject@simulation.result)]], "method") <- "lasso"
                
                # simulation for kknn
                td <- mapply(kknnImpute, local.selectedGenes[1:3],
                             MoreArgs = list(as.data.frame(t(local.simData[, -which(colnames(local.simData) %in% local.selectedCells)])),
                                             as.data.frame(t(local.simData[, local.selectedCells]))))
                theObject@simulation.result[[length(theObject@simulation.result)+1]] <- td
                attr(theObject@simulation.result[[length(theObject@simulation.result)]], "drop-percentage") <- p
                attr(theObject@simulation.result[[length(theObject@simulation.result)]], "method") <- "kknn"
              }
            }
            return(theObject)
          })