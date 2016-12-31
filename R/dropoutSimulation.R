#' The DropoutSimulation class
#' 
#' This class has slots to store results of different simulations
#' using Lasso and kknn methods.
#' 
#' @slot dropout.percentage Vector containing values between 0 and 1
#'  representing the percentages of dropout to simulate over repitions
#' @slot n an integer specifying the number of repititions for
#'  a particular dropout percentage
#' @slot simulation.result.genes a list storing results of all simulations
#'  where genes are imputed as function of cells
#' @slot simulation.result.cells a list storing results of all simulations
#'  where cells are imputed as functions of genes
#'  
#' @importFrom methods new
#' @exportClass DropoutSimulation

DropoutSimulation <- setClass("DropoutSimulation",
                              slots = 
                                c(
                                  dropout.percentage = "vector",
                                  n = "numeric",
                                  simulation.result.genes = "list",
                                  simulation.result.cells = "list"
                                ),
                              prototype = list(
                                simulation.result.genes = list(),
                                simulation.result.cells = list()
                              ))


# method to run simulations on genes with dependency on cells
setGeneric(name = "simulateDropoutGene",
           def = function(theObject, expressionData, dropoutPercentage, n) {
             standardGeneric("simulateDropoutGene")
           })

setMethod(f = "simulateDropoutGene",
          signature = "DropoutSimulation",
          definition = function(theObject, expressionData, dropoutPercentage, n) {
            theObject@dropout.percentage <- dropoutPercentage
            theObject@n <- n
            
            # select 500 most abundantly expressed genes
            local.means <- data.frame(rowMeans(expressionData))
            names(local.means) <- c("rowMeans")
            rownames(local.means)[order(local.means$rowMeans, decreasing = TRUE)[1:500]] -> local.selectedGenes
            expressionData[local.selectedGenes, ] -> local.simData
            for (p in theObject@dropout.percentage) {
              for (simID in 1:n) {
                sample(colnames(expressionData), p*ncol(expressionData)) -> local.selectedCells
                
                # simulation for lasso
                td <- mapply(lassoImpute, local.selectedGenes,
                             MoreArgs = list(as.matrix(t(local.simData[, -which(colnames(local.simData) %in% local.selectedCells)])), 
                                             as.matrix(t(local.simData[, local.selectedCells]))))
                theObject@simulation.result.genes[[length(theObject@simulation.result.genes)+1]] <- td
                attr(theObject@simulation.result.genes[[length(theObject@simulation.result.genes)]], "drop-percentage") <- p
                attr(theObject@simulation.result.genes[[length(theObject@simulation.result.genes)]], "method") <- "lasso"
                
                # simulation for kknn
                td <- mapply(kknnImpute, local.selectedGenes,
                             MoreArgs = list(as.data.frame(t(local.simData[, -which(colnames(local.simData) %in% local.selectedCells)])),
                                             as.data.frame(t(local.simData[, local.selectedCells]))))
                theObject@simulation.result.genes[[length(theObject@simulation.result.genes)+1]] <- td
                attr(theObject@simulation.result.genes[[length(theObject@simulation.result.genes)]], "drop-percentage") <- p
                attr(theObject@simulation.result.genes[[length(theObject@simulation.result.genes)]], "method") <- "kknn"
              }
            }
            return(theObject)
          })

# method to run simulations on cells with dependency on genes
setGeneric(name = "simulateDropoutCells",
           def = function(theObject, expressionData, dropoutPercentage, n) {
             standardGeneric("simulateDropoutCells")
           })

setMethod(f = "simulateDropoutCells",
          signature = "DropoutSimulation",
          definition = function(theObject, expressionData, dropoutPercentage, n) {
            theObject@dropout.percentage <- dropoutPercentage
            theObject@n <- n
            
            # select 500 most abundantly expressed genes
            local.means <- data.frame(rowMeans(expressionData))
            names(local.means) <- c("rowMeans")
            rownames(local.means)[order(local.means$rowMeans, decreasing = TRUE)[1:500]] -> local.selectedGenes
            expressionData[local.selectedGenes, ] -> local.simData
            colnames(local.simData) <- gsub("[_?]", "", gsub("^[0-9]", "X", colnames(local.simData), perl = TRUE), perl = TRUE)
            for (p in theObject@dropout.percentage) {
              for (simID in 1:n) {
                sample(rownames(local.simData), p*nrow(local.simData)) -> local.dropoutGenes
                # simulation for lasso
                td <- mapply(lassoImpute, colnames(local.simData),
                             MoreArgs = list(as.matrix(local.simData[-which(rownames(local.simData) %in% local.dropoutGenes),]),
                                             as.matrix(local.simData[local.dropoutGenes,]) ))
                theObject@simulation.result.cells[[length(theObject@simulation.result.cells)+1]] <- td
                attr(theObject@simulation.result.cells[[length(theObject@simulation.result.cells)]], "drop-percentage") <- p
                attr(theObject@simulation.result.cells[[length(theObject@simulation.result.cells)]], "method") <- "lasso"
                
                # simulation for kknn
                td <- mapply(kknnImpute, colnames(local.simData),
                             MoreArgs = list(as.data.frame(local.simData[-which(rownames(local.simData) %in% local.dropoutGenes),]),
                                             as.data.frame(local.simData[local.dropoutGenes,]) ))
                theObject@simulation.result.cells[[length(theObject@simulation.result.cells)+1]] <- td
                attr(theObject@simulation.result.cells[[length(theObject@simulation.result.cells)]], "drop-percentage") <- p
                attr(theObject@simulation.result.cells[[length(theObject@simulation.result.cells)]], "method") <- "kknn"
              }
            }
            return(theObject)
          })