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
#' @name DropoutSimulation
#' @rdname DropoutSimulation
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


#' Simulations on genes with dependency on cells
#' 
#' @param theObject an object of the DropoutSimulation class
#' @param expressionData the data matrix containing gene expression values,
#'  genes as rows and cells as columns
#' @param dropoutPercentage a vector specifying dropout percentages to simulate
#' @param n the number of simulations for each dropout
#' 
#' @rdname DropoutSimulation
#' @docType methods
#' @exportMethod simulateDropoutGene
setGeneric(name = "simulateDropoutGene",
           def = function(theObject, expressionData, dropoutPercentage, n) {
             standardGeneric("simulateDropoutGene")
           })

#' @rdname DropoutSimulation
#' @docType methods
#' @export
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

#' Simulations on cells with dependency on genes
#' 
#' @param theObject an object of the DropoutSimulation class
#' @param expressionData the data matrix containing gene expression values,
#'  genes as rows and cells as columns
#' @param dropoutPercentage a vector specifying dropout percentages to simulate
#' @param n the number of simulations for each dropout
#' 
#' @rdname DropoutSimulation
#' @docType methods
#' @exportMethod simulateDropoutCells
setGeneric(name = "simulateDropoutCells",
           def = function(theObject, expressionData, dropoutPercentage, n) {
             standardGeneric("simulateDropoutCells")
           })

#' @rdname DropoutSimulation
#' @docType methods
#' @export
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

#' Analysis of simulation results
#' 
#' This function plots graphs to help analyse between lasso and knn approaches.
#' 
#' @param theObject an object of the DropoutSimulation class containing simulation results
#' @param type the type of graph to plot, "compare" will plot difference of mse(also the default value)
#' @param p the dropout percentage to plot, default: 0 will consider the whole range
#' 
#' @rdname DropoutSimulation
#' @docType methods
#' @importFrom ggplot2 ggplot aes labs geom_density
#' @exportMethod plot
setGeneric(name = "plot",
           def = function(theObject, type = "compare", p = 0) {
             standardGeneric("plot")
           })

#' @rdname DropoutSimulation
#' @docType methods
#' @export
setMethod(f = "plot",
          signature = "DropoutSimulation",
          definition = function(theObject, type = "compare", p = 0){
            
            # plot type 1
            # comparison over all dropout percentages
            if(type == "compare"){
              plotData <- vector()
              for(i in seq(1, length(theObject@simulation.result.cells)-1, by=2)) {
                if(attr(theObject@simulation.result.cells[[i]], "drop-percentage") == p || 0 == p) {
                  plotData <- c(plotData,
                                mapply('-',
                                       theObject@simulation.result.cells[[i]]["mse",],
                                       theObject@simulation.result.cells[[i+1]]["mse",]))
                }
              }
              pd <- as.data.frame(plotData)
              plotObject <- ggplot(pd, aes(plotData)) + geom_density()
              if(0 == p){
                plotObject <- plotObject + labs(title = "MSE differences (lasso-knn) over all simulations", x = "mse(lasso)-mse(knn)")
              } else {
                plotObject <- plotObject + labs(title = paste("MSE differences (lasso-knn) over", p, "dropout"), x = "mse(lasso)-mse(knn)")
              }
              return(plotObject)
            } else {
                stop("Invalid parameters specified")
            }
          })