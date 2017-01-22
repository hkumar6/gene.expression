#' The mixedSpeciesImput class
#' 
#' This class has slots to store results of different simulations
#' using Lasso and kknn methods.
#' 
#' @slot simulation.result.genes a list storing results of all simulations
#'  where genes are imputed as function of cells
#' @slot simulation.result.cells a list storing results of all simulations
#'  where cells are imputed as functions of genes
#'  
#' @importFrom methods new
#' @exportClass mixedSpeciesImput
#' 
#' @examples
#' t <- mixedSpeciesImput()

mixedSpeciesImput <- setClass("mixedSpeciesImput",
                           slots = 
                             c(
                               simulation.result.genes = "list",
                               simulation.result.cells = "list"
                             ),
                           prototype = list(
                             simulation.result.genes = list(),
                             simulation.result.cells = list()
                           ))


# method to run simulations on genes with dependency on cells
setGeneric(name = "mixedSpeciesGene",
           def = function(theObject, mixedSpeciesdata) {
             standardGeneric("mixedSpeciesGene")
           })

setMethod(f = "mixedSpeciesGene",
          signature = "mixedSpeciesImput",
          definition = function(theObject, mixedSpeciesdata) {
            
            # reduce gene dimension to 100 mouse and 100 human genes with the highest gene expression
            human.simData <- mixedSpeciesdata[grep("HUMAN",rownames(mixedSpeciesdata)),]
            mouse.simData <- mixedSpeciesdata[grep("MOUSE",rownames(mixedSpeciesdata)),]
            
            local.means.human <- data.frame(rowMeans(human.simData))
            names(local.means.human) <- c("rowMeans")
            rownames(local.means.human)[order(local.means.human$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.human
            local.means.mouse <- data.frame(rowMeans(mouse.simData))
            names(local.means.mouse) <- c("rowMeans")
            rownames(local.means.mouse)[order(local.means.mouse$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.mouse
            
            simData.human <- mixedSpeciesdata[selectedGenes.human,]
            simData.mouse <- mixedSpeciesdata[selectedGenes.mouse,]
            
            simData <- rbind(simData.human,simData.mouse)
            
            human.cell <- colnames(simData)[which(colSums(simData.human) > colSums(simData.mouse))]
            human.cell.simData <- simData[,which(colSums(simData.human) > colSums(simData.mouse))]
            human.cell <- unlist(lapply(human.cell, function(x){ toString(c("HUMAN",x))}))
            colnames(human.cell.simData) <- human.cell
            
            mouse.cell <- colnames(simData)[which(colSums(simData.human) < colSums(simData.mouse))]
            mouse.cell.simData <- simData[,which(colSums(simData.human) < colSums(simData.mouse))]
            mouse.cell <- unlist(lapply(mouse.cell, function(x){ toString(c("MOUSE",x))}))
            colnames(mouse.cell.simData) <- mouse.cell
            
            simData <- cbind(human.cell.simData, mouse.cell.simData)

            # simulation for lasso
            td <- mapply(lasso.mixed.data, rownames(simData),
                         MoreArgs = list(t(simData)))
            theObject@simulation.result.genes[[length(theObject@simulation.result.genes)+1]] <- td
            attr(theObject@simulation.result.genes[[length(theObject@simulation.result.genes)]], "method") <- "lasso"
            
            # simulation for lasso Negative Binomial
            td <- mapply(lasso.negbin.mixed.data, rownames(simData),
                         MoreArgs = list(t(simData)))
            theObject@simulation.result.genes[[length(theObject@simulation.result.genes)+1]] <- td
            attr(theObject@simulation.result.genes[[length(theObject@simulation.result.genes)]], "method") <- "lassoNegBin"
            
            # simulation for kknn
            rownames(simData) <- gsub("[_:-]", "", rownames(simData), perl = TRUE)
            td <- mapply(kknn.mixed.data, rownames(simData),
                         MoreArgs = list(t(as.data.frame(simData))))
            theObject@simulation.result.genes[[length(theObject@simulation.result.genes)+1]] <- td
            attr(theObject@simulation.result.genes[[length(theObject@simulation.result.genes)]], "method") <- "kknn"
            
            return(theObject)
          })

# method to run simulations on cells with dependency on genes
setGeneric(name = "mixedSpeciesCells",
           def = function(theObject, mixedSpeciesdata) {
             standardGeneric("mixedSpeciesCells")
           })

setMethod(f = "mixedSpeciesCells",
          signature = "mixedSpeciesImput",
          definition = function(theObject, mixedSpeciesdata) {
            
            # reduce gene dimension to 100 mouse and 100 human genes with the highest gene expression
            human.simData <- mixedSpeciesdata[grep("HUMAN",rownames(mixedSpeciesdata)),]
            mouse.simData <- mixedSpeciesdata[grep("MOUSE",rownames(mixedSpeciesdata)),]
            
            local.means.human <- data.frame(rowMeans(human.simData))
            names(local.means.human) <- c("rowMeans")
            rownames(local.means.human)[order(local.means.human$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.human
            local.means.mouse <- data.frame(rowMeans(mouse.simData))
            names(local.means.mouse) <- c("rowMeans")
            rownames(local.means.mouse)[order(local.means.mouse$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.mouse
            
            simData.human <- mixedSpeciesdata[selectedGenes.human,]
            simData.mouse <- mixedSpeciesdata[selectedGenes.mouse,]
            
            simData <- rbind(simData.human,simData.mouse)
            
            human.cell <- colnames(simData)[which(colSums(simData.human) > colSums(simData.mouse))]
            human.cell.simData <- simData[,which(colSums(simData.human) > colSums(simData.mouse))]
            human.cell <- unlist(lapply(human.cell, function(x){ toString(c("HUMAN",x))}))
            colnames(human.cell.simData) <- human.cell
            
            mouse.cell <- colnames(simData)[which(colSums(simData.human) < colSums(simData.mouse))]
            mouse.cell.simData <- simData[,which(colSums(simData.human) < colSums(simData.mouse))]
            mouse.cell <- unlist(lapply(mouse.cell, function(x){ toString(c("MOUSE",x))}))
            colnames(mouse.cell.simData) <- mouse.cell
            
            simData <- cbind(human.cell.simData, mouse.cell.simData)

            # simulation for lasso
            td <- mapply(lasso.mixed.data, colnames(simData),
                         MoreArgs = list(simData))
            theObject@simulation.result.cells[[length(theObject@simulation.result.cells)+1]] <- td
            attr(theObject@simulation.result.cells[[length(theObject@simulation.result.cells)]], "method") <- "lasso"
            
            # simulation for lasso Negative Binomial
            td <- mapply(lasso.negbin.mixed.data, colnames(simData),
                         MoreArgs = list(simData))
            theObject@simulation.result.cells[[length(theObject@simulation.result.cells)+1]] <- td
            attr(theObject@simulation.result.cells[[length(theObject@simulation.result.cells)]], "method") <- "lassoNegBin"
            
            
            # simulation for kknn
            colnames(simData) <- gsub("[ ,]", "", colnames(simData), perl = TRUE)
            td <- mapply(kknn.mixed.data, colnames(simData),
                       MoreArgs = list(as.data.frame(simData)))
            theObject@simulation.result.cells[[length(theObject@simulation.result.cells)+1]] <- td
            attr(theObject@simulation.result.cells[[length(theObject@simulation.result.cells)]], "method") <- "kknn"
            
            return(theObject)
          })