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
            
            # reduce randomly gene dimension by 50%
            sample(rownames(mixedSpeciesdata), 0.5*nrow(mixedSpeciesdata)) -> selectedGenes
            simData <- mixedSpeciesdata[selectedGenes,]
            
            # simulation for lasso
            td <- mapply(lasso.mixed.data, rownames(simData)[1:3],
                         MoreArgs = list(t(simData)))
            theObject@simulation.result.genes[[length(theObject@simulation.result.genes)+1]] <- td
            attr(theObject@simulation.result.genes[[length(theObject@simulation.result.genes)]], "method") <- "lasso"
            
            # simulation for kknn
            #TODO
            
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
            
            # reduce randomly gene dimension by 50%
            sample(rownames(mixedSpeciesdata), 0.5*nrow(mixedSpeciesdata)) -> selectedGenes
            simData <- mixedSpeciesdata[selectedGenes,]
            selectedGenes <- rownames(simData)
            # simulation for lasso
            td <- mapply(lasso.mixed.data, colnames(simData)[1:3],
                         MoreArgs = list(simData))
            theObject@simulation.result.cells[[length(theObject@simulation.result.cells)+1]] <- td
            attr(theObject@simulation.result.cells[[length(theObject@simulation.result.cells)]], "method") <- "lasso"
            
            # simulation for kknn
            # TODO
            
            return(theObject)
          })