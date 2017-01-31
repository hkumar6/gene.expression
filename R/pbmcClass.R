#' The PBMC class
#' 
#' This class has methods for the purpose of studying the pbmc dataset for
#' gene correlation.
#' @exportClass PBMCGenes

PBMCGenes <- setClass("PBMCGenes",
                      slots = 
                        c(
                          result.genes.kknn = "matrix",
                          gene.cor = "matrix"
                        ))


setGeneric(name = "imputePBMC",
           def = function(pbmcData) {
             standardGeneric("imputePBMC")
           })
setMethod(f = "imputePBMC",
          signature = "ANY",
          definition = function(pbmcData) {
            selectedGenes <- colnames(pbmc)[order(colMeans(pbmc), decreasing = TRUE)[1:1000]]
            simData <- pbmcData[,selectedGenes]
            print(dim(simData))
            # kknn.pbmc(selectedCols[1],  as.data.frame(simData))
            theObject <- PBMCGenes()
            theObject@result.genes.kknn <- mapply(kknn.pbmc, selectedGenes, MoreArgs = list(as.data.frame(simData)))
            return(theObject)
          })