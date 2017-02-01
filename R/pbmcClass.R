#' The PBMC class
#' 
#' This class has methods for the purpose of studying the pbmc dataset for
#' gene correlation.
#' @exportClass PBMCGenes

PBMCGenes <- setClass("PBMCGenes",
                      slots = 
                        c(
                          result.genes = "list",
                          gene.cor = "matrix"
                        ))


setGeneric(name = "imputePBMC",
           def = function(pbmcData) {
             standardGeneric("imputePBMC")
           })
setMethod(f = "imputePBMC",
          signature = "ANY",
          definition = function(pbmcData) {
            p <- 0.5
            selectedRows <- sample(1:nrow(pbmcData), p*nrow(pbmcData))
            selectedCols <- sample(colnames(pbmcData), p*ncol(pbmcData))
            simData <- pbmcData[selectedRows, selectedCols]
            print(dim(simData))
            theObject <- PBMCGenes()
            theObject@result.genes <- mapply(kknn.pbmc, selectedCols[1:3], MoreArgs = list(as.data.frame(simData)))
            return(theObject)
          })