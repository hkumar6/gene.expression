#' The PBMC class
#' 
#' This class has methods for the purpose of studying the pbmc dataset for
#' gene correlation.
#' @exportClass PBMCGenes

PBMCGenes <- setClass("PBMCGenes",
                      slots = 
                        c(
                          result.genes = "matrix",
                          gene.cor = "matrix"
                        ))


setGeneric(name = "imputePBMC",
           def = function(pbmcData) {
             standardGeneric("imputePBMC")
           })
setMethod(f = "imputePBMC",
          signature = "ANY",
          definition = function(pbmcData) {
            
          })