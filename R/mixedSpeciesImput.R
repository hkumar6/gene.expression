#' The mixedSpeciesImput class
#'
#' This class has slots to store results of different simulations
#' using Lasso, Random Forest and kknn methods.
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
    
    #Reduce Number of Cells
    #Seperate between human and mouse genes
    human.simData <- mixedSpeciesdata[grep("HUMAN",rownames(mixedSpeciesdata)),]
    mouse.simData <- mixedSpeciesdata[grep("MOUSE",rownames(mixedSpeciesdata)),]
    simData <- rbind(human.simData,mouse.simData)
    
    #count number of gene expression along cells seperately for human and mouse genes
    human.counts <- colSums(human.simData != 0)
    mouse.counts <- colSums(mouse.simData != 0)
    
    #exclude cells that have a very similar count of gene expression in both species
    mixed.cells <- sort(abs(human.counts - mouse.counts))[1:20]
    mixed.cells <- names(mixed.cells)
    simData <- simData[,!(names(simData) %in% mixed.cells)]
    human.counts <- human.counts[!(names(human.counts) %in% mixed.cells)]
    mouse.counts <- mouse.counts[!(names(mouse.counts) %in% mixed.cells)]
    
    #seperate the data according to the counts
    human.cells <- colnames(simData)[which(human.counts > mouse.counts)]
    mouse.cells <- colnames(simData)[which(human.counts < mouse.counts)]
    human.cell.simData <- simData[,which(human.counts > mouse.counts)]
    mouse.cell.simData <- simData[,which(human.counts < mouse.counts)]
    
    #rename the cells starting with HUMAN or MOUSE
    human.cells <- unlist(lapply(human.cells, function(x){ toString(c("HUMAN",x))}))
    mouse.cells <- unlist(lapply(mouse.cells, function(x){ toString(c("MOUSE",x))}))
    colnames(mouse.cell.simData) <- mouse.cells
    colnames(human.cell.simData) <- human.cells
    
    #select 100 most expressed cells from human and 100 from mouse
    #local.means.human <- data.frame(colMeans(human.cell.simData))
    #names(local.means.human) <- c("colMeans")
    #rownames(local.means.human)[order(local.means.human$colMeans, decreasing = TRUE)[1:100]] -> selectedCells.human
    #local.means.mouse <- data.frame(colMeans(mouse.cell.simData))
    #names(local.means.mouse) <- c("colMeans")
    #rownames(local.means.mouse)[order(local.means.mouse$colMeans, decreasing = TRUE)[1:100]] -> selectedCells.mouse
    
    #simData.human <- human.cell.simData[,selectedCells.human]
    #simData.mouse <- mouse.cell.simData[,selectedCells.mouse]
    
    #unite both data sets
    simData <- cbind(human.cell.simData, mouse.cell.simData)
    
    #Reduce Number of Genes (100 mouse and 100 human most expressed genes)
    local.means.human <- data.frame(rowMeans(human.simData))
    names(local.means.human) <- c("rowMeans")
    rownames(local.means.human)[order(local.means.human$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.human
    local.means.mouse <- data.frame(rowMeans(mouse.simData))
    names(local.means.mouse) <- c("rowMeans")
    rownames(local.means.mouse)[order(local.means.mouse$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.mouse
    
    simData.human <- simData[selectedGenes.human,]
    simData.mouse <- simData[selectedGenes.mouse,]
    
    simData <- rbind(simData.human,simData.mouse)
    
    # simulation for random forest
    rownames(simData) <- gsub("[_:-]", "", rownames(simData), perl = TRUE)
    td <- mapply(randomForest.mixed.data, rownames(simData),
    MoreArgs = list(t(as.data.frame(simData))))
    theObject@simulation.result.genes[[length(theObject@simulation.result.genes)+1]] <- td
    attr(theObject@simulation.result.genes[[length(theObject@simulation.result.genes)]], "method") <- "randomForest"
    
    
    simData <- cbind(human.cell.simData, mouse.cell.simData)
    
    # simulation for random forest
    rownames(simData) <- gsub("[_:-]", "", rownames(simData), perl = TRUE)
    td <- mapply(randomForest.mixed.data, rownames(simData),
    MoreArgs = list(t(as.data.frame(simData))))
    theObject@simulation.result.genes[[length(theObject@simulation.result.genes)+1]] <- td
    attr(theObject@simulation.result.genes[[length(theObject@simulation.result.genes)]], "method") <- "randomForest"
    
    
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
    
    #Reduce Number of Cells
    #Seperate between human and mouse genes
    human.simData <- mixedSpeciesdata[grep("HUMAN",rownames(mixedSpeciesdata)),]
    mouse.simData <- mixedSpeciesdata[grep("MOUSE",rownames(mixedSpeciesdata)),]
    simData <- rbind(human.simData,mouse.simData)
    
    #count number of gene expression along cells seperately for human and mouse genes
    human.counts <- colSums(human.simData != 0)
    mouse.counts <- colSums(mouse.simData != 0)
    
    #exclude cells that have a very similar count of gene expression in both species
    mixed.cells <- sort(abs(human.counts - mouse.counts))[1:20]
    mixed.cells <- names(mixed.cells)
    simData <- simData[,!(names(simData) %in% mixed.cells)]
    human.counts <- human.counts[!(names(human.counts) %in% mixed.cells)]
    mouse.counts <- mouse.counts[!(names(mouse.counts) %in% mixed.cells)]
    
    #seperate the data according to the counts
    human.cells <- colnames(simData)[which(human.counts > mouse.counts)]
    mouse.cells <- colnames(simData)[which(human.counts < mouse.counts)]
    human.cell.simData <- simData[,which(human.counts > mouse.counts)]
    mouse.cell.simData <- simData[,which(human.counts < mouse.counts)]
    
    #rename the cells starting with HUMAN or MOUSE
    human.cells <- unlist(lapply(human.cells, function(x){ toString(c("HUMAN",x))}))
    mouse.cells <- unlist(lapply(mouse.cells, function(x){ toString(c("MOUSE",x))}))
    colnames(mouse.cell.simData) <- mouse.cells
    colnames(human.cell.simData) <- human.cells
    
    #select 100 most expressed cells from human and 100 from mouse
    #local.means.human <- data.frame(colMeans(human.cell.simData))
    #names(local.means.human) <- c("colMeans")
    #rownames(local.means.human)[order(local.means.human$colMeans, decreasing = TRUE)[1:100]] -> selectedCells.human
    #local.means.mouse <- data.frame(colMeans(mouse.cell.simData))
    #names(local.means.mouse) <- c("colMeans")
    #rownames(local.means.mouse)[order(local.means.mouse$colMeans, decreasing = TRUE)[1:100]] -> selectedCells.mouse
    
    #simData.human <- human.cell.simData[,selectedCells.human]
    #simData.mouse <- mouse.cell.simData[,selectedCells.mouse]
    
    #unite both data sets
    simData <- cbind(human.cell.simData, mouse.cell.simData)
    
    #Reduce Number of Genes (100 mouse and 100 human most expressed genes)
    local.means.human <- data.frame(rowMeans(human.simData))
    names(local.means.human) <- c("rowMeans")
    rownames(local.means.human)[order(local.means.human$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.human
    local.means.mouse <- data.frame(rowMeans(mouse.simData))
    names(local.means.mouse) <- c("rowMeans")
    rownames(local.means.mouse)[order(local.means.mouse$rowMeans, decreasing = TRUE)[1:100]] -> selectedGenes.mouse
    
    simData.human <- simData[selectedGenes.human,]
    simData.mouse <- simData[selectedGenes.mouse,]
    
    simData <- rbind(simData.human,simData.mouse)
    # simulation for random forest
    td <- mapply(randomForest.mixed.data, colnames(simData),
    MoreArgs = list(simData))
    theObject@simulation.result.cells[[length(theObject@simulation.result.cells)+1]] <- td
    attr(theObject@simulation.result.cells[[length(theObject@simulation.result.cells)]], "method") <- "randomForest"
    
    
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