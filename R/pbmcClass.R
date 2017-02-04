#' The PBMC class
#' 
#' This class has methods for the purpose of studying the pbmc dataset for
#' gene correlation.
#' @exportClass PBMCGenes

PBMCGenes <- setClass("PBMCGenes",
                      slots = 
                        c(
                          result.genes.kknn = "matrix",
                          result.genes.lasso = "matrix",
                          result.genes.rf = "matrix",
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
            # kknn.pbmc(selectedCols[1],  as.data.frame(simData))
            theObject <- PBMCGenes()
            # theObject@result.genes.kknn <- mapply(kknn.pbmc, selectedGenes[60:70], MoreArgs = list(as.data.frame(simData)))
            # theObject@result.genes.rf <- mapply(randomForest.pbmc, selectedGenes, MoreArgs = list(as.data.frame(simData)))
            theObject@result.genes.lasso <- mapply(lasso.pbmc, selectedGenes, MoreArgs = list(as.matrix(simData)))
            return(theObject)
          })





selectedGenes <- intersect(colnames(pbmc)[order(colMeans(pbmc), decreasing = TRUE)], rownames(go.human.adj))[1:500]

# lasso analysis
pbmc.lasso.predict <- mapply(lasso.pbmc, selectedGenes, MoreArgs = list(as.matrix(pbmc[,selectedGenes])))
pbmc.lasso.cor <- cor(pbmc.lasso.predict)

ok.genes <- intersect(colnames(pbmc.lasso.cor), rownames(go.human.adj))

pbmc.lasso.cor.ok <- pbmc.lasso.cor[ok.genes, ok.genes]

pbmc.lasso.pathway.cor <- unlist(lapply(1:ncol(go.human.adj), function(x){
  pathway.genes <- intersect(names(which(go.human.adj[,x] == 1)), ok.genes)
  return(mean(upperTriangle(abs(pbmc.lasso.cor.ok[pathway.genes, pathway.genes]))))
}))

pathways.ok <- which(!is.na(pbmc.lasso.pathway.cor))

pbmc.lasso.pathway.cor.ok <- pbmc.lasso.pathway.cor[pathways.ok]

pbmc.lasso.random.pathway.cor <- unlist(lapply(pathways.ok, function(x){
  # pathway.genes <- intersect(names(which(go.human.adj[,x] == 1)), ok.genes)
  # random.gene <- sample(intersect(names(which(go.human.adj[,x] != 1)), ok.genes), 1)
  # random.genes <- c(sample(pathway.genes, length(pathway.genes)-1), random.gene)
  random.genes <- sample(ok.genes, length(intersect(names(which(go.human.adj[,x] == 1)), ok.genes)))
  return(mean(upperTriangle(abs(pbmc.lasso.cor.ok[random.genes, random.genes]))))
}))

plotData <- data.frame(pbmc.lasso.pathway.cor.ok, pbmc.lasso.random.pathway.cor)

ggplot(plotData, aes(x=pbmc.lasso.pathway.cor.ok, y=pbmc.lasso.random.pathway.cor)) + geom_point() + geom_line(aes(y=pbmc.lasso.pathway.cor.ok))

ggplot(plotData, aes(x=pbmc.lasso.random.pathway.cor)) + geom_density() + labs(x="Correlation between random set of genes") + geom_vline(xintercept = median(pbmc.lasso.pathway.cor.ok), linetype="dashed")

sum(pbmc.lasso.pathway.cor.ok>pbmc.lasso.random.pathway.cor)/length(pbmc.lasso.pathway.cor.ok)


# kknn analysis
pbmc.kknn.predict <- mcmapply(kknn.pbmc, selectedGenes, MoreArgs = list(as.data.frame(pbmc[,selectedGenes])), mc.cores = 4)



# random forest
pbmc.rf.predict <- mcmapply(randomForest.pbmc, selectedGenes, MoreArgs = list(as.data.frame(pbmc[,selectedGenes])), mc.cores = 4)
pbmc.rf.cor <- cor(pbmc.rf.predict)

ok.genes <- intersect(colnames(pbmc.rf.cor), rownames(go.human.adj))

pbmc.rf.cor.ok <- pbmc.rf.cor[ok.genes, ok.genes]

pbmc.rf.pathway.cor <- unlist(lapply(1:ncol(go.human.adj), function(x){
  pathway.genes <- intersect(names(which(go.human.adj[,x] == 1)), ok.genes)
  return(mean(upperTriangle(abs(pbmc.rf.cor.ok[pathway.genes, pathway.genes]))))
}))

pathways.ok <- which(!is.na(pbmc.rf.pathway.cor))

pbmc.rf.pathway.cor.ok <- pbmc.rf.pathway.cor[pathways.ok]

pbmc.rf.random.pathway.cor <- unlist(lapply(pathways.ok, function(x){
  # pathway.genes <- intersect(names(which(go.human.adj[,x] == 1)), ok.genes)
  # random.gene <- sample(intersect(names(which(go.human.adj[,x] != 1)), ok.genes), 1)
  # random.genes <- c(sample(pathway.genes, length(pathway.genes)-1), random.gene)
  random.genes <- sample(ok.genes, length(intersect(names(which(go.human.adj[,x] == 1)), ok.genes)))
  return(mean(upperTriangle(abs(pbmc.rf.cor.ok[random.genes, random.genes]))))
}))

plotData <- data.frame(pbmc.rf.pathway.cor.ok, pbmc.rf.random.pathway.cor)

ggplot(plotData, aes(x=pbmc.lasso.pathway.cor.ok, y=pbmc.lasso.random.pathway.cor)) + geom_point() + geom_line(aes(y=pbmc.lasso.pathway.cor.ok))

ggplot(plotData, aes(x=pbmc.rf.random.pathway.cor)) + geom_density() + labs(x="Correlation between random set of genes") + geom_vline(xintercept = median(pbmc.rf.pathway.cor.ok), linetype="dashed")

sum(pbmc.lasso.pathway.cor.ok>pbmc.lasso.random.pathway.cor)/length(pbmc.lasso.pathway.cor.ok)
