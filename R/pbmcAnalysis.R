original.cor <- cor(pbmc)

ok.genes <- intersect(colnames(original.cor), rownames(go.human.adj))

original.cor.ok <- original.cor[ok.genes, ok.genes]

original.pathway.cor <- unlist(lapply(1:ncol(go.human.adj), function(x){
  pathway.genes <- intersect(names(which(go.human.adj[,x] == 1)), ok.genes)
  return(mean(upperTriangle(abs(original.cor.ok[pathway.genes, pathway.genes]))))
}))

pathways.ok <- which(!is.na(original.pathway.cor))

original.pathway.cor.ok <- original.pathway.cor[pathways.ok]

random.pathway.cor <- unlist(lapply(pathways.ok, function(x){
  pathway.genes <- intersect(names(which(go.human.adj[,x] == 1)), ok.genes)
  random.gene <- sample(intersect(names(which(go.human.adj[,x] != 1)), ok.genes), 1)
  random.genes <- c(sample(pathway.genes, length(pathway.genes)-1), random.gene)
  # random.genes <- sample(ok.genes, length(intersect(names(which(go.human.adj[,x] == 1)), ok.genes)))
  return(mean(upperTriangle(abs(original.cor.ok[random.genes, random.genes]))))
}))

plotData <- data.frame(original.pathway.cor.ok, random.pathway.cor)

ggplot(plotData, aes(x=original.pathway.cor.ok, y=random.pathway.cor)) + geom_point() + geom_line(aes(y=original.pathway.cor.ok))

ggplot(plotData, aes(x=random.pathway.cor)) + geom_density() + labs(x="Correlation between random set of genes") + geom_vline(xintercept = median(original.pathway.cor.ok), linetype = "dashed")

sum(original.pathway.cor.ok>random.pathway.cor)/length(original.pathway.cor.ok)


# analysis based on common ontologies - no use
n <- 3
results = matrix(nrow=100, ncol=2)
for(i in 1:100) {
  test1 <- sample(ok.genes, n)
  results[i,1] <- sum(which(colSums(go.human.adj[test1,]) == n))
  results[i,2] <- mean(upperTriangle(abs(original.cor.ok[test1, test1])))
}


# pbmc kknn analysis

pbmc.kknn.cor <- cor(pbmc.kknn.predict)

ok.genes <- intersect(colnames(pbmc.kknn.cor), rownames(go.human.adj))

pbmc.kknn.cor.ok <- pbmc.kknn.cor[ok.genes, ok.genes]

pbmc.kknn.pathway.cor <- unlist(lapply(1:ncol(go.human.adj), function(x){
  pathway.genes <- intersect(names(which(go.human.adj[,x] == 1)), ok.genes)
  return(mean(upperTriangle(abs(pbmc.kknn.cor.ok[pathway.genes, pathway.genes]))))
}))

pathways.ok <- which(!is.na(pbmc.kknn.pathway.cor))

pbmc.kknn.pathway.cor.ok <- original.pathway.cor[pathways.ok]

pbmc.kknn.random.pathway.cor <- unlist(lapply(pathways.ok, function(x){
  # pathway.genes <- intersect(names(which(go.human.adj[,x] == 1)), ok.genes)
  # random.gene <- sample(intersect(names(which(go.human.adj[,x] != 1)), ok.genes), 1)
  # random.genes <- c(sample(pathway.genes, length(pathway.genes)-1), random.gene)
  random.genes <- sample(ok.genes, length(intersect(names(which(go.human.adj[,x] == 1)), ok.genes)))
  return(mean(upperTriangle(abs(pbmc.kknn.cor.ok[random.genes, random.genes]))))
}))

plotData <- data.frame(pbmc.kknn.pathway.cor.ok, pbmc.kknn.random.pathway.cor)

ggplot(plotData, aes(x=original.pathway.cor.ok, y=random.pathway.cor)) + geom_point() + geom_line(aes(y=original.pathway.cor.ok))

ggplot(plotData, aes(x=pbmc.kknn.random.pathway.cor)) + geom_density() + labs(x="Correlation between random set of genes") + geom_vline(xintercept = median(pbmc.kknn.pathway.cor.ok), linetype = "dashed")

sum(original.pathway.cor.ok>random.pathway.cor)/length(original.pathway.cor.ok)
