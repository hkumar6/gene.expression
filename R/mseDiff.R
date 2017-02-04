n <- floor(runif(1)*100)
n
lasso10 <- Filter(function(x) {return((attr(x, "method") == "lasso") && (attr(x, "drop-percentage") == 0.1))}, cca_genes@simulation.result.genes)[[n]]
lasso40 <- Filter(function(x) {return((attr(x, "method") == "lasso") && (attr(x, "drop-percentage") == 0.4))}, cca_genes@simulation.result.genes)[[n]]
lasso70 <- Filter(function(x) {return((attr(x, "method") == "lasso") && (attr(x, "drop-percentage") == 0.7))}, cca_genes@simulation.result.genes)[[n]]

kknn10 <- Filter(function(x) {return((attr(x, "method") == "kknn") && (attr(x, "drop-percentage") == 0.1))}, cca_genes@simulation.result.genes)[[n]]
kknn40 <- Filter(function(x) {return((attr(x, "method") == "kknn") && (attr(x, "drop-percentage") == 0.4))}, cca_genes@simulation.result.genes)[[n]]
kknn70 <- Filter(function(x) {return((attr(x, "method") == "kknn") && (attr(x, "drop-percentage") == 0.7))}, cca_genes@simulation.result.genes)[[n]]

rf10 <- Filter(function(x) {return((attr(x, "method") == "randomForest") && (attr(x, "drop-percentage") == 0.1))}, cca_genes@simulation.result.genes)[[n]]
rf40 <- Filter(function(x) {return((attr(x, "method") == "randomForest") && (attr(x, "drop-percentage") == 0.4))}, cca_genes@simulation.result.genes)[[n]]
rf70 <- Filter(function(x) {return((attr(x, "method") == "randomForest") && (attr(x, "drop-percentage") == 0.7))}, cca_genes@simulation.result.genes)[[n]]

lassoMse10 <- unlist(lasso10["mse",])
lassoMse40 <- unlist(lasso40["mse",])
lassoMse70 <- unlist(lasso70["mse",])

kknnMse10 <- unlist(kknn10["mse",])
kknnMse40 <- unlist(kknn40["mse",])
kknnMse70 <- unlist(kknn70["mse",])

rfMse10 <- unlist(rf10["mse",])
rfMse40 <- unlist(rf40["mse",])
rfMse70 <- unlist(rf70["mse",])
# rowIndex <- 1:length(lassoMse)

lk <- lassoMse70 - kknnMse70
lr <- lassoMse70 - rfMse70
kr <- kknnMse70 - rfMse70


# dat <- data.frame(cond=factor(rep(c("lasso10", "kknn10", "rf10", "lasso40", "kknn40", "rf40", "lasso70", "kknn70", "rf70"))), values = c(lassoMse10, kknnMse10, rfMse10, lassoMse40, kknnMse40, rfMse40, lassoMse70, kknnMse70, rfMse70))
# ggplot(dat, aes(x=cond, y=values, fill=cond)) + geom_boxplot(alpha = 0.7)

dat <- data.frame(cond=factor(rep(c("lasso-rf", "kknn-rf", "lasso-kknn"))), values = c(lr, kr, lk))
ggplot(dat, aes(x=cond, y=values, fill=cond)) + geom_boxplot(alpha = 0.7) + coord_flip()



# plotData <- data.frame(lk, lr, kr)
# ggplot(plotData, aes(x=lk, color="lasso-kknn")) + geom_density() + geom_density(aes(x=lr, color="lasso-rf")) + geom_density(aes(x=kr, color="kknn-rf")) + labs(x="Differences")
