Averagexp <- apply(expr,1,mean)
expr500 <- sort(Averagexp, decreasing=TRUE)[1:500]
Nexpr <- expr[which(Averagexp >= min(expr500)),]
Ntexpr <- t(Nexpr)
##################################################
require(glmnet)
library(Hmisc) # for the function rcorr for the spearman coefficient
k <- round((0.1)*405)
vek <- round(runif(k,1,405))

LASSO <- function(gene){
  Ntexpr.gen <- Ntexpr[,-c(which(colSums(Ntexpr) == sum(gene)))] #checked for non-equal sums (not good but works)
  gene.missing <- gene[-vek]
  gene.fill <- gene[vek]
  Ntexpr.train <- Ntexpr.gen[-vek,]
  Ntexpr.test <- Ntexpr.gen[vek,]
  fit.expr <- glmnet(Ntexpr.train,gene.missing,family = "gaussian",standardize = TRUE)
  CV.expr <- cv.glmnet(Ntexpr.train,gene.missing, type.measure="mse")
  gene.vek.predicted <- predict(fit.expr, newx = Ntexpr.test, s = c(CV.expr$lambda.min))
  mse <- sum((gene.vek.predicted-gene.fill)^2)/k
  c <- rcorr(gene.vek.predicted,gene.fill, type = "spearman")
  Spear_corr <- c$r
  cov <- sum(coef(fit.expr, CV.expr$lambda.min)>0)
  result <- list(mse = mse, Spear_corr = Spear_corr, real = gene.fill, predicted = gene.vek.predicted,
                 number.cov = cov)
  return(result)
}

result <- apply(Ntexpr,2, LASSO)

mse <- unlist(lapply(result, function(x) x$mse))
unlist(lapply(result, function(x) x$number.cov))

Pred_gene <- do.call(cbind, lapply(result, function(x) x$predicted))
colnames(Pred_gene) <- colnames(Ntexpr)

Real_gene <- do.call(cbind, lapply(result, function(x) x$real))
colnames(Real_gene) <- colnames(Ntexpr)

cov <- unlist(lapply(result, function(x) x$number.cov))

