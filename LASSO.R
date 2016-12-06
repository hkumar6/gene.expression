################ Initialize
Averagexp <- apply(expr,1,mean)
expr500 <- sort(Averagexp, decreasing=TRUE)[1:500]
Nexpr <- expr[which(Averagexp >= min(expr500)),]
################# define drop-out genes
Actb <- Nexpr[which(rownames(Nexpr)=="Actb"),]
Sod1 <- Nexpr[which(rownames(Nexpr)=="Sod1"),]
Gapdh <- Nexpr[which(rownames(Nexpr)=="Sod1"),]
################# let them drop-out
k <- round(0.05*length(Actb))
vek <- round(runif(k,1,405)) #randomly set 5% to zero

Actb.missing <- Actb[-vek]
Actb.fill <- Actb[vek]
Sod1.missing <- Sod1[-vek]
Sod1.fill <- Sod1[vek]
Gapdh.missing <- Gapdh[-vek]
Gapdh.fill <- Gapdh[vek]

##########################    BEGINNING LASSO

#THEORY
#Elastic net
#the penality is defined for 0 <= alpha <= 1 with
#(1-alpha)/2 ||beta||_2^2 + alpha ||beta||_1

#########################   LASSO approach & Prediction Actb 
Nexpr <- Nexpr[-(which(rownames(Nexpr)=="Actb")),]
Ntexpr <- t(Nexpr)

Nexpr.train <- Nexpr.Gapdh[,-vek]
Nexpr.test <- Nexpr.Gapdh[,vek]
Ntexpr.train <- t(Nexpr.train)
Ntexpr.test <- t(Nexpr.test)

fit.expr <- glmnet(Ntexpr.train,Actb.missing,family = "gaussian",standardize = TRUE, alpha = 1)
plot.glmnet(fit.expr,xvar = "lambda")
CV.expr <- cv.glmnet(Ntexpr.train,Actb.missing, type.measure="mse")

Actb.vek.predicted <- predict(fit.expr, newx = Ntexpr.test, s = c(CV.expr$lambda.min))

res <- sum((Actb.fill-as.vector(Actb.vek.predicted))^2)/length(Actb.vek.predicted)

Actb.predicted <- c(Actb[-vek],Actb.vek.predicted)
Actb.new <- c(Actb[-vek],Actb[vek])
plot(Actb.vek.predicted)
points(Actb[vek], col = "red")

#########################   LASSO approach & Prediction Sod1 
Nexpr <- Nexpr[-(which(rownames(Nexpr)=="Sod1")),]
Ntexpr <- t(Nexpr)

Nexpr.train <- Nexpr.Gapdh[,-vek]
Nexpr.test <- Nexpr.Gapdh[,vek]
Ntexpr.train <- t(Nexpr.train)
Ntexpr.test <- t(Nexpr.test)

fit.expr <- glmnet(Ntexpr.train, Sod1.missing, family = "gaussian",standardize = TRUE, alpha = 1)
plot.glmnet(fit.expr, xvar = "lambda")
CV.expr <- cv.glmnet(Ntexpr.train, Sod1.missing,family = "gaussian", type.measure="mse")

Sod1.vek.predicted <- predict(fit.expr, newx = Ntexpr.test, s = c(CV.expr$lambda.min))

sum(coef(fit.expr, CV.expr$lambda.min)>0)
res <- sum((Sod1.fill-as.vector(Sod1.vek.predicted))^2)

plot(Sod1.vek.predicted)
points(Sod1[vek], col = "red")

#########################   LASSO approach & Prediction Gapdh
Nexpr.Gapdh <- Nexpr[-(which(rownames(Nexpr)=="Gapdh")),]

Nexpr.train <- Nexpr.Gapdh[,-vek]
Nexpr.test <- Nexpr.Gapdh[,vek]
Ntexpr.train <- t(Nexpr.train)
Ntexpr.test <- t(Nexpr.test)

fit.expr <- glmnet(Ntexpr.train, Gapdh.missing, family = "gaussian",standardize = TRUE, alpha = 1)
plot.glmnet(fit.expr, xvar = "lambda")
CV.expr <- cv.glmnet(Ntexpr.train, Gapdh.missing,family = "gaussian", type.measure="mse")

Gapdh.vek.predicted <- predict(fit.expr, newx = Ntexpr.test, s = c(CV.expr$lambda.min))

sum(coef(fit.expr, CV.expr$lambda.min)>0)
res <- sum((Gapdh.fill-as.vector(Gapdh.vek.predicted))^2)

plot(Gapdh.vek.predicted)
points(Gapdh[vek], col = "red")

####################    mrse with different samples
gene <- Gapdh

cov_sample <- c()
mrse_sample <- c()
for (i in 1:40){
  vek <- round(runif(k,1,405))
  gene.missing <- gene[-vek]
  gene.fill <- gene[vek]
  Nexpr.train <- Nexpr[,-vek]
  Nexpr.test <- Nexpr[,vek]
  Ntexpr.train <- t(Nexpr.train)
  Ntexpr.test <- t(Nexpr.test)
  fit.expr <- glmnet(Ntexpr.train,gene.missing,family = "gaussian",standardize = TRUE)
  CV.expr <- cv.glmnet(Ntexpr.train,gene.missing,family = "gaussian", type.measure="mse")
  gene.vek.predicted <- predict(fit.expr, newx = Ntexpr.test, s = c(CV.expr$lambda.min))
  cov <- sum(coef(fit.expr, CV.expr$lambda.min)>0)
  cov_sample <- c(cov_sample, cov)
  res <- sqrt(sum((gene.fill-as.vector(gene.vek.predicted))^2)*(1/20))
  mrse_sample <- c(mrse_sample, res)
}
plot(sort(mrse_sample))

##################    mrse for different drop-out percentage

cov_perc <- c()
mrse_perc <- c()
for (i in 1:20){
  k <- round((i/100)*length(gene))
  vek <- round(runif(k,1,405)) 
  gene.missing <- gene[-vek]
  gene.fill <- gene[vek]
  Nexpr.train <- Nexpr[,-vek]
  Nexpr.test <- Nexpr[,vek]
  Ntexpr.train <- t(Nexpr.train)
  Ntexpr.test <- t(Nexpr.test)
  fit.expr <- glmnet(Ntexpr.train,gene.missing,family = "gaussian",standardize = TRUE)
  CV.expr <- cv.glmnet(Ntexpr.train,gene.missing,family = "gaussian", type.measure="mse")
  gene.vek.predicted <- predict(fit.expr, newx = Ntexpr.test, s = c(CV.expr$lambda.min))
  cov <- sum(coef(fit.expr, CV.expr$lambda.min)>0)
  cov_perc <- c(cov_perc, cov)
  res <- sqrt(sum((gene.fill-as.vector(gene.vek.predicted))^2)*(1/20))
  mrse_perc <- c(mrse_perc, res)
}

plot(mrse_perc)
########################################   plot genes

dev.off()
par(oma=c(1,1,3,1),mar=c(2,2,3,2),mfrow=c(2,2))
ind <- unlist(lapply(vek, function(x){which(sort(Gapdh)==Gapdh[x])}))
plot(sort(Gapdh), main = "sorted Gapdh at 5%")
points(ind,Gapdh.vek.predicted,col="red")
plot(Gapdh.vek.predicted,Gapdh.fill, main = "predicted vs. real values at 5%", xlab = "predicted",
     ylab = "real", xlim = c(min(Gapdh.fill,Gapdh.vek.predicted),max(Gapdh.fill,Gapdh.vek.predicted)), 
     ylim = c(min(Gapdh.fill,Gapdh.vek.predicted),max(Gapdh.fill,Gapdh.vek.predicted)))
lines(1:15,1:15,col ="red")
ind <- unlist(lapply(mrse_sample,function (x){which(sort(mrse_sample)==x)}))
plot(sort(mrse_sample), main = "MRSE - different samples at 5%", col = "blue"
     ,ylim = c(min(mrse_sample,(cov_sample[ind])*(max(mrse_sample)/max(cov_sample))),max(mrse_sample)))
par(new=TRUE)
plot(cov_sample[ind],type = "l", axes = FALSE)
axis(4,ylim=c(min(cov_sample),max(cov_sample)))
par(new = FALSE)
plot(mrse_perc, main = "MRSE - drop-out percentage", col = "blue")
par(new=TRUE)
plot(cov_perc, axes = FALSE, type = "l")
axis(4,ylim = c(min(cov_perc),max(cov_perc)))
mtext(text="Gapdh",side=3,line=0,outer=TRUE)




