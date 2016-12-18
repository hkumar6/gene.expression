#that when you use the data.frame() function, character variables are imported 
#as factors or categorical variables. Use the str() function to get to know more about your data frame.
#str(writers_df)

# folder
setwd("~/documents/progetto_biocomputazionale")
# upload dataset
load("~/documents/progetto_biocomputazionale/treutlein2016.rdata")
#instal randomForest with install.packages("randomForest")
library("randomForest")
# Load the party package. install.packages("party") It will automatically load other required packages.
library(party)
#install ggplot2 instal kknn with install.packages("ggplot2")
library("ggplot2")
# prepare data for simulation
# prepare data for simulation
dataMatrix=data.frame(matrix(expr))
names(dataMatrix) = c("matrix")
selectedGenes = rownames(dataMatrix)[order(dataMatrix$matrix, decreasing = TRUE)[1:22524]] 
selectedCells = sample(colnames(expr), 0.10*ncol(expr)) #0.10 then, 0.20, ....,0.70

# prepare data for simulation
expr[selectedGenes ] -> simData
simData[ -which(colnames(simData) %in% selectedCells)] -> simData.learn
simData[ selectedCells] -> simData.test



# transpose the data for randomForest library
data.frame(t(simData.test)) -> simData.test
data.frame(t(simData.learn)) -> simData.learn

#built predictive model
fit <- randomForest(expr ~ ., formula=as.formula(Gadph), data=expr, ntree=2000, ntry=10)
#fit <- randomForest(Gapdh ~ ., data=expr,importance=TRUE,proximity=TRUE, ntree=2000, ntry=10)










model <- randomForest(taste ~ . - quality, data = train)



# transpose the data for randomForest library





#predict.randomForest predict method for random forest objects
#Description
#Prediction of test data using random forest.
#Usage
## S3 method for class  randomForest 
predict(object, newdata, type="response",
        norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE,
        cutoff, ...)