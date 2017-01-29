#' kknn imputation - single column
#' 
#' Impute results for a specific gene or cell using kknn
#' 
#' @param id string for the cell or gene to impute, must be a column name in the data
#' @param simData.learn training data as data.frame
#' @param simData.test testing data, for predictions, as data.frame
#' @param kmaxParam integer specifying the maximum number of neighbors to consider for
#'  training data
#' @param mixedSpeciesData TRUE/FALSE depending on whether the data belongs to the
#'  mixed species experiment. If TRUE, additional measures are included for analysis
#' 
#' @return data.frame containing results of the simulation
#'      mean squared error, Spearman correlation, optimal parameters for kknn
#'
#' @importFrom stats as.formula
#' @importFrom kknn train.kknn kknn
#' @importFrom Hmisc rcorr
#' @export

kknnImpute <- function(id, simData.learn, simData.test, kmaxParam=40, mixedSpeciesData = FALSE, predicted = FALSE) {
  trainingInfo <- train.kknn(as.formula(paste(id, ".", sep = "~")), simData.learn, kmaxParam, kernel = c("triangular", "rectangular", "epanechnikov", "optimal"), distance = 2)
  r <- kknn(as.formula(paste(id, ".", sep = "~")), simData.learn, simData.test, distance = 2, kernel = trainingInfo$best.parameters$kernel, k = trainingInfo$best.parameters$k)
  
  if(predicted) {
    return(r$fitted.values)
  }
  
  mse <- sum((simData.test[,id]-r$fitted.values)^2)/(dim(simData.test[id])[1])
  x <- simData.test[id]
  x$predicted = r$fitted.values
  if (dim(x)[1]>4) {
    c <- rcorr(as.matrix(x), type = "spearman")
    sc <- c$r[1,2]
  } else {
    sc <- NA
  }
  outputList <- data.frame(mse = mse, Spear_corr = sc,
                     optimalK = trainingInfo$best.parameters$k,
                     optimalKernel = trainingInfo$best.parameter$kernel)
  if(TRUE == mixedSpeciesData) {
    outputList$perc.human <- length(grep("HUMAN", rownames(simData.learn)[as.vector(r$C)]))/sum(dim(simData.learn))
    if(length(grep("HUMAN", id)) > 0) {
      zeroPredictions = sum(x$predicted[grep("MOUSE", rownames(simData.test))] == 0)
      outputList$mse = mean(x$predicted[grep("MOUSE", rownames(simData.test))])
    } else {
      zeroPredictions = sum(x$predicted[grep("HUMAN", rownames(simData.test))] == 0)
      outputList$mse = mean(x$predicted[grep("HUMAN", rownames(simData.test))])
    }
    outputList$zero.predicted = zeroPredictions
  }
  return(outputList)
}

#' kknn imputation - multiple columns
#' 
#' Impute results for multiple genes or cells in testing data using kknn
#' 
#' @param selectedList vector of strings for the cells or genes to impute, must be a column name in the data
#' @param simData.learn training data as data.frame
#' @param simData.test testing data, for predictions, as data.frame
#' 
#' @return data.frame containing results of the simulation
#'      mean squared error, Spearman correlation, optimal parameters for kknn
#'
#' @importFrom stats as.formula
#' @importFrom kknn train.kknn kknn
#' @importFrom Hmisc rcorr
#' @export
kknnImputeAll <- function(selectedList, simData.learn, simData.test) {
  simulationResult <- mapply(kknnImpute, selectedList,
                             MoreArgs = list(
                               simData.learn,
                               simData.test
                             ))
  return(simulationResult)
}