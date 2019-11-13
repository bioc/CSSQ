#' Identify possible combinations
#'
#' This function creates a data frame of all possbile combinations of sample 
#' labels.
#' This information is utilized by \code{\link{calculatePvalue}} to calculate 
#' the P-value using column permutation method.
#' 
#' @param trueLabel The true labels for the samples.
#' @param comparison A vector containing the comparison to be made. Names here 
#' need to correspond to the sample groups in the sample file (Eg. c("G1",G2") 
#' means the comparison G1/G2).
#' @param numSamples Number of samples in the dataset.
#' 
#' @return A data frame with possible combinations of samples other the true 
#' intended comparison.
#' @seealso \code{\link{DBAnalyze}} which calls this function and 
#' \code{\link{getNewLabels}} which this function calls
#' 
#' @importFrom utils combn

getComparisons <- function(trueLabel,comparison,numSamples){
    combns <- combn(seq_len(numSamples),numSamples/2);
    allLabels <- data.frame(vapply(seq_len(ncol(combns)/2),function(x) getNewLabels(trueLabel,comparison,numSamples,combns,x),character(numSamples)));
    returnLabels <- allLabels[,vapply(seq_len(ncol(allLabels)), function(x) !all(allLabels[,x] == trueLabel),logical(1))];
    return(returnLabels);
}