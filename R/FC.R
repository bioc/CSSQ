#' Calculates Fold change values for the given label and comparison
#'
#' This calculates the Fold change for the given comparison.
#' 
#' @param preprocessedData A 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object from \code{\link{preprocessData}}.
#' @param comparison A vector containing the comparison to be made. 
#' Names here need to correspond to the sample groups in the sample file 
#' (Eg. c("G1",G2") means the comparison G1/G2).
#' @param label A vector containing the labels to use for the samples in 
#' preprocessedData.
#' @param numSamples Number of samples in the dataset.
#' 
#' @return A vector that is the Fold change for the comparison and labels given.
#' @seealso \code{\link{DBAnalyze}} which calls this function 
#' 
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @keywords internal

calculateFC <- function(preprocessedData,label,comparison,numSamples){
    norm_data <- assays(preprocessedData)$normCount
    colnames(norm_data) <- label
    restrucNormData <- cbind(norm_data[,names(norm_data) %in% (comparison[1])],norm_data[,names(norm_data) %in% (comparison[2])])
    FC <- data.frame(apply(restrucNormData,1,function(x) ifelse(mean(x[seq_len(numSamples/2)]) > mean(x[seq(((numSamples/2)+1),numSamples)]), mean(x[seq_len(numSamples/2)]) / mean(x[seq(((numSamples/2)+1),numSamples)]), -1*(mean(x[seq(((numSamples/2)+1),numSamples)])/mean(x[seq_len(numSamples/2)])))))
    colnames(FC)<- "FC"
    return(FC)
}


