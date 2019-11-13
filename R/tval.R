#' Calculates modified T-statistics values for the given label and comparison
#'
#' This calculates the modified T-statistics for the given comparison.
#' 
#' @param preprocessedData A 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object from \code{\link{preprocessData}}.
#' @param comparison A vector containing the comparison to be made. Names here 
#' need to correspond to the sample groups in the sample file (Eg. c("G1",G2") 
#' means the comparison G1/G2).
#' @param label A vector containing the labels to use for the samples in 
#' preprocessedData.
#' @param numSamples Number of samples in the dataset.
#' 
#' @return A vector that is the modified T-statistics for the comparison and 
#' labels given.
#' @seealso \code{\link{DBAnalyze}} which calls this function 
#' 
#' @import GenomicRanges
#' @import SummarizedExperiment


calculateTvalue <- function(preprocessedData,label,comparison,numSamples){
    norm_data <- assays(preprocessedData)$normCount;
    var_data <- assays(preprocessedData)$varData;
    colnames(norm_data) <- label;
    colnames(var_data) <- label;
    restrucNormData <- cbind(norm_data[,names(norm_data) %in% (comparison[1])],norm_data[,names(norm_data) %in% (comparison[2])]);
    restrucVarMatrix <- cbind(var_data[,names(var_data) %in% (comparison[1])],var_data[,names(var_data) %in% (comparison[2])]);
    group1MeanVar <- rowMeans(restrucVarMatrix[,seq_len(numSamples/2)]);
    group2MeanVar <- rowMeans(restrucVarMatrix[,((numSamples/2)+1):ncol(restrucVarMatrix)]);
    denom <- sqrt(group1MeanVar/(numSamples/2) + group2MeanVar/(numSamples/2));
    combDT <- cbind(restrucNormData,denom);
    tstat <- apply(combDT,1,function(x) ((mean(x[seq_len(numSamples/2)])- mean(x[seq(((numSamples/2)+1),numSamples)]))/x[length(x)]));
    return(tstat);
}
