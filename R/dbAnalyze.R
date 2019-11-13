#' Performs differential binding analysis
#'
#' This is a wrapper function that performs the different parts of
#' differential binding analysis. Returns a 
#' \code{\link[GenomicRanges]{GRanges-class}} with a calculated P-value and 
#' Fold change for each region.
#' 
#' @param preprocessedData A 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object from \code{\link{preprocessData}}.
#' @param comparison A vector containing the comparison to be made. Names here 
#' need to correspond to the sample groups in the sample file (Eg. c("G1",G2") 
#' means the comparison G1/G2).
#' 
#' @return A \code{\link[GenomicRanges]{GRanges-class}} object 
#' containing the regions along with their P-values and Fold change for the 
#' comparison.
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' exRange <- GRanges(seqnames=c("chr1","chr2","chr3","chr4"),
#' ranges=IRanges(start=c(1000,2000,3000,4000),end=c(1500,2500,3500,4500)));
#' sampleInfo <- read.table(system.file("extdata", "sample_info.txt", 
#' package="CSSQ",mustWork = TRUE),sep="\t",header=TRUE)
#' exCount <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),nrow=4,ncol=4)
#' exData <- SummarizedExperiment(assays = list(list(ansCount=exCount)),
#' rowRanges=exRange,colData=sampleInfo)
#' normExData <- normalizeData(exData,numClusters=2)
#' res <- DBAnalyze(normExData,comparison=c("HSMM","HESC"))
#' res

DBAnalyze <- function(preprocessedData,comparison=c()){
    numSamples <- nrow(colData(preprocessedData));
    otherComparisons <- getComparisons(colData(preprocessedData)[,2],comparison,numSamples);
    trueTstat <- calculateTvalue(preprocessedData,label = colData(preprocessedData)[,2],comparison,numSamples);
    trueFC <- calculateFC(preprocessedData,label = colData(preprocessedData)[,2],comparison,numSamples);
    compare_tstats <- vapply(seq_len(ncol(otherComparisons)), function(x) calculateTvalue(preprocessedData,label=otherComparisons[,x],comparison,numSamples),double(length(rowRanges(preprocessedData))));
    adjPval <- calculatePvalue(trueTstat,compare_tstats);
    regionRange <- rowRanges(preprocessedData);
    values(regionRange) <- cbind(values(regionRange),DataFrame(adj.pval = adjPval),DataFrame(trueFC));
    return(regionRange);
}