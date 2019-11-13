#' Quantify region level count data
#'
#' The input is the set of regions and the sample information. It will 
#' calculate the number of reads falling in each region for each sample.
#' Returns a 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} 
#' object with regions, sample informationa and counts for all samples. 
#' 
#' @param regionBed A bed file containing the list of regions that are being 
#' analyzed.
#' @param sampleInfo Object from \code{\link{preprocessData}} containing 
#' sample information.
#' @param sampleDir Location of the input sample files in `sampleInfo` file. 
#' @param backgroundSubtract Logical indicating if background correction 
#' should be performed .
#' @param ... Additional arguments passed on to \code{\link{getBgSubVal}}.
#' 
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} 
#' containing the regions, 
#'    sample information and counts for all samples.
#' @seealso \code{\link{getBgSubVal}} which this function calls. 
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import GenomicAlignments
#' @importFrom utils read.table
#' @import IRanges
#' @import Rsamtools
#' @export
#' @examples
#' sampleInfo <- read.table(system.file("extdata", "sample_info.txt", 
#' package="CSSQ",mustWork = TRUE),sep="\t",header=TRUE)
#' countData <- getRegionCounts(system.file("extdata", "chr19_regions.bed", 
#' package="CSSQ"),sampleInfo,
#' sampleDir = system.file("extdata", package="CSSQ"))
#' countData
#' head(assays(countData)$countData)
#' colData(countData)
#' rowRanges(countData)

getRegionCounts <- function(regionBed,sampleInfo,sampleDir = ".",backgroundSubtract=TRUE,...){
    regionList <- read.table(regionBed);
    regionRange <- GRanges(seqnames=regionList$V1,ranges=IRanges(start=regionList$V2,end=regionList$V3));
    sampleInfo[,3] <- vapply(sampleInfo[,3], function(x) paste(sampleDir,x,sep="/"),character(1))
    sampleInfo[,5] <- vapply(sampleInfo[,5], function(x) paste(sampleDir,x,sep="/"),character(1))
    analysisInfo <- SummarizedExperiment(rowRanges=regionRange,colData=sampleInfo);
    if (ncol(sampleInfo) < 5){
        print ("No Input/Control sample information provided. Defaulting to no background correction.");
        NormbgSubCounts <- data.frame(vapply(seq_len(nrow(colData(analysisInfo))),function(x) getBgSubVal(analysisInfo,sampleIndex = x,backgroundSubtract=FALSE,...),double(length(regionRange))));
    }
    else{
        NormbgSubCounts <- data.frame(vapply(seq_len(nrow(colData(analysisInfo))),function(x) getBgSubVal(analysisInfo,sampleIndex = x,backgroundSubtract=backgroundSubtract,...),double(length(regionRange))));
    }
    colnames(NormbgSubCounts) <- colData(analysisInfo)[,1];
    countData <- SummarizedExperiment(assays = list(list(countData=NormbgSubCounts)),rowRanges=rowRanges(analysisInfo),colData=colData(analysisInfo));
    return(countData);
}