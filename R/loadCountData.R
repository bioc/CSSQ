#' Load count data from input file. 
#'
#' It converts input count file and a bed file with regions 
#' into a  \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object. 
#' 
#' @param countFile A path to file containing the count data for the dataset. 
#' This should be a tab separated file with sample names as header. 
#' @param regionBed A bed file containing the list of regions that are being 
#' analyzed.
#' @param sampleInfo Object from \code{\link{preprocessData}} containing sample
#' information.
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object containing
#' the region information,
#' sample information and the count data. 
#' 
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom utils read.table
#' @import IRanges
#' @export
#' @examples
#' countData <- loadCountData(system.file("extdata", "sample_count_data.txt", 
#' package="CSSQ",mustWork = TRUE),system.file("extdata", "chr19_regions.bed",
#' package="CSSQ"),
#' read.table(system.file("extdata", "sample_info.txt", package="CSSQ",
#' mustWork = TRUE),
#' sep="\t",header=TRUE))
#' countData

loadCountData <- function(countFile,regionBed,sampleInfo) {
    regionList <- read.table(regionBed);
    regionRange <- GRanges(seqnames=regionList$V1,ranges=IRanges(start=regionList$V2,end=regionList$V3));
    countData <- read.table(countFile,sep="\t",header=TRUE);
    analysisInfo <- SummarizedExperiment(assays = list(countData=countData),rowRanges=regionRange,colData=sampleInfo);
    return(analysisInfo);
}