#' Perform quantification and normalization of count data
#'
#' This function quantifies each each region for a sample and performs
#' background correction and normalization as instructed. 
#' Returns a vector of count information for the input regions.
#' 
#' @param analysisInfo A 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} 
#' object with regions and sample information from within 
#' \code{\link{getRegionCounts}}.
#' @param sampleIndex Index of the sample to process.
#' @param normalizeReadDepth Logical indicating if count data should be 
#' normalized for library sequencing depth. When TRUE (default), counts will 
#' be normalized for sequencing depth for each library. When FALSE, no such 
#' normalization will be performed and raw counts will be used. 
#' @param normalizeLength Logical indicating if count data should be 
#' normalized to the length of the regions. When TRUE, count data will be 
#' normalized for the length of the region being analyzed. 
#' When FALSE (default), no such normalization will be performed. 
#' @param backgroundSubtract Logical indicating if background correction 
#' should be performed.
#' When TRUE (default), background subtraction will be performed after 
#' length and depth normalization if applicable. When FALSE, no background 
#' subtraction will be performed.  
#' @param countMode Count method passed on to 
#' \code{\link[GenomicAlignments]{summarizeOverlaps}}
#' (\code{"Union"}, \code{"IntersectionNotEmpty"}, \code{"IntersectionStrict"} 
#' or \code{"user supplied function"}).
#' @param ignore.strand A logical indicating if strand should be considered 
#' when matching. 
#' Passed on to \code{\link[GenomicAlignments]{summarizeOverlaps}}.
#' @param inter.feature A logical indicating if the `r countMode` should be 
#' aware of overlapping features. When TRUE, reads mapping to multiple 
#' features are dropped (i.e., not counted). When FALSE (default), these reads
#' are retained and a count is assigned to each feature they map to. Passed on 
#' to \code{\link[GenomicAlignments]{summarizeOverlaps}}.
#' 
#' @return A vector containing the counts for all the regions. 
#' @seealso \code{\link{getRegionCounts}} which calls this function 
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import SummarizedExperiment
#' @import Rsamtools
#' @export
#' @examples
#' regionBed <- read.table(system.file("extdata", "chr19_regions.bed", 
#' package="CSSQ",mustWork = TRUE))
#' sampleInfo <- read.table(system.file("extdata", "sample_info.txt", 
#' package="CSSQ",mustWork = TRUE),sep="\t",header=TRUE)
#' sampleInfo[,3] <- sapply(sampleInfo[,3], 
#' function(x) system.file("extdata", x, package="CSSQ"))
#' sampleInfo[,5] <- sapply(sampleInfo[,5], 
#' function(x) system.file("extdata", x, package="CSSQ"))
#' regionRange <- GRanges(seqnames=regionBed$V1,
#' ranges=IRanges(start=regionBed$V2,end=regionBed$V3));
#' analysisInfo <- SummarizedExperiment(rowRanges=regionRange,
#' colData=sampleInfo)
#' NormbgSubCounts <- data.frame(sapply(c(1:nrow(colData(analysisInfo))),
#' function(x) getBgSubVal(analysisInfo,sampleIndex = x,backgroundSubtract=TRUE,
#' normalizeReadDepth=TRUE,normalizeLength=FALSE,countMode="Union", 
#' ignore.strand=TRUE,inter.feature=FALSE)))
#' NormbgSubCounts

getBgSubVal <- function(analysisInfo,sampleIndex,normalizeReadDepth=TRUE,normalizeLength=FALSE,backgroundSubtract=TRUE,countMode="Union", ignore.strand=TRUE,inter.feature=FALSE){
    params <- ScanBamParam(what = scanBamWhat(), flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE, isDuplicate = FALSE));
    if (backgroundSubtract == TRUE){
        in_val <- summarizeOverlaps(features=rowRanges(analysisInfo),reads=as.character(colData(analysisInfo)[,5][sampleIndex]),mode=countMode,ignore.strand=TRUE,param=params,inter.feature=inter.feature);
        if (normalizeReadDepth == TRUE){
            in_val <- assay(in_val)*(1000000/colData(analysisInfo)[,6][sampleIndex]);
        }
        else{
            in_val <- assay(in_val)[1];
        }
    }
    else{
        in_val <- 0;
    }
    ip_val <- summarizeOverlaps(features=rowRanges(analysisInfo),reads=as.character(colData(analysisInfo)[,3][sampleIndex]),mode=countMode,ignore.strand=ignore.strand,param=params,inter.feature=inter.feature);
    if (normalizeReadDepth == TRUE){
        ip_val <- assay(ip_val)*(1000000/colData(analysisInfo)[,4][sampleIndex]);
    }
    else{
        ip_val <- assay(ip_val)[1];
    }
    returnVal <- as.numeric(ip_val-in_val);
    if (normalizeLength == TRUE){
        width = width(ranges(rowRanges(analysisInfo)));
        returnVal <- returnVal / width;
    }
    return(returnVal);
}