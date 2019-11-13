#' Perform quantification and normalization of count data
#'
#' This function quantifies each each region for a sample and performs
#' background correction and normalization as instructed. 
#' Returns a vector of count information for the input regions.
#' 
#' @param countData A 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object from \code{\link{getRegionCounts}} with count data.
#' @param noNeg A Logical parameter indicating how to deal with negative 
#' values. 
#' When TRUE (default), all negative values will be moved to 0 before 
#' transforming. 
#' When FALSE, the signs will be maintained while the transformation will be 
#' applied to the absolute value. 
#' @param plotData A logical parameter indicating whether to make plots of the
#' data distribution.
#' When TRUE, a histogram will be plotted for the data before and after 
#' transformation. 
#' When FALSE (default), no plots will be made. 
#' 
#' @return A 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} object 
#' containing the anscombe transformed count data as the assay. 
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @export
#' @examples
#' exRange <- GRanges(seqnames=c("chr1","chr2","chr3","chr4"),
#' ranges=IRanges(start=c(1000,2000,3000,4000),end=c(1500,2500,3500,4500)));
#' sampleInfo <- read.table(system.file("extdata", "sample_info.txt", 
#' package="CSSQ",mustWork = TRUE),sep="\t",header=TRUE)
#' exCount <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),nrow=4,ncol=4)
#' exData <- SummarizedExperiment(assays = list(list(countData=exCount)),
#' rowRanges=exRange,colData=sampleInfo)
#' ansExData <- ansTransform(exData)
#' assays(ansExData)$ansCount

ansTransform <- function(countData,noNeg = TRUE,plotData=FALSE) {
    rawData <- assays(countData)$countData;
    if (noNeg == TRUE){
        rawData[rawData <0] <- 0;
        y <- 2*(sqrt(rawData+(3/8)));
    }
    else{
        y <- (2*(sqrt(abs(rawData)+(0))))*(abs(rawData)/rawData);
    }
    ansCount <- SummarizedExperiment(assays = list(list(ansCount=y)),rowRanges=rowRanges(countData),colData=colData(countData));
    if (plotData == TRUE){
    tmp <- vapply(seq_len(nrow(colData(countData))), function(x) plotDist(assays(countData)$countData[,x],assays(ansCount)$ansCount[,x],as.character(colData(countData)[,1][x])),integer(1));
    }
    return(ansCount);
}
