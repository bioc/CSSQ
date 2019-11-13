#' Normalize anscombe transformed data
#'
#' This function iterates over \code{\link{kmeansNormalize}} to perform 
#' normalization for all samples in the dataset. It returns an 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} 
#' object with normalized counts, cluster information and the variance of that 
#' cluster for that sample. 
#' 
#' @param ansData 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} object 
#' from \code{\link{ansTransform}}.
#' @param numClusters A number indicating the number of clusters to use for 
#' k-means clustering.
#'
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' containing the 
#' normalized counts, cluster information and the variance of the cluster in 
#' the sample. 
#' @seealso \code{\link{kmeansNormalize}} which this function calls. 
#' @import GenomicRanges
#' @import SummarizedExperiment
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
#' assays(normExData)$normCount

normalizeData <- function(ansData,numClusters=4) {
    info <- vapply (seq_len(ncol(assays(ansData)$ansCount)), function(x) kmeansNormalize(assays(ansData)$ansCount[,x],numClusters=numClusters),list(character,double,double));
    clusterData <- data.frame(vapply(seq_len(ncol(info)),function(x) unlist(info[,x][1],use.names=FALSE),character(length(rowRanges(ansData)))));
    normCount <- data.frame(vapply(seq_len(ncol(info)),function(x) unlist(info[,x][2],use.names=FALSE),double(length(rowRanges(ansData)))));
    geneVars <- data.frame(vapply(seq_len(ncol(info)),function(x) unlist(info[,x][3],use.names=FALSE),double(length(rowRanges(ansData)))));
    colnames(clusterData) <- colData(ansData)[,1];
    colnames(normCount) <- colData(ansData)[,1];
    colnames(geneVars) <- colData(ansData)[,1];
    normInfo <- SummarizedExperiment(assays = list(list(normCount=normCount,clusterData=clusterData,varData=geneVars)),rowRanges=rowRanges(ansData),colData=colData(ansData));
    return(normInfo);
}