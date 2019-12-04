#' Wrapper function to preprocess the data
#'
#' This is a wrapper function that calls the functions to preprocess the data.
#' It results in a 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} object
#'  normalized counts and meta data that can be used by  
#' \code{\link{DBAnalyze}}.
#'
#' @param inputRegions A bed file  the regions to analyze.
#' @param sampleInfoFile A tab separated file  all sample information. 
#' The following are the columns that are present in the file. 
#' * Sample Name : Names for the samples. 
#' * Group : The group the sample belongs.
#' * IP : The name of the sample bam file. 
#' * IP_aligned_reads : The number of aligned reads in the sample. This is used
#' in depth normalization process.
#' * IN : The name of the sample's control bam file.
#' * IN_aligned_reads : The number of aligned reads in the control file. This 
#' is used in depth normalization process.
#' @param sampleDir Location of the input sample files in `sampleInfoFile` file.
#' (default: ".")
#' Name,Group/Label,IP bam location,IP number of reads,IN bam location,
#' IN number of reads).
#' @param inputCountData The path to the file  count data. This parameter 
#' is used when directly loading count data from a file. This should be a tab 
#' separated file  sample names as header. 
#' @param numClusters A numerical parameter indicating the number of clusters 
#' to use in the normalization step. Passed on to \code{\link{normalizeData}}.
#' (default: 4)
#' @param noNeg A logical parameter indicating how to deal  negative 
#' values. It is passed to \code{\link{ansTransform}}. (default: TRUE)
#' @param plotDataToPDF A logical parameter indicating whether to make plots of the 
#' data distribution to a separate PDF file for each sample.
#' It is passed on to passed to \code{\link{ansTransform}}.
#' (default: FALSE)
#' @param ... Additional arguments passed on to \code{\link{getRegionCounts}}.
#' 
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' containing the normalized counts, cluster information, the variance of the 
#' cluster in the sample and metadata. 
#' @seealso \code{\link{getRegionCounts}}, \code{\link{ansTransform}} and 
#' \code{\link{normalizeData}} which this function calls
#' 
#' @importFrom utils read.table
#' @export
#' @examples
#' processedData <- preprocessData(system.file("extdata", "chr19_regions.bed",
#' package="CSSQ"),system.file("extdata", "sample_info.txt", package="CSSQ"),
#' sampleDir = system.file("extdata", package="CSSQ"),
#' numClusters=4,noNeg=TRUE,plotDataToPDF=FALSE)
#' processedData

preprocessData <- function(inputRegions,sampleInfoFile,sampleDir = ".",inputCountData,numClusters=4,noNeg=TRUE,plotDataToPDF=FALSE,...){
    sampleInfo <- read.table(sampleInfoFile,header=TRUE,sep="\t")
    if (missing(inputCountData)){
    countData <- getRegionCounts(inputRegions,sampleInfo,sampleDir=sampleDir,...)
    }
    else{
    countData <- loadCountData(inputCountData,inputRegions,sampleInfo)
    }
    ansCountData <- ansTransform(countData,noNeg=noNeg,plotDataToPDF=plotDataToPDF)
    normInfo <- normalizeData(ansCountData,numClusters=numClusters)
    return (normInfo)
}