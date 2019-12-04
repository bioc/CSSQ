#' Plot data distribution histograms
#'
#' This function is to plot data distribution histogram before
#' and after anscombe transformation.
#' 
#' @param countData A 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object from \code{\link{getRegionCounts}} with count data.
#' @param ansCount A 
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#' object from \code{\link{ansTransform}} with anscombe transformed data.
#' @param sampleName Name of the sample being plotted.
#' @param plotDataToPDF A logical parameter indicating whether to make plots of the
#' data distribution to a separate PDF file for each sample.
#' When TRUE, a histogram will be plotted for the data before and after 
#' transformation. 
#' When FALSE, no plots will be made. (default: FALSE)
#' 
#' @return A list of the histogram of the count data before and after
#' anscombe transformation if plotDataToPDF == FALSE. 
#' None if plotDataToPDF == TRUE.
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @export
#' @examples
#' exRange <- GRanges(seqnames=c("chr1","chr2","chr3","chr4"),
#' ranges=IRanges(start=c(1000,2000,3000,4000),end=c(1500,2500,3500,4500)))
#' sampleInfo <- read.table(system.file("extdata", "sample_info.txt", 
#' package="CSSQ",mustWork = TRUE),sep="\t",header=TRUE)
#' exCount <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),nrow=4,ncol=4)
#' exData <- SummarizedExperiment(assays = list(countData=exCount),
#' rowRanges=exRange,colData=sampleInfo)
#' ansExData <- ansTransform(exData)
#' plotEx <- plotDist(exData,ansExData,"HESC_R1")
#' plotEx[[1]]

plotDist <- function(countData,ansCount,sampleName,plotDataToPDF=FALSE) {
    sampeIndex <- try(which(as.character(colData(countData)[,1]) == sampleName)[[1]],silent=TRUE)
    if (inherits(sampeIndex, 'try-error')){stop("Invalid sample name")}
    countVec <- assays(countData)$countData[,sampeIndex]
    ansVec <- assays(ansCount)$ansCount[,sampeIndex]
    plot1 <- ggplot(data.frame(countVec), aes(countVec)) + geom_histogram(colour="black",fill="white",bins=100) + labs(title="Non transformed data",x="Signal level", y = "Frequency")+ theme_classic() + theme(plot.title = element_text(hjust = 0.5)) 
    plot2 <- ggplot(data.frame(ansVec), aes(ansVec)) + geom_histogram(colour="black",fill="white",bins=100) + labs(title="Anscombe transformed data",x="Anscombe transformed signal level", y = "Frequency")+ theme_classic() + theme(plot.title = element_text(hjust = 0.5))
    if (plotDataToPDF == FALSE){
    return(list(plot1,plot2))
    }
    else{
    pdf (paste(sampleName,"plots.pdf",sep="_"))
    print(plot1)
    print(plot2)
    dev.off()
    }
}