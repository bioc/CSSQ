#' Plot data distribution histograms
#'
#' This function is to plot data distribution histogram before
#' and after anscombe transformation.
#' 
#' @param countVec Vector containing the count data for a sample.
#' @param ansVec Vector containing the anscombe transformed count data for a 
#' sample.
#' @param sampleName Name of the sample being plotted.
#' 
#' @return None
#' @import ggplot2
#' @importFrom grDevices dev.off pdf

plotDist <- function(countVec,ansVec,sampleName) {
    pdf (paste(sampleName,"plots.pdf",sep="_"));
    plot1 <- ggplot(data.frame(countVec), aes(countVec)) + geom_histogram(colour="black",fill="white",bins=100) + labs(title="Count data",x="Signal level", y = "Frequency")+ theme_classic() + theme(plot.title = element_text(hjust = 0.5)) ;
    plot2 <- ggplot(data.frame(ansVec), aes(ansVec)) + geom_histogram(colour="black",fill="white",bins=100) + labs(title="Anscombe transformed data",x="Anscombe transformd signal level", y = "Frequency")+ theme_classic() + theme(plot.title = element_text(hjust = 0.5));
    print (plot1);print(plot2);
    dev.off();
}