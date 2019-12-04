#' Perform k-means clustering, normalize anscombe data and calculate cluster 
#' variances for a sample. 
#'
#' This function performs normalization on the anscombe transformed data by
#' clustering them using k-means algorithmn and utilizing the information from
#' clusters. It returns an DataFrame object  normalized counts,
#' cluster information and the variance of that cluster for that sample. 
#' 
#' @param ansDataVec Anscombe transformed count data for a sample.
#' @param numClusters A number indicating the number of clusters to use for 
#' k-means clustering. (default: 4)
#'
#' @return DataFrame containing the normalized counts, 
#'    cluster information and the variance of the cluster in the sample. 
#' @seealso \code{\link{normalizeData}} which iterates over this function.
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom stats kmeans
#' @export
#' @examples
#' exCount <- c(1,2,3,4,5,6,7,8,9,10)
#' kmeansEx <- kmeansNormalize(exCount,numClusters=2)
#' kmeansEx

kmeansNormalize <- function(ansDataVec,numClusters=4) {
        y <- kmeans(ansDataVec,numClusters,nstart=20)
        b <- data.frame(y$centers)[order(data.frame(y$centers)$y.centers),,drop=FALSE]
        clusterLabels <- c(letters[seq_len(numClusters)])
        b$id <- clusterLabels
        b$id_o <- row.names(b)
        z <- data.frame(cbind(ansDataVec,y$cluster))
        colnames(z) <- c("vals","cluster")
        for(j in seq_len(numClusters)){
            z[z$cluster==b$id_o[j],][2] <- b$id[j]
        }
        cluster_data <- z$cluster
        calc_max <- (mean(z[z$cluster == clusterLabels[numClusters],]$vals)) + (3 * sd(z[z$cluster == clusterLabels[numClusters],]$vals))
        z$norm <- z$vals/calc_max
        varValues <- data.frame(vapply(seq_len(numClusters), function(x) var(z[z$cluster == clusterLabels[x],]$norm),double(1)))
        row.names(varValues) <- clusterLabels
        variance <- vapply(cluster_data,function(x) as.numeric(varValues[row.names(varValues)==x,][1]),double(1))
        names(variance) <- NULL
        norm_data <- z$norm
        return(list(clus=cluster_data,norm=norm_data,var=variance))
}