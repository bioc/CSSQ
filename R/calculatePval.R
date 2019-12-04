#' Calculates P-value for the regions
#'
#' This calculates the adjusted P-values for the regions using column 
#' permutation and Benjamini Hochberg correction methods. 
#' 
#' @param trueTstat The T-statistics value calculated using 
#' \code{\link{calculateTvalue}} for the intended comparison.
#' @param compare_tstats The T-statistics value calculated using 
#' \code{\link{calculateTvalue}} for all other comparisons possible.
#' 
#' @return A vector that is the adjusted P-value for the intended comparison.
#' @seealso \code{\link{DBAnalyze}} which calls this function 
#' 
#' @importFrom stats p.adjust


calculatePvalue <- function(trueTstat,compare_tstats){
    abs_compare_tstats <- abs(compare_tstats)
    num_higher <- vapply(abs(trueTstat),function(x) (sum(abs_compare_tstats > x)),integer(1))
    numCompare <- ncol(compare_tstats)
    pval <- num_higher/(length(trueTstat) * numCompare)
    adjpval <- p.adjust(pval,method="BH",length(pval))
    return(adjpval)
}