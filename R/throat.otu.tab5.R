#' Filtered OTU count table from 16S sequencing of the throat microbiome samples
#'
#' This table was derived from throat.otu.tab, after filtering out OTUs that are present (having non-zero counts) 
#' in less than five samples. This filter reduced the number of OTUs from 856 to 233. 
#'
#' @docType data
#'
#' @usage data("throat.otu.tab5")
#'
#' @format A data frame with 60 observations on 233 variables.
#'
#' @keywords datasets
#'
#' @examples
#' data(throat.otu.tab5)
"throat.otu.tab5"