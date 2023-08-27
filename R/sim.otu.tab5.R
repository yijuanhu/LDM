#' Filtered OTU count table of the simulated microbiome samples
#'
#' This table was derived from sim.otu.tab, after filtering out OTUs that are present (having non-zero counts) 
#' in less than five samples. This filter reduced the number of OTUs from 813 to 593.
#'
#' @docType data
#'
#' @usage data("sim.otu.tab5")
#'
#' @format A data frame with 100 observations on 593 variables.
#'
#' @keywords datasets
#'
#' @examples
#' data(sim.otu.tab5)
"sim.otu.tab5"