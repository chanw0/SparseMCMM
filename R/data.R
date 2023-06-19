#' Simulated Data
#' A simulated datalist for test and illustration. There are 100 subjects (50 cases and 50 controls) and 10 taxa.
#'
#' @format A data list with three components:
#' \describe{
#'   \item{Treatment}{A treatment vector for 100 subjects. 0 represents the control group and 1 represents the case group.}
#'   \item{otu.com}{A 100 x 10 numeric matrix containing compositional microbiome data. Each row represents a subject, and each column represents a
#'   taxon. The row sum equals 1.}
#'   \item{outcome}{A outcome vector for 100 subjects.}
#' }
"SimulatedData"
