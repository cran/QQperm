#' An example collapsing data matrix and the associated case/control status
#'
#' @format A list containing both the genotype data matrix and the case/control statuses:
#'
#'  \describe{
#'     \item{data}{a named data matrix, with a row for each gene and a column for each sample. Each cell in the matrix then represents an indicator function that reflects whether the given sample had none (0) or at least one (1) qualifying variant in the corresponding gene. It also accepts counts of qualifying variants for a given (gene,sample) pair, but will treat these as a (0/1) indicator.}
#'    \item{is.case}{case/control indicator. The order of the indicators here must match that of the order of the samples defined in the data matrix.}
#' }
"example.data"


#' An example collapsing data matrix and the associated case/control status
#'
#' @format A list containing both the genotype data matrix and the case/control statuses:
#'
#'  \describe{
#'     \item{data}{a named data matrix, with a row for each gene and a column for each sample. Each cell in the matrix then represents an indicator function that reflects whether the given sample had none (0) or at least one (1) qualifying variant in the corresponding gene. It also accepts counts of qualifying variants for a given (gene,sample) pair, but will treat these as a (0/1) indicator.}
#'    \item{is.case}{case/control indicator. The order of the indicators here must match that of the order of the samples defined in the data matrix.}
#' }
"toy.data"

#' Distributions of expected and observed P-values from the igm.data dataset.
#'
#' @format A list containing both the observed and the permutation-based expected P-Values.
#'
#'  \describe{
#'    \item{observed}{A vector for observed P-values from the true case/control assignments.}
#'    \item{perm}{A vector for expected P-values (permutation-based expected p-value distribution).}
#' }
"example.Ps"

