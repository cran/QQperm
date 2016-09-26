#' Read sample file and genotype collapsing matrix in IGM (Insitute for Genomic Medicine, Columbia University) format.
#'
#' @author Slave Petrovski and Quanli Wang
#' @param sample.file The input sample file.
#' @param matrix.file The input matrix file.
#' @param filter.list An optional list of genes to exclude from the analysis.
#'
#' @return Returns a list containing the data matrix (matrix) and case(1)/control(0) phenotype indicators (is.case).
#'
#' @examples
#' #igm.data <- igm.read.data(sample.file, matrix.file)
#'
igm.read.data <- function(sample.file, matrix.file, filter.list = NULL) {
  #get sample and phenotype label information
  ped <- read.table(sample.file,header = FALSE, sep = '\t', as.is = TRUE)
  samples <- ped[,2]
  is.case <- as.matrix(ped[,6] ==2)
  rownames(is.case) <- samples

  #get count matrix
  data <- read.table(matrix.file,header = TRUE, sep = '\t', as.is = TRUE)
  genes <- data[,1]
  data <- data[,samples]
  data[data > 0] <- 1

  #convert data to matrix
  matrix <- as.matrix(data[])
  rownames(matrix)<-genes

  if (!is.null(filter.list)) {
    genes.common <- intersect(genes,filter.list)
    matrix <- matrix[genes.common,]
  }
  out <- list()
  out$data <- matrix
  out$is.case <- is.case
  out
}

#' The permutation-based empirical NULL distribution of P-values is generated through label switching and permutation of the true case/control assignment. To achieve this, for a given matrix it randomly permutes the case and control labels of the original configuration and then recomputes the two-tail Fisher's Exact test for all genes. This is repeated n.permutation (e.g., 1000) times. For each of the n.permutations the p-values are ordered and then the mean of each rank-ordered estimate is taken across the n.permutations, i.e., the average 1st order statistic, the average 2nd order statistic, etc. These then represent the empirical estimates of the expected ordered p-values (expected -log10(p-values)). This empirical-based expected p-value distribution no longer depends on an assumption that the p-values are uniformly distributed under the null.
#' @author Slave Petrovski and Quanli Wang
#' @param matrix The input genotype matrix, with rows for genes and columns for samples.
#' @param is.case The case(1)/control(0) indicator.
#' @param n.permutations The number of label permutations requested. Recommended number is 1,000.
#'
#' @return Returns a list containing the observed P-values based on the correct case/control assignment (observed) and the permuted P-Values reflecting the appropriate NULL distribution for the given case/control configuration (perm).
#'
#' @examples
#' #Ps  <- igm.get.pvalues(matrix, is.case)
#'
igm.get.pvalues <-function(matrix, is.case, n.permutations = 1000) {
  #Number of cases and controls
  n.samples <- length(is.case)
  n.cases <- sum(is.case)
  n.controls <- n.samples - n.cases

  #pre-compute all possible contingency.table for performance
  #this will create a look up table
  contingency.table = matrix(0,2,2)
  Fisher.precompute <- matrix(0, n.cases + 1, max(rowSums(matrix)) + 1)
  for (i in 1: dim(Fisher.precompute)[1]) {
    for (j in 1:dim(Fisher.precompute)[2]) {
      contingency.table[1,1] <- i-1
      contingency.table[1,2] <- n.cases - contingency.table[1,1]
      contingency.table[2,1] <- j -1
      contingency.table[2,2] <- n.controls - contingency.table[2,1]
      Fisher.precompute[i,j] <- fisher.test(contingency.table)$p.value
    }
  }

  #permutation, save all p-values just in case median will be needed later on
  P.Values <- matrix(1,dim(matrix)[1],n.permutations)
  total.1 <- rowSums(matrix)
  for (i in 1: n.permutations) {
    K <- sample.int(n.samples, size = n.cases, replace = FALSE)
    Labels.1.1 <- rowSums(matrix[,K])
    Labels.0.1 <- total.1 - Labels.1.1
    P.Values[,i] <- sort(Fisher.precompute[cbind(Labels.1.1+1,Labels.0.1+1)])
  }
  P.perm <- rowMeans(P.Values)


  #compute observed (true case control configration) p-values
  K <- which(is.case)
  Labels.1.1 <- rowSums(matrix[,K])
  Labels.0.1 <- total.1 - Labels.1.1
  P.observed <- sort(Fisher.precompute[cbind(Labels.1.1+1,Labels.0.1+1)])

  out <- list()
  out$perm <- P.perm
  out$observed <- P.observed
  out
}
