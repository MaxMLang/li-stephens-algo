#' Efficient Forward Algorithm for Haplotype Analysis
#'
#' Computes the forward probabilities for a given haplotype sequence with an
#' optimized approach, reducing computational complexity by using a modified
#' induction step that leverages a shared component of the probabilities.
#'
#' @param haps A matrix of observed haplotypes with individuals in rows and
#'        loci in columns.
#' @param hap A vector representing the haplotype for which the forward
#'        probabilities are to be computed.
#' @param error The error rate used in the emission probability calculations.
#'        Defaults to 0.1.
#'
#' @return A matrix of forward probabilities with the same dimensions as
#'         `haps`. Each element \code{alpha[k, t]} represents the probability
#'         of observing the sequence from the first position up to position `t`,
#'         given haplotype `k`, calculated more efficiently than the original
#'         forward algorithm.
#'
#' @details This version of the forward algorithm introduces an efficiency
#'          improvement by modifying the induction step to compute a shared
#'          component (`phi`) across all haplotypes at each locus. This reduces
#'          the overall number of computations required. The function assumes
#'          a constant recombination rate (`rho`), set internally to 0.001.
#'          The efficiency gain is particularly significant for large datasets.
#'
#' @examples
#' # Example usage:
#' haps <- matrix(c(1, 1, 0, 0, 1, 0, 1, 0), nrow = 2, byrow = TRUE)
#' hap <- c(1, 0)
#' alpha <- forward2(haps, hap)
#' print(alpha)
#'
#' @export
forward2 <- function(haps, hap, error = 0.1) {
  K <- nrow(haps)
  T <- ncol(haps)
  
  # Initialize the alpha matrix
  alpha <- matrix(0, nrow = K, ncol = T)
  
  # Transition probability (assuming constant recombination rate)
  rho <- 0.001  # Example recombination rate, adjust as necessary
  
  # Initialization
  for (k in 1:K) {
    bk1 <- ifelse(hap[1] == haps[k, 1], 1 - error, error)
    alpha[k, 1] <- (1 / K) * bk1  # Uniform initial probability
  }
  
  # Induction with the modified approach
  for (t in 2:T) {
    # Compute phi, which is constant for all k
    phi <- (1 - (1 - rho)) / K * sum(alpha[, t - 1])
    
    for (k in 1:K) {
      emission_prob <- ifelse(hap[t] == haps[k, t], 1 - error, error)
      alpha[k, t] <-
        (phi + ((1 - rho) * alpha[k, t - 1])) * emission_prob
      
    }
  }
  
  return(alpha)
}
