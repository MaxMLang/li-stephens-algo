#' Backward Algorithm for Haplotype Analysis
#'
#' Computes the backward probabilities for a given haplotype sequence, using
#' a specified error rate and a fixed recombination rate.
#'
#' @param haps A matrix of observed haplotypes with individuals in rows and
#'        loci in columns.
#' @param hap A vector representing the haplotype for which the backward
#'        probabilities are to be computed.
#' @param error The error rate used in the emission probability calculations.
#'        Defaults to 0.1.
#'
#' @return A matrix of backward probabilities with the same dimensions as
#'         `haps`. Each element \code{beta[k, t]} represents the probability
#'         of observing the sequence from position `t` to the end, given
#'         haplotype `k`.
#'
#' @details The function calculates backward probabilities using a Markov
#'          model with a specified error rate for the emission probabilities
#'          and a fixed recombination rate (`rho`) for the transition
#'          probabilities. The backward algorithm is used in the context of
#'          hidden Markov models and is essential for computing the
#'          probabilities of haplotype sequences. The last column of the
#'          returned matrix is initialized with ones, and the probabilities
#'          are computed backwards from the last position to the first.
#'
#'          The recombination rate `rho` is set internally to 0.001. This
#'          represents the probability of a recombination event between
#'          adjacent loci.
#'
#' @examples
#' # Example usage:
#' haps <- matrix(c(1, 1, 0, 0, 1, 0, 1, 0), nrow = 2, byrow = TRUE)
#' hap <- c(1, 0)
#' beta <- backward(haps, hap)
#' print(beta)
#'
#' @export
#' 
backward <- function(haps, hap, error = 0.1) {
  K <- nrow(haps)
  T <- ncol(haps)
  
  # Initialize the beta matrix with ones in the last column
  beta <- matrix(0, nrow = K, ncol = T)
  beta[, T] <- 1
  
  # Transition probability 
  rho <- 0.001  # Example recombination rate
  
  # Induction step
  for (t in (T - 1):1) {
    for (k in 1:K) {
      # Calculate the sum for each k
      sum_beta <- 0
      for (i in 1:K) {
        # Transition probability Aki
        Aki <- ifelse(k == i, rho / K + 1 - rho, rho / (K))
        
        # Emission probability bit+1
        bit_plus_1 <-
          ifelse(hap[t + 1] == haps[i, t + 1], 1 - error, error)
        
        # Update sum_beta
        sum_beta <- sum_beta + Aki * bit_plus_1 * beta[i, t + 1]
      }
      # Update beta matrix
      beta[k, t] <- sum_beta
    }
  }
  
  return(beta)
}
