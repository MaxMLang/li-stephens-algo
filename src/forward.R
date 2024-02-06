#' Forward Algorithm for Haplotype Analysis
#'
#' Computes the forward probabilities for a given haplotype sequence using
#' a specified error rate and a constant recombination rate.
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
#'         given haplotype `k`.
#'
#' @details The function calculates forward probabilities using a Markov
#'          model with a specified error rate for the emission probabilities
#'          and a constant recombination rate (`rho`), set internally to 0.001,
#'          for the transition probabilities. The initialization step assigns
#'          uniform initial probabilities to each haplotype, and the induction
#'          step updates these probabilities across loci.
#'
#' @examples
#' # Example usage:
#' haps <- matrix(c(1, 1, 0, 0, 1, 0, 1, 0), nrow = 2, byrow = TRUE)
#' hap <- c(1, 0)
#' alpha <- forward(haps, hap)
#' print(alpha)
#'
#' @export
forward <- function(haps, hap, error = 0.1) {
  K <- nrow(haps)
  T <- ncol(haps)
  
  # Initialize the alpha matrix
  alpha <- matrix(0, nrow = K, ncol = T)
  
  # Transition probability 
  rho <- 0.001  # Example recombination rate given in the task
  
  # Initialization
  for (k in 1:K) {
    bk1 <- ifelse(hap[1] == haps[k, 1], 1 - error, error)
    alpha[k, 1] <- (1 / K) * bk1  # Uniform initial probability
  }
  
  # Induction
  for (t in 2:T) {
    for (k in 1:K) {
      sum_alpha <- 0
      for (i in 1:K) {
        Aik <- ifelse(i == k, (rho) / K + 1 - rho, rho / K)
        emission_prob <-
          ifelse(hap[t] == haps[k, t], 1 - error, error)
        sum_alpha <-
          sum_alpha + alpha[i, t - 1] * Aik * emission_prob
      }
      alpha[k, t] <- sum_alpha
    }
  }
  
  return(alpha)
}
