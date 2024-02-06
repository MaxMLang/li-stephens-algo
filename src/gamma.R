#' Calculation of Gamma Probabilities for Haplotype Analysis
#'
#' Computes the gamma probabilities, which represent the posterior probabilities
#' of being in a particular state at a particular time, given the entire haplotype
#' sequence. This function uses the original `forward` and `backward` algorithms.
#'
#' @param haps A matrix of observed haplotypes with individuals in rows and
#'        loci in columns.
#' @param hap A vector representing the haplotype for which the gamma
#'        probabilities are to be computed.
#' @param error The error rate used in the emission probability calculations.
#'        Defaults to 0.1.
#'
#' @return A matrix of gamma probabilities with the same dimensions as `haps`.
#'         Each element \code{gamma_matrix[k, t]} represents the posterior probability
#'         of haplotype `k` at position `t`, given the entire sequence.
#'
#' @details The gamma probabilities are calculated by multiplying the forward
#'          probabilities (alpha) by the backward probabilities (beta) and
#'          normalizing by the probability of the observed sequence. This
#'          approach provides insights into the state distribution across the
#'          sequence given the observations.
#'
#' @examples
#' # Example usage:
#' haps <- matrix(c(1, 1, 0, 0, 1, 0, 1, 0), nrow = 2, byrow = TRUE)
#' hap <- c(1, 0)
#' gamma_matrix <- gamma(haps, hap)
#' print(gamma_matrix)
#'
#' @export
#' 
gamma <- function(haps, hap, error = 0.1) {
  # Calculate the forward (alpha) and backward (beta) matrices
  alpha <- forward(haps, hap, error)
  beta <- backward(haps, hap, error)
  
  # Initialize the gamma matrix
  K <- nrow(haps)
  T <- ncol(haps)
  gamma_matrix <- matrix(0, nrow = K, ncol = T)
  
  # Calculate P(O = o) using the sum of the last column of alpha
  P_O_o <- sum(alpha[, T])
  
  # Calculate gamma values with normalization
  for (t in 1:T) {
    for (k in 1:K) {
      gamma_matrix[k, t] <- alpha[k, t] * beta[k, t] / P_O_o
    }
  }
  
  return(gamma_matrix)
}
