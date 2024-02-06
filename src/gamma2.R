#' Efficient Calculation of Gamma Probabilities for Haplotype Analysis
#'
#' Computes the gamma probabilities with an efficiency improvement by utilizing
#' the optimized `forward2` and `backward2` functions for calculating the
#' forward and backward probabilities, respectively.
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
#'         of haplotype `k` at position `t`, given the entire sequence, calculated
#'         more efficiently than with the original gamma function.
#'
#' @details This version of the gamma calculation introduces efficiency improvements
#'          by utilizing the `forward2` and `backward2` functions, which are optimized
#'          versions of the forward and backward algorithms. These optimizations
#'          reduce the computational load without compromising the accuracy of the
#'          calculated probabilities.
#'
#' @examples
#' # Example usage:
#' haps <- matrix(c(1, 1, 0, 0, 1, 0, 1, 0), nrow = 2, byrow = TRUE)
#' hap <- c(1, 0)
#' gamma_matrix <- gamma2(haps, hap)
#' print(gamma_matrix)
#'
#' @export

gamma2 <- function(haps, hap, error = 0.1) {
  # Calculate the forward2 (alpha) and backward2 (beta) matrices
  alpha <- forward2(haps, hap, error)
  beta <- backward2(haps, hap, error)
  
  # Initialize the gamma matrix
  K <- nrow(haps)
  T <- ncol(haps)
  gamma_matrix <- matrix(0, nrow = K, ncol = T)
  
  # Calculate P(O = o) using the sum of the last column of alpha
  normalization_factor <- sum(alpha[, T])
  
  # Calculate gamma values with normalization
  for (t in 1:T) {
    for (k in 1:K) {
      gamma_matrix[k, t] <-
        alpha[k, t] * beta[k, t] / normalization_factor
    }
  }
  return(gamma_matrix)
}
