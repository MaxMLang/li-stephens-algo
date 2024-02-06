#' Efficient Backward Algorithm for Haplotype Analysis
#'
#' Computes the backward probabilities for a given haplotype sequence with
#' an optimized approach, reducing the computational complexity by
#' pre-calculating a shared component of the probabilities.
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
#' @details This version of the backward algorithm introduces an efficiency
#'          improvement by pre-calculating the sum of the product of beta and
#'          the emission probabilities for all haplotypes at each locus (\code{phi}).
#'          This sum is then used to update the backward probabilities for each
#'          haplotype, significantly reducing the number of computations required.
#'          The function assumes a constant recombination rate (`rho`), set
#'          internally to 0.001, across all loci. This rate represents the
#'          probability of a recombination event between adjacent loci. The
#'          modification enhances performance, especially for large datasets,
#'          without compromising the accuracy of the backward probabilities.
#'
#' @examples
#' # Example usage:
#' haps <- matrix(c(1, 1, 0, 0, 1, 0, 1, 0), nrow = 2, byrow = TRUE)
#' hap <- c(1, 0)
#' beta <- backward2(haps, hap)
#' print(beta)
#'
#' @export

backward2 <- function(haps, hap, error = 0.1) {
  K <- nrow(haps)
  T <- ncol(haps)
  
  # Initialize the beta matrix
  beta <- matrix(0, nrow = K, ncol = T)
  beta[, T] <- 1
  
  # Transition probability (assuming constant recombination rate)
  rho <- 0.001  # Example recombination rate
  
  # Induction step with the modified approach
  for (t in (T - 1):1) {
    # Calculate phi, which is a sum over all i for beta 
    # and emission probabilities
    phi <- 0
    for (i in 1:K) {
      bit_plus_1 <- ifelse(hap[t + 1] == haps[i, t + 1], 1 - error, error)
      phi <- phi + beta[i, t + 1] * bit_plus_1
    }
    phi <- phi * (rho / K)
    
    for (k in 1:K) {
      bit_plus_1 <- ifelse(hap[t + 1] == haps[k, t + 1], 1 - error, error)
      beta[k, t] <- phi + (1 - rho) * beta[k, t + 1] * bit_plus_1
    }
  }
  
  return(beta)
}
