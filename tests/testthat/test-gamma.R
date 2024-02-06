test_that("Column sums of the gamma matrix should sum to 1", {
  # Set the parameters
  K <-
    sample(5:20, 1)  # Number of haplotypes, randomly chosen between 5 and 20
  T <-
    sample(5:20, 1)  # Number of loci, randomly chosen between 5 and 20
  error <- 0.1  # Error rate
  
  # Create haplotypes and hap with random 0 or 1 entries
  haps <-
    matrix(sample(0:1, K * T, replace = TRUE),
           nrow = K,
           ncol = T)
  hap <- sample(0:1, T, replace = TRUE)
  
  
  # Call the gamma function
  gamma_matrix <- gamma(haps, hap, error)
  
  # Check that the sum of each column is 1
  for (t in 1:T) {
    expect_equal(sum(gamma_matrix[, t]), 1)
  }
})
