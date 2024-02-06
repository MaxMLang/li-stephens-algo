test_that("backward and backward2 give the same output", {
  K <- 10  # Number of haplotypes
  T <- 15  # Number of loci
  error <- 0.1  # Error rate
  
  # Create random haplotypes and hap
  haps <-
    matrix(sample(0:1, K * T, replace = TRUE),
           nrow = K,
           ncol = T)
  hap <- sample(0:1, T, replace = TRUE)
  
  # Run both functions
  result_backward <- backward(haps, hap, error)
  result_backward2 <- backward2(haps, hap, error)
  
  # Compare the results
  expect_equal(result_backward, result_backward2)
})