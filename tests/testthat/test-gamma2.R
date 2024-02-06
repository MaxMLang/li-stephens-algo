test_that("gamma and gamma2 give the same output", {
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
  result_gamma <- gamma(haps, hap, error)
  result_gamma2 <- gamma2(haps, hap, error)
  
  # Compare the results
  expect_equal(result_gamma, result_gamma2)
})