test_that("Forward function behaves correctly when haplotypes always match",
          {
            K <- sample(5:20, 1) 
            T <- sample(5:20, 1) 
            error <- 0.1      
            
            # Create haplotypes and hap with all zeros
            haps <- matrix(0, nrow = K, ncol = T)
            hap <- rep(0, T)
            
            alpha_matrix <- forward(haps, hap, error)
            
            # Check the first column
            expect_equal(alpha_matrix[, 1], rep((1 - error) / K, K))
            
            # Check each successive column
            for (t in 2:T) {
              expect_equal(alpha_matrix[, t], 
                           alpha_matrix[, t - 1] * (1 - error))
            }
          })