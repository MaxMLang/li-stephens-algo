K_values <- seq(5, 25, 1)  
T_values <- seq(5, 300, 25)  

# Function to create a random haplotype matrix
generate_haps <- function(K, T) {
  matrix(sample(0:1, K * T, replace = TRUE),
         nrow = K,
         ncol = T)
}

# Measure execution time for varying K (keeping T constant)
results_K <- data.frame(K = integer(), Time = numeric())
for (K in K_values) {
  haps <- generate_haps(K, 100)  # Keeping T constant at 100
  hap <- sample(0:1, 100, replace = TRUE)
  time_taken <- microbenchmark(forward(haps, hap), times = 10)$time
  results_K <- rbind(results_K, data.frame(K, Time = mean(time_taken)))
}

# Measure execution time for varying T (keeping K constant)
results_T <- data.frame(T = integer(), Time = numeric())
for (T in T_values) {
  haps <- generate_haps(10, T)  # Keeping K constant at 10
  hap <- sample(0:1, T, replace = TRUE)
  time_taken <- microbenchmark(forward(haps, hap), times = 10)$time
  results_T <- rbind(results_T, data.frame(T, Time = mean(time_taken)))
}


# Measure execution time for varying K (keeping T constant)
results_K_back <- data.frame(K = integer(), Time = numeric())
for (K in K_values) {
  haps <- generate_haps(K, 100)  # Keeping T constant at 100
  hap <- sample(0:1, 100, replace = TRUE)
  time_taken <- microbenchmark(backward(haps, hap), times = 10)$time
  results_K_back <- rbind(results_K_back, data.frame(K, Time = mean(time_taken)))
}

# # Measure execution time for varying T (keeping K constant)
results_T_back <- data.frame(T = integer(), Time = numeric())
for (T in T_values) {
  haps <- generate_haps(10, T)  # Keeping K constant at 10
  hap <- sample(0:1, T, replace = TRUE)
  time_taken <- microbenchmark(backward(haps, hap), times = 10)$time
  results_T_back <- rbind(results_T_back, data.frame(T, Time = mean(time_taken)))
}




ggplot(data = results_K, aes(x = K, y = Time)) +
  geom_line() +
  geom_point() +
  ggtitle("Execution Time vs. Number of Haplotypes (K)", 
          subtitle = "Forward Function") +
  xlab("Number of Haplotypes (K)") +
  ylab("Execution Time (microseconds)")

# Plot for varying T (keeping K constant)
ggplot(data = results_T, aes(x = T, y = Time)) +
  geom_line() +
  geom_point() +
  ggtitle("Execution Time vs. Number of Loci (T)", 
          subtitle = "Forward Function") +
  xlab("Number of Loci (T)") +
  ylab("Execution Time (microseconds)")

# Print tables for results_K and results_T
print(xtable(results_K, caption = "Tabular Results for Forward Function varying K"), type = "latex", comment = FALSE)
print(xtable(results_T, caption = "Tabular Results for Forward Function varying T"), type = "latex", comment = FALSE)


# Plot for varying K (keeping T constant)
ggplot(data = results_K_back, aes(x = K, y = Time)) +
  geom_line() +
  geom_point() +
  ggtitle("Execution Time vs. Number of Haplotypes (K)", 
          subtitle = "Backward Function") +
  xlab("Number of Haplotypes (K)") +
  ylab("Execution Time (microseconds)")

# Plot for varying T (keeping K constant)
ggplot(data = results_T_back, aes(x = T, y = Time)) +
  geom_line() +
  geom_point() +
  ggtitle("Execution Time vs. Number of Loci (T)", 
          subtitle = "Backward Function") +
  xlab("Number of Loci (T)") +
  ylab("Execution Time (microseconds)")


# Print tables for results_K_back and results_T_back
print(xtable(results_K_back, caption = "Tabular Results for Backward Function varying K"), type = "latex", comment = FALSE)
print(xtable(results_T_back, caption = "Tabular Results for Backward Function varying T"), type = "latex",  comment = F)


K_values <- seq(5, 50, 10)  # Varying K values
T <- 10  # Fixed T value for simplicity

# Function to create a random haplotype matrix
generate_haps <- function(K, T) {
  matrix(sample(0:1, K * T, replace = TRUE),
         nrow = K,
         ncol = T)
}

# Measure execution time
results <-
  data.frame(K = integer(),
             Time = numeric(),
             Function = character())

for (K in K_values) {
  haps <- generate_haps(K, T)
  hap <- sample(0:1, T, replace = TRUE)
  
  time_gamma <- microbenchmark(gamma(haps, hap), times = 10)$time
  results <-
    rbind(results, data.frame(K, Time = mean(time_gamma), Function = "gamma"))
  
  time_gamma2 <- microbenchmark(gamma2(haps, hap), times = 10)$time
  results <-
    rbind(results, data.frame(K, Time = mean(time_gamma2), Function = "gamma2"))
}


# Create a plot
ggplot(results, aes(x = K, y = Time, color = Function)) +
  geom_line() +
  geom_point() +
  ggtitle("Execution Time of Gamma and Gamma2") +
  xlab("Number of Haplotypes (K)") +
  ylab("Execution Time (microseconds)") +
  scale_color_viridis_d() +
  theme_minimal()

# Print the results table
print(results)

print(xtable(results, by = "Region_Name"), type = "latex", comment = FALSE)
