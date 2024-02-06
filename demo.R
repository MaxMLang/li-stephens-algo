# Read the reference and target haplotype files
refpanel <- read.table("data/refpanel.txt", header = TRUE, row.names = 1)
samples <- read.table("data/samples.txt", header = TRUE, row.names = 1)

# Assuming gamma2 function is already defined

# Get the population labels for the reference panel
refpop_labels <- rownames(refpanel)
YRI_indices <-
  which(grepl("YRI", refpop_labels))  # Indices for YRI haplotypes
CEU_indices <-
  which(grepl("CEU", refpop_labels))  # Indices for CEU haplotypes

# Process the first 5 target haplotypes
for (i in 1:5) {
  hap <- as.numeric(samples[i,])  # Convert to numeric values
  
  # Calculate the gamma matrix
  gamma_matrix <- gamma2(as.matrix(refpanel), hap)
  
  # Sum contributions from YRI and CEU haplotypes separately
  YRI_contribution <- rowSums(gamma_matrix[YRI_indices,])
  CEU_contribution <- rowSums(gamma_matrix[CEU_indices,])
  
  # Print results
  cat(paste0("Sample ", rownames(samples)[i], ":\n"))
  cat("YRI Contribution:", sum(YRI_contribution), "\n")
  cat("CEU Contribution:", sum(CEU_contribution), "\n")
  
  # Check if the sample has an entirely African genetic background
  if (sum(YRI_contribution) > 0.99 * length(hap)) {
    cat("This sample suggests an entirely African genetic background.\n\n")
  } else {
    cat("This sample does not suggest an entirely African genetic background.\n\n")
  }
}
