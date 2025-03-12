# Load necessary libraries
library(stats)

Files = list.files(path="Training/")
labels <- as.numeric(gsub("[^0-9]", "", Files))

Files = list.files(path="Training/")
# Files
X = matrix(NA, nrow=150, ncol = 108000)

for(i in seq_along(Files)){
  Im = readImage(paste0("Training/",Files[i]))
  ri=as.vector(Im[,,1])
  gi=as.vector(Im[,,2])
  bi=as.vector(Im[,,3])
  X[i,] = cbind(ri,gi,bi)
}
dim(X)

# Define parameter grid for tuning
n_comp_thresholds <- c(0.80, 0.85, 0.90, 0.95)  # Cumulative variance thresholds for PCA
k_values <- c(1, 3, 5)  # Number of neighbors for KNN
percent_thresholds <- c(0.1, 0.2, 0.3)  # Thresholds for novel person detection
num_splits <- 5  # Number of cross-validation splits
num_persons_out <- 2  # Number of persons to leave out for testing

# Run the tuning function
tuning_results <- lppo_tuning(
  data = X, 
  labels = labels,
  num_persons_out = num_persons_out,
  n_comp_thresholds = n_comp_thresholds,
  k_values = k_values,
  percent_thresholds = percent_thresholds,
  num_splits = num_splits
)


# Access results simply
cat("Best k:", tuning_results$k, "\n")
cat("Best percent threshold:", tuning_results$percent_threshold, "\n")
cat("Best distance function:", tuning_results$distance_function, "\n")
cat("Best accuracy:", tuning_results$accuracy, "\n")

