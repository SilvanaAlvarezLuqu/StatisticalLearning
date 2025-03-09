#   TUNNING
#---------------
Files = list.files(path="../Training/")
labels <- as.numeric(gsub("[^0-9]", "", Files))

# Range of general parameters:
k_values <- 1:5  # number of neighbors
n_comp_thresholds <- seq(0.8,0.95,0.05)

# Specific parameters
threshold_euclidean <- seq(42000, 56000, by = 2000) # 8 values

# Tunning
euc<-knn_classifier(train_data, train_labels, test_data, k=5, threshold=50000,
                           pca_model=NULL, n_comp=NULL, 
                    distance_func=euclidean_distance)
  
tunning_euclidean <- lppo_tuning(data = X, labels = labels, num_persons_out = 2,
                                 n_comp_thresholds, k_values, dist_threshold = threshold_euclidean,
                                 num_splits = 5, distance_func = euclidean_distance)


# Perform LPPO with parameter tuning
tuning_mahalanobis <- lppo_tuning(data = X, labels = labels, num_persons_out = 2,
                                  n_comp_thresholds, k_values, dist_threshold = threshold_euclidean,
                                  num_splits = 5, distance_func = mahalanobis_distance)



# Display best k, threshold, and accuracy
print(paste("Best k:", tuning_results$best_k))
print(paste("Best threshold:", tuning_results$best_threshold))
print(paste("Best accuracy:", round(tuning_results$best_accuracy * 100, 2), "%"))

# Optionally, you can print the full results
print(tuning_results$results)

write.table(tuning_results$results, file = "results_15comp.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

