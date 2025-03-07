#--------------------------------------------
# Define PCA function
#--------------------------------------------

PCA <- function(X_train, n_comp){
  M <- colMeans(X_train)    # Mean
  G <- X_train - M        # Centered Data
  G_ <- t(G)                # Transpose centered data
  n_train <- nrow(X_train)        # Number of train observations 
    
  Sigma_s <- (G%*%G_)/(n_train-1) # Compute Sigma small
  eig=eigen(Sigma_s)        # Obtain Sigma_small eigen vectors/values
  e_vec_s=eig$vectors       # select eigen vector of sigma small
  
  P = G_%*%e_vec_s          # Eigen vectors of Sigma_large 
  L = diag(eig$values)      # Diagonal matrix of eigenvalues
  D =  round(eig$values / sum(eig$values),3) # Prop of variance explained
  
  # Select top principal components
  top_comp <- P[, 1:n_comp]
  
  # Return PCA components and mean for projection
  list(components = top_comp, mean = M, vectors = P,
       values = eig$values, var_exp = D)
}

# Function to project data onto PCA space
project_pca <- function(data, pca_model) {
  centered_data <- data - pca_model$mean # Center the data
  projected_data <- as.matrix(centered_data) %*% pca_model$components
  return(projected_data)
}

#--------------------------------------------
# Define Distances Functions
#--------------------------------------------

#   MAHALANOBIS
#-----------------
mahalanobis_distance <- function(test_point, train_data, pca_model, n_comp) {
  
  inv_cov <- diag(1/pca_model$values[1:n_comp]) # Compute inverse covariance matrix
  
  # Compute Mahalanobis distance for each training sample
  distances <- apply(train_data, 1, function(row) {
    diff <- row - test_point  # Element-wise subtraction
    sqrt(t(diff) %*% inv_cov %*% diff)  # Mahalanobis distance formula
  })
  
  return(distances)
}

#   EUCLIDEAN
#-----------------

euclidean_distance <- function(test_point, train_data) {
  distances <- apply(train_data, 1, function(row) {
    sqrt(sum((row - test_point)^2))  # Euclidean distance formula
  })
  return(distances)
}


#   SSE MODIFIED
#------------------

sse_mod_distance <- function(test_point, train_data) {
  distances <- apply(train_data, 1, function(row) {
    sum((row - test_point)^2)/(sum((row)^2)*sum((test_point)^2))  # Euclidean distance formula
  })
  return(distances)
}

#   WEIGHTED ANGLE-BASED DISTANCE
#----------------------------------

w_angle_distance <- function(test_point, train_data, pca_model, n_comp) {
  lambda <- pca_model$values[1:n_comp]
  z <- (1/lambda) # Compute inverse covariance matrix
  
  # Compute Mahalanobis distance for each training sample
  distances <- apply(train_data, 1, function(row) {
    numerator <- sum(z*row*test_point)
    denominator <- sqrt(sum(row^2)*sum(test_point^2))
    - numerator/denominator
  })
  
  return(distances)
}


#--------------------------------------------
# Define function for k-NN classification
#--------------------------------------------
knn_classifier <- function(train_data, train_labels, test_data, k, threshold, pca_model, n_comp) {
  
  predictions <- c()
  
  for (i in 1:nrow(test_data)) {
    
    test_point <- test_data[i, ]
    
    # Compute Mahalanobis distances from the test point to all training points
    distances <- mahalanobis_distance(test_point, train_data, pca_model, n_comp)
    
    # Get indices of k nearest neighbors
    neighbors <- order(distances)[1:k]
    
    # Extract their labels
    neighbor_labels <- train_labels[neighbors]
    
    # Majority voting: most common label among k neighbors
    predicted_label <- names(which.max(table(neighbor_labels)))
    
    # If average distance of k neighbors is too high, classify as 0 (unknown)
    if (mean(distances[neighbors]) > threshold) {
      predicted_label <- 0
    }
    
    # Store prediction
    predictions <- c(predictions, predicted_label)
  }
  
  return(predictions)
}

#--------------------------------------------
#           PARAMETER  TUNNING
#--------------------------------------------

#   FUNCTION
#---------------
lppo_tuning <- function(data, labels, num_persons_out, n_comp, 
                        k_values, threshold_values, num_splits) {
  unique_persons <- unique(labels)
  best_accuracy <- 0
  best_k <- NULL
  best_threshold <- NULL
  results <- list()
  
  for (k in k_values) {
    for (threshold in threshold_values) {
      split_accuracies <- c()  # Store accuracies for each split
      
      for (split in 1:num_splits) {
        # Randomly select 2 people to put all their images in the test set
        test_persons <- sample(unique_persons, num_persons_out)  # Randomly pick 2 people
        
        # For the remaining 23 people, 5 images in train, 1 image in test
        remaining_persons <- setdiff(unique_persons, test_persons)
        
        train_data <- data[(labels %in% remaining_persons), ]
        train_labels <- labels[(labels %in% remaining_persons)]
        
        # Test set: 1 image per person from the 23 left
        test_data <- data[labels %in% remaining_persons, ]
        test_labels <- labels[labels %in% remaining_persons]
        
        # Ensure each of the remaining 23 people has 1 image in the test set
        selected_test_indices <- sample(1:nrow(test_data), 23)
        test_data <- test_data[selected_test_indices, ]
        test_labels <- test_labels[selected_test_indices]
        
        # Add all 6 images of the 2 fully left-out persons in the test set
        full_test_data <- data[labels %in% test_persons, ]
        full_test_labels <- labels[labels %in% test_persons]
        
        test_data <- rbind(test_data, full_test_data)
        test_labels <- c(test_labels, full_test_labels)
        
        # Compute PCA on the training set
        pca_model <- PCA(train_data, n_comp)
        
        # Project train and test data onto PCA space
        train_pca <- project_pca(train_data, pca_model)
        test_pca <- project_pca(test_data, pca_model)
        
        # Run k-NN classification with Mahalanobis distance
        predictions <- knn_classifier(train_pca, train_labels, test_pca, k,
                                      threshold, pca_model, n_comp)
        
        # Calculate accuracy for this split
        accuracy <- mean(predictions == test_labels)
        split_accuracies <- c(split_accuracies, accuracy)
      }
      
      # Calculate average accuracy for this combination of k and threshold
      avg_accuracy <- mean(split_accuracies)
      
      # Update best k and threshold if this combination is better
      if (avg_accuracy > best_accuracy) {
        best_accuracy <- avg_accuracy
        best_k <- k
        best_threshold <- threshold
      }
      
      # Store results for each combination of k and threshold
      results[[paste("k =", k, "threshold =", threshold)]] <- avg_accuracy
      
      # Update progress
      current_combination <- current_combination + 1
      cat(sprintf("\rProgress: %d/%d combinations (%.2f%%)", 
                  current_combination, total_combinations, 
                  (current_combination / total_combinations) * 100))
    }
  }
  
  # Move to the next line after the progress updates
  cat("\n")
  
  # Return the best parameters and accuracy
  list(best_k = best_k, best_threshold = best_threshold,
       best_accuracy = best_accuracy, results = results)
}

#   TUNNING
#---------------
Files = list.files(path="../Training/")
labels <- as.numeric(gsub("[^0-9]", "", Files))

# Define the k and threshold ranges
k_values <- 1:5  # Example range for k (number of neighbors)

#  quantile(m_dist, seq(0.6,1,0.05))
#     60%      65%      70%      75%      80%      85%      90%      95%     100% 
#1250.423 1281.545 1308.429 1339.554 1374.580 1413.649 1463.031 1524.398 1751.476 
threshold_values <- seq(1450, 1800, by = 50)  # Example range for threshold (distance)
# threshold_values <- c(1200)
# threshold based on the quantile distribution
# ncomp: 15 bc acumulate the 90% 

total_combinations <- length(k_values) * length(threshold_values)
current_combination <- 0

# Perform LPPO with parameter tuning
tuning_results <- lppo_tuning(X, labels, num_persons_out = 2, n_comp = 15, 
                              k_values = k_values, threshold_values = threshold_values, num_splits = 10)

# Display best k, threshold, and accuracy
print(paste("Best k:", tuning_results$best_k))
print(paste("Best threshold:", tuning_results$best_threshold))
print(paste("Best accuracy:", round(tuning_results$best_accuracy * 100, 2), "%"))

# Optionally, you can print the full results
print(tuning_results$results)

write.table(tuning_results$results, file = "results_15comp.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

