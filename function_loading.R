#--------------------------------------------
# Define PCA function
#--------------------------------------------

PCA <- function(X_train, n_comp){
  M <- colMeans(X_train)    # Mean
  G <- X_train - M          # Centered Data
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
# Define Mahalanobis distance function
#--------------------------------------------

mahalanobis_distance <- function(test_point, train_data, pca_model, n_comp) {
  
  inv_cov <- diag(1/pca_model$values[1:n_comp]) # Compute inverse covariance matrix
  
  # Compute Mahalanobis distance for each training sample
  distances <- apply(train_data, 1, function(row) {
    diff <- row - test_point  # Element-wise subtraction
    sqrt(t(diff) %*% inv_cov %*% diff)  # Mahalanobis distance formula
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

## FISHER VER
lppo_tuning_FISHER <- function(data, labels, num_persons_out, n_comp, 
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
        
        # Ensure each of the remaining 23 people has 1 image in the test set
        # it randomly samples 
        selected_test_indices <- sample(1:nrow(train_data), 23)
        test_data <- train_data[selected_test_indices, ]
        test_labels <- train_labels[selected_test_indices]
        
        # this should all still line up, i think, since the mod is the same
        selected_train_indices = setdiff(1:length(train_labels),selected_test_indices) #this looks good i think
        train_data <- data[selected_train_indices, ]
        train_labels <- labels[selected_train_indices]
        
        # we want to do this AFTER removing the test sample
        # this should be the number of classes remaining, so 23
        train_class_means = vector(length=length(1:max(unique(train_labels))))
        for (label in unique(train_labels)){
          # this should give us the mean of each class that remains in the train_labels
            train_class_means[label] = colMeans(train_data[(train_labels %in% label),])
        }
        
        # Add all 6 images of the 2 fully left-out persons in the test set
        full_test_data <- data[labels %in% test_persons, ]
        full_test_labels <- labels[labels %in% test_persons]
        
        test_data <- rbind(test_data, full_test_data)
        test_labels <- c(test_labels, full_test_labels)
        
        # Compute PCA on the training set
        pca_model <- PCA(train_data, n_comp)
        # Project train and test data onto PCA space
        train_pca <- project_pca(train_data, pca_model)
        # get class means here, from train_pca
        test_pca <- project_pca(test_data, pca_model)
        # produce fisher samples from these
        
        fisher_model = Fisher(train_pca, n_comp, train_class_means, train_labels)
        
        train_fisher = project_fisher(train_pca, n_comp,fisher_model)
        test_fisher = project_fisher(test_pca, n_comp,fisher_model)
        
        print(fisher_model$values[1:n_comp])
        # Run k-NN classification with Mahalanobis distance
        predictions <- knn_classifier(train_fisher, train_labels, test_fisher, k,
                                      threshold, fisher_model, n_comp)
        print("got here")
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

