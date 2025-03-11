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
# Define Distances Functions
#--------------------------------------------

#   MAHALANOBIS
#-----------------
mahalanobis_distance <- function(test_point, train_data, pca_model, n_comp) {
  
  inv_cov <- diag(1/pca_model$values[1:n_comp]) # Compute inverse covariance matrix
  
  # Compute Mahalanobis distance for each training sample
  distances <- apply(train_data, 1, function(row) {
    diff <- row - test_point  # Element-wise subtraction
    sqrt(abs(t(diff) %*% inv_cov %*% diff))  # Mahalanobis distance formula
    # throwing in an absolute value to skip complex number errors
  })
  
  return(distances)
}

#   EUCLIDEAN
#-----------------

euclidean_distance <- function(test_point, train_data) {
  # Ensure test_point is a vector, not a matrix or data frame
  if (!is.null(dim(test_point))) {
    test_point <- as.vector(test_point)
  }
  
  # Calculate distances using vectorized operations when possible
  distances <- apply(train_data, 1, function(row) {
    sqrt(abs(sum((row - test_point)^2)))  # Euclidean distance formula
    # same skip here
  })
  
  return(distances)
}



#   SSE MODIFIED
#------------------

sse_mod_distance <- function(test_point, train_data, 
                             pca_model=NULL, n_comp=NULL) {
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
    denominator <- sqrt(abs(sum(row^2)*sum(test_point^2))) #same abs() complex number skip
    - numerator/denominator
  })
  
  return(distances)
}


#--------------------------------------------
# Define function for k-NN classification
#--------------------------------------------
knn_classifier <- function(train_data, train_labels, test_data, k, threshold,
                           pca_model=NULL, n_comp=NULL,
                           distance_func=mahalanobis_distance) {
  
  predictions <- c()
  
  for (i in 1:nrow(test_data)) {
    
    test_point <- test_data[i, ]
    
    # Compute the selected distances from the test point to all training points
    distances <- distance_func(test_point, train_data, pca_model, n_comp)
    
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
lppo_tuning <- function(data, labels, num_persons_out, n_comp_thresholds, 
                        k_values, dist_threshold, num_splits,
                        distance_funcs = list(
                          "mahalanobis" = mahalanobis_distance,
                          "euclidean" = euclidean_distance,
                          "sse_mod" = sse_mod_distance,
                          "w_angle" = w_angle_distance
                        )) {
  unique_persons <- unique(labels)
  best_accuracy <- 0
  best_k <- NULL
  best_threshold <- NULL
  best_n_comp <- NULL
  best_distance <- NULL
  results_df <- data.frame()
  
  # Cast labels to character to ensure consistent comparison
  labels <- as.character(labels)
  
  total_combinations <- length(n_comp_thresholds) * length(k_values) * 
    length(dist_threshold) * length(distance_funcs) * num_splits
  current_combination <- 0
  
  for (split in 1:num_splits) {
    cat(sprintf("Processing split %d/%d\n", split, num_splits))
    
    set.seed(782 + split)
    test_persons <- sample(unique_persons, num_persons_out)
    remaining_persons <- setdiff(unique_persons, test_persons)
    
    test_indices_p1 <- which(labels %in% test_persons)
    remaining_indices <- which(labels %in% remaining_persons)
    
    set.seed(782 + split)  
    test_indices_p2 <- sample(remaining_indices, round(length(remaining_indices) * 0.2))
    
    test_indices <- c(test_indices_p1, test_indices_p2)
    
    test_data <- data[test_indices, ]
    test_labels <- labels[test_indices]
    
    train_data <- data[-test_indices, ]
    train_labels <- labels[-test_indices]
    
    pca_model <- PCA(train_data, min(nrow(train_data)-1, ncol(train_data)))
    cumulative_variance <- cumsum(pca_model$var_exp)/sum(pca_model$var_exp)
    
    n_comp_values <- sapply(n_comp_thresholds, function(th) {
      min(which(cumulative_variance >= th))
    })
    
    n_comp_values <- unique(n_comp_values)
    
    for (n_comp in n_comp_values) {
      pca_model_n <- list(components = pca_model$vectors[, 1:n_comp], 
                          mean = pca_model$mean, 
                          values = pca_model$values)
      
      train_pca <- project_pca(train_data, pca_model_n)
      test_pca <- project_pca(test_data, pca_model_n)
      
      # Calculate appropriate threshold values for each distance function if needed
      # For example, euclidean distances might need different thresholds than Mahalanobis
      
      for (dist_name in names(distance_funcs)) {
        dist_func <- distance_funcs[[dist_name]]
        
        # Adjust thresholds for this distance function if needed
        local_thresholds <- dist_threshold
        if (dist_name == "euclidean") {
          # For euclidean, it might need much larger thresholds
          # Calculate an appropriate scale based on the data
          sample_distances <- c()
          for (i in 1:min(100, nrow(train_pca))) {
            sample_distances <- c(sample_distances, 
                                  mean(euclidean_distance(train_pca[i,], train_pca[-i,])))
          }
          median_dist <- median(sample_distances)
          local_thresholds <- median_dist * dist_threshold  # Scale thresholds appropriately
        }
        
        for (k in k_values) {
          for (threshold_idx in 1:length(local_thresholds)) {
            threshold <- local_thresholds[threshold_idx]
            original_threshold <- dist_threshold[threshold_idx]  # For reporting
            
            # Custom KNN implementation for this distance function
            predictions <- c()
            
            for (i in 1:nrow(test_pca)) {
              test_point <- test_pca[i, ]
              
              # Call distance function with proper parameters
              if (dist_name %in% c("mahalanobis", "w_angle")) {
                distances <- dist_func(test_point, train_pca, pca_model_n, n_comp)
              } else {
                distances <- dist_func(test_point, train_pca)
              }
              
              # Make sure k is not larger than available data
              k_actual <- min(k, length(distances))
              
              # Get indices of k nearest neighbors
              k_indices <- order(distances)[1:k_actual]
              
              # Extract their labels
              k_labels <- train_labels[k_indices]
              
              # Calculate frequencies of labels
              label_counts <- table(k_labels)
              
              if (length(label_counts) > 0) {
                # Find the most frequent label
                predicted_label <- names(which.max(label_counts))
                
                # Check if the average distance exceeds the threshold
                if (mean(distances[k_indices]) > threshold) {
                  predicted_label <- "0"  # Unknown category
                }
              } else {
                predicted_label <- "0"  # Default if no neighbors found
              }
              
              predictions <- c(predictions, predicted_label)
            }
            
            # Calculate accuracy
            accuracy <- mean(predictions == test_labels)
            
            # Include debugging info
            if (accuracy < 0.01) {
              cat(sprintf("\nVery low accuracy (%f): dist=%s, k=%d, threshold=%.2f\n", 
                          accuracy, dist_name, k, threshold))
              
              # Print a sample of predictions
              sample_indices <- sample(1:length(predictions), 
                                       min(5, length(predictions)))
              
              cat("Sample predictions vs true labels:\n")
              for (idx in sample_indices) {
                cat(sprintf("Pred: %s, True: %s\n", 
                            predictions[idx], test_labels[idx]))
              }
              
              # Check if all predictions are "0" (unknown)
              if (all(predictions == "0")) {
                cat("All predictions are 0 (unknown) - threshold may be too low\n")
              }
            }
            
            result_row <- data.frame(
              split = split, 
              n_comp = n_comp, 
              k = k, 
              threshold = original_threshold,  # Use original threshold for consistency
              scaled_threshold = threshold,    # Record actual threshold used
              distance = dist_name,
              accuracy = accuracy,
              var_explained = cumulative_variance[n_comp]
            )
            
            results_df <- rbind(results_df, result_row)
            
            current_combination <- current_combination + 1
            cat(sprintf("\rProgress: %d/%d combinations (%.2f%%)", 
                        current_combination, total_combinations, 
                        (current_combination / total_combinations) * 100))
          }
        }
      }
    }
  }
  
  avg_results <- aggregate(
    accuracy ~ n_comp + k + threshold + distance, 
    data = results_df, 
    FUN = mean
  )
  
  best_idx <- which.max(avg_results$accuracy)
  best_params <- avg_results[best_idx, ]
  
  best_var_exp <- mean(results_df$var_explained[results_df$n_comp == best_params$n_comp])
  
  # Additional analysis by distance function
  dist_performance <- aggregate(
    accuracy ~ distance, 
    data = avg_results, 
    FUN = function(x) c(max = max(x), mean = mean(x))
  )
  
  # Get details about the best configuration for each distance function
  best_by_distance <- lapply(names(distance_funcs), function(d) {
    dist_results <- avg_results[avg_results$distance == d, ]
    if (nrow(dist_results) == 0) return(NULL)
    
    best_idx <- which.max(dist_results$accuracy)
    if (length(best_idx) == 0) return(NULL)
    
    best_row <- dist_results[best_idx, ]
    return(list(
      distance = d,
      accuracy = best_row$accuracy,
      k = best_row$k,
      n_comp = best_row$n_comp,
      threshold = best_row$threshold
    ))
  })
  names(best_by_distance) <- names(distance_funcs)
  
  return(list(
    best_k = best_params$k,
    best_threshold = best_params$threshold,
    best_n_comp = best_params$n_comp,
    best_distance = best_params$distance,
    best_accuracy = best_params$accuracy,
    var_explained = best_var_exp,
    detailed_results = results_df,
    avg_results = avg_results,
    distance_performance = dist_performance,
    best_by_distance = best_by_distance
  ))
}

####################################################################################

new_lppo_tuning_fisher <- function(data, labels, num_persons_out, n_comp_thresholds,f_thresh, 
                        k_values, dist_threshold, num_splits,
                        distance_funcs = list(
                          "mahalanobis" = mahalanobis_distance,
                          "euclidean" = euclidean_distance,
                          "sse_mod" = sse_mod_distance,
                          "w_angle" = w_angle_distance
                        )) {
  unique_persons <- unique(labels)
  best_accuracy <- 0
  best_k <- NULL
  best_threshold <- NULL
  best_n_comp <- NULL
  best_distance <- NULL
  results_df <- data.frame()
  
  # Cast labels to character to ensure consistent comparison
  labels <- as.character(labels)
  
  total_combinations <- length(n_comp_thresholds) * length(k_values) * 
    length(dist_threshold) * length(distance_funcs) * num_splits
  current_combination <- 0
  
  for (split in 1:num_splits) {
    cat(sprintf("Processing split %d/%d\n", split, num_splits))
    
    set.seed(782 + split)
    test_persons <- sample(unique_persons, num_persons_out)
    remaining_persons <- setdiff(unique_persons, test_persons)
    
    test_indices_p1 <- which(labels %in% test_persons)
    remaining_indices <- which(labels %in% remaining_persons)
    
    set.seed(782 + split)  
    test_indices_p2 <- sample(remaining_indices, round(length(remaining_indices) * 0.2))
    
    test_indices <- c(test_indices_p1, test_indices_p2)
    
    test_data <- data[test_indices, ]
    test_labels <- labels[test_indices]
    
    train_data <- data[-test_indices, ]
    train_labels <- labels[-test_indices]
    
    
    # we want to do this AFTER removing the test sample
    # this should be the number of classes remaining, so 23
    train_class_means = vector(length=length(1:max(unique(train_labels))))
    for (label in unique(train_labels)){
      # this should give us the mean of each class that remains in the train_labels
      train_class_means[label] = colMeans(train_data[(train_labels %in% label),])
    }
    
    
    pca_model <- PCA(train_data, min(nrow(train_data)-1, ncol(train_data)))
    cumulative_variance <- cumsum(pca_model$var_exp)/sum(pca_model$var_exp)
    
    n_comp_values <- sapply(n_comp_thresholds, function(th) {
      min(which(cumulative_variance >= th))
    })
    
    n_comp_values <- unique(n_comp_values)
    
    for (n_comp in n_comp_values) {
      pca_model_n <- list(components = pca_model$vectors[, 1:n_comp], 
                          mean = pca_model$mean, 
                          values = pca_model$values)
      
      train_pca <- project_pca(train_data, pca_model_n)
      test_pca <- project_pca(test_data, pca_model_n)
      #swap in Fisher distances, see what crashes
      # n_comp = round(n_comp/2) #fuck it, let's see what happens
      for (f_comp in f_thresh){
        fisher_model = Fisher(train_pca[, 1:f_comp], f_comp, train_class_means, train_labels)
        
        
        train_fisher = project_fisher(train_pca[, 1:f_comp], f_comp,fisher_model)
        test_fisher = project_fisher(test_pca[, 1:f_comp], f_comp,fisher_model)
        
        
        # Calculate appropriate threshold values for each distance function if needed
        # For example, euclidean distances might need different thresholds than Mahalanobis
        for (dist_name in names(distance_funcs)) {
          dist_func <- distance_funcs[[dist_name]]
          # Adjust thresholds for this distance function if needed
          local_thresholds <- dist_threshold
          if (dist_name == "euclidean") {
            # For euclidean, it might need much larger thresholds
            # Calculate an appropriate scale based on the data
            sample_distances <- c()
            for (i in 1:min(100, nrow(train_fisher))) {
              sample_distances <- c(sample_distances, 
                                    mean(euclidean_distance(train_fisher[i,], train_fisher[-i,])))
            }
            median_dist <- median(sample_distances)
            local_thresholds <- median_dist * dist_threshold  # Scale thresholds appropriately
          }
          
          for (k in k_values) {
            for (threshold_idx in 1:length(local_thresholds)) {
              threshold <- local_thresholds[threshold_idx]
              original_threshold <- dist_threshold[threshold_idx]  # For reporting
              
              # Custom KNN implementation for this distance function
              predictions <- c()
              
              for (i in 1:nrow(test_fisher)) {
                test_point <- test_fisher[i, ]
                
                # Call distance function with proper parameters
                if (dist_name %in% c("mahalanobis", "w_angle")) {
                  distances <- dist_func(test_point, train_fisher, fisher_model, f_comp)
                } else {
                  distances <- dist_func(test_point, train_fisher)
                }
                
                # Make sure k is not larger than available data
                k_actual <- min(k, length(distances))
                
                # Get indices of k nearest neighbors
                k_indices <- order(distances)[1:k_actual]
                
                # Extract their labels
                k_labels <- train_labels[k_indices]
                
                # Calculate frequencies of labels
                label_counts <- table(k_labels)
                
                if (length(label_counts) > 0) {
                  # Find the most frequent label
                  
                  predicted_label <- names(which.max(label_counts))
                  # print(mean(distances[k_indices]))
                  # print(threshold)
                  
                  # Check if the average distance exceeds the threshold
                  if (mean(distances[k_indices]) > threshold) {
                    predicted_label <- "0"  # Unknown category
                  }
                } else {
                  predicted_label <- "0"  # Default if no neighbors found
                }
                
                predictions <- c(predictions, predicted_label)
              }
              
              # Calculate accuracy
              accuracy <- mean(predictions == test_labels)
              
              # Include debugging info
              if (accuracy < 0.01) {
                cat(sprintf("\nVery low accuracy (%f): dist=%s, k=%d, threshold=%.2f\n", 
                            accuracy, dist_name, k, threshold))
                
                # Print a sample of predictions
                sample_indices <- sample(1:length(predictions), 
                                         min(5, length(predictions)))
                
                cat("Sample predictions vs true labels:\n")
                for (idx in sample_indices) {
                  cat(sprintf("Pred: %s, True: %s\n", 
                              predictions[idx], test_labels[idx]))
                }
                
                # Check if all predictions are "0" (unknown)
                if (all(predictions == "0")) {
                  cat("All predictions are 0 (unknown) - threshold may be too low\n")
                }
              }
              
              result_row <- data.frame(
                split = split, 
                n_comp = n_comp, 
                k = k, 
                threshold = original_threshold,  # Use original threshold for consistency
                scaled_threshold = threshold,    # Record actual threshold used
                distance = dist_name,
                accuracy = accuracy,
                var_explained = cumulative_variance[f_comp]
              )
              
              results_df <- rbind(results_df, result_row)
              
              current_combination <- current_combination + 1
              cat(sprintf("\rProgress: %d/%d combinations (%.2f%%)", 
                          current_combination, total_combinations, 
                          (current_combination / total_combinations) * 100))
            }
          }
        }
      }
    }
  }
  
  avg_results <- aggregate(
    accuracy ~ n_comp + k + threshold + distance, 
    data = results_df, 
    FUN = mean
  )
  
  best_idx <- which.max(avg_results$accuracy)
  best_params <- avg_results[best_idx, ]
  
  best_var_exp <- mean(results_df$var_explained[results_df$n_comp == best_params$n_comp])
  
  # Additional analysis by distance function
  dist_performance <- aggregate(
    accuracy ~ distance, 
    data = avg_results, 
    FUN = function(x) c(max = max(x), mean = mean(x))
  )
  
  # Get details about the best configuration for each distance function
  best_by_distance <- lapply(names(distance_funcs), function(d) {
    dist_results <- avg_results[avg_results$distance == d, ]
    if (nrow(dist_results) == 0) return(NULL)
    
    best_idx <- which.max(dist_results$accuracy)
    if (length(best_idx) == 0) return(NULL)
    
    best_row <- dist_results[best_idx, ]
    return(list(
      distance = d,
      accuracy = best_row$accuracy,
      k = best_row$k,
      n_comp = best_row$n_comp,
      threshold = best_row$threshold
    ))
  })
  names(best_by_distance) <- names(distance_funcs)
  
  return(list(
    best_k = best_params$k,
    best_threshold = best_params$threshold,
    best_n_comp = best_params$n_comp,
    best_distance = best_params$distance,
    best_accuracy = best_params$accuracy,
    var_explained = best_var_exp,
    detailed_results = results_df,
    avg_results = avg_results,
    distance_performance = dist_performance,
    best_by_distance = best_by_distance
  ))
}

####################################################################################












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
        accuracy <- round(mean(predictions == test_labels),3)
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

