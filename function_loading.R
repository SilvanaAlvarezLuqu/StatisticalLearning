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
mahalanobis_distance <- function(test_point, train_data, model, n_comp) {
  train_data=train_data[,1:n_comp]
  test_point=test_point[1:n_comp]
  
  inv_cov <- diag(1/model$values[1:n_comp]) # Compute inverse covariance matrix
  # Compute Mahalanobis distance for each training sample
  distances <- apply(train_data, 1, function(row) { #for each row
    diff <- row - test_point  # Element-wise subtraction
    sqrt(t(diff) %*% inv_cov %*% diff)  # Mahalanobis distance formula
  })
  
  return(distances)
}

#   EUCLIDEAN
#-----------------

euclidean_distance <- function(test_point, train_data, pca_model=NULL, n_comp=NULL) {
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
knn_classifier <- function(train_data, train_labels, test_data, k, 
                           percent_threshold, pca_model=NULL, n_comp=NULL,
                           distance_func=mahalanobis_distance) {
  
  # 1. Calculate within-class and between-class distances
  within_distances <- c()
  between_distances <- c()
  distances <- matrix(NA, nrow(train_data), nrow(train_data))
  
  # Calculate pairwise distances for training data
  for (i in 1:nrow(train_data)) {
    test_point <- train_data[i,]  
    distances[i,] <- distance_func(test_point, train_data, pca_model, n_comp)
  }
  
  # Separate within and between class distances
  n <- nrow(train_data)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (train_labels[i] == train_labels[j]) {
        within_distances <- c(within_distances, distances[i, j])
      } else {
        between_distances <- c(between_distances, distances[i, j])
      }
    }
  }
  
  max_within <- max(within_distances)
  q10_between <- quantile(between_distances, 0.1)
  
  # 2. Set threshold
  if (max_within < q10_between) {
    # There's a clear gap - use it
    threshold <- max_within + percent_threshold*(q10_between - max_within)
  } else {
    threshold <- quantile(within_distances, 1-percent_threshold)
  }
  # threshold <- quantile(within_distances, percent_threshold)
  # threshold <- 400
  cat("Threshold value:", threshold, "\n")
  
  # 3. Evaluate test samples - SIMPLIFIED APPROACH
  predictions <- character(nrow(test_data))  # Pre-allocate with correct type
  min_distances <- numeric(nrow(test_data))  # Store minimum distances
  
  for (i in 1:nrow(test_data)) {
    test_point <- test_data[i, ]
    
    # Get distances to all training points
    distances_test <- distance_func(test_point, train_data, pca_model, n_comp)
    min_dist <- min(distances_test)
    min_distances[i] <- min_dist
    
    # FIRST apply threshold check (completely separate from KNN)
    if (min_dist > threshold) {
      predictions[i] <- "0"  # Explicitly as character
      cat("Sample", i, "distance", min_dist, "> threshold", threshold, "→ label 0\n")
    } else {
      # Only do KNN for samples that pass threshold
      neighbors <- order(distances_test)[1:k]
      neighbor_labels <- train_labels[neighbors]
      
      # Handle potential issues with which.max
      label_counts <- table(neighbor_labels)
      most_common_idx <- which.max(label_counts)
      predicted_label <- names(label_counts)[most_common_idx]
      
      # Explicit conversion to character
      predictions[i] <- as.character(predicted_label)
      cat("Sample", i, "distance", min_dist, "<= threshold", threshold, "→ KNN label", predicted_label, "\n")
    }
  }
  
  # Return results
  return(list(
    within = within_distances, 
    between = between_distances,
    thres = threshold,
    pred = predictions, 
    min_distances = min_distances
  ))
}

#--------------------------------------------
#           PARAMETER  TUNNING
#--------------------------------------------

#   FUNCTION
#---------------
lppo_tuning <- function(data, labels, num_persons_out, n_comp_thresholds, 
                        k_values, percent_thresholds, num_splits,
                        distance_funcs = list(
                          "mahalanobis" = mahalanobis_distance,
                          "euclidean" = euclidean_distance,
                          "sse_mod" = sse_mod_distance,
                          "w_angle" = w_angle_distance
                        )) {
  unique_persons <- unique(labels)
  best_accuracy <- 0
  best_k <- NULL
  best_percent <- NULL  
  best_n_comp <- NULL
  best_distance <- NULL
  results_df <- data.frame()
  
  # Cast labels to character to ensure consistent comparison
  labels <- as.character(labels)
  
  total_combinations <- length(n_comp_thresholds) * length(k_values) * 
    length(percent_thresholds) * length(distance_funcs) * num_splits
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
    cumulative_variance <- cumsum(pca_model$var_exp) #/sum(pca_model$var_exp)
    
    n_comp_values <- sapply(n_comp_thresholds, function(th) {
      min(which(cumulative_variance >= th))
    })
    
    n_comp_values <- unique(n_comp_values)
    
    for (n_comp in n_comp_values) {
      pca_model_n <- list(components = pca_model$vectors[, 1:n_comp], 
                          mean = pca_model$mean, 
                          values = pca_model$values[1:n_comp])
      
      train_pca <- project_pca(train_data, pca_model_n)
      test_pca <- project_pca(test_data, pca_model_n)
      
      for (dist_name in names(distance_funcs)) {
        dist_func <- distance_funcs[[dist_name]]
        
        for (k in k_values) {
          for (percent_threshold in percent_thresholds) {
            # Add error handling
            tryCatch({
              result <- knn_classifier(
                train_data = train_pca,
                train_labels = train_labels,
                test_data = test_pca,
                k = k,
                percent_threshold = percent_threshold,
                pca_model = pca_model_n,
                n_comp = n_comp,
                distance_func = dist_func
              )
              predictions <- result$pred
              
              # Calculate accuracy
              accuracy <- mean(predictions == test_labels)
              
              # Check if threshold exists
              if (is.null(result$thres) || length(result$thres) == 0) {
                threshold_value <- NA
              } else {
                threshold_value <- result$thres
              }
              
              result_row <- data.frame(
                split = split, 
                n_comp = n_comp, 
                k = k, 
                percent = percent_threshold,
                distance = dist_name,
                accuracy = accuracy,
                var_explained = cumulative_variance[n_comp],
                threshold_value = threshold_value
              )
              
              results_df <- rbind(results_df, result_row)
            }, error = function(e) {
              cat(sprintf("\nError with dist=%s, k=%d, percent=%.2f: %s\n", 
                          dist_name, k, percent_threshold, e$message))
            })
            
            current_combination <- current_combination + 1
            cat(sprintf("\rProgress: %d/%d combinations (%.2f%%)", 
                        current_combination, total_combinations, 
                        (current_combination / total_combinations) * 100))
          }
        }
      }
    }
  }
  
  # Calculate average results across splits
  avg_results <- aggregate(
    accuracy ~ n_comp + k + percent + distance,
    data = results_df, 
    FUN = mean
  )
  
  # Find best parameters
  best_idx <- which.max(avg_results$accuracy)
  best_params <- avg_results[best_idx, ]
  
  # Get the average threshold value for the best configuration
  best_config_results <- results_df[
    results_df$n_comp == best_params$n_comp & 
      results_df$k == best_params$k & 
      results_df$percent == best_params$percent & 
      results_df$distance == best_params$distance,
  ]
  best_threshold <- mean(best_config_results$threshold_value, na.rm = TRUE)
  
  # Calculate variance explained
  best_var_exp <- mean(best_config_results$var_explained)
  
  # Return only the essential information
  return(list(
    k = best_params$k,
    percent_threshold = best_params$percent,
    n_components = best_params$n_comp,
    distance_function = best_params$distance,
    accuracy = best_params$accuracy,
    threshold_value = best_threshold,
    variance_explained = best_var_exp
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

