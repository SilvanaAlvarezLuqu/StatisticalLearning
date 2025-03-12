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
# Define FDA function
#--------------------------------------------

# Produce model
Fisher <- function(PCA, labels){
  overallmeans = colMeans(PCA)
  labels = as.numeric(labels)
  n_features = ncol(PCA)  # Number of features in PCA
  n_classes = length(labels) #i think this lets us skip feeding it n_comp
  
  # Get class means
  cmeans = matrix(NA, nrow=max(labels), ncol = n_features) #max because we want to have empty mean vecs
  for (label in unique(labels)){
    cmeans[label,] = colMeans(PCA[(labels == label),])
  }
  
  # Get Sb
  Sb = matrix(0, nrow=n_features, ncol = n_features)
  
  for (label in unique(labels)){
    base_matrix = cmeans[label,] - overallmeans
    unit_m = sum(labels == label) * base_matrix %*% t(base_matrix)
    Sb = Sb + unit_m
  }
  
  # Get Sw
  Sw = matrix(0, nrow=n_features, ncol = n_features)
  for (label in unique(labels)){
    class_matrix = PCA[(labels == label),] #our mini matrix for each class
    # B_matrix = class_matrix - cmeans[label,] #remove the class mean from all observations
    B_matrix = sweep(class_matrix, 2, cmeans[label,], "-")  # Remove the class mean from all observations
    D_matrix = matrix(0, nrow=n_features, ncol = n_features) #for storing each class' sums
    
    for (j in 1:sum(labels == label)){ # horrific code but whatever
      C_matrix = B_matrix[j,] %*% t(B_matrix[j,]) #each observations' matrix
      D_matrix = D_matrix + C_matrix #add em up
    }
    Sw = Sw + D_matrix
  }
  
  # Get W and form our projection
  w=eigen(solve(Sw)%*%Sb); #number of eigenvals that arent 0 should be num of classes - 1
  variance_matrix =  round(w$values / sum(w$values),3)
  
  
  
  # Return PCA components and mean for projection
  # Get only the real parts of eigenvals and eigenvecs
  list(mean = overallmeans, vectors = Re(w$vectors), var_exp = variance_matrix,
       values=Re(w$values))
}

# Project our data using the model w we got
project_fisher = function(PCA,w){
  proj=as.matrix(PCA)%*%w$vectors #I THINK this is the projection we want
  proj
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

#   SSE MODIFIED
#------------------

sse_mod_distance <- function(test_point, train_data, pca_model=NULL, n_comp=NULL) {
  # Apply the same normalization if PCA model is provided
  if (!is.null(pca_model) && !is.null(n_comp)) {
    lambda <- pca_model$values[1:n_comp]
    norm_factor <- sqrt(1/lambda)
    
    train_data_norm <- sweep(train_data, 2, norm_factor, "*")
    test_point_norm <- test_point * norm_factor
    
    distances <- apply(train_data_norm, 1, function(row) {
      sum((row - test_point_norm)^2)/(sum((row)^2)*sum((test_point_norm)^2))
    })
  } else {
    # Original calculation
    distances <- apply(train_data, 1, function(row) {
      sum((row - test_point)^2)/(sum((row)^2)*sum((test_point)^2))
    })
  }
  return(distances)
}

#   WEIGHTED ANGLE-BASED DISTANCE
#----------------------------------

w_angle_distance <- function(test_point, train_data, pca_model, n_comp) {
  lambda <- pca_model$values[1:n_comp]
  z <- (1/lambda) 

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

#    TRAIN-TEST SPLIT
#-------------------------
create_train_test_split <- function(data, labels, num_persons_out, split_seed = NULL) {
  # Cast labels to character to ensure consistent comparison
  labels <- as.character(labels)
  
  # Get unique persons
  unique_persons <- unique(labels)
  
  # Set seed if provided
  if (!is.null(split_seed)) {
    set.seed(split_seed)
  }
  
  # Sample test persons
  test_persons <- sample(unique_persons, num_persons_out)
  remaining_persons <- setdiff(unique_persons, test_persons)
  
  # Get indices for test set part 1 (all data from test_persons)
  test_indices_p1 <- which(labels %in% test_persons)
  
  # Get indices for remaining persons
  remaining_indices <- which(labels %in% remaining_persons)
  
  # Set seed again if provided (to ensure reproducibility)
  if (!is.null(split_seed)) {
    set.seed(split_seed)
  }
  
  # Sample additional test indices from remaining persons
  test_indices_p2 <- sample(remaining_indices, round(length(remaining_indices) * 0.2))
  
  # Combine all test indices
  test_indices <- c(test_indices_p1, test_indices_p2)
  
  # Create test and train datasets
  test_data <- data[test_indices, ]
  test_labels <- labels[test_indices]
  
  train_data <- data[-test_indices, ]
  train_labels <- labels[-test_indices]
  
  # Return the result as a list
  return(list(
    train_data = train_data,
    train_labels = train_labels,
    test_data = test_data,
    test_labels = test_labels,
    test_indices = test_indices
  ))
}

#   FUNCTION
#---------------
lppo_tuning <- function(data, labels, num_persons_out, n_comp_thresholds, 
                        k_values, percent_thresholds, num_splits,
                        distance_funcs = list(
                          "mahalanobis" = mahalanobis_distance,
                          "sse_mod" = sse_mod_distance,
                          "w_angle" = w_angle_distance
                        )) {
  # Create objects
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
    
    # train/test split
    split_seed <- 782 + split
    split_data <- create_train_test_split(
      data = data,
      labels = labels,
      num_persons_out = num_persons_out,
      split_seed = split_seed
    )
    
    train_data <- split_data$train_data
    train_labels <- split_data$train_labels
    test_data <- split_data$test_data
    test_labels <- split_data$test_labels
    
    # PCA and projection
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
      
      # K_NN tuning
      
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
              
              # correct images labeled
              correct_predictions <- sum(predictions == test_labels)
              
              # Identify impostors
              impostors <- names(table(test_labels)[table(test_labels)==6]) # all the images in the test set
              correct_impostors <- sum(predictions== 0 & test_labels %in% impostors)
              
              accuracy <- (correct_predictions + correct_impostors) / length(predictions)
              
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
    variance_explained = best_var_exp,
    avg_results = avg_results
  ))
}

####################################################################################


lppo_tuning_FISHER <- function(data, labels, num_persons_out, n_comp_thresholds, 
                        k_values, percent_thresholds, num_splits,
                        distance_funcs = list(
                          "mahalanobis" = mahalanobis_distance,
                          "sse_mod" = sse_mod_distance,
                          "w_angle" = w_angle_distance
                        )) {
  # Create objects
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
    
    # train/test split
    split_seed <- 782 + split
    split_data <- create_train_test_split(
      data = data,
      labels = labels,
      num_persons_out = num_persons_out,
      split_seed = split_seed
    )
    
    train_data <- split_data$train_data
    train_labels <- split_data$train_labels
    test_data <- split_data$test_data
    test_labels <- split_data$test_labels
    
    # FDA and projection
    f_model <- Fisher(train_data, train_labels)
    cumulative_variance <- cumsum(f_model$var_exp) #/sum(pca_model$var_exp)
    
    n_comp_values <- sapply(n_comp_thresholds, function(th) {
      min(which(cumulative_variance >= th))
    })
    
    n_comp_values <- unique(n_comp_values)
    
    for (n_comp in n_comp_values) {
      f_model_n <- list(vectors = f_model$vectors[, 1:n_comp], 
                          mean = f_model$mean, 
                          values = f_model$values[1:n_comp])
      
      train_f <- project_fisher(train_data, f_model_n)
      test_f <- project_fisher(test_data, f_model_n)
      
      # K_NN tuning
      
      for (dist_name in names(distance_funcs)) {
        dist_func <- distance_funcs[[dist_name]]
        
        for (k in k_values) {
          for (percent_threshold in percent_thresholds) {
            # Add error handling
            tryCatch({
              result <- knn_classifier(
                train_data = train_f,
                train_labels = train_labels,
                test_data = test_f,
                k = k,
                percent_threshold = percent_threshold,
                pca_model = f_model_n,
                n_comp = n_comp,
                distance_func = dist_func
              )
              predictions <- result$pred
              
              #Calculate accuracy
              # correct images labeled
              correct_predictions <- sum(predictions == test_labels)
              
              # Identify impostors
              impostors <- names(table(test_labels)[table(test_labels)==6]) # all the images in the test set
              correct_impostors <- sum(predictions== 0 & test_labels %in% impostors)
              
              accuracy <- (correct_predictions + correct_impostors) / length(predictions)
              
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








