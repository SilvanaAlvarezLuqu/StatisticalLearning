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
  # threshold <- quantile(within_distances, percentile_threshold)
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

# PRUEBA
#---------------
n_comp<-10
pca_model <- PCA(train_data, min(nrow(train_data)-1, ncol(train_data)))
pca_model_n <- list(components = pca_model$vectors[, 1:n_comp], 
                    mean = pca_model$mean, 
                    values = pca_model$values[1:n_comp])

train_pca <- project_pca(train_data, pca_model_n)
test_pca <- project_pca(test_data, pca_model_n)
length(train_labels)

predictions <- knn_classifier(
  train_data = train_pca,
  train_labels = train_labels,
  test_data = test_pca,
  k = 3,
  percent_threshold = 0.2,
  pca_model = pca_model_n,
  n_comp = n_comp
)

hist(predictions$within) 
hist(predictions$between)

dist<- c(predictions$within,predictions$between)
hist(dist)
abline(v=predictions$thres)

predictions$thres
quantile(dist, seq(0,1,0.1))  
  quantile(predictions$within, seq(0,1,0.1))  
  quantile(predictions$between, seq(0,1,0.1))  
  quantile(predictions$min_distances, seq(0,1,0.1))  
predictions$pred
hist(predictions$distances_test)

(n<- length(predictions$pred))
(sum(predictions$pred == test_labels))/n # + sum(predictions$pred==0 & test_labels %in% c("14", "6"))

#####################################################################################################


#   EUCLIDEAN
#-----------------

euclidean_distance <- function(test_point, train_data, pca_model=NULL, n_comp=NULL) {
  # If PCA model is provided, use it to normalize dimensions
  if (!is.null(pca_model) && !is.null(n_comp)) {
    
    lambda <- pca_model$values[1:n_comp]
    norm_factor <- sqrt(1/lambda)
    
    # Apply normalization to both train_data and test_point
    train_data_norm <- sweep(train_data, 2, norm_factor, "*")
    test_point_norm <- test_point * norm_factor
    
    # Calculate normalized distances
    distances <- apply(train_data_norm, 1, function(row) {
      sqrt(sum((row - test_point_norm)^2))
    })
  } else {
    # Original calculation if no normalization needed
    distances <- apply(train_data, 1, function(row) {
      sqrt(sum((row - test_point)^2))
    })
  }
  
  return(distances)
}




