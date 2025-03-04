# Define Mahalanobis distance function
mahalanobis_dist <- function(data, query, inv_cov_matrix) {
  diff <- t(data) - query  # Subtract query image from all others
  dist <- apply(diff, 2, function(x) sqrt(t(x) %*% inv_cov_matrix %*% x))  # Calculate distance
  return(dist)
}

knn_classifier <- function(query_image, train_data, train_labels, k, threshold, S_inv) {
  # Calculate distances using Mahalanobis metric
  distances <- mahalanobis_dist(train_data, query_image, S_inv)
  
  # Get k nearest neighbors
  nearest_neighbors <- order(distances)[1:k]
  neighbor_labels <- train_labels[nearest_neighbors]
  
  # Minimum distance among the k neighbors
  min_distance <- min(distances[nearest_neighbors])
  
  # If minimum distance is greater than threshold, classify as unknown (0)
  if (min_distance > threshold) {
    return(0)  # Unknown person
  } else {
    # Predict the most common label among the k neighbors
    return(as.numeric(names(sort(table(neighbor_labels), decreasing = TRUE))[1]))
  }
}


# Loop over k and threshold values
results <- data.frame(k = integer(), threshold = numeric(), accuracy = numeric())

# Define parameter grid for k and threshold
k_values <- c(1, 3, 5, 7)         # Candidate k values
threshold_values <- seq(1500, 1800, by = 50)  # Candidate threshold values (for Mahalanobis distance)

# Cross-validation setup: LOI-CV (leave-one-image-out)
person_ids <- unique(labels)  # Unique person identifiers
num_people <- length(person_ids)

# Precompute covariance matrix and its inverse
cov_matrix <- cov(pcas[,1:30])
S_inv <- solve(cov_matrix)

# Loop over k and threshold values for tuning
for (k in k_values) {
  for (threshold in threshold_values) {
    
    accuracy <- 0  # Initialize accuracy
    
    # Perform Leave-One-Image-Out CV
    for (person_id in person_ids) {
      # Find the indices of all images for the current person
      person_indices <- which(labels == person_id)
      
      # Perform Leave-One-Image-Out for this person (use 5 images for training, 1 for testing)
      for (i in 1:length(person_indices)) {
        # Choose the i-th image of this person as the test image and leave one person "out of the DB"
        test_idx <- person_indices[i]
        train_idx <- setdiff(person_indices, c(test_idx,25))  # Use all other images of this person for training
        
        # Training data and labels (include all images of other persons)
        X_train <- rbind(pcas[train_idx, 1:30], pcas[which(labels != person_id), 1:30])  # Include images of other people
        y_train <- c(labels[train_idx], labels[which(labels != person_id)])
        
        # Test data (1 image for testing)
        query_image <- pcas[test_idx, 1:30]
        true_label <- labels[test_idx]
        
        # Classify the test image using k-NN
        prediction <- knn_classifier(query_image, X_train, y_train, k, threshold, S_inv)
        
        # Calculate accuracy for this fold
        accuracy <- accuracy + (prediction == true_label)
      }
    }
    
    # Average accuracy for this combination of k and threshold
    accuracy <- accuracy / (num_people * length(person_indices))  # Normalize by the total number of test images
    
    # Store the results
    results <- rbind(results, data.frame(k = k, threshold = threshold, accuracy = accuracy))
  }
}

# Display the results
print(results)

# Find the best parameters
best_params <- results[which.max(results$accuracy), ]
print(best_params)
