# Create a data frame to store results
results <- data.frame(k = integer(), threshold = numeric(), accuracy = numeric())

# Define your parameter grid for k and threshold
Files = list.files(path="Training/")
labels <- gsub("[^0-9]", "", Files)
labels<- as.numeric(labels)

k_values <- c(1, 3, 5, 7)         # Candidate k values


distances <- matrix(NA, 150, 150)
for (i in 1:150){
  query_image <- pcas[i, 1:30]  # Example: using the 10th image as a query
  distances[i,] <- mahalanobis_dist(pcas[,1:30], query_image, S_inv)
  
}

mean(distances); max(distances)
# [1] 1162.67
# [1] 1751.476
threshold_values <- seq(1200, 1800,by = 50)  # Candidate threshold values (for Mahalanobis distance)

cov_matrix <- cov(pcas[,1:30])  # Compute the covariance matrix for PCA-transformed data
cov_inv <- solve(cov_matrix)  # Inverse of the covariance matrix


# Iterate over all combinations of k and threshold
for (k in k_values) {
  for (threshold in threshold_values) {
    
    accuracies <- c()  # Store accuracies for each fold
    
    for (person_id in 1:25) {  # 25 persons, leave one out for testing
      # Split the data into training and test sets
      test_indices <- which(labels == person_id)  # Find images of the person to be tested
      train_indices <- setdiff(1:nrow(X), test_indices)  # Use the rest for training
      
      X_train_fold <- X[train_indices, ]
      labels_train_fold <- labels[train_indices]
      X_test_fold <- X[test_indices, ]
      labels_test_fold <- labels[test_indices]
      
      # Compute inverse covariance matrix for PCA (use your PCA transformation here)
      S_inv <- diag(1/pca_info[["EigVal"]][1:30])  # Replace with your own S_inv computation
      
      # Predict using k-NN classifier
      predictions <- sapply(1:nrow(X_test_fold), function(i) {
        as.numeric(knn_classifier(X_train_fold, labels_train_fold, X_test_fold[i, ], k, threshold, S_inv))
      })
      
      # Calculate accuracy for this fold
      accuracy <- mean(predictions == labels_test_fold)
      accuracies <- c(accuracies, accuracy)
    }
    
    # Store average accuracy for the current parameter combination
    results <- rbind(results, data.frame(k = k, threshold = threshold, accuracy = mean(accuracies)))
  }
}

# Get the best parameter combination
best_params <- results[which.max(results$accuracy), ]
print(best_params)

########################################################################################
# Define candidate parameters
k_values <- c(1, 3, 5, 7)         # Candidate k values
mean(distances); max(distances)
threshold_values <- seq(1200, 1800,by = 50)  # Candidate threshold values (for Mahalanobis distance)
n_repeats <- 10                   # Number of repeats for inner CV (to stabilize k tuning)


X<- pcas[,1:30]

# Data frame to store outer loop results
outer_results <- data.frame(threshold = numeric(),
                            avg_unknown_acc = numeric(),
                            best_k_overall = numeric())

# Outer loop: Leave-One-Person-Out (simulate unknown detection)
for (th in threshold_values) {
  outer_unknown_accs <- c()  # Store unknown accuracy for each outer fold
  best_k_list <- c()         # Store the best k from each outer fold
  
  for (person_id in unique(labels)) {
    # Outer fold: leave all images of one person out as "unknown"
    outer_test_idx <- which(labels == person_id)         # Test set: all images of this person
    outer_train_idx <- setdiff(1:nrow(X), outer_test_idx)   # Training set: images of all other persons
    
    X_train_outer <- X[outer_train_idx, ]
    labels_train_outer <- labels[outer_train_idx]
    X_test_outer <- X[outer_test_idx, ]
    
    # ---------------------------
    # Inner loop: Tune k on outer training set
    # ---------------------------
    inner_results <- data.frame(k = numeric(), acc = numeric())
    
    for (k in k_values) {
      inner_accs <- c()  # Accuracies for different splits inside the outer training set
      
      # For stability, repeat random splits (each person in the outer training set is split into training/validation)
      for (rep in 1:n_repeats) {
        inner_train_idx <- c()
        inner_val_idx <- c()
        
        # For each person in the outer training set, pick one image for validation and the rest for training
        for (pid in unique(labels_train_outer)) {
          person_indices <- which(labels_train_outer == pid)
          # Only perform split if there are at least 2 images for that person
          if (length(person_indices) > 1) {
            val_idx <- sample(person_indices, 1)
            train_idx <- setdiff(person_indices, val_idx)
          } else {
            # If only one image is available, use it for both training and validation (degenerate case)
            val_idx <- person_indices
            train_idx <- person_indices
          }
          inner_train_idx <- c(inner_train_idx, train_idx)
          inner_val_idx <- c(inner_val_idx, val_idx)
        }
        
        X_train_inner <- X_train_outer[inner_train_idx, ]
        labels_train_inner <- labels_train_outer[inner_train_idx]
        X_val_inner <- X_train_outer[inner_val_idx, ]
        labels_val_inner <- labels_train_outer[inner_val_idx]
        
        # Use a very high threshold to avoid rejecting known faces in inner CV
        high_threshold <- 1e6  
        S_inv <- diag(1/pca_info[["EigVal"]][1:30])
        
        # Predict on the inner validation set using our k-NN classifier
        inner_preds <- sapply(1:nrow(X_val_inner), function(i) {
          knn_classifier(X_train_inner, labels_train_inner, X_val_inner[i, ], k, high_threshold, S_inv)
        })
        inner_accs <- c(inner_accs, mean(inner_preds == labels_val_inner))
      }
      
      inner_results <- rbind(inner_results,
                             data.frame(k = k, acc = mean(inner_accs)))
    }
    
    # Best k from the inner loop for this outer fold
    best_k_inner <- inner_results$k[which.max(inner_results$acc)]
    best_k_list <- c(best_k_list, best_k_inner)
    
    # ---------------------------
    # Evaluate on the outer test set using the candidate threshold (th) and best k from inner loop
    # ---------------------------
    S_inv <- diag(1/pca_info[["EigVal"]][1:30])
    outer_preds <- sapply(1:nrow(X_test_outer), function(i) {
      knn_classifier(X_train_outer, labels_train_outer, X_test_outer[i, ], best_k_inner, th, S_inv)
    })
    # For unknown detection, we expect the classifier to return 0 for these outer test images
    unknown_acc <- mean(outer_preds == 0)
    outer_unknown_accs <- c(outer_unknown_accs, unknown_acc)
  }
  
  # Average unknown accuracy across all outer folds for this threshold candidate
  avg_unknown_acc <- mean(outer_unknown_accs)
  # For reporting, we take the median best k across outer folds
  best_k_overall <- median(best_k_list)
  
  outer_results <- rbind(outer_results,
                         data.frame(threshold = th,
                                    avg_unknown_acc = avg_unknown_acc,
                                    best_k_overall = best_k_overall))
}

print(outer_results)

# Select the best threshold based on the highest unknown accuracy
best_overall <- outer_results[which.max(outer_results$avg_unknown_acc), ]
cat("Best threshold:", best_overall$threshold, "\n")
cat("Best k (median across outer folds):", best_overall$best_k_overall, "\n")

# Optionally, you can further test the final classifier on external photos:
S_inv <- diag(1 / pca$EigVal)
final_external_preds <- sapply(1:nrow(external_photos), function(i) {
  knn_classifier(X, labels, external_photos[i, ], best_overall$best_k_overall, best_overall$threshold, S_inv)
})
external_acc <- mean(final_external_preds == 0)
cat("External unknown accuracy:", external_acc, "\n")

