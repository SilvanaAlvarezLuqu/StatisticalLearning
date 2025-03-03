# Apply Mahalanobis distance
#---------------------------------
S_inv <- diag(1/pca_info[["EigVal"]][1:30])

mahalanobis_dist <- function(X_train, x_query, S_inv) {
  diffs <- sweep(X_train, 2, x_query)  # Subtract query image from all rows
  dists <- sqrt(rowSums((diffs %*% S_inv) * diffs))  # Compute distance
  return(dists)
}

query_image <- pcas[10, 1:30]  # Example: using the 10th image as a query
distance <- mahalanobis_dist(pcas[,1:30], query_image, S_inv)

distances <- matrix(NA, 150, 150)
for (i in 1:150){
  query_image <- pcas[i, 1:30]  # Example: using the 10th image as a query
  distances[i,] <- mahalanobis_dist(pcas[,1:30], query_image, S_inv)
  
}

hist(distances)
quantile(distances, probs = seq(0, 1, 0.1))


# Define that the photo does not belong to the data set if it exceeds
# the 95% quantile of the distances (1525 in mahalanobis distance)


# KNN CLASSIFIER
#-------------------
knn_classifier <- function(X_train, labels, x_query, k, threshold, S_inv) {
  distances <- mahalanobis_dist(X_train, x_query, S_inv)
  neighbors <- order(distances)[1:k]  # Get k nearest neighbors
  nearest_labels <- labels[neighbors]
  
  # Majority vote for classification
  best_match <- names(which.max(table(nearest_labels)))
  best_distance <- min(distances[neighbors])
  
  # Apply threshold: if distance is too large, return 0 (unknown)
  if (best_distance > threshold) {
    return(0)
  } else {
    return(best_match)
  }
}

