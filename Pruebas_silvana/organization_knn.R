#--------------------------------------------
# Define Mahalanobis distance function
#--------------------------------------------

mahalanobis_dist <- function(data, query, inv_cov_matrix) {
  diff <- t(data) - query  # Subtract query image from all others
  dist <- apply(diff, 2, function(x) sqrt(t(x) %*% inv_cov_matrix %*% x))  # Calculate distance
  return(dist)
}

#--------------------------------------------
# Define function for k-NN classification
#--------------------------------------------
knn_classifier <- function(query_image, train_data, train_labels, k, threshold, S_inv) {
  
  # Calculate distances
  distances <- mahalanobis_dist(train_data, query_image, S_inv)
  
  # Get the k nearest neighbors
  nearest_neighbors <- order(distances)[1:k]
  
  # Get the labels of the k nearest neighbors
  neighbor_labels <- train_labels[nearest_neighbors]
  
  # Predict the class (majority vote)
  predicted_class <- as.numeric(names(sort(table(neighbor_labels), decreasing = TRUE))[1])
  
  # Apply the threshold based on the minimum distance
  min_distance <- min(distances[nearest_neighbors])  # Get the minimum distance to the neighbors
  if (min_distance <= threshold) {
    return(predicted_class)  # Classify as belonging to the database
  } else {
    return(0)  # Classify as "unknown" (0)
  }
}

#--------------------------------------------
#           PARAMETER  TUNNING
#--------------------------------------------
# Table of results
results <- data.frame(k = integer(), threshold = numeric(), accuracy = numeric())

# Grid searching
k_values <- c(1, 3, 5, 7)         # Candidate k values
threshold_values <- seq(1500, 1800, by = 50)  # Candidate threshold values (for Mahalanobis distance)

# Labels
Files = list.files(path="Training/")
labels <- gsub("[^0-9]", "", Files)
labels<- as.numeric(labels)

person_ids <- unique(labels)  # Unique person identifiers
num_people <- length(person_ids)


S_inv <- diag(1 / pca_info[["EigVal"]][1:30])
