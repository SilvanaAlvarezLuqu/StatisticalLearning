
library(OpenImageR)
Im=OpenImageR::readImage("Training/1AT.jpg")
dim(Im)
red=Im[,,1]
green=Im[,,2]
blue=Im[,,3]
imageShow(red)
blue
# concatenate them in a long vector
as.vector(red)
youhavetodo=c(as.vector(red),as.vector(green), as.vector(blue))

# size of var covar matrix= 3x3, for principal components
m=cbind(as.vector(red),as.vector(green),as.vector(blue));m
dim(m)
out=princomp(m)
pc1=matrix(out$scores[,1],nrow = nrow(red),ncol = ncol(green))
imageShow(pc1)
pc2=matrix(out$scores[,2],nrow = nrow(red),ncol = ncol(green))
imageShow(pc2)
pc3=matrix(out$scores[,3],nrow = nrow(red),ncol = ncol(green))
imageShow(pc3)

out$sdev/sum(out$sdev)


###ACTUAL THING
library(dplyr)
Files = list.files(path="Training/");Files
#m = vector("list", length(Files))
m = matrix(NA, nrow=150, ncol = 108000)

for(i in seq_along(Files)){
  Im = readImage(paste0("Training/",Files[i]))
  ri=as.vector(Im[,,1])
  gi=as.vector(Im[,,2])
  bi=as.vector(Im[,,3])
  m[i,] = cbind(ri,gi,bi)
}
dim(m) ## m == X in formula
m[1,1]
obs_mean = colMeans(m)
m2 = m-colMeans(m) #subtracting every mean from every row
n=nrow(m2)
m_ =t(m2)
Sigma_s= (m2%*%m_)/(n-1)
dim(Sigma_s)
eig=eigen(Sigma_s)
e_vec_s=eig$vectors
e_vec_l = m_%*%e_vec_s #== P in formula
dim(e_vec_l)
e_vals=diag(eig$values)
round(eig$values/sum(eig$values),3)

PCAs=e_vec_l[1:10,]%*%m2
PCAs=m2%*%e_vec_l

dim(PCAs)
# n = 36000
num<-c(10:19,1,20:25,2:9)
label=c()
for (i in num){
  temp = rep(i,6)
  label=c(label,temp)
}
M=cbind(as.factor(label),m2)





# KNN func
# Custom kNN function
knn <- function(train, test, labels, k) {
  # Ensure inputs are matrices
  train <- as.matrix(train)
  test <- as.matrix(test)
  
  # Number of test samples
  n_test <- nrow(test)
  
  # Initialize predictions
  predictions <- vector("list", n_test)
  
  # Loop through each test sample
  for (i in 1:n_test) {
    # Compute Euclidean distances between the test point and all training points
    distances <- sqrt(rowSums((train - matrix(test[i, ], nrow = nrow(train), ncol = ncol(train), byrow = TRUE))^2))
    
    # Find the indices of the k nearest neighbors
    nearest_indices <- order(distances)[1:k]
    
    # Get the labels of the k nearest neighbors
    nearest_labels <- labels[nearest_indices]
    
    # For classification: Use majority voting
    if (is.factor(labels)) {
      predictions[[i]] <- names(sort(table(nearest_labels), decreasing = TRUE))[1]
    }
    # For regression: Use the average
    else if (is.numeric(labels)) {
      predictions[[i]] <- mean(nearest_labels)
    }
  }
  
  # Return predictions as a vector
  return(unlist(predictions))
}

# Example usage
set.seed(123)

# Generate some example data
train_data <- matrix(rnorm(100), ncol = 2)  # 50 training samples, 2 features
train_labels <- factor(sample(c("A", "B"), 50, replace = TRUE))  # Binary classification labels
test_data <- matrix(rnorm(10), ncol = 2)  # 5 test samples, 2 features

# Set k (number of neighbors)
k <- 3

# Make predictions
predictions <- knn(train_data, test_data, train_labels, k)

# Print predictions
print(predictions)








