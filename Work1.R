
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
library(OpenImageR)
library(dplyr)
Files = list.files(path="Training/");Files
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
m2 = m-colMeans(m) #subtracting every colmean from every col
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


### FISHER
# YOU HAVE TO USE THE PCA DATA!!!!!!!!
# in summary, you get Sb, Sw and then multiply them to get w, then with w you can get the fisher value
# n == 150, p == 108000
# HERE'S THE TRICK, we compare each one of our classes to the OVERALL means and such
# or rather, to the ALTERNATIVE MEAN, so all the classes EXCEPT the one we're looking at
# we have classes c1-c25
num<-c(10:19,1,20:25,2:9)
ncol = 25
PCASSS = PCAs[,1:ncol]
overallmeans = colMeans(PCASSS)
cmeans = matrix(NA, nrow=25, ncol = ncol)
# get class means
pcopy = PCASSS
for (i in num){
  cmeans[i,] = colMeans(pcopy[1:6,])
  pcopy = pcopy[-c(1:6),]
}
cmeans[1,]

# Get Sb
# in plaintext, for each matrix, you subtract the overall mean of each column, like before
# from each class mean matrix, multiply the transposed vectors, times the number of observations
# (n) in each class, in our case, 6, then you get the sum of all those outputed and scaled matrices?

ni=6 #we can fortunately hardcode it since we know theres 6 observations per class
Sb = matrix(0, nrow=ncol, ncol = ncol)
for (i in num){
  base_matrix = cmeans[i,] - overallmeans
  unit_m = ni * base_matrix %*% t(base_matrix)
  Sb = Sb + unit_m
}

# I KNOW HOW TO GET SW
# OK SO ITS A 6 ROW BY 150 MATRIX,
# you get it by subtracting THE CLASS MEAN from each CLASS OBSERVATION
# then you multiply the transpose matrix
pcopy = PCASSS
Sw = matrix(0, nrow=ncol, ncol = ncol)
for (i in num){
  A_matrix = pcopy[c(1:6),] #our mini matrix for each class
  B_matrix = A_matrix - cmeans[i,] #remove the class mean from all observations
  D_matrix = matrix(0, nrow=ncol, ncol = ncol) #for storing each class' sums
  for (j in 1:6){
      C_matrix = B_matrix[j,] %*% t(B_matrix[j,]) #each observations' matrix
      D_matrix = D_matrix + C_matrix #add em up
    }
  Sw = Sw + D_matrix
  pcopy = pcopy[-c(1:6),]
}
Sw
# Get W, finally

w=eigen(solve(Sw)%*%Sb);w #number of eigenvals that arent 0 should be num of classes - 1

proj=as.matrix(PCASSS[,1:ncol])%*%w$vectors #I THINK this is the projection we want, of 25 vars
dim(proj)
proj
round(w$values / sum(w$values),3)
# make a class mean getter function



# Functionalize
Fisher <- function(PCA, n_comp, cmeans, labels){
  overallmeans = colMeans(PCA)
  
  # Get Sb
  Sb = matrix(0, nrow=n_comp, ncol = n_comp)

  for (label in unique(labels)){
    base_matrix = cmeans[label] - overallmeans
    # gross way to get n value, but should work
    unit_m = nrow(PCA[(labels %in% label),]) * base_matrix %*% t(base_matrix)
    Sb = Sb + unit_m
  }
  # Get Sw
  Sw = matrix(0, nrow=n_comp, ncol = n_comp)
  for (label in unique(labels)){
    class_matrix = PCA[(labels %in% label),] #our mini matrix for each class
    B_matrix = class_matrix - cmeans[label] #remove the class mean from all observations
    D_matrix = matrix(0, nrow=n_comp, ncol = n_comp) #for storing each class' sums

    for (j in 1:nrow(PCA[(labels %in% label),])){ # horrific code but whatever
      C_matrix = B_matrix[j,] %*% t(B_matrix[j,]) #each observations' matrix
      D_matrix = D_matrix + C_matrix #add em up
    }
    Sw = Sw + D_matrix
  }
  
  # Get W and form our projection
  w=eigen(solve(Sw)%*%Sb);w #number of eigenvals that arent 0 should be num of classes - 1
  variance_matrix =  round(w$values / sum(w$values),3)
  
  # Return PCA components and mean for projection
  list(mean = obs_mean, vectors = w$vectors, var_exp = variance_matrix,
       values=w$values)
}

project_fisher = function(PCA,n_comp,w){
  proj=as.matrix(PCA[,1:n_comp])%*%w$vectors #I THINK this is the projection we want, of 25 vars
  proj
  }

fish=Fisher(PCAs,25)
out=project
out
# FISHER KNN CLASSIFIER
source("function_loading.R")

Files = list.files(path="Training/")
labels <- as.numeric(gsub("[^0-9]", "", Files))
labels
# Define the k and threshold ranges
k_values <- 1:5  # Example range for k (number of neighbors)
threshold_values <- seq(1450, 1800, by = 50)  # Example range for threshold (distance)
total_combinations <- length(k_values) * length(threshold_values)
current_combination <- 0

# Perform LPPO with parameter tuning
tuning_results <- lppo_tuning_FISHER(m, labels, num_persons_out = 2, n_comp = 5, 
                              k_values = k_values, threshold_values = threshold_values, num_splits = 10)

# Display best k, threshold, and accuracy
print(paste("Best k:", tuning_results$best_k))
print(paste("Best threshold:", tuning_results$best_threshold))
print(paste("Best accuracy:", round(tuning_results$best_accuracy * 100, 2), "%"))

# Optionally, you can print the full results
print(tuning_results$results)

write.table(tuning_results$results, file = "results_15comp.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")








