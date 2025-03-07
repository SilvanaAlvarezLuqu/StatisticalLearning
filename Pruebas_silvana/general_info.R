Files = list.files(path="../Training/")
# Files
X = matrix(NA, nrow=150, ncol = 108000)

for(i in seq_along(Files)){
  Im = readImage(paste0("../Training/",Files[i]))
  ri=as.vector(Im[,,1])
  gi=as.vector(Im[,,2])
  bi=as.vector(Im[,,3])
  X[i,] = cbind(ri,gi,bi)
}
dim(X)


pca <- PCA(X, n_comp=150)
proj <- project_pca(X, pca)

#   Mahalanobis distance
#--------------------------
m_dist <- matrix(NA, nrow(X), nrow(X))
for (i in 1:nrow(X)) {
  test_point <- proj[i,]  # Example: using the 10th image as a query
  m_dist[i,] <-mahalanobis_distance(test_point, proj, pca, n_comp=150)
}

summary(m_dist)
hist(m_dist)
quantile(m_dist, seq(0.6,1,0.05))

threshold_mahalanobis <- seq(1450, 1800, by = 50)

#   Euclidean distance
#--------------------------
euclidean_dist <- matrix(NA, nrow(X), nrow(X))
for (i in 1:nrow(X)) {
  test_point <- proj[i,]  # Example: using the 10th image as a query
  euclidean_dist[i,] <-euclidean_distance(test_point, proj)
}

summary(euclidean_dist)
hist(euclidean_dist)
quantile(euclidean_dist, seq(0.6,1,0.05))

threshold_euclidean <- seq(42000, 56000, by = 1000)

#   SSE Modified distance
#--------------------------
sse_mod_dist <- matrix(NA, nrow(X), nrow(X))
for (i in 1:nrow(X)) {
  test_point <- proj[i,]  # Example: using the 10th image as a query
  sse_mod_dist[i,] <-sse_mod_distance(test_point, proj)
}

summary(sse_mod_dist)
hist(sse_mod_dist)
quantile(sse_mod_dist, seq(0.6,1,0.05))

threshold_sse_mod <- seq(3.4e-10, 2.2e-09, by = 2e-010)

#  Weighted Angle-based distance
#----------------------------------
w_angle_dist <- matrix(NA, nrow(X), nrow(X))
for (i in 1:nrow(X)) {
  test_point <- proj[i,]  # Example: using the 10th image as a query
  w_angle_dist[i,] <-w_angle_distance(test_point, proj, pca, n_comp=150)
}

summary(w_angle_dist)
hist(w_angle_dist)
quantile(w_angle_dist, seq(0.6,1,0.05))

threshold_w_angle <- seq(-1.5e-04, 2e-04, by = .5e-04)