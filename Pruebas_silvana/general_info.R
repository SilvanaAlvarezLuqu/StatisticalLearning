pca <- PCA(X,30)
proj <- project_pca(X, pca)

m_dist <- matrix(NA, nrow(X), nrow(X))
for (i in 1:nrow(X)) {
  test_point <- proj[i,]  # Example: using the 10th image as a query
  m_dist[i,] <-mahalanobis_distance(test_point, proj, pca, n_comp=30)
}

summary(m_dist)
hist(m_dist)
quantile(m_dist, seq(0.6,1,0.05))
