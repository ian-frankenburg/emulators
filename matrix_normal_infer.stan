#include dmatrix_normal.stan
data{
  int n;
  int p;
  matrix[n, p] Y;
  matrix[n, p] X;
}
parameters{
  vector[p] beta;
  corr_matrix[n] U; 
  corr_matrix[p] V; 
  vector<lower=0>[n] sigma_u; 
  vector<lower=0>[p] sigma_v; 
}
transformed parameters{
  cov_matrix[n] Sigma_u = quad_form_diag(U, sigma_u); 
  cov_matrix[n] Sigma_v = quad_form_diag(V, sigma_v); 
  matrix[p, p] B = diag_matrix(beta);
  matrix[n, p] M = X * B;
  // for spatial model the covariance will be either a tensor normal or a 
  // matrix normal with one covariance structure being a kronecker product to 
  // capture the spatial grid
}
model{
  beta ~ normal(0, 1);
  sigma_u ~ cauchy(0, 5); // prior on the standard deviations
  sigma_v ~ cauchy(0, 5); // prior on the standard deviations
  U ~ lkj_corr(1); // LKJ prior on the correlation matrix 
  V ~ lkj_corr(1); // LKJ prior on the correlation matrix 
  target += dmatrix_normal(Y, M, U, V);
}
