data {
  int<lower=0> N; // number of observations
  int p; // dimension of calibration parameters
  vector[N] y; // observed outcomes
  matrix[N, p] x; // input covariate matrix
  matrix[N, p] theta; // simulation parameters at design points
  row_vector[p] theta_star; // predict response at these parameter values
  row_vector[p] x_star; // predict response at these covariate values
}

parameters {
  vector[p] lambda;
  row_vector<lower=0,upper=1>[p] rho;
}

transformed parameters{
  matrix[rows(theta), rows(theta)] C; // covariance matrix
  matrix[rows(theta), rows(theta)] L; // covariance matrix
  vector[N] sigma=rep_vector(1,N);
  row_vector[p] beta;
  beta = -5.0 * log(rho);
  C = diag_matrix(rep_vector(1,rows(theta)));
  
  for (i in 1:(rows(theta)-1)) {
    for (j in (i+1):rows(theta)) {
      C[i, j] = exp(-dot_product(beta .* (theta[i, ] - theta[j, ]),(theta[i, ] - theta[j, ])));
      C[j, i] = C[i, j];
    }
  }
  L = cholesky_decompose(C);
}

model {
  lambda ~ normal(0,10);
  rho ~ beta(1,.5);
  y ~ multi_normal_cholesky(x*lambda, L);
}

generated quantities{
  vector[rows(theta)] gamma_pred;
  vector[N] mu_pred;
  vector[N] sigma_pred;
  vector[N] y_pred=rep_vector(0,N);
  matrix[rows(theta), rows(theta)] Sigma_inv = inverse(C); // covariance matrix
  for(i in 1:rows(theta)){
    gamma_pred[i] = exp(-dot_product(beta .* (theta_star - theta[i, ]),(theta_star - theta[i, ])));
  }
  for(i in (1):N){
    mu_pred[i] = dot_product(lambda,x_star)+dot_product(gamma_pred,Sigma_inv*(y -  x*lambda));
    sigma_pred[i] = (1-dot_product(gamma_pred,Sigma_inv*gamma_pred));
    y_pred[i] = normal_rng(mu_pred[i],sigma_pred[i]);
  }
}
