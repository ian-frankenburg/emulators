functions{
}

data {
  int nobs; // number field nobservations
  int nsims; // number computer simulations
  int p; // number of non-calibrating inputs to simulation; for my exampe p=1
  int q; // number of calibration parameters; for my example q=1
  matrix[p, nobs] x_obs; // field observations inputs
  vector[nobs] y_obs; // noisy outcomes
  matrix[p, nobs*nsims] x_sims; // simulation locs
  matrix[q, nobs*nsims] t_sims; // simulation locs
  vector[nsims*nobs] y_sims; // outcome of simulations
}

transformed data{
  vector[nobs+nsims*nobs] y = append_row(y_obs,y_sims);
  vector[nobs+nsims*nobs] ones = rep_vector(1,nobs+nsims*nobs);
  matrix[p,nobs*nsims+nobs] x = append_col(x_obs,x_sims);
}

parameters {
  vector<lower=0,upper=1>[2] theta;
  vector<lower=0>[p+q] beta_eta;
  vector<lower=0>[p] beta_delta;
  real<lower=0> lambda_delta; 
  real<lower=0> lambda_eta; 
  real<lower=0> sigma2; 
  real mu;
}

transformed parameters{
  // beta_delta: correlation parameter for bias term
  // beta_e: correlation parameter of observation error
  matrix[q,nobs] theta_mat;
  matrix[q,nobs*nsims+nobs] t;
  for(i in 1:nobs){
    theta_mat[,i] = theta;
  }
  t = append_col(theta_mat,t_sims);
}

model {
  matrix[p+q, nsims*nobs+nobs] xt; // join simulator inputs, design points, and calibration parameters
  matrix[nsims*nobs+nobs, nsims*nobs+nobs] Sigma_eta; // simulator covariance
  matrix[nobs, nobs] Sigma_delta; // simulator covariance
  matrix[nsims*nobs+nobs, nsims*nobs+nobs] Sigma_z; // covariance matrix
  matrix[nsims*nobs+nobs, nsims*nobs+nobs] L; // cholesky decomposition of covariance matrix 
  vector[p+q] temp_eta;
  vector[p] temp_delta;
  
  xt = append_row(x,t);
  Sigma_delta = diag_matrix(rep_vector(1/lambda_delta+1e-10, nobs));

  for (i in 1:(nobs-1)) {
    for (j in (i+1):(nobs-1)) {
      temp_delta = x[,i] - x[,j];
      // Sigma_delta[i, j] = beta_delta .* temp_delta * temp_delta';
      Sigma_delta[i, j] = exp(-quad_form(diag_matrix(beta_delta),temp_delta))/lambda_delta;
      // Sigma_delta[i, j] = exp(-Sigma_delta[i, j]) / lambda_delta;
      Sigma_delta[j, i] = Sigma_delta[i, j];
    }
  }

  Sigma_eta = diag_matrix(rep_vector(1/lambda_eta+1e-10, nsims*nobs+nobs));
  
  for (i in 1:(nsims*nobs+nobs-1)) {
    for (j in (i+1):(nsims*nobs+nobs-1)) {
      temp_eta = xt[,i] - xt[,j];
      Sigma_eta[i, j] = exp(-quad_form(diag_matrix(beta_eta),temp_eta))/lambda_eta;
      Sigma_eta[j, i] = Sigma_eta[i, j];
    }
  }
  
  
  
  Sigma_z = Sigma_eta;
  Sigma_z[1:nobs,1:nobs] = Sigma_z[1:nobs,1:nobs] + diag_matrix(rep_vector(sigma2,nobs)) + Sigma_delta;
  
  theta[2] ~ uniform(.5,1);
  theta[1] ~ uniform(0,.5);
  beta_eta[1:(p+q)] ~ normal(0,1);
  beta_delta[1:(p)] ~ normal(0,1);
  sigma2 ~ normal(0,1); // gamma (shape, rate)
  lambda_eta ~ gamma(1,1);
  lambda_delta ~ gamma(1,1);
  L = cholesky_decompose(Sigma_z); // cholesky decomposition
  y ~ multi_normal_cholesky(mu*ones, L);
  mu ~ normal(0,1);
}

generated quantities{
}

