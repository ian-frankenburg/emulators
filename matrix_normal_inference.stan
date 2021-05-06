#include matrix_normal_lpdf.stan
data {
  int T; // number of spatial locations
  int S; // length of time series
  matrix[T,S] Y; // input tensor of time series
  int p; // dimension of calibration parameters
  matrix[num_series, p] theta; // simulated parameter locations
  row_vector[p] theta_pred; // predict curve at these parameter values
  vector[lags] inits; // initial values for prediction
}
transformed data{
  int lags=1; //for simplicity, assume AR(1)
}

transformed data{
  real F[num_series,lags,T-lags];
  for(n in (lags+1):T){
    for(k in 1:num_series){
      for(j in 1:lags){
        F[k,j,n-lags] = y[k,n-j]; 
      }
    }
  }
}

parameters {
  real<lower=0> tau;
  vector[lags] phi[T];
  row_vector<lower=0,upper=1>[p] rho;
  matrix[num_series,lags] F_init;
}

transformed parameters{
  matrix[rows(theta), rows(theta)] C; // covariance matrix
  matrix[rows(theta), rows(theta)] L; // covariance matrix
  row_vector[p] beta;
  vector<lower=1>[T] sigma=rep_vector(1,T);
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
  to_vector(F_init) ~ multi_normal(rep_vector(0,num_series), 10*diag_matrix(rep_vector(1,num_series)));
  phi[1,] ~ normal(0,1);
  for(i in 2:(T)){
    phi[i,] ~ multi_normal(phi[i-1,],diag_matrix(rep_vector(tau,lags)));
  }
  tau ~ std_normal();
  rho ~ beta(1,.3);
  for (n in 1:T){
    if(n==1){
      y[,n] ~ multi_normal_cholesky(F_init*phi[1,], sqrt(sigma[1])*L);
    }else{
      y[,n] ~ multi_normal_cholesky(to_matrix(F[,,n-lags])*phi[n,], sqrt(sigma[n])*L);
    }
  }
}

generated quantities{
  matrix[T-lags,lags] Fpred;
  vector[rows(theta)] gamma_z;
  vector[T-lags] mu_z;
  vector[T-lags] sigma_z;
  vector[T] ypred;
  matrix[rows(theta), rows(theta)] Sigma_inv = inverse(C); // covariance matrix
  matrix[T-lags, rows(theta)] errors;
  ypred[1:lags] = inits;
  for(j in 1:rows(theta)){
    gamma_z[j] = exp(-dot_product(beta .* (theta_pred - theta[j, ]),(theta_pred - theta[j, ])));
  }
  for (n in (lags+1):T){
    errors[n-lags,] = to_row_vector(y[,n] - (to_matrix(F[,,n-lags])*phi[n-lags,]));
  }
  for(i in (lags+1):T){
    for(j in 1:lags){
      Fpred[i-lags,j] = ypred[i-j];
    }
    mu_z[i-lags] = Fpred[i-lags,]*phi[i-lags,]+sigma[i-lags]*dot_product(gamma_z,Sigma_inv*errors[i-lags,]');
    sigma_z[i-lags] = 1/sigma[i-lags]*(1-dot_product(gamma_z,Sigma_inv*gamma_z));
    ypred[i] = normal_rng(mu_z[i-lags],sqrt(sigma_z[i-lags]));
  }
}

