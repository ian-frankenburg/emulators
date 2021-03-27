data {
  int T; // length of time series
  int num_series; // number of time series simulations
  matrix[T, num_series] y; // input matrix of time series
  int p; // dimension of calibration parameters
  matrix[num_series, p] theta; // simulated parameter locations
  int num_pred;
  matrix[num_pred, p] theta_pred; // predict curve at these parameter values
  int lags; //lag for AR(lags) model
  vector[num_pred] inits; // initial values for prediction
}

transformed data{
  real F[num_series,lags,T-lags];
  for(n in (lags+1):T){
    for(j in 1:lags){
      for(k in 1:num_series){
        F[k,j,n-lags] = y[n-j,k]; 
      }
    }
  }
}

parameters {
  vector<lower=0>[1] tau;
  vector[lags] phi[T-lags];
  row_vector<lower=0,upper=1>[p] rho;
  //row_vector[p] alpha;
  //real<lower=0> sig;
  vector<lower=1>[T-lags] sigma;
  vector[T-lags] mu;
}

transformed parameters{
  matrix[rows(theta), rows(theta)] C; // covariance matrix
  matrix[rows(theta), rows(theta)] L; // covariance matrix
  //vector[T-lags] sigma=rep_vector(1,T-lags);
  row_vector[p] beta;
  row_vector[p] alpha=rep_row_vector(0,p);
  //vector[T-lags] mu=rep_vector(0,T-lags);

  beta = -5.0 * log(rho);
  C = diag_matrix(rep_vector(1,rows(theta)));
  for (i in 1:(rows(theta)-1)) {
    for (j in (i+1):rows(theta)) {
      C[i, j] = exp(-dot_product(beta .* (theta[i, ] - theta[j, ]),(theta[i, ] - theta[j, ])));
      C[j, i] = C[i, j];
    }
  }
  L = cholesky_decompose(C);
  // phi[1,] = tau[1] * phi_raw[1,];
  // for(t in 2:(T-lags)){
  //   phi[t,] = phi[t-1,] + tau[1] * phi_raw[t,];
  // }
}

model {
  matrix[rows(theta),p] alphaAlltemp;
  vector[rows(theta)] alphaAll;
  alphaAlltemp = diag_post_multiply(theta,alpha);
  for(j in 1:rows(theta)){
    alphaAll[j] = dot_product(alphaAlltemp[j,],rep_vector(1,p));
  }
  phi[1] ~ normal(0,1);
  sigma ~ normal(1,1);
  mu[1] ~ normal(0,1);
  alpha ~ normal(0,5);
  for(i in 2:(T-lags)){
    phi[i] ~ normal(phi[i-1],tau[1]);
    mu[i] ~ normal(mu[i-1],tau[1]);
  }
  tau[1] ~ std_normal();
  rho ~ beta(1,.3);
  for (n in (lags+1):T) {
    y[n,] ~ multi_normal_cholesky(mu[n-lags]+alphaAll+to_matrix(F[,,n-lags])*phi[n-lags,], sqrt(sigma[n-lags])*L);
  }
}

generated quantities{
  matrix[T-lags,lags] Fpred;
  vector[rows(theta)] gamma_z;
  vector[T-lags] mu_z;
  vector[T-lags] sigma_z;
  vector[num_pred] ypred[T];
  matrix[rows(theta), rows(theta)] Sigma_inv = inverse(C); // covariance matrix
  matrix[T-lags, rows(theta)] errors;
  matrix[rows(theta),p] alphaAlltemp;
  vector[rows(theta)] alphaAll;
  alphaAlltemp = diag_post_multiply(theta,alpha);
  ypred[1,] = inits;
  for(k in 1:rows(theta_pred)){
    for(j in 1:rows(theta)){
      alphaAll[j] = dot_product(alphaAlltemp[j,],rep_vector(1,p));
      gamma_z[j] = exp(-dot_product(beta .* (theta_pred[k,] - theta[j, ]),(theta_pred[k,] - theta[j, ])));
    }
    for (n in (lags+1):T){
      errors[n-lags,] = y[n,] - (mu[n-lags]+alphaAll + to_matrix(F[,,n-lags])*phi[n-lags,])';
    }
    for(i in (lags+1):T){
      for(j in 1:lags){
        Fpred[i-lags,j] = ypred[i-j,k];
      }
      mu_z[i-lags] = mu[i-lags]+dot_product(alpha,theta_pred[k,])+Fpred[i-lags,]*phi[i-lags,]+sigma[i-lags]*dot_product(gamma_z,Sigma_inv*errors[i-lags,]');
      sigma_z[i-lags] = 1/sigma[i-lags]*(1-dot_product(gamma_z,Sigma_inv*gamma_z));
      ypred[i,k] = normal_rng(mu_z[i-lags],sigma_z[i-lags]);
    }
  }
}
