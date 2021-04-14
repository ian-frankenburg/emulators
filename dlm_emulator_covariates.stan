data {
  int T; // length of time series
  int num_series; // number of time series simulations
  matrix[T, num_series] y; // input matrix of time series
  int p; // dimension of calibration parameters
  matrix[num_series, p] theta; // simulated parameter locations
  int q; // dimension of covariates
  matrix[T, q] x; // input matrix of time series
  row_vector[p] theta_pred; // predict curve at these parameter values
  row_vector[q] x_pred; // predict curve at these covariate values
  int lags; //lag for AR(lags) model
  vector[lags] inits; // initial values for prediction
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
  row_vector[p] alpha;
  row_vector[q] lambda;
  row_vector<lower=0,upper=1>[p] rho;
}

transformed parameters{
  matrix[rows(theta), rows(theta)] C; // covariance matrix
  matrix[rows(theta), rows(theta)] L; // covariance matrix
  real sig=1;
  vector[T-lags] sigma=rep_vector(sig,T-lags);
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
  matrix[rows(theta),p] alphaAlltemp;
  vector[rows(theta)] alphaAll;
  matrix[T,q] lambdaAlltemp;
  vector[T] lambdaAll;
  alphaAlltemp = diag_post_multiply(theta,alpha);
  for(j in 1:rows(theta)){
    alphaAll[j] = dot_product(alphaAlltemp[j,],rep_vector(1,p));
  }
  lambdaAlltemp = diag_post_multiply(x,lambda);
  for(j in 1:T){
    lambdaAll[j] = dot_product(lambdaAlltemp[j,],rep_vector(0,q));
  }
  phi[1] ~ normal(0,1);
  for(i in 2:(T-lags)){
    phi[i] ~ normal(phi[i-1],tau);
  }
  tau[1] ~ std_normal();
  rho ~ beta(1,.3);
  for (n in (lags+1):T) {
    y[n,] ~ multi_normal_cholesky(alphaAll + lambdaAll[n-lags]+to_matrix(F[,,n-lags])*phi[n-lags,], sqrt(sigma[n-lags])*L);
  }
}

generated quantities{
  matrix[T-lags,lags] Fpred;
  vector[rows(theta)] gamma_z;
  vector[T-lags] mu_z;
  vector[T-lags] sigma_z;
  vector[T] ypred=rep_vector(0,T);
  matrix[rows(theta), rows(theta)] Sigma_inv = inverse(C); // covariance matrix
  matrix[T-lags, rows(theta)] errors;
  vector[rows(theta)] temp;
  matrix[rows(theta),p] alphaAlltemp;
  vector[rows(theta)] alphaAll;
  matrix[rows(theta),q] lambdaAlltemp;
  vector[rows(theta)] lambdaAll;
  alphaAlltemp = diag_post_multiply(theta,alpha);
  lambdaAlltemp = diag_post_multiply(x,lambda);
  for(j in 1:rows(theta)){
    alphaAll[j] = dot_product(alphaAlltemp[j,],rep_vector(1,p));
    lambdaAll[j] = dot_product(lambdaAlltemp[j,],rep_vector(0,q));
    gamma_z[j] = exp(-dot_product(beta .* (theta_pred - theta[j, ]),(theta_pred - theta[j, ])));
  }
  for (n in (lags+1):T){
    errors[n-lags,] = y[n,] - (alphaAll + lambdaAll[n-lags] + to_matrix(F[,,n-lags])*phi[n-lags,])';
  }
  ypred[1:lags] = inits;
  for(i in (lags+1):T){
    for(j in 1:lags){
      Fpred[i-lags,j] = ypred[i-j];
    }
    mu_z[i-lags] = dot_product(alpha,theta_pred)+dot_product(lambda,x_pred)+Fpred[i-lags,]*phi[i-lags,]+sigma[i-lags]*dot_product(gamma_z,Sigma_inv*errors[i-lags,]');
    sigma_z[i-lags] = 1/sigma[i-lags]*(1-dot_product(gamma_z,Sigma_inv*gamma_z));
    ypred[i] = normal_rng(mu_z[i-lags],sigma_z[i-lags]);
  }
}
