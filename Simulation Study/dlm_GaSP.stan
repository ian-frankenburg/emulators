data {
  int T; // length of time series
  int num_series; // number of time series simulations
  matrix[num_series,T] y; // input matrix of time series
  int p; // dimension of calibration parameters
  matrix[num_series, p] theta; // simulated parameter locations
  int lags; //lag for AR(lags) model
  row_vector[T] z;
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
  vector<lower=0>[T] tau;
  real<lower=0> sig_z;
  vector[lags] phi[T];
  row_vector<lower=0,upper=1>[p] rho;
  vector<lower=0,upper=1>[p] calibrate;
}

transformed parameters{
  matrix[rows(theta), rows(theta)] C; // covariance matrix
  matrix[rows(theta), rows(theta)] L; // covariance matrix
  row_vector[p] beta;
  vector<lower=1>[T] sigma=rep_vector(1,T);
  vector<lower=0,upper=1>[rows(theta)] r;
  beta = -5.0 * log(rho);
  C = diag_matrix(rep_vector(1,rows(theta)));
  for (i in 1:(rows(theta)-1)) {
    for (j in (i+1):rows(theta)) {
      C[i, j] = exp(-dot_product(beta .* (theta[i, ] - theta[j, ]),(theta[i, ] - theta[j, ])));
      C[j, i] = C[i, j];
    }
  }
  L = inverse(C);
  for(j in 1:rows(theta)){
    r[j] = exp(-dot_product(beta .* (calibrate' - theta[j, ]),(calibrate' - theta[j, ])));
  }
}

model {
  phi[1,] ~ normal(0,1);
  for(i in 2:(T)){
    phi[i,] ~ multi_normal(phi[i-1,],diag_matrix(rep_vector(tau[i],lags)));
  }
  tau ~ std_normal();
  rho ~ beta(1,.3);
  calibrate[1] ~ uniform(0,.5);
  calibrate[2] ~ uniform(.5,1);
  sig_z ~ gamma(1,1);
  for (n in 2:T){
    y[,n] ~ multi_normal(to_matrix(F[,,n-lags])*phi[n,], sqrt(sigma[n])*C);
    z[n]~normal(phi[n,]*z[n-1]+r'*L*(y[,n]-to_matrix(F[,,n-lags]) * phi[n,]), sqrt(sig_z+1-r'*L*r));
  }
}

// generated quantities{
  //   matrix[T-lags,lags] Fpred;
  //   vector[rows(theta)] gamma_z;
  //   vector[T-lags] mu_z;
  //   vector[T-lags] sigma_z;
  //   vector[T] ypred;
  //   matrix[rows(theta), rows(theta)] Sigma_inv = inverse(C); // covariance matrix
  //   matrix[T-lags, rows(theta)] errors;
  //   ypred[1:lags] = inits;
  //   for(j in 1:rows(theta)){
    //     gamma_z[j] = exp(-dot_product(beta .* (theta_pred - theta[j, ]),(theta_pred - theta[j, ])));
    //   }
  //   for (n in (lags+1):T){
    //     errors[n-lags,] = to_row_vector(y[,n] - (to_matrix(F[,,n-lags])*phi[n-lags,]));
    //   }
  //   for(i in (lags+1):T){
    //     for(j in 1:lags){
      //       Fpred[i-lags,j] = ypred[i-j];
      //     }
    //     mu_z[i-lags] = Fpred[i-lags,]*phi[i-lags,]+sigma[i-lags]*dot_product(gamma_z,Sigma_inv*errors[i-lags,]');
//     sigma_z[i-lags] = 1/sigma[i-lags]*(1-dot_product(gamma_z,Sigma_inv*gamma_z));
//     ypred[i] = normal_rng(mu_z[i-lags],sqrt(sigma_z[i-lags]));
//   }
// }