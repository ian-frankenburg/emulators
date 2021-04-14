data {
  int T; // length of time series
  int num_series; // number of time series simulations
  matrix[num_series,T] y; // input matrix of time series
  int p; // dimension of calibration parameters
  matrix[num_series-1, p] theta; // simulated parameter locations
  int lags; //lag for AR(lags) model
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
  vector[lags] phi[T];
  row_vector<lower=0,upper=1>[p] rho;
  row_vector<lower=0,upper=1>[p] calibrate;
  //vector<lower=0>[T-lags] sigma;
  real<lower=0> tau;
}

transformed parameters{
  matrix[num_series, p] theta_joint = append_row(theta,calibrate);
  matrix[num_series, num_series] C; // covariance matrix
  matrix[num_series, num_series] L; // covariance matrix
  vector[T] sigma=rep_vector(1,T);
  row_vector[p] beta;
  beta = -5.0 * log(rho);
  C = diag_matrix(rep_vector(1,num_series));
  for (i in 1:(num_series-1)) {
    for (j in (i+1):num_series) {
      C[i, j] = exp(-dot_product(beta .* (theta_joint[i, ] - theta_joint[j, ]),(theta_joint[i, ] - theta_joint[j, ])));
      C[j, i] = C[i, j];
    }
  }
  L = cholesky_decompose(C);
}

model {
  phi[1,] ~ normal(0,1);
  for(i in 2:(T)){
    phi[i,] ~ multi_normal(phi[i-1,], diag_matrix(rep_vector(tau,lags)));
  }
  tau ~ std_normal();
  rho ~ beta(1,.3);
  calibrate[1] ~ uniform(0,1);
  calibrate[2] ~ uniform(0,1);
  for (n in 2:T){
    y[,n] ~ multi_normal_cholesky(to_matrix(F[,,n-lags])*phi[n,], sqrt(sigma[n])*L);
  }
}
