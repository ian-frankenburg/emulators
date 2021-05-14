data {
  int T; // length of time series
  int num_series; // number of time series simulations
  matrix[num_series,T] yold; // input matrix of time series
  int p; // dimension of calibration parameters
  matrix[num_series, p] theta; // simulated parameter locations
  row_vector[T] y;
  int lags; //lag for AR(lags) model
  vector[lags] phi[T];
  row_vector<lower=0,upper=1>[p] rho;
  real<lower=0> tau;
}

transformed data{
  real F[num_series,lags,T-lags];
  row_vector[p] beta = -5.0 * log(rho);
  vector[T] sigma=rep_vector(1,T);
  for(n in (lags+1):T){
    for(k in 1:num_series){
      for(j in 1:lags){
        F[k,j,n-lags] = yold[k,n-j]; 
      }
    }
  }
}

parameters {
  real<lower=0> sig_z;
  row_vector<lower=0,upper=1>[p] calibrate;
}

transformed parameters{
  matrix[num_series, num_series] C; // covariance matrix
  matrix[num_series, num_series] L; // covariance matrix
  vector<lower=0,upper=1>[rows(theta)] r;
  for(j in 1:rows(theta)){
    r[j] = exp(-dot_product(beta .* (calibrate - theta[j, ]),(calibrate - theta[j, ])));
  }
  C = diag_matrix(rep_vector(1,num_series));
  for (i in 1:(num_series-1)) {
    for (j in (i+1):num_series) {
      C[i, j] = exp(-dot_product(beta .* (theta[i, ] - theta[j, ]),(theta[i, ] - theta[j, ])));
      C[j, i] = C[i, j];
    }
  }
  L = inverse(C);
}

model {
  calibrate[1] ~ uniform(0,1);
  calibrate[2] ~ uniform(0,1);
  sig_z ~ gamma(1,1);
  // predictive mean comprised of TVAR parameter + component from the training emulation runs
  for(t in 2:T){
    y[t]~normal(phi[t]*y[t-1]+r'*L*(yold[,t]-to_matrix(F[,,t-lags]) * phi[t,]), sqrt(sig_z+1-r'*L*r));
  }
}
