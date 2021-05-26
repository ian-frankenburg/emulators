functions{
  matrix kron(matrix A, matrix B) {
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  // for (i in 1:m) {
  //   for (j in 1:n) {
  //     int row_start;
  //     int row_end;
  //     int col_start;
  //     int col_end;
  //     row_start = (i - 1) * p + 1;
  //     row_end = (i - 1) * p + p;
  //     col_start = (j - 1) * q + 1;
  //     col_end = (j - 1) * q + 1;
  //     C[row_start:row_end, col_start:col_end] = A[i, j] * B;
  //   }
  // }
  for (i in 1:m)
    for (j in 1:n)
      for (k in 1:p)
        for (l in 1:q)
          C[p*(i-1)+k,q*(j-1)+l] = A[i,j]*B[k,l];
          
  return C;
}
}
data {
  int T; // length of time series
  int num_series; // number of time series simulations
  int p; // dimension of calibration parameters
  matrix[num_series, p] theta; // simulated parameter locations
  row_vector[p] theta_pred; // predict curve at these parameter values
  int lags; //lag for AR(lags) model
  vector[lags] inits; // initial values for prediction
  int S;// n^2 spatial grid points
  matrix[num_series*S,T] y; // input matrix of time series
}

transformed data{
  real F[num_series*S,lags,T-lags];
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
  //matrix[num_series*S,lags] F_init;
}

transformed parameters{
  matrix[rows(theta), rows(theta)] C; // covariance matrix
  matrix[rows(theta)*S, rows(theta)*S] L; // covariance matrix
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
  // matrix size compatibilty error!
  L = cholesky_decompose(kron(diag_matrix(rep_vector(1,S)),C));
}

model {
  //to_vector(F_init) ~ multi_normal(rep_vector(0,num_series*S), 10*diag_matrix(rep_vector(1,num_series*S)));
  phi[1,] ~ normal(0,1);
  for(i in 2:(T)){
    phi[i,] ~ multi_normal(phi[i-1,],diag_matrix(rep_vector(tau,lags)));
  }
  tau ~ std_normal();
  rho ~ beta(1,.3);
  for (n in 2:T){
    y[,n] ~ multi_normal_cholesky(to_matrix(F[,,n-lags])*phi[n,], sqrt(sigma[n])*L);
  }
}

generated quantities{
  matrix[T-lags,lags] Fpred;
  vector[rows(theta)] gamma_z;
  vector[T-lags] mu_z;
  vector[T-lags] sigma_z;
  vector[T] ypred;
  matrix[rows(theta), rows(theta)] Sigma_inv = inverse(L); // covariance matrix
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
