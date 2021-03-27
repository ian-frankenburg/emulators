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
  matrix[lags,lags] P0 = diag_matrix(rep_vector(1,lags)); 
  for(n in (lags+1):T){
    for(j in 1:lags){
      for(k in 1:num_series){
        F[k,j,n-lags] = y[n-j,k]; 
      }
    }
  }
}

parameters {
  real<lower=0>tau;
  vector[lags] m0; 
  row_vector<lower=0,upper=1>[p] rho;
  //real<lower=0> sig;
  //vector<lower=1>[T-lags] sigma;
  // vector[T-lags] mu;
}

transformed parameters{
  matrix[rows(theta), rows(theta)] C; // covariance matrix
  row_vector<lower=0>[p] beta;
  vector[T-lags] sigma=rep_vector(1,T-lags);
  beta = -5.0 * log(rho);
  C = diag_matrix(rep_vector(1+1e-8,rows(theta)));
  for (i in 1:(rows(theta)-1)) {
    for (j in (i+1):rows(theta)) {
      C[i, j] = exp(-dot_product(beta .* (theta[i, ] - theta[j, ]),(theta[i, ] - theta[j, ])));
      C[j, i] = C[i, j];
    }
  }
  C = C;
}

model {
  matrix[rows(theta), rows(theta)] L; // covariance matrix
  //Filtering 
  // prediction equations
  vector[lags] m_pred[T-lags]; //filtered mean x_t|y_1:t
  real P_pred[lags,lags, T-lags]; //predicted var x_t|y_1:t-1
  // update equations
  vector[rows(theta)] v; //filtered mean x_t|y_1:t
  real S[rows(theta),rows(theta), T-lags]; //filtered var x_t|y_1:t
  matrix[lags,rows(theta)] K; //filtered var x_t|y_1:t
  vector[lags] m[T-lags]; //filtered mean x_t|y_1:t
  real P[lags,lags, T-lags]; //filtered var x_t|y_1:t
  
  m_pred[1,] = m0;
  P_pred[,,1] = to_array_2d(P0);

  for (t in 1:(T-lags)) {
    //Prediction step
    if (t>1) {
        m_pred[t,] = m[t-1,];
        P_pred[,,t] = to_array_2d(to_matrix(P[,,t-1]) + diag_matrix(rep_vector(tau,lags)));
    }

    //Update step
    v = y[t,]' - to_matrix(F[,,t])*m_pred[t,];
    S[,,t] = to_array_2d(to_matrix(F[,,t])*to_matrix(P_pred[,,t])*to_matrix(F[,,t])' + sigma[t]*C);
    K = to_matrix(P_pred[,,t]) * to_matrix(F[,,t])' * inverse(to_matrix(S[,,t]));
    m[t,] = m_pred[t,] + K * v; 
    P[,,t] = to_array_2d(to_matrix(P_pred[,,t]) - K * to_matrix(S[,,t]) * K');
  }
  m0 ~ normal(0,5);
  tau ~ gamma(1,1);
  rho ~ beta(1,.3);
  //sigma ~ inv_gamma(1,1);
  for (t in (lags+1):T) {
    L = cholesky_decompose(to_matrix(S[,,t-lags]));
    y[t,] ~ multi_normal_cholesky(to_matrix(F[,,t-lags])*m_pred[t-lags,], L);
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
  vector[lags] phi[T-lags];
  
  //Filtering 
  // prediction equations
  vector[lags] m_pred[T-lags]; //filtered mean x_t|y_1:t
  real P_pred[lags,lags, T-lags]; //predicted var x_t|y_1:t-1
  // update equations
  vector[rows(theta)] v; //filtered mean x_t|y_1:t
  real S[rows(theta),rows(theta), T-lags]; //filtered var x_t|y_1:t
  matrix[lags,rows(theta)] K; //filtered var x_t|y_1:t
  vector[lags] m[T-lags]; //filtered mean x_t|y_1:t
  real P[lags,lags, T-lags]; //filtered var x_t|y_1:t
  
  m_pred[1,] = m0;
  P_pred[,,1] = to_array_2d(P0);

  for (t in 1:(T-lags)) {
    //Prediction step
    if (t>1) {
        m_pred[t,] = m[t-1,];
        P_pred[,,t] = to_array_2d(to_matrix(P[,,t-1]) + diag_matrix(rep_vector(tau,lags)));
    }

    //Update step
    v = y[t,]' - to_matrix(F[,,t])*m_pred[t,];
    S[,,t] = to_array_2d(to_matrix(F[,,t])*to_matrix(P_pred[,,t])*to_matrix(F[,,t])' + sigma[t]*C);
    K = to_matrix(P_pred[,,t]) * to_matrix(F[,,t])' * inverse(to_matrix(S[,,t]));
    m[t,] = m_pred[t,] + K * v; 
    P[,,t] = to_array_2d(to_matrix(P_pred[,,t]) - K * to_matrix(S[,,t]) * K');
  }
  
  for(t in 1:(T-lags)){
    phi[t,] = multi_normal_rng(m_pred[t,],to_matrix(P_pred[,,t]));
  }
  ypred[1,] = inits;
  for(k in 1:rows(theta_pred)){
    for(j in 1:rows(theta)){
      gamma_z[j] = exp(-dot_product(beta .* (theta_pred[k,] - theta[j, ]),(theta_pred[k,] - theta[j, ])));
    }
    for (n in (lags+1):T){
      errors[n-lags,] = y[n,] - (to_matrix(F[,,n-lags])*phi[n-lags,])';
    }
    for(i in (lags+1):T){
      for(j in 1:lags){
        Fpred[i-lags,j] = ypred[i-j,k];
      }
      mu_z[i-lags] = Fpred[i-lags,]*phi[i-lags,]+sigma[i-lags]*dot_product(gamma_z,Sigma_inv*errors[i-lags,]');
      sigma_z[i-lags] = 1/sigma[i-lags]*(1-dot_product(gamma_z,Sigma_inv*gamma_z));
      ypred[i,k] = normal_rng(mu_z[i-lags],sigma_z[i-lags]);
    }
  }
}
