data {
  int<lower=0> N; // length of time series
  int p; //  dimension of latent state e.g. lag order if AR dynamic regression matrix
  int m;
  matrix[N,m] ymat;
  int AR;
}

transformed data{
  real F[p*m,1];
  real G[p,p,N];
  vector[m*N] y = to_vector(ymat);
  
  
  
  
  
  
  
  
  for(k in 1:N){
    for(i in 1:p){
      if(k<=p){
        if(AR==1){
          F[i,1,k] = y[k];
        }else{
          F[i,1,k] = 1;
        }
      }else{
        if(AR==1){
          F[i,1,k] = y[k-i];
        }else{
          F[i,1,k] = 1;
        }
      }
    }
    G[,,k] = to_array_2d(diag_matrix(rep_vector(1, p)));
  }
}

parameters{
  vector[p] m_init;
  vector<lower=0>[p] c_init;
  real<lower=0,upper=1> delta;
}
transformed parameters{
  vector[p] a[N];
  matrix[p, p] R;
  real C[p,p,N];
  matrix[p, p] P;
  vector[N] Q;
  vector[p] m[N];
  vector[N] f;
  vector[N] error;
  // induce correlation here
  matrix[p, p] C_init = diag_matrix(c_init);
  real V = 1;
  vector[N] ll = rep_vector(0, N);
  for(t in 1:N){
    if(t==1){
      // predict
      a[t] = to_matrix(G[,,t])*m_init;
      C[,,t] = to_array_2d(C_init);
      P = to_matrix(G[,,t]) * to_matrix(C[,,t]) * to_matrix(G[,,t])';
      R = to_matrix(G[,,t]) * to_matrix(C[,,t]) * to_matrix(G[,,t])' + (1-delta)/delta*P;
    }else{
      // predict
      a[t] = to_matrix(G[,,t])*m[t-1];
      P = to_matrix(G[,,t]) * to_matrix(C[,,t-1]) * to_matrix(G[,,t])';
      R = to_matrix(G[,,t]) * to_matrix(C[,,t-1]) * to_matrix(G[,,t])' + (1-delta)/delta*P;
    }
    // marginalize
    f[t] = to_vector(F[,1,t])' * a[t];
    Q[t] = to_vector(F[,1,t])' * R * to_vector(F[,1,t]) + V;
    ll[t] = normal_lpdf(y[t] | f[t], Q[t]);
    // filter
    error[t] = y[t] - f[t];
    m[t] = a[t] + R * to_vector(F[,1,t]) * 1/Q[t] * error[t];
    C[,,t] = to_array_2d(R - R * to_vector(F[,1,t]) * 1/Q[t] * to_vector(F[,1,t])' * R');
  }
}

model {
  delta ~ beta(5,1);
  c_init ~ std_normal();
  m_init ~ normal(0,1);
  target += sum(ll);
}
generated quantities{
  vector[N] ypred;
  for(t in 1:N)
    ypred[t] = normal_rng(f[t], Q[t]);
}

