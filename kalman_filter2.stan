data {
  int<lower=0> N; // length of time series
  int q; // number of time series
  int p; //  dimension of latent state e.g. lag order if AR dynamic regression matrix
  vector[q] y[N];
}

transformed data{
  // matrix[N-p,p] F;
  // for(n in (p+1):N){
  //   for(j in 1:p){
  //     F[n-p,j] = y[n-j];
  //   }
  // }
  real F[p,q,N];
  real G[p,p,N];
  for(k in 1:N){
    for(i in 1:p){
      for(j in 1:q){
        F[i,j,k] = 1;
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
  real Q[q,q,N];
  vector[p] m[N];
  vector[q] f[N];
  vector[q] error[N];
  matrix[p, p] C_init = diag_matrix(c_init);
  matrix[q,q] V = diag_matrix(rep_vector(1, q));
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
    f[t] = to_matrix(F[,,t])' * a[t];
    Q[,,t] = to_array_2d(to_matrix(F[,,t])' * R * to_matrix(F[,,t]) + V);
    ll[t] = multi_normal_lpdf(y[t] | f[t], to_matrix(Q[,,t]));
    // filter
    error[t] = y[t] - f[t];
    m[t] = a[t] + R * to_matrix(F[,,t]) * inverse(to_matrix(Q[,,t])) * error[t];
    C[,,t] = to_array_2d(R - R * to_matrix(F[,,t]) * inverse(to_matrix(Q[,,t])) * 
    to_matrix(R * to_matrix(F[,,t]))');
  }
}

model {
  delta ~ beta(1,1);
  c_init ~ std_normal();
  m_init ~ normal(y[1],.1);
  target += sum(ll);
}
generated quantities{
  vector[q] ypred[N];
  for(t in 1:N)
    ypred[t] = multi_normal_rng(f[t], to_matrix(Q[,,t]));
}

