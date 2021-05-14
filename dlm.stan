data {
  int<lower=0> N;
  vector[N] y;
  int npred;
  int p; // lag length
}

transformed data{
  matrix[N-p,p] F;
  for(n in (p+1):N){
    for(j in 1:p){
      F[n-p,j] = y[n-j];
    }
  }
}

parameters {
  real<lower=0> tau;
  real<lower=0> sigma;
  real alpha;
  vector<lower=-1,upper=1>[p] phi_raw[N-p];
}

transformed parameters{
  vector[p] phi[N-p];
  phi[1,] = tau*phi_raw[1,];
  for(i in 2:(N-p)){
    phi[i,] = phi[i-1,] + tau*phi_raw[i,];
  }
}

model {
  int kk=1;
  //vector[N-p] mu = alpha;
  alpha ~ normal(0,10);
  for(i in 1:(N-p)){
    phi_raw[i,] ~ std_normal();
  }
  sigma ~ normal(0,1);
  for (n in (p+1):N) {
    y[n] ~ normal(alpha+F[n-p,]*phi[n-p,], sigma);
  }
}

generated quantities{
  vector[N] ypred;
  row_vector[p] Fpred;
  ypred[1:p] = y[1:p];
  for (n in (p+1):N) {
    for(j in 1:p){
      Fpred[j] = ypred[n-j];
    }
    ypred[n] = normal_rng(alpha+Fpred*phi[n-p,],sigma);
  }
}
