# run vanilla GaSP file first
model_data = list('T'=length(scaley),num_series=sims,
                  y=t(matrix(y_design,ncol=sims,nrow=length(scaley))),p=2,z=c(scaley),
                  theta=theta_design,num_pred=1,
                  theta_pred=c(1,1),lags=1,inits=as.array(ytest[1:1]))
file <- file.path("dlm_calibrator.stan") 
model <- cmdstan_model(file)
mcmc <- model$sample(data = model_data,
                     max_treedepth = 10,adapt_delta = .8,
                     chains = 4,
                     iter_warmup = 1000,
                     iter_sampling =1000,
                     refresh = 25, parallel_chains = getOption("mc.cores", 4))
hist(mcmc$draws("calibrate[1]"),breaks=100)
hist(mcmc$draws("calibrate[2]"),breaks=100)
mcmc$summary("calibrate")

infer = rev(mcmc$summary(variables="calibrate")$median)
q05 = mcmc$summary(variables="calibrate")$q5
q95 = mcmc$summary(variables="calibrate")$q95

plot(x=sim_pars[,1],y=sim_pars[,2],xlab="beta",ylab="gamma",pch=19,cex=2,xlim=c(.5,1),ylim=c(0,.5))
points(x=sim_test[1],y=sim_test[2],col="red",pch=19,cex=2)
points(x=infer[1],y=infer[2],col="blue",pch=19,cex=2)
abline(v=q05[1])
abline(v=q95[1])
abline(h=q05[2])
abline(h=q95[2])

mcmc$summary("calibrate")


# int T; // length of time series
# int num_series; // number of time series simulations
# matrix[num_series,T] y; // input matrix of time series
# int p; // dimension of calibration parameters
# matrix[num_series-1, p] theta; // simulated parameter locations
# int lags; //lag for AR(lags) model
# vector[lags] phi[T];
# row_vector<lower=0,upper=1>[p] rho;
# real<lower=0> tau;
phi_draw=mcmc$draws("phi")
phi=matrix(0,nrow=nobs,ncol=1)
for(i in 1:dim(phi_draw)[3]){
  phi[i] = mean(phi_draw[,,i])
}
plot(c(phi))
rho_draw=mcmc$draws("rho")
rho=c(0,0)
for(i in 1:dim(rho_draw)[3]){
  rho[i] = mean(rho_draw[,,i])
}
tau=mean(mcmc$draws("tau"))
model_data = list('T'=nrow(ysims),num_series=ncol(ysims),
                  yold=t(ysims),y=ytest,p=2,
                  theta=matrix(sim_pars,ncol=2),phi=phi,rho=rho,tau=tau,
                  theta_pred=sim_test,lags=1,inits=as.array(ytest[1:1]))
file <- file.path("dlm_calibrator_modularization.stan")
model <- cmdstan_model(file)
mcmc2 <- model$sample(data = model_data,
                      max_treedepth = 10,adapt_delta = .9,
                      chains = 4,
                      iter_warmup = 500,
                      iter_sampling =500,
                      refresh = 500, parallel_chains = getOption("mc.cores", 4))
mcmc2$summary(variables="calibrate")
plot(x=sim_pars[,1],y=sim_pars[,2],xlab="beta",ylab="gamma",pch=19,cex=2)
points(x=sim_test[1],y=sim_test[2],col="red",pch=19,cex=2)
points(x=mean(mcmc2$draws("calibrate[1]")),y=mean(mcmc2$draws("calibrate[2]")),col="blue",pch=19,cex=2)
