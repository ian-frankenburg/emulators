library(lhs)
library(geostats)
library(deSolve)
library(cmdstanr)
library(rstan)
sir <- function(t, y, params) {
  with(as.list(c(params, y)), {
    S = y[1]
    I = y[2]
    dS = -beta*S*I/N
    dI = beta*S*I/N-gamma*I
    list(c(dS,dI))
  })
}
N=1000
nobs = 20
x = seq(1,nobs,length.out = nobs)
simulation<- function(beta,gamma,x){
  out = ode(c(N-1,1), x, sir, list(beta=beta,gamma=gamma, method="ode45"))[,3]
  out
}
y = simulation(.7,.25,x) + rnorm(length(x),1,1)
plot(x,y,xlim=c(min(x),max(x)),pch=19)

############ FOR DYNAMIC OUTPUT ############
# nsims=10
# x_design = rep(x,nsims)
# cube = randomLHS(nsims,k = 1)
# theta_design = qunif(cube,0,10)
# y_design=rep(0,nsims*nobs)
# c_design_rep=c()
# for(i in 1:nsims){
#   ##y_design[(1+(i-1)*length(x)):(i*length(x))] = simulation(theta_design[i],x)
#   #y_design[i] = simulation(theta_design[i],x[i])
#   #c_design_rep = c(rep(theta_design[i],nobs),c_design_rep)
# }
# plot(x_design,y_design,col="red",cex=2)
# points(x,y,xlim=c(0,max(x)),pch=19)
########################
## FOR UNIVARIATE OUTPUT
sims=6
nsims = sims*nobs
x_design = rep(x,nsims/length(x))
cube = randomLHS(nsims/length(x),k = 2)
theta_design = matrix(0,nrow=nsims/length(x),ncol=2)
theta_design[,2] = qunif(cube[,2],.2,.3)
theta_design[,1] = qunif(cube[,1],theta_design[,2],1)
y_design=rep(0,nsims)
k=1
for(i in 1:nrow(theta_design)){
  y_design[(1+(k-1)*length(x)):(k*length(x))] = simulation(theta_design[i,1],theta_design[i,2],x)
  k=k+1
}
plot(x,simulation(.7,.25,x),type="l",col=rgb(0, 0, 0, alpha = 1),lwd=6,lty=1,ylim=c(0,max(y_design)),xlab="Days",ylab="Case Counts",main="Emulation & Calibration Simulation")
#points(simulation(.68,0.32,x),type="l",col=rgb(.698, .133, .133, alpha = .9),lwd=6)
points(x_design,y_design, col = rgb(39/255, 116/255, 174/255, alpha = 0.7),
       pch = 16, cex = 2)
points(x,y, col = rgb(1, 184/255, 28/255, 0, alpha = 0.7),
       pch = 16, cex = 2.5)
legend(1, 400, legend=c("Truth", "Observed","Simulated","Estimated"),
       col=c("black", rgb(1, 184/255, 28/255, 0, alpha = 0.7),rgb(39/255, 116/255, 174/255),rgb(.698, .133, .133)), lty=1, cex=0.8,seg.len=1,lwd=8)
# int nobs; // number field nobservations
# int nsims; // number computer simulations
# int p; // number of non-calibrating inputs to simulation; for my exampe p=1
# int c; // number of calibration parameters; for my example c=1
# matrix[nobs, p] x; // field observations inputs
# vector[nobs] y; // noisy outcomes
# matrix[nsims,p] x_design; // simulation locs
# matrix[nsims,c] c_design; // simulation locs
# vector[nsims] y_design; // outcome of simulations
#real sigma2;
c_design = rep(0)
k=1
for(i in 1:length(theta_design)){
  c_design[(1+(k-1)*length(x)):(k*length(x))] = rep(theta_design[i],length(x))
  k=k+1
}
model_data = list(nobs=length(y),
                  nsims=sims,
                  q=2,
                  x_obs=matrix(x,ncol=1),
                  y_obs=y,
                  x_sims=matrix(x_design,ncol=sims),
                  t_sims=t(theta_design),
                  y_sims=matrix(y_design,ncol=sims))
file <- file.path("multivariate_w_univariate_2.stan")
model <- cmdstan_model(file)
mcmc <- model$sample(data = model_data,
                     max_treedepth = 10,adapt_delta = .99,
                     chains = 3,
                     iter_warmup = 250,
                     iter_sampling = 250,
                     refresh = 1, parallel_chains = getOption("mc.cores", 4))
stanfit <- rstan::read_stan_csv(mcmc$output_files())
summary(stanfit,pars=c("theta[1]","theta[2]","sigma2","mu"))$summary
a=traceplot(stanfit,pars=c("theta"))
traceplot(stanfit,pars=c("theta","sigma2"))
arr <- as.array(stanfit)
matplot(arr[,c(1,2,3),1])





# sir <- function(t, y, params) {
#   with(as.list(c(params, y)), {
#     S = y[1]
#     I = y[2]
#     dS = -beta*S*I/N
#     dI = beta*S*I/N-gamma*I
#     list(c(dS,dI))
#   })
# }
# N=1000
# x = seq(1,20,length.out = nobs)
# pts = 10
# design = expand.grid(beta=seq(0,1,length.out = pts),gamma=seq(.05,1,length.out = pts))
# Y = matrix(0,nrow=nrow(design),ncol=nobs)
# for(i in 1:nrow(design)){
#   Y[i,] = ode(c(N-1,1), seq(0,nobs-1), sir, list(beta=design[i,1],gamma=design[i,2]), method="ode45")[,3]
# }
# matplot(t(Y),type="l")
# cases = ode(c(N-1,1), seq(0,nobs-1), sir, list(beta=.8,gamma=.2), method="ode45")[,3]+rnorm(nobs,1)
# points(cases)
# 
# y = cases + rnorm(length(x),0,1)
# plot(x,y,xlim=c(0,max(x)),ylim=c(0,max(y)),pch=19)
# xx = sample(seq(1,20,length.out = 10),size = nobs,replace = F) # field observations
# plot(y[sort(xx)])
# x=sort(xx)
# nsims=10
# x_design = rep(x,nsims)
# cube = randomLHS(length(x_design),k = 1)
# theta_design = cube#qunif(cube,0,1)
# y_design=c()
# for(i in 1:length(theta_design)){
#   y_design[i] =  ode(c(N-1,1), seq(0,nobs-1), sir, list(beta=theta_design[i],gamma=.2), method="ode45")[x_design[i],3]
# }
# points(matrix(x_design,nrow=1),matrix(y_design,nrow=1),col="red")
# 
# # int nobs; // number field nobservations
# # int nsims; // number computer simulations
# # int p; // number of inputs to simulation; for my exampe p=1
# # int c; // number of calibration parameters; for my example c=1
# # matrix[nobs,p] x; // field observations inputs
# # vector[nobs] y; // noisy outcomes
# # matrix[nsims,p] x_design; // simulation locs
# # matrix[nsims,c] c_design; // simulation locs
# # vector[nsims] y_design; // outcome of simulations
# model_data = list(nobs=length(y),
#                   nsims=length(x_design),
#                   p=1,
#                   c=1,
#                   x=matrix(x,ncol=1),
#                   y=y,
#                   x_design=matrix(x_design,ncol=1),
#                   c_design=matrix(theta_design,ncol=1),
#                   y_design=y_design)
# file <- file.path("simple_emulation.stan")
# model <- cmdstan_model(file)
# mcmc <- model$sample(data = model_data,
#                      max_treedepth = 10,adapt_delta = .8,
#                      chains = 4,
#                      iter_warmup = 50,
#                      iter_sampling = 50,
#                      refresh = 1, parallel_chains = getOption("mc.cores", 4))
# stanfit <- rstan::read_stan_csv(mcmc$output_files())
# summary(stanfit,pars=c("theta","sigma2"))$summary
# traceplot(stanfit,pars=c("theta","sigma2"))
