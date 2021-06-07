library(lhs)
library(geostats)
library(deSolve)
library(cmdstanr)
library(rstan)
library(truncnorm)
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
y = yraw= simulation(.7,.2,x)
for(i in 1:length(y)){
  #yraw[i]=(rnbinom(1,mu = y[i], size=100))+1
  yraw[i]=(rpois(1,y[i]))
}
plot(x,(yraw),xlim=c(min(x),max(x)),pch=19)
for(i in 1:length(y)){
  y[i]=(yraw[i])
}
#x = seq(1,nobs,b=2)
#y = y[x]

sims=7
nsims=sims*length(y)
# x_design = sample(0:10,size = nsims,replace = T) # field observations
# x_design = sample(0:10,size = nsims,replace = T) # field observations
x_design = rep(x,sims)
cube = randomLHS(sims,k = 2)
theta_design = matrix(0,nrow=sims,ncol=2)
theta_design[,1] = qunif(cube[,1],0,.5)
theta_design[,2] = qunif(cube[,2],.5,1)
y_design=rep(0,sims*length(y))
k=1
for(i in 1:nrow(theta_design)){
  y_design[(1+(k-1)*length(x)):(k*length(x))] = (sqrt(simulation(theta_design[i,2],theta_design[i,1],seq(1,nobs,length.out = nobs))))
  k=k+1
}
scaley = (sqrt(y))
#x=(x-min(x))/(max(x)-min(x))
#plot(x,scaley,type="l",col=rgb(0, 0, 0, alpha = 1),lwd=6,lty=1,ylim=c(min(y_design),max(y_design)),xlab="Days",ylab="Case Counts",main="Emulation & Calibration Simulation")
#points(x,(simulation(.626,.2,1:nobs)),type="l",col=rgb(.698, .133, .133, alpha = .9),lwd=6)
#x_design=(x_design-min(x_design))/(max(x_design)-min(x_design))
# points(x_design,y_design, col = rgb(39/255, 116/255, 174/255, alpha = 0.7),
#        pch = 16, cex = 2)
plot(x,scaley, col = rgb(1, 184/255, 28/255, 0, alpha = 0.7),
     pch = 16, cex = 2.5)
k=1
lines(y_design[(1+(k-1)*length(x)):(k*length(x))],lwd=5, col = rgb(39/255, 116/255, 174/255, alpha = 0.7))
for(k in 2:nrow(theta_design)){
  lines(y_design[(1+(k-1)*length(x)):(k*length(x))],lwd=5, col = rgb(39/255, 116/255, 174/255, alpha = 0.7))
}
points(x,scaley, pch=19,cex=2,ylim=c(0,max(y_design)))
legend(1, 20, legend=c("Observed","Simulated","Estimated"),
       col=c("black", rgb(39/255, 116/255, 174/255),rgb(.698, .133, .133)), lty=1, cex=0.8,seg.len=1,lwd=8)
c_design = matrix(0,nrow=2,ncol=length(y)*sims)
k=1
for(i in 1:nrow(theta_design)){
  c_design[,(1+(k-1)*length(x)):(k*length(x))] = theta_design[i,]
  k=k+1
}
model_data = list(nobs=length(y),
                  nsims=sims,init=1,
                  p=1,
                  q=2,
                  x_obs=matrix(x,ncol=length(x)),
                  y_obs=c(scaley),
                  x_sims=matrix(x_design,ncol=length(x)*sims),
                  t_sims=(matrix(c_design,ncol=sims*length(y),nrow=2)),
                  y_sims=(y_design))
file <- file.path("GaSP.stan")
model <- cmdstan_model(file)
mcmc <- model$sample(data = model_data,
                     max_treedepth = 10,adapt_delta = .8,
                     chains = 4,
                     iter_warmup = 500,
                     iter_sampling = 500,
                     refresh = 1, parallel_chains = getOption("mc.cores", 4))
mcmc$summary("theta")
stanfit <- rstan::read_stan_csv(mcmc$output_files())
hist(mcmc$draws("theta[1]"),breaks=100)
hist(mcmc$draws("theta[2]"),breaks=100)


summary(stanfit,pars=c("theta[1]","theta[2]","sigma2","lambda_eta"))$summary
a=traceplot(stanfit,pars=c("theta"))
traceplot(stanfit,pars=c("theta"))
hist(mcmc$draws("theta[1]"),breaks=100)
hist(mcmc$draws("theta[2]"),breaks=100)


betaDraws <- as.array(stanfit)[,c(2,3,4),1]
plot(betaDraws[,1],type="l",lwd=4)
lines(betaDraws[,2],type="l",lwd=4,col="blue")
lines(betaDraws[,3],type="l",lwd=4, col="red")

gammaDraws <- as.array(stanfit)[,c(2,3,4),2]
plot(gammaDraws[,1],type="l",lwd=4)
lines(gammaDraws[,2],type="l",lwd=4,col="blue")
lines(gammaDraws[,3],type="l",lwd=4, col="red")

betaDraws = c(mcmc$draws("theta[1]"))
gammaDraws = c(mcmc$draws("theta[2]"))
for(i in 1:200){
  lines(x,sqrt(simulation(sample(betaDraws,1),sample(gammaDraws,1),1:nobs)),type="l",col=rgb(.698, .133, .133, alpha = .05),lwd=6)
}

p = c()
for(t in 1:30){
  p[t] = 1/sqrt(4*pi*alpha*t)*exp(-1/(4*alpha*t)*)
}
plot(p)


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