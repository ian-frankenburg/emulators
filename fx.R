library(quantmod)
library(cmdstanr)
library(tfplot)
getFX("USD/NZD")
getFX("USD/AUD")
y1 = c((ts(USDNZD)))
y = c((ts(USDAUD)))
getSymbols("GOOGL", from = '2019-01-01',
           to = "2021-03-01",warnings = FALSE,
           auto.assign = TRUE)
y1=(as.numeric(GOOGL$GOOGL.Close))
y2=y1+rnorm(length(y1),0,.05)
z=c(scale(log(y1),scale=T))
matplot(z,type="l")
#y=cbind(y1,y2)
matplot(z,type="l")
# int<lower=0> N; // length of time series
# int q; // number of time series
# int p; //  lag order if AR dynamic regression matrix
# vector[q] y[N];
matplot(z,type="l")
compiled = cmdstan_model(stan_file="dlm.stan")
n_obs = 100
y=z[1:n_obs]
model_dat = list("y"=y, p=1, q=1, N=length(y), AR=0)
onestep=c()
for(i in 1:20){
  sample <- compiled$sample(
    data = model_dat,
    chains = 2,
    parallel_chains = 8,refresh=50,
    iter_warmup= 250,iter_sampling = 250,adapt_delta = .8,
    max_treedepth = 10
  )
  onestep[i] = sample$summary("onestep")$mean
  y = z[1:(n_obs+i)]
  model_dat = list("y"=y, p=1, q=1, N=length(y), AR=1)
}
plot(c(z[1:100],onestep))
lines(z[1:120])
lines(ypred$mean,pch=19)
sample$summary("delta")
ypred=sample$summary("f")
matplot(ypred$mean,type="l", lwd = 2, lty=1, col="darkblue")
matplot(y,add=T,col="darkred",type="l",lwd=3)
matplot(ypred$q5,type="l",col="darkblue",add=T,lwd=2)
matplot(ypred$q95,type="l",col="darkblue",add=T,lwd=2)
