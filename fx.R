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
y1=log(as.numeric(GOOGL$GOOGL.Close))
y2=y1+rnorm(length(y1),0,.05)
y=cbind(y1,y2)
matplot(cbind(y1,y2),type="l")
#y=cbind(y1,y2)
matplot(y,type="l")
# int<lower=0> N; // length of time series
# int q; // number of time series
# int p; //  lag order if AR dynamic regression matrix
# vector[q] y[N];
model_dat = list("y"=matrix(y,ncol=2),"p"=1,q=2,N=nrow(y))
compiled = cmdstan_model(stan_file="kalman_filter2.stan")
sample <- compiled$sample(
  data = model_dat,
  chains = 2,
  parallel_chains = 8,refresh=5,
  iter_warmup= 250,iter_sampling = 250,adapt_delta = .9,
  max_treedepth = 10
)

ypred=matrix(sample$summary("ypred")$mean,ncol=2)
matplot(ypred,type="l")
matplot(ypred[,2],type="l",col="darkblue",add=T)
