kappa = 1
theta = 5
N = 1000

#Question 1
gammaSim <- function(N, kappa, theta){
  x <- rgamma(N, shape = kappa, scale = theta)
  hist(x, freq = FALSE, breaks = 30, col = "grey")
  
  #Plotting lines
  xfit = seq(min(gam), max(gam), length = 2000)
  yfit = dgamma(xfit, shape = kappa, scale = theta)
  
  lines(xfit, yfit)
  
  q=quantile(x,c(.1,.3,.5,.7,.9))
  abline(v=q)
  
  return (q)
}
#par(mfrow=(1, 1))
q = gammaSim(N, kappa, theta)
q

#Question 2
gammaPred <- function(par, q){

  #Parameters
  kappa = par[1]
  theta = par[2]
  
  #Cuts for quantiles
  cuts = c(0, q, Inf)
  p = c()
  
  for (i in (1:length(cuts)-1)){
    p[i] = pgamma(cuts[i+1], kappa, scale = theta) - pgamma(cuts[i], kappa, scale = theta)
  }
  return(p)
}
gammaPred(c(kappa, theta), q)

#Question 3
N = 500
kappa = 4
theta = 3
q = gammaSim(N, kappa, theta)
freq = c()

freq = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*N

#Question 4
par = c(kappa, theta)
gsqfun <- function(par, q, freq){
  p = gammaPred(par, q)
  p = p/sum(p)
  #print (p)
  
  pf = p*sum(freq)
  pf
  
  gslev=1:length(freq)
  gslev=gslev[freq>0]
  
  gsq=2*sum(freq[gslev]*log(freq[gslev]/pf[gslev]))
  return (gsq)
}

#Question 5 and 6
pars = data.frame(c(0, 0))

for (i in c(1:10)){
  par = c(2, 3)
  q = gammaSim(N, kappa, theta)
  comp = gsqfun(par, q, freq)
  
  while(TRUE){
    fit=optim(par,gsqfun,freq=freq,q=q)
    par=fit$par
    if((comp-fit$value)<=.001) break
    comp=fit$value
  }
  pars[i] = par  
}

library(matrixStats)
#First row of the dataframe is kappas
kappa_mean = apply(pars, 1, mean)[1]
kappa_sd = apply(pars, 1, sd)[1] 


#Second row is theta
theta_mean = apply(pars, 1, mean)[2]
theta_sd = apply(pars, 1, sd)[2] 



#Question 7
N = 50
pars = data.frame(c(0, 0))

for (i in c(1:10)){
  par = c(2, 3)
  q = gammaSim(N, kappa, theta)
  comp = gsqfun(par, q, freq)
  
  while(TRUE){
    fit=optim(par,gsqfun,freq=freq,q=q)
    par=fit$par
    if((comp-fit$value)<=.001) break
    comp=fit$value
  }
  pars[i] = par  
}

#First row of the dataframe is kappas
kappa_mean2 = apply(pars, 1, mean)[1]
kappa_sd2 = apply(pars, 1, sd)[1] 


#Second row is theta
theta_mean2 = apply(pars, 1, mean)[2]
theta_sd2 = apply(pars, 1, sd)[2] 


#Reading values for 6 and 7
kappa_mean 
kappa_mean2

kappa_sd
kappa_sd2


theta_mean
theta_mean2


theta_sd
theta_sd2


#The means remain around the same but the standard deviation for the runs with 10 samples is high
#This could be because low sample size causes non-consistent fitting of the curve; similar to 
#how probability of coin tosses starting from 1 or 0 would stabilize in the limit.