exgsim=function(mu,sig,tau,N){
  norm=rnorm(N,mu,sig)
  exp=rexp(N,1/tau)
  exg=norm+exp
  if(min(exg)<0) print("RTs < 0 ; Try different par values")
  q=quantile(exg,c(.1,.3,.5,.7,.9))
  freq_bin=c(.1,.2,.2,.2,.2,.1)*N
  return(list(q, freq_bin))
}

#Simulate data from gamma and returns quantiles of the simulated distribution

gammaSim <- function(N, kappa, theta){
  x <- rgamma(N, shape = kappa, scale = theta)
  #hist(x, freq = FALSE, breaks = 30, col = "grey")
  
  #Plotting lines
  #xfit = seq(min(gam), max(gam), length = 2000)
  #yfit = dgamma(xfit, shape = kappa, scale = theta)
  
  #lines(xfit, yfit)
  
  q=quantile(x,c(.1,.3,.5,.7,.9))
  #abline(v=q)
  freq_bin=c(.1,.2,.2,.2,.2,.1)*N
  return (list(q, freq_bin))
}


pexg=function(q,mu,sig,tau){
  pnorm( ((1/tau)*(q-mu)),0,((1/tau)*sig) ) - 
    exp(  -(1/tau)*(q-mu) + ((1/tau)*sig)^2/2 + 
            log(pnorm(((1/tau)*(q-mu)),((1/tau)*sig)^2,
                      ((1/tau)*sig))) )
}

#this is something called the "error function"
# that the ExGaussian pdf function needs
erfc <- function(q) 2*pnorm(q*sqrt(2),lower=FALSE)

#pdf for an ExGaussian
# (You don't have to understand all this math)
dexg=function(q,mu,sig,tau) ((1/tau)/2) * 
  exp( ((1/tau)/2)*(2*mu+(1/tau)*sig^2-2*q) ) * 
  erfc( (mu+(1/tau)*sig^2-q)/(sqrt(2)*sig) )



#Returns predicted proportion of scores
# in each RT bin separated by a set of
# quantiles.
# par is a parameter vector with mu, sig,
# and tau in that order
# q is a set of quantiles
pred=function(par,q){
  mu=par[1]
  sig=par[2]
  tau=par[3]
  p=c()
  bcuts=c(0,q,Inf) #cutoffs for RT bins
  for(bi in 1:(length(bcuts)-1)){
    p[bi]=pexg(bcuts[bi+1],mu,sig,tau)-
      pexg(bcuts[bi],mu,sig,tau)
  }
  
  return(p)
}


#computes g-squared
# par is a parameter vector with mu, sig,
# and tau in that order
# dat is an observed set of frequencies
# q is a set of observed quantiles
gsqfun=function(par,f,q){
  p=pred(par,q) #get predicted proportions
  
  #If a cell is 0 in the predictions but not in
  # the data, fill in the prediction with a small
  # number. This is a trick to keep the function
  # from returning infinity.
  p[p==0 & f>0] = .0001
  p=p/sum(p) #renormalize so p adds up to 1 again
  
  pf=c() #predicted frequencies in each RT bin
  pf=p*sum(f)
  
  #select RT bins that will be used to
  # calculate G^2. Bins with 0 counts
  # will be skipped to avoid values of -Inf. 
  # (This could only happen if you used fixed cutoffs
  # for the RT bins instead of taking quantiles
  # from the observed RTs)
  gslev=1:length(f)
  gslev=gslev[f>0]
  
  gsq=2*sum(f[gslev]*log(f[gslev]/pf[gslev]))
}

exgauss_fit = function(q, f){
  #starting parameter values for optimization
  par=c(.5,.3,.5)
  
  #optimization routine like above but
  # for quantile-based likelihood
  comp=gsqfun(par,f,q)
  out=0 #just a variable to control when the loop stops
  while(out==0){
    fit=optim(par,fn=gsqfun,f=f,q=q)
    par=fit$par
    if((comp-fit$value)<=.001) out=1
    comp=fit$value
  }
  par3=par
  
  gsq=gsqfun(par3,f,q)
  
  ls = list(par3, gsq)
  return (ls)
}


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

gsqfun_gamma <- function(par, q, freq){
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

gamma_fit = function(q, f){
  par = c(2, 3)
  comp = gsqfun_gamma(par, q, f)
  
  while(TRUE){
    fit=optim(par,gsqfun_gamma,freq=f,q=q)
    par=fit$par
    if((comp-fit$value)<=.001) break
    comp=fit$value
  }
  gamma_par = par
  gsq = gsqfun_gamma(par, q, f)
  
  ls = list(gamma_par, gsq)
  return (ls)
}


cons=function(cpar){
  penf=0
  pmn=.05
  if(cpar[1]<pmn){
    penf=penf+(pmn-cpar[1]) 
    cpar[1]=pmn
  }
  pmn=.1
  if(cpar[2]<pmn){
    penf=penf+(pmn-cpar[2]) 
    cpar[2]=pmn
  }
  #pmn=.01
  #if(cpar[3]<pmn){
  #  penf=penf+(pmn-cpar[3]) 
  #  cpar[3]=pmn
  #}
  ls=list(cpar,penf)
  return(ls)
}



waldPred <- function(par, q){
  
  #Parameters
  a = par[1]
  v = par[2]
  
  #Cuts for quantiles
  cuts = c(0, q, Inf)
  p = c()
  
  for (i in (1:(length(cuts)-1))){
    p[i] = pwald(cuts[i+1], a, v) - pwald(cuts[i], a, v)
  }
  return(p)
}

gsqfun_wald <- function(par, q, freq){
  
  ls = cons(par)
  cpar = ls[[1]]
  pen = ls[[2]]
  p = waldPred(cpar, q)
  p = p/sum(p)
  #print (p)
  
  pf = p*sum(freq)
  pf
  
  gslev=1:length(freq)
  gslev=gslev[freq>0]
  
  gsq=2*sum(freq[gslev]*log(freq[gslev]/pf[gslev]))
  gsq = gsq + pen*gsq
  return (gsq)
}

wald_fit = function(q, f){
  par = c(0.15, 0.3)
  comp = gsqfun_wald(par, q, f)
  
  while(TRUE){
    fit=optim(par,gsqfun_wald,freq=f,q=q)
    par=fit$par
    if((comp-fit$value)<=.001) break
    comp=fit$value
  }
  wald_par = par
  gsq = gsqfun_wald(par, q, f)
  
  ls = list(wald_par, gsq)
  return (ls)
}


cons=function(cpar){
  penf=0
  pmn=.05
  if(cpar[1]<pmn){
    penf=penf+(pmn-cpar[1]) 
    cpar[1]=pmn
  }
  pmn=.1
  if(cpar[2]<pmn){
    penf=penf+(pmn-cpar[2]) 
    cpar[2]=pmn
  }
  #pmn=.01
  #if(cpar[3]<pmn){
  #  penf=penf+(pmn-cpar[3]) 
  #  cpar[3]=pmn
  #}
  ls=list(cpar,penf)
  return(ls)
}



waldPred <- function(par, q){
  
  #Parameters
  a = par[1]
  v = par[2]
  
  #Cuts for quantiles
  cuts = c(0, q, Inf)
  p = c()
  
  for (i in (1:(length(cuts)-1))){
    p[i] = pwald(cuts[i+1], a, v) - pwald(cuts[i], a, v)
  }
  return(p)
}

gsqfun_wald_2con <- function(par, q1, freq1, q2, freq2){
  
  ls = cons(par[1:2])
  cpar = ls[[1]]
  pen = ls[[2]]
  p = waldPred(cpar, q1)
  p = p/sum(p)
  #print (p)
  
  pf = p*sum(freq1)
  pf
  
  gslev=1:length(freq1)
  gslev=gslev[freq1>0]
  
  gsq1=2*sum(freq1[gslev]*log(freq1[gslev]/pf[gslev]))
  gsq1 = gsq1 + pen*gsq1
  
  #ls = cons(c(par[1], par[3]))
  #For fixed condition uncomment below. For v-free parameter, uncomment above
  ls = cons(par[1:2])
  cpar = ls[[1]]
  pen = ls[[2]]
  p = waldPred(cpar, q2)
  p = p/sum(p)
  #print (p)
  
  pf = p*sum(freq2)
  pf
  
  gslev=1:length(freq2)
  gslev=gslev[freq2>0]
  
  gsq2=2*sum(freq2[gslev]*log(freq2[gslev]/pf[gslev]))
  gsq2 = gsq2 + pen*gsq2
  
  gsq = gsq1 + gsq2
  
  
  return (gsq)
}

wald_fit_2con = function(q1, f1, q2, f2){
  par = c(0.15, 0.3, 0.3)
  comp = gsqfun_wald(par, q1, f1, q2, f2)
  #print (comp)
  while(TRUE){
    fit=optim(par,gsqfun_wald,freq1=f1,q1=q1, freq2=f2, q2 = q2)
    # print (fit)
    par=fit$par
    if((comp-fit$value)<=.001) break
    comp=fit$value
  }
  wald_par = par
  gsq = gsqfun_wald(par, q1, f1, q2, f2)
  
  ls = list(wald_par, gsq)
  return (ls)
}
