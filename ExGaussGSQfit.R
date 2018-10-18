# before running this, you need a "dat"
#  vector with a set of observed (or simulated)
#  RTs.


# returns quantiles of a distribution
getqs=function(dat){
  q=quantile(dat,c(.1,.3,.5,.7,.9))
  return(q)
}

#q=getqs(dat)
q = q_exg
#Creates a vector of observed frequency
#  counts in each RT bin separated by
#  quantiles from the exgsim function
makefreqs=function(q,N){
  freqs=c(.1,.2,.2,.2,.2,.1)*N
  return(freqs)
}

#runs the function above
f=makefreqs(q,N)

#cdf for an ExGaussian
# (You don't have to understand all this math)
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
