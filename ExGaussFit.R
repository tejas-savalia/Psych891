
#ExGaussian parameters
mu=.7
sig=.15
tau=.3

#number of simulated data trials
N=300

#Just to get end points for
# nice, consistent, axes on graphs
xmn=0
xmx=qexp(.999,1/tau)+qnorm(.999,mu,sig)
ymn=0
ymx=dnorm(mu,mu,sig)
  
#Simulates (randomly samples) data from an
#  ExGaussian and returns quantiles of the
#  simulated distribution
exgsim=function(mu,sig,tau,N){
  norm=rnorm(N,mu,sig)
  exp=rexp(N,1/tau)
  exg=norm+exp
  hist(exg,30,col="grey",border="grey",freq=F,
       xlim=c(xmn,xmx),ylim=c(ymn,ymx))
  q=quantile(exg,c(.1,.3,.5,.7,.9))
  abline(v=q)
  return(q)
}

#runs the sim function above
q=exgsim(mu,sig,tau,N)
q
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

#Adds a line to a histogram of sampled data
#  showing the true ExGaussian function from
#  which the data were sampled. You have to
#  have a historgram plotted already for this to
#  work.
addtrue=function(mu,sig,tau){
  vq=seq(xmn,xmx,length.out = 500)
  points(vq,dexg(vq,mu,sig,tau),type="l",lwd=2)
}

#runs function above
addtrue(mu,sig,tau)

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

#pred(c(mu, sigma, tau), q)
#computes g-squared
# par is a parameter vector with mu, sig,
# and tau in that order
# dat is an observed set of frequencies
# q is a set of observed quantiles
gsqfun=function(par,dat,q){
  p=pred(par,q) #get predicted proportions
  print(p)
  #If a cell is 0 in the predictions but not in
  # the data, fill in the prediction with a small
  # number. This is a trick to keep the function
  # from returning infinity.
  p[p==0 & dat>0] = .0001
  p=p/sum(p) #renormalize so p adds up to 1 again
  
  pf=c() #predicted frequencies in each RT bin
  pf=p*sum(dat)
  
  #select RT bins that will be used to
  # calculate G^2. Bins with 0 counts
  # will be skipped to avoid values of -Inf. 
  # (This could only happen if you used fixed cutoffs
  # for the RT bins instead of taking quantiles
  # from the observed RTs)
  gslev=1:length(dat)
  gslev=gslev[dat>0]
  
  gsq=2*sum(dat[gslev]*log(dat[gslev]/pf[gslev]))
}

#starting parameter values for optimization
par=c(.5,.25,.5)

#optimization routine.
# It is repeated in a loop bc
# it doesn't always get to the
# global minimum the first time.
# Keeps running until it isn't
# improving the fit anymore.
comp=gsqfun(par,f,q)
out=0 #just a variable to control when the loop stops
while(out==0){
  fit=optim(par,gsqfun,dat=f,q=q)
  par=fit$par
  if((comp-fit$value)<=.001) out=1
  comp=fit$value
}


#Adds a line to a histogram of sampled data
#  showing the best-fitting ExGaussian function.
#  You have to have a historgram plotted already 
#  for this to work.
addfit=function(par){
  mu=par[1]
  sig=par[2]
  tau=par[3]
  vq=seq(xmn,xmx,length.out = 500)
  points(vq,dexg(vq,mu,sig,tau),type="l",col="red",lty=2,lwd=2)
}

#runs function above
addfit(par)

