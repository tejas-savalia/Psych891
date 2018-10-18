

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



#f1 and f2 are frequencies between RT bins from
# 2 conditions
#q1 and q2 are quantiles from 2 conditions

AllFixedFit=function(f1,q1,f2,q2){

  # function to apply constraints to parameter
  #  values during fit to avoid illegal values.
  #  Also computes and returns penalty factor (penf)
  #  that is larger for attempted par values
  #  that go farther into the illegal range
  #  Works for a single condition.
  cons=function(cpar){
    penf=0
    pmn=.01
    if(cpar[1]<pmn){
      penf=penf+(pmn-cpar[1]) 
      cpar[1]=pmn
    }
    pmn=.001
    if(cpar[2]<pmn){
      penf=penf+(pmn-cpar[2]) 
      cpar[2]=pmn
    }
    pmn=.01
    if(cpar[3]<pmn){
      penf=penf+(pmn-cpar[3]) 
      cpar[3]=pmn
    }
    ls=list(cpar,penf)
    return(ls)
  }
  
  
  
  #Returns predicted proportion of scores
  # in each RT bin separated by a set of
  # quantiles.
  # Works for a single condition
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
  gsqfun=function(par,f,q){
    #CONDITION 1
    #f=f1
    #q=q1
    ls=cons(par[1:3])
    cpar=ls[[1]]
    penf=ls[[2]]
    p=pred(cpar,q) #get predicted proportions
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
    #add in penalty for illegal parameter values
    gsq = gsq + penf*gsq
    
    #CONDITION 2
    #f=f2
    #q=q2
    #ls=cons(par[1:3])
    #cpar=ls[[1]]
    #penf=ls[[2]]
    #p=pred(cpar,q) #get predicted proportions
    #If a cell is 0 in the predictions but not in
    # the data, fill in the prediction with a small
    # number. This is a trick to keep the function
    # from returning infinity.
    #p[p==0 & f>0] = .0001
    #p=p/sum(p) #renormalize so p adds up to 1 again
    #pf=c() #predicted frequencies in each RT bin
    #pf=p*sum(f)
    #select RT bins that will be used to
    # calculate G^2. Bins with 0 counts
    # will be skipped to avoid values of -Inf. 
    # (This could only happen if you used fixed cutoffs
    # for the RT bins instead of taking quantiles
    # from the observed RTs)
    #gslev=1:length(f)
    #gslev=gslev[f>0]
    #gsq2=2*sum(f[gslev]*log(f[gslev]/pf[gslev]))
    #add in penalty for illegal parameter values
    #gsq2 = gsq2 + penf*gsq2
    
    #gsq=gsq1+gsq2
    return(gsq)
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
    fit=optim(par,gsqfun,f1=f,q1=q)
    ls=cons(fit$par)
    par=ls[[1]]
    if((comp-fit$value)<=.01) out=1
    comp=fit$value
  }
  
  ls=list(par,comp)
  return(ls)
}
