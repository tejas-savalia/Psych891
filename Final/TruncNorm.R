
mu=6
sig=5
lo=3
hi=30

dtnorm=function(x,mu=0,sig=1,lo=-Inf,hi=Inf){
  den=rep(0,length(x))
  den[x>=lo & x<=hi] = dnorm(x,mu,sig)/(pnorm(hi,mu,sig)-pnorm(lo,mu,sig))
  return(den)
}

ptnorm=function(x,mu=0,sig=1,lo=-Inf,hi=Inf){
  cden=rep(0,length(x))
  cden[x>=lo & x<=hi] = (pnorm(x,mu,sig)-pnorm(lo,mu,sig))/(pnorm(hi,mu,sig)-pnorm(lo,mu,sig))
  cden[x>hi]=1
  return(cden)
}

qtnorm=function(p,mu=0,sig=1,lo=-Inf,hi=Inf){
  q=rep(0,length(p))
  #corresponding probability on standard normal
    #proportion of standard normal perserved after truncation
    pkeep=pnorm(hi,mu,sig)-pnorm(lo,mu,sig)
    #proportion between lo cut and q on a nontruncated nornal
    ps=p*pkeep
    #total proportion below q on a nontruncated normal
    pstot=ps+pnorm(lo,mu,sig)
    #q
    q=qnorm(pstot,mu,sig)
  return(q)
}

rtnorm=function(N,mu=0,sig=1,lo=-Inf,hi=Inf){
  p=runif(N)
  x=qtnorm(p,mu,sig,lo,hi)
  return(x)
}

sam=rtnorm(10000,mu,sig,lo,hi)
hist(sam,40,freq=F)

xx=seq(lo,hi,length.out=400)
points(xx,dtnorm(xx,mu,sig,lo,hi),type='l')


