
dc=.1 #diffusion coefficient
a=.4
v=.3
N=1000

#Simulates data from Wald dist. by
# running a 1-bound random walk on 
# each simulated trial
# This one is slow to run, but much
# more interpretable in terms of how
# the RTs are being generated for each
# trial.
waldSim2=function(a,v,N){
  #time steps for random walk (in seconds)
  ts=.001
  #scale average evidence step
  # and SD in evidence steps based
  # on time steps
  dav=v*ts #average change in position of evidence accumulation process
  dsd=sqrt(dc^2*ts) #standard deviation in position changes
  
  rt=c()
  
  for(trl in 1:N){
    out=0
    pos=0 #position of evidence acc. process
    cnt=0 #counter for time steps
    while(out==0){
      pos=pos+rnorm(1,dav,dsd)
      cnt=cnt+1
      if(pos>a) out=1
    }
    rt[trl]=cnt*ts
  }
  return(rt)
}


#Simulates data from Wald dist. 
# This one uses a sampling trick
# that is much faster but hard to
# interpret. 
waldSim=function(a,v,N){
  mu=a/v
  lambda=a^2/dc^2
  y=rnorm(N)^2
  x=mu + ((mu^2*y)/(2*lambda)) - ((mu/(2*lambda))*sqrt(4*mu*lambda*y+mu^2*y^2))
  z=runif(N)
  rt=c()
  rt[z<=(mu/(mu+x))] = x[z<=(mu/(mu+x))]
  rt[z>(mu/(mu+x))] = mu^2/x[z>(mu/(mu+x))]
  return(rt)
}  

#Simulates data from Wald dist. 
# Returns quantiles isntead of all RTs 
waldSimq=function(a,v,N){
  mu=a/v
  lambda=a^2/dc^2
  y=rnorm(N)^2
  x=mu + ((mu^2*y)/(2*lambda)) - ((mu/(2*lambda))*sqrt(4*mu*lambda*y+mu^2*y^2))
  z=runif(N)
  rt=c()
  rt[z<=(mu/(mu+x))] = x[z<=(mu/(mu+x))]
  rt[z>(mu/(mu+x))] = mu^2/x[z>(mu/(mu+x))]
  q=quantile(rt,c(.1,.3,.5,.7,.9))
  return(q)
}  


dwald2=function(t,mu,lambda) sqrt(lambda/(2*pi*t^3)) * exp((-lambda*(t-mu)^2)/(2*mu^2*t))

dwald=function(t,a,v){
  mu=a/v
  lambda=a^2/dc^2
  d=sqrt(lambda/(2*pi*t^3)) * exp((-lambda*(t-mu)^2)/(2*mu^2*t))
  return(d)
}


pwald=function(t,a,v){
  dwald1=function(t){
    mu=a/v
    lambda=a^2/dc^2
    d=sqrt(lambda/(2*pi*t^3)) * exp((-lambda*(t-mu)^2)/(2*mu^2*t))
    return(d)
  }
  p=rep(0,length(t))
  for(i in 1:length(t)){
    if(t[i]>0) p[i]=integrate(dwald1,0,t[i])$val
  }
  return(p)
}



vt=seq(0,10,length.out=300)


mu=a/v
lambda=a^2/dc^2

#plot(vt,dwald2(vt,mu,lambda),type='l',lwd=5)
#hist(waldSim2(a,v,N),30,freq=F)
hist(waldSim(a,v,N),30,freq=F)
points(vt,dwald(vt,a,v),type='l',lty=3,col="red",lwd=3)


