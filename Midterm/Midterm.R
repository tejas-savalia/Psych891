setwd('C:\\Users\\cemnlcmaplab\\Documents\\Psych 891 Modelling cog dynamics\\Codes\\Midterm')
source('wald.R')
source('TruncNorm.R')

nsim = 2000 
N = 250 #Number of Trials


#True parameters
#Boundary height from a truncated normal
a = rtnorm(nsim, mu=0.15, sig = 0.05, lo = 0.05, hi = 0.30) 
#Drift rate from a truncated normal
v = rtnorm(nsim, mu = 0.3, sig = 0.1, lo = 0.1, hi = 0.5)

q_wald = c()
#Expected frequencies
f = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*N

#Fit parameters to save
f_a = c()
f_v = c()

f_wald_gsq = c()


cnt = 0
upfreq = 50


################################################################################

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

######################################################################################################

#Running them loops


si= 1
cnt = 0
upfreq = 50
for(si in 1:nsim){
  
  
  #runs the sim function above
  q_wald=waldSimq(a[si], v[si], N)
  f_wald = f
  
  #Creates a vector of observed frequency
  #  counts in each RT bin separated by
  #  quantiles from the exgsim function
  
  ls=wald_fit(q_wald, f_wald)
  #return fit parameters and gsq values. 
  
  fpar=ls[[1]]
  gval=ls[[2]]
  
  f_a[si]=fpar[1]
  f_v[si]=fpar[2]
  
  f_wald_gsq[si]=gval
  
  
  
  
  #show progress updates
  
  if(cnt%%upfreq==0) print(cnt)
  cnt=cnt+1
}



##########################################################################################
#Plots
par(mfrow = c(2, 2))
hist(f_wald_gsq, 100, freq = F, main = "G_sq fits")
vg = seq(0, max(f_wald_gsq), length.out = 2000)
points(vg, dchisq(vg, 3), type = 'l', lwd = 2)

plot(a, f_a)
abline(a = 0, b = 1, col = 'blue', lwd = 2)

plot(v, f_v)
abline(a = 0, b = 1, col = 'blue', lwd = 2)


#########################################################################
#clear environment for the next question
rm(list = ls())

#########################################################################



setwd('C:\\Users\\cemnlcmaplab\\Documents\\Psych 891 Modelling cog dynamics\\Codes\\Midterm')
source('wald.R')
source('TruncNorm.R')
nsim = 2000 
N = 250 #Number of Trials

#True parameters
#Boundary height from a truncated normal
a = rtnorm(nsim, mu=0.15, sig = 0.05, lo = 0.05, hi = 0.30) 
#Drift rate from a truncated normal
v1 = rtnorm(nsim, mu = 0.3, sig = 0.1, lo = 0.1, hi = 0.5)
v_ef = rtnorm(nsim, mu = 0.05, sig = 0.025, lo = 0.01, hi = 0.1)
v2 = v1 + v_ef

q1_wald = c()
q2_wald = c()

#Expected frequencies
f1 = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*N
f2 = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*N

#Fit parameters to save
f_a = c()
f_v1 = c()
f_v2 = c()

f_wald_gsq = c()


cnt = 0
upfreq = 50


################################################################################

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

gsqfun_wald <- function(par, q1, freq1, q2, freq2){
  
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

wald_fit = function(q1, f1, q2, f2){
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

######################################################################################################

#Running them loops


si= 1
cnt = 0
upfreq = 50
for(si in 1:nsim){
  
  
  #runs the sim function above
  q1_wald=waldSimq(a[si], v1[si], N)
  q2_wald = waldSimq(a[si], v2[si], N)
  f1_wald = f1
  f2_wald = f2
  #Creates a vector of observed frequency
  #  counts in each RT bin separated by
  #  quantiles from the exgsim function
  
  ls=wald_fit(q1 = q1_wald, f1 = f1_wald, q2 = q2_wald, f2 = f2_wald)
  #return fit parameters and gsq values. 
  
  fpar=ls[[1]]
  gval=ls[[2]]
  
  f_a[si]=fpar[1]
  f_v1[si]=fpar[2]
  f_v2[si] = fpar[3]
  f_wald_gsq[si]=gval
  
  
  
  
  #show progress updates
  
  if(cnt%%upfreq==0) print(cnt)
  cnt=cnt+1
}



##########################################################################################
#Plots
par(mfrow = c(2, 2))
hist(f_wald_gsq, 100, freq = F, main = "G_sq fits")
vg = seq(0, max(f_wald_gsq), length.out = 2000)
points(vg, dchisq(vg, 7), type = 'l', lwd = 2)

plot(a, f_a)
abline(a = 0, b = 1, col = 'blue', lwd = 2)

plot(v1, f_v1, pch = 1, col = 'red')
abline(a = 0, b = 1, col = 'blue', lwd = 2)

points(v2, f_v2, pch = 3, col = 'green')
abline(a = 0, b = 1, col = 'blue', lwd = 2)

###########################################################################################
###########################################################################################

load('Exp1.RData')
N = length(e1dat[1,]) #Number of trials
n_participants = length(e1dat[, 1]) #Number of participants
f = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*N

f_wald_a_participant_1 = c()
f_wald_v_participant_1 = c()

f_wald_gsq_participant_1 = c()

for (i in 1:n_participants){
  q = quantile(e1dat[i,], c(0.1, 0.3, 0.5, 0.7, 0.9))
  ls = wald_fit(q, f)
  fpar = ls[[1]]
  gval = ls[[2]]
  
  f_wald_a_participant_1[i] = fpar[1]
  f_wald_v_participant_1[i] = fpar[2]
  
  f_wald_gsq_participant_1[i] = gval

}


load('Exp2.RData')
N = length(e2dat[,,1 ])/4 #Number of trials

f1 = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*N
f2 = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*N

f_wald_a_participant_21 = c()
f_wald_v1_participant_21 = c()
f_wald_v2_participant_21 = c()

f_wald_gsq_participant_21 = c()


n_participants = length(e2dat[1, 1, ])
#Number of participants. Ignoring the conditions for now

#f = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*N

#f_wald_a_participant_21 = c()
#f_wald_v_participant_21 = c()

#f_wald_gsq_participant_21 = c()
i = 1
for (i in 1:n_participants){
  
  q1 = quantile(e2dat[,,i][e2dat[,,i][,1] == 1, 2], c(0.1, 0.3, 0.5, 0.7, 0.9))
  #Quantiles for condition 2
  q2 = quantile(e2dat[,,i][e2dat[,,i][,1] == 2, 2], c(0.1, 0.3, 0.5, 0.7, 0.9))
  
  #q = quantile(e2dat[,,i][,2], c(0.1, 0.3, 0.5, 0.7, 0.9))
  ls = wald_fit(q1, f1, q2, f2)
  fpar = ls[[1]]
  gval = ls[[2]]
  
  f_wald_a_participant_21[i] = fpar[1]
  f_wald_v_participant_21[i] = fpar[2]
  
  f_wald_gsq_participant_21[i] = gval
  
}

########################################################################################

load('Exp2.RData')
N = length(e2dat[,,1])/4 #Number of trials; divided by num conditions(2) and number of columns (2)
n_participants = length(e2dat[1, 1, ])
#Number of participants. Dividing by two for 2 conditions

f1 = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*N
f2 = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*N

f_wald_a_participant_22 = c()
f_wald_v1_participant_22 = c()
f_wald_v2_participant_22 = c()

f_wald_gsq_participant_22 = c()
i = 1
for (i in 1:n_participants){
  #Quantiles for condition 1
  q1 = quantile(e2dat[,,i][e2dat[,,i][,1] == 1, 2], c(0.1, 0.3, 0.5, 0.7, 0.9))
  #Quantiles for condition 2
  q2 = quantile(e2dat[,,i][e2dat[,,i][,1] == 2, 2], c(0.1, 0.3, 0.5, 0.7, 0.9))
  ls = wald_fit(q1, f1, q2, f2)
  #ls = wald_fit(q, f)
  fpar = ls[[1]]
  gval = ls[[2]]
  
  f_wald_a_participant_22[i] = fpar[1]
  f_wald_v1_participant_22[i] = fpar[2]
  f_wald_v2_participant_22[i] = fpar[3]
  
  f_wald_gsq_participant_22[i] = gval
  
}


