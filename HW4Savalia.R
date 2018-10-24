
#setwd("C:/Users/Jeff/Documents/RT seminar/")

nsim=1000 #number of simulated data sets to fit


#True ExGaussian parameters
#mu=.6
#sig=.1
#tau=.5
#mu2=.61


true_mu = runif(nsim, 0.8, 1.2)
true_sig = runif(nsim, 0.1, 0.2)
true_tau = runif(nsim, 0.4, 0.8)

true_shape = runif(nsim, 5, 10)
true_scale = runif(nsim, 0.2, 0.35)

#number of trials per condition in each simulated data set
N=500 #condition 1
#N2=1000 #condition 2

#Simulates (randomly samples) data from an
#  ExGaussian and returns quantiles of the
#  simulated distribution
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

#source("HW2Savalia.R")
#pars
#runs the sim function above
exg = exgsim(true_mu[1],true_sig[1],true_tau[1],N)
gam = gammaSim(N, true_shape[1], true_scale[1])
q_exg = c()
f_exg = c()

q_gam = gam[1]
f_gam = gam[2]

fmu=c() #vectors of fit values from simulation (all pars from mu free model)
fsig=c()
ftau=c()
fmu2=c()
fgsq=c() # for the morel with mu free across conds
fgsq2=c() # for the model with all pars fixed across conds


f_exg_gsq = c()
f_gam_gsq = c()

fkap = c() #kappa 
fth = c()  #theta  

cnt=0
upfreq=20 #frequency of updates when it is running


#######################################################################################################
#Functions to fit exgauss. Till the next hashes.


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

#########################################################################################################
#Functions to fit gamma from hw 2
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

######################################################################################################

#Running them loops


#si= 1

upfreq = 40
for(si in 1:nsim){

  
  #runs the sim function above
  exg=exgsim(true_mu[si],true_sig[si],true_tau[si],N)
  q_exg[si]=exg[1]
  f_exg[si] = exg[2]
  
  #Creates a vector of observed frequency
  #  counts in each RT bin separated by
  #  quantiles from the exgsim function
  
  ls=exgauss_fit(q_exg[[si]], f_exg[[si]])
  #return fit parameters and gsq values. Modify the function in exgaussgsqfit.r
  
  fpar=ls[[1]]
  gval=ls[[2]]
  
  fmu[si]=fpar[1]
  fsig[si]=fpar[2]
  ftau[si]=fpar[3]

  f_exg_gsq[si]=gval
  
  
  
  ls=gamma_fit(q_exg[[si]], f_exg[[si]])
  #return fit parameters and gsq values. Modify the function in hw2
  fpar=ls[[1]]
  gval=ls[[2]]
  
  fkap[si] = fpar[1]
  fth[si] = fpar[2]
  
  f_gam_gsq[si]=gval
  
  #show progress updates
  cnt=cnt+1
  if(cnt%%upfreq==0) print(cnt)
  
}



par(mfrow=c(2,2),mar=c(5,5,2,2))
hist(f_exg_gsq,50,freq=F, main = "G_squared of exg fits")
 vg=seq(0,max(f_exg_gsq),length.out=1000)
 points(vg,dchisq(vg,2),type='l',lwd=2)
hist(f_gam_gsq,50,freq=F, main = "G_squared of gamma fits")
 #vg=seq(0,max(f_gam_gsq),length.out=1000)
 #points(vg,dchisq(vg,2),type='l',lwd=2)
hist(f_gam_gsq-f_exg_gsq,50,freq=F, main = "G_squared of difference of fits")
 #vg=seq(0,max(fgsq2),length.out=1000)
 #points(vg,dchisq(vg,1),type='l',lwd=2)
 
plot(true_mu, fmu, pch = 1, col = 'black', xlim = c(-1, 2), ylim = c(-1, 2))
points(true_sig, fsig, pch = 4, col = 'blue')
points(true_tau, ftau, pch = 2, col = 'red')



#Q3
#The first plot shows g_squared fits being similar to a chi_squared with 2 degrees of freedom; 
#this is expected because ex_gauss has 3 free parameters and the data has 5.
#The g_squared in the second one is fairly high; this is also expected because gamma is not the distribution
#that generated the data
#The difference of fit plot shows high values too; owing to low g_squared for exgauss and high for
#gamma distributions
#The scatter plot shows high correlation between the true and fit parameters. This is expected because 
# the fit parameters are for the same distribution that generated the data.


#Q4

q_exg = c()
f_exg = c()

q_gam = gam[1]
f_gam = gam[2]

fmu=c() #vectors of fit values from simulation (all pars from mu free model)
fsig=c()
ftau=c()
fmu2=c()
fgsq=c() # for the morel with mu free across conds
fgsq2=c() # for the model with all pars fixed across conds


f_exg_gsq = c()
f_gam_gsq = c()

fkap = c() #kappa 
fth = c()  #theta  

cnt=0


upfreq = 40
for(si in 1:nsim){
  
  
  #runs the sim function above
  gam = gammaSim(N, true_shape[si], true_scale[si])
  q_gam[si]=gam[1]
  f_gam[si] = gam[2]
  

  ls=exgauss_fit(q_gam[[si]], f_gam[[si]])
  #return fit parameters and gsq values. Modify the function in exgaussgsqfit.r
  
  fpar=ls[[1]]
  gval=ls[[2]]
  
  fmu[si]=fpar[1]
  fsig[si]=fpar[2]
  ftau[si]=fpar[3]
  
  f_exg_gsq[si]=gval
  
  
  
  ls=gamma_fit(q_gam[[si]], f_gam[[si]])
  #return fit parameters and gsq values. Modify the function in hw2
  fpar=ls[[1]]
  gval=ls[[2]]
  
  fkap[si] = fpar[1]
  fth[si] = fpar[2]
  
  f_gam_gsq[si]=gval
  
  #show progress updates
  cnt=cnt+1
  if(cnt%%upfreq==0) {
    print(cnt)
    #print(q_gam[si])
  }
}

par(mfrow=c(2,2),mar=c(5,5,2,2))
hist(f_exg_gsq,50,freq=F, main = "G_squared of exg fits")
#vg=seq(0,max(f_exg_gsq),length.out=1000)
#points(vg,dchisq(vg,3),type='l',lwd=2)

hist(f_gam_gsq,50,freq=F, main = "G_squared of gamma fits")
vg=seq(0,max(f_gam_gsq),length.out=1000)
points(vg,dchisq(vg,3),type='l',lwd=2)

hist(f_exg_gsq-f_gam_gsq,50,freq=F, main = "G_squared of difference of fits")
#vg=seq(0,max(fgsq2),length.out=1000)
#points(vg,dchisq(vg,1),type='l',lwd=2)
plot(true_shape, fkap, pch = 1, col = 'black', xlim = c(0, max(true_shape)), ylim = c(0, max(fkap)))

points(true_scale,fth, pch = 4, col = 'blue')



#Q5
#The ex-gaussian fits decent here as evident by high density at lower g-squared values. 
#Expected because ex-gaussian has more free parameters than the data generated from gamma distribution (3 vs 2)

#The gamma distribution seems to fit well to its own data too. The density aligns with the bold
#line for a chi-squared with 3 degrees of freedom indicating gamma's 2 parameters subtracted from data's
#5 degrees of freedom

#The difference of fit is fairly low, indicative of ex-gaussian fitting just as well as gamma.
#Eye-balling it, there *does* seem to be a higher density on positive values, indicating gamma fits 
#a little better; not sure if significantly better, though.

#Since Gammafit fits well, the parameters are highly correlated as well. 

#Q6
#The ex-gaussian fits the data generated from gamma better than the other way around. 
#It's indicated by the shape of the g_squared distributions; exg is much dense near low g_sq values
#for data generated from gamma distribution than the other way around. The difference of fit is also
#centered around zero.

#Q7
devA = 550
devB = 545
devC = 552
devD = 541

fA = 3
fB = 4
fC = 3
fD = 5

N = 1000

#bic = dev + vlogN

bicA = devA + fA * log(N)
bicB = devB + fB * log(N)
bicC = devC + fC * log(N)
bicD = devD + fD * log(N)

wAB = (devA/devB)*(N^((fB-fA)/2))
wBC = (devB/devC)*(N^((fC-fB)/2))
wCD = (devC/devD)*(N^((fD-fC)/2))
wDA = (devD/devA)*(N^((fA-fD)/2))

#BIC Values are 570.72, 572.63, 572.72, 575.54 respectively
#Relative BIC weights are 31.91, 0.031, 1020.33, 0.0009 respectively. 
#(Did not compute for all possible 10 combinations)


##########################################################################################

setwd('C:\\Users\\cemnlcmaplab\\Documents\\Psych 891 Modelling cog dynamics\\Codes\\Midterm')


load('Exp1.RData')
N = length(e1dat[1,]) #Number of trials
n_participants = length(e1dat[, 1]) #Number of participants
f = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*N

f_exg_mu_participant_1 = c()
f_exg_sig_participant_1 = c()
f_exg_tau_participant_1 = c()
f_exg_gsq_participant_1 = c()

for (i in 1:n_participants){
  q = quantile(e1dat[i,], c(0.1, 0.3, 0.5, 0.7, 0.9))
  ls = exgauss_fit(q, f)
  fpar = ls[[1]]
  gval = ls[[2]]
  
  f_exg_mu_participant_1[i] = fpar[1]
  f_exg_sig_participant_1[i] = fpar[2]
  f_exg_tau_participant_1[i] = fpar[3]
  
  f_exg_gsq_participant_1[i] = gval
  
}


f_gam_kap_participant_1 = c()
f_gam_th_participant_1 = c()

f_gam_gsq_participant_1 = c()

for (i in 1:n_participants){
  q = quantile(e1dat[i,], c(0.1, 0.3, 0.5, 0.7, 0.9))
  ls = gamma_fit(q, f)
  fpar = ls[[1]]
  gval = ls[[2]]
  
  f_gam_kap_participant_1[i] = fpar[1]
  f_gam_th_participant_1[i] = fpar[2]
  
  f_gam_gsq_participant_1[i] = gval
  
}




