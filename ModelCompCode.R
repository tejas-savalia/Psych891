
#quantiles from two conditions (you can change these
#  to see what happens if you like)
q1=c(0.6417067,0.8118283,0.967615,1.235573,1.807976)
q2=c(0.6579347,0.8144716,0.9833921,1.270491,1.811078)

#frequencies in RT ins from two conditions (you can change the N's
#  to see what happens if you like)
N1 = 1000
N2 = 1000
f1=c(.1,.2,.2,.2,.2,.1)*N1
f2=c(.1,.2,.2,.2,.2,.1)*N2


#Run the fit for the more flexible model (mu can differ between conditions)
source("Exg2conMufree.R")
ls=MuFreeFit(f1,q1,f2,q2)
par1=ls[[1]]
gsq1=ls[[2]]

#Run the fit for the less flexible model (mu cannot differ between conditions)
source("Exg2conAllFixed.R")
ls=AllFixedFit(f1,q1,f2,q2)
par2=ls[[1]]
gsq2=ls[[2]]


#g-squared test
gdiff=gsq2-gsq1
p=1-pchisq(gdiff,1)

# deviance of the "full" model with as many
#  df as the data. Sum of the deviance from
#  each condition
fpred=c(.1,.2,.2,.2,.2,.1)
devf = -2*(dmultinom(fpred*N1,N1,fpred,log=T)+
             dmultinom(fpred*N2,N2,fpred,log=T))

#gsq equals model deviance minus deviance of the full model
#  so model deviance equals full model deviance plus gsq
dev1=gsq1+devf
dev2=gsq2+devf

ddiff=dev2-dev1
#difference in deviance is the same as the difference
# in g-squared, so you can do the same hypothesis test
# as above using the deviance

# log likelihoods
ll1 = dev1/-2
ll2 = dev2/-2

# likelihood ratio
lr = exp(ll1-ll2)
# or
#lr = exp((-.5*dev1)-(-.5*dev2))

# likelihood weight
lw1 = lr/(lr+1)
lw2 = 1/(lr+1)


aic1 = dev1 + 2*4
aic2 = dev2 + 2*3

aicc1 = dev1 + 2*4 + ((2*4*(4+1))/((N1+N2)-4-1))
aicc2 = dev2 + 2*3 + ((2*3*(3+1))/((N1+N2)-3-1))


#AIC ratio
aicr = exp((-.5*aic1)-(-.5*aic2))
#or
# aicr = exp((-.5*(aic1-aic2)))


#AIC weights
aicw1 = aicr/(1+aicr)
aicw2 = 1/(1+aicr)

# same using W & F (2004) Equation 4
 # Difference from minimum aic
 daic1 = aic1 - min(aic1,aic2)
 daic2 = aic2 - min(aic1,aic2)

 w1 = exp(-.5*daic1) / (exp(-.5*daic1) + exp(-.5*daic2))
 w2 = exp(-.5*daic2) / (exp(-.5*daic1) + exp(-.5*daic2))
