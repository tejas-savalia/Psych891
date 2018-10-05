#Q 1
#Function to return density of a wald distribution for specific RT value

dwald <- function(rt, mu, lambda){
  sqrt(lambda/(2*pi*rt^3))*exp(-lambda*(rt-mu)/(2*mu^2*rt))
}

#dwald(1, 2, 1)

#library(SuppDists)

#Q2
f = c(45, 67, 23, 88)
n = sum(f)
fprop = f/n
devf =  -2*(dmultinom(f,n,fprop,log=T))
gsq_model = 3.2
dev_model = devf+gsq_model
#18.94523

#Q3
dev_model = 140
ll = dev_model/-2
#70

#Q4
p_a = 9/(9+1) 
p_b = 1/(9+1)

#Probability that the data came from model A is 0.9. 
#Probability that the data came from model B is 0.1.

#Q5

MA_AIC = 450
MB_AIC = 452
MC_AIC = 446

min_aic = 446

del_AIC_A = 4
del_AIC_B = 6
del_AIC_C = 0

denom = exp(-del_AIC_A/2) + exp(-del_AIC_B/2) + exp(-del_AIC_C/2)

w_AIC_A = exp((-del_AIC_A/2))/denom
w_AIC_B = exp((-del_AIC_B/2))/denom
w_AIC_C = exp((-del_AIC_C/2))/denom

#AIC Weights of models A, B and C are 0.114, 0.04, 0.843 respectively

#Q6

gdiff = 18.1-12.3
p = pchisq(gdiff, df=3, lower.tail = F)
